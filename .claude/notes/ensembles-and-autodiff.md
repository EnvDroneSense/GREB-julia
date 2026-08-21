# Ensembles / multi-run and differentiability

**Status:** Investigated — costs measured, differentiability partially blocked. Not started.

What was actually tried for ensemble/multi-run support and for making the model
differentiable, with measured costs and the specific blockers found. Read this
before re-opening either topic.

---

This document is a planning/audit reference, not an implementation record —
no code in `src/` changed as part of producing it. It answers, with actual
measurements taken against this repo (commit `a3d443d`, Julia 1.12.6, Windows,
14 logical CPUs), the two questions raised for future work: (1) what would a
multi-run/ensemble mechanism cost and where would it strain, and (2) what,
concretely, blocks differentiating GREB, and how far does today's code already
get you for free. All numbers below were produced by small scratch scripts
run against a `Pkg.develop`-linked copy of this package in an isolated
scratch environment (not committed to this repo); the scripts are reproducible
from the snippets included here.

## Part 1 — Multi-run/ensemble mechanism: measured costs

### Per-instance memory footprint (measured via `Base.summarysize`)

| Struct | Measured size | Dominant contributor |
|---|---|---|
| `ClimateFields` | **744.81 MiB** | 29 `Array{Float64,3}` fields at `(xdim, ydim, nstep_yr)` = (96, 48, 730) = 3,363,840 cells each |
| `CirculationWorkspace` | 1373.27 KiB (≈1.34 MiB) | 38 `Matrix{Float64}` fields at `(xdim, ydim)` = 4,608 cells each |
| `ModelState` | 396.61 KiB | 10 `Matrix{Float64}` `(xdim, ydim)` accumulators |
| `MonthlyAccumulator` | 468.72 KiB | 13 `Matrix{Float64}` `(xdim, ydim)` accumulators |

**Implication for Part 1:** ensemble memory scales almost entirely with how
many `ClimateFields` instances are simultaneously live, not with the small
per-timestep workspaces. Each concurrent ensemble member needs its own
`ClimateFields` (required for isolation — see below), so N-way concurrency
costs roughly `N × 745 MiB` just for that one struct. On a machine with, say,
16 GiB free, that caps naive full-concurrency at ~20 simultaneous members
before other memory pressure (OS, JLD2-loaded climatology if used, GC
headroom) becomes the real constraint — worth knowing before assuming
"spawn one task per config" scales to hundreds of members without also
capping in-flight concurrency (e.g. a semaphore/worker-pool bound rather than
unconditional `Threads.@spawn` per element).

### Measured `Threads.@spawn`-per-config speedup — and a real gotcha found

Setup: 6 configs (`co2_concentration` = 300..400 by 20), each running
`greb_model!(RunSpec(flux=0, ctrl=1, scnr=0), cfg; jld2_dir="", fields=ClimateFields())`
(no real dataset — `jld2_dir=""` uses all-zero climatology, so this measures
mechanism overhead, not real per-year science-run cost). Measured with
`julia -t 4`:

```
sequential (6 runs): 13.009 s   (2.168 s/run)
threaded (Threads.@spawn, 4 threads): 6.361 s
speedup = 2.05x
results identical (sequential vs threaded): true
```

**Correctness confirmed**: sequential and threaded runs produced bit-identical
results, validating that giving each task its own `PhysicsConfig`/
`ClimateFields` (as the mechanism must) achieves real isolation under
concurrency — no cross-talk observed.

**Speedup is sub-linear (2.05x on 4 threads), and the reason is structural,
not noise**: `greb_model!` already allocates three distinct
`CirculationWorkspace` instances internally (`ws`, `ws_a`, `ws_q`) specifically
so `tendencies!` can run its own 2-way `Threads.@spawn` split
(`src/tendencies.jl:19-28`) whenever `Threads.nthreads() > 1`. That gate is
unconditional — it doesn't know or care whether the calling context is
already inside another spawned task. So an ensemble-level `Threads.@spawn`
per config, run under `-t 4`, produces **nested** parallelism: 6 outer tasks
each trying to spawn 2 more inner tasks, all competing for the same 4 OS
threads. Julia's scheduler multiplexes this correctly (no deadlock, no
corruption — confirmed by the identical results above) but the
oversubscription eats into the naive 4x-on-4-threads expectation. **This is
worth designing around explicitly** in any real Part 1 implementation — e.g.
by having the ensemble runner pass `ws_a=ws_q` (forcing `tendencies!`'s
sequential branch) for the inner call when outer-level concurrency is already
in play, or by capping outer concurrency well below `Threads.nthreads()`.
Neither is implemented here; this is a scoping note.

### A real (pre-existing) hazard found while measuring this: `redirect_stdout` is not concurrency-safe

Both the existing test suite (`test/runtests.jl`) and the benchmark harness
use `redirect_stdout(devnull) do ... end` to silence `greb_model!`'s
`println` calls — safe today because tests run sequentially. While measuring
the timing above, wrapping each ensemble task's `greb_model!` call in
`redirect_stdout(devnull)` and running it under `Threads.@spawn` caused
**output from some tasks to silently vanish**, including — critically — the
script's own subsequent `println` calls for the timing results themselves,
which never appeared even though the process exited normally (exit code 0).
`redirect_stdout` mutates a single global/task-shared stream target; two
tasks concurrently entering/exiting their `do` blocks race on save-and-restore
and can leave the "real" stdout redirected to `devnull` even after both
`do`-blocks have returned. **This isn't an AD or physics issue — it's a
concurrency hazard in the existing print-suppression idiom that any Part 1
implementation must not reuse verbatim inside spawned tasks.** (Removing the
`redirect_stdout` wrapper — accepting console noise, or capturing per-task
output into a private `IOBuffer` instead — reproduced clean, reliable timing
immediately.) Flagging this now because it's exactly the kind of thing that
would otherwise surface as a flaky/nondeterministic ensemble test later.

### The `cfg.co2_scenario` mutation hazard, confirmed by reading the source

`src/model.jl:354` and `:359` do `cfg.co2_scenario = load_co2_scenario_jld2(...)`
/ `load_custom_co2_scenario(...)` for RCP/SSP/historical/custom-CO2
experiments — i.e. `greb_model!` writes back into whatever `PhysicsConfig`
object the caller passed in. Confirmed via direct read of both call sites.
Any ensemble mechanism must guarantee each config gets its own `PhysicsConfig`
instance (not a shared template mutated per member) or two members running
the same CO2-table experiment concurrently could clobber each other's table.

## Part 2 — Differentiability: audit, with things actually tried

### Struct genericity audit — exact field counts (all hardcoded to `Float64`)

| Struct | Total fields | Array-typed (`Vector{Float64}`/`Matrix{Float64}`/`Array{Float64,3}`) fields | Non-array fields |
|---|---|---|---|
| `CirculationWorkspace` (`src/state.jl:8-71`) | 42 | 42 (4 `Vector`, 38 `Matrix`) | 0 |
| `ClimateFields` (`src/state.jl:240-292`) | 38 | 38 (9 `Matrix`, 29 `Array{·,3}`) | 0 |
| `MonthlyAccumulator` (`src/state.jl:152-167`) | 14 | 13 `Matrix` | 1 (`count::Int`) |
| `SurfaceState` (`src/state.jl:138-143`) | 4 | 4 `Matrix` | 0 |
| `ModelState` (`src/state.jl:335-350`) | 11 | 10 `Matrix` | 1 (`sw_solar_forcing::Float64`) |
| `MonthlyRecord` (`src/state.jl:364`, `NamedTuple`) | 13 | 13 `Matrix{Float64}` | 0 |
| `PhysicsConfig` (`src/config.jl:10-73`) | 27 | 0 arrays; 6 `Float64` scalars (`co2_concentration`, `earth_sun_distance_pct`, `c_q`, `c_rq`, `c_omega`, `c_omegastd`) + `co2_scenario::Dict{Int,Float64}` | 20 `Bool`/`Int`/`Symbol`/`String` |

**Total: 120 array-typed fields across 6 structs**, every one hardcoded to
`Float64`, zero of them generic. Every constructor (`CirculationWorkspace()`,
`ClimateFields()`, etc.) calls `zeros(Float64, ...)`/`ones(Float64, ...)`
literally.

Every function that takes one of these structs as an argument does so with a
concrete type annotation. Counted by grep across `src/`:

```
src/circulation.jl:8   src/config.jl:3      src/io.jl:5
src/model.jl:7         src/output.jl:7      src/state.jl:2
src/tendencies.jl:4    src/physics/hydrology.jl:2
src/physics/ocean.jl:4 src/physics/radiation.jl:4
```
**46 concrete-type annotations across 10 files** would need to become
type-parameter-friendly (`ws::CirculationWorkspace` → `ws::CirculationWorkspace{T}`,
or loosened to accept any `T`) for a value to flow through as anything other
than `Float64`. This is the scope of the "genericize the structs" prerequisite
named in the original audit request — now with an exact count instead of an
estimate.

### Empirically confirmed: exactly where a `Dual` gets rejected today

Ran `ForwardDiff.Dual` values against the actual structs and functions (not
hypothetically — actual `try`/`catch` probes against this codebase):

| Probe | Result |
|---|---|
| `forcing(1, ForwardDiff.Dual(20.0,1.0), cfg, fields, nothing)` where `cfg.experiment = :co2_sine_wave` | **Succeeds.** `year` has no type annotation in `forcing`'s signature, so a `Dual` flows through its arithmetic (`cos`, `*`, `+`) unmodified. |
| `forcing(1, ForwardDiff.Dual(2000.0,1.0), cfg, fields, nothing)` where `cfg.experiment = :rcp85` | **Succeeds** (returns constant `CO2=340.0`) — this branch never actually touches `year` numerically, so nothing to differentiate, but no error either. |
| `ws.dTa_crcl[1,1] = ForwardDiff.Dual(1.0,1.0)` (assigning into a `CirculationWorkspace` buffer) | **Fails**: `MethodError: no method matching Float64(::ForwardDiff.Dual{...})` — confirms the `Matrix{Float64}` field type is a hard wall, not just a style choice. |
| `cfg.co2_concentration = ForwardDiff.Dual(340.0,1.0)` (assigning into `PhysicsConfig`) | **Fails**, same `MethodError`. Confirms `PhysicsConfig`'s numeric fields are equally hardcoded — differentiating with respect to *any* physical parameter that lives in `PhysicsConfig` (not just the array structs) hits this same wall. |
| `PhysicsConfig(co2_concentration=ForwardDiff.Dual(340.0,1.0))` (via the `@kwdef` keyword constructor) | **Fails**, same error — the keyword constructor doesn't provide an escape hatch. |
| `ForwardDiff.derivative(co2 -> forcing(1, 2000.0, PhysicsConfig(experiment=:full_model, co2_concentration=co2), fields, nothing).CO2, 340.0)` | **Fails** for the same reason, one level removed — even though `forcing`'s fast path for `:full_model` is just `CO2 = cfg.co2_concentration` (a trivial identity), the value has to pass through the `PhysicsConfig` struct's `Float64` field first. |

**Practical reading:** today, you can only differentiate GREB with respect to
a bare scalar argument that a function reads directly (like `forcing`'s
`year`) — the moment a value has to pass through `PhysicsConfig` or any of
the state structs, even trivially, AD breaks. This makes the audited
prerequisite very concrete: **most physically interesting differentiation
targets (any `PhysicsConfig` field: CO2 level, hydrology parameters, orbital
parameters) are blocked today, not just full-field/full-grid gradients.**

### Empirically confirmed: `@turbo` doesn't error on `Dual`, it silently drops SIMD

Contrary to the original assumption that `@turbo` would hard-fail on
non-`Float64`/`Float32` arrays, testing it directly showed:

```
┌ Warning: LoopVectorization.check_args on your inputs failed; running
│ fallback @inbounds @fastmath loop instead.
└ ...
Dual @turbo sum OK: [Dual{Nothing}(0.32,2.0), Dual{Nothing}(0.72,2.0), Dual{Nothing}(0.83,2.0)]
```
`@turbo` type-checks its inputs (`LoopVectorization.check_args`), and when
given `Dual`-typed arrays it **falls back to a plain `@inbounds @fastmath`
loop** rather than throwing — meaning a `@turbo` loop, if it were reached
with `Dual` arrays, would produce a *correct* result, just with a runtime
warning and full loss of SIMD vectorization (the entire performance reason
`@turbo` exists in this codebase). This changes the framing of the original
audit: `@turbo` is not a hard correctness blocker for `ForwardDiff` the way it
would be a hard blocker for Zygote's source-to-source tracing — it's a
**silent performance cliff**. Since GREB relies on `@turbo` across ~35 hot-loop
sites for its actual runtime performance (see `benchmark/` results), any
future full-model AD path through those loops would need to either accept
that performance cliff for the differentiated code path, or provide an
explicit non-`@turbo` fallback kernel that's exercised only when the element
type isn't a hardware-native float. Exact `@turbo` site count, verified by
grep: **35** (`src/circulation.jl`: 15, `src/physics/hydrology.jl`: 9,
`src/physics/radiation.jl`: 4, `src/physics/ocean.jl`: 3, `src/output.jl`: 2,
`src/model.jl`: 1, `src/state.jl`: 1).

### Memory cost of genericizing, quantified

```
sizeof(Float64)                              =  8 bytes
sizeof(ForwardDiff.Dual{Nothing,Float64,1})  = 16 bytes  (2.0x)
sizeof(ForwardDiff.Dual{Nothing,Float64,4})  = 40 bytes  (5.0x)
```
Forward-mode AD's per-element cost scales linearly with the number of
simultaneously-differentiated scalar inputs (the `Dual`'s "tag width"). Applied
to `ClimateFields`'s measured 744.81 MiB:
- Differentiating w.r.t. **1** scalar parameter (`Dual{...,1}`, what the
  minimal example below uses): **≈1.46 GiB** per `ClimateFields` instance.
- Differentiating w.r.t. **4** scalar parameters simultaneously (a small
  forward-mode Jacobian, e.g. `co2_concentration` + 3 hydrology
  coefficients at once): **≈3.64 GiB** per instance.

This is why forward-mode AD (`ForwardDiff`) is a reasonable choice for "a
handful of scalar parameters" but scales badly toward "gradient with respect
to every grid cell of a spatial field" (thousands of tag dimensions) — that
regime is what reverse-mode tools (`Enzyme`) exist for, at the cost of the
mutation/`@turbo` compatibility questions above being harder to answer without
a dedicated investigation.

### Minimal validated gradient example — actually run, with numbers

Target: `forcing(it, year, cfg, fields, icmn_ctrl)` (`src/tendencies.jl:83`),
`:co2_sine_wave` branch (`tendencies.jl:148-149`):
`CO2 = 340.0 + 170.0 + 170.0 * cos(2π * (year - 13.0) / 30.0)`. Chosen because
it's the one path confirmed above to accept a `Dual` with zero source changes.

```julia
cfg = create_experiment_config(:full_model); cfg.experiment = :co2_sine_wave
fields = ClimateFields()
f(year) = forcing(1, year, cfg, fields, nothing).CO2
ad = ForwardDiff.derivative(f, 20.0)
```

Actual measured results at `year = 20.0`:

| Method | Value | \|diff from AD\| |
|---|---|---|
| `ForwardDiff.derivative` | **−35.40967037699586** | — |
| Central finite difference, h=1e-2 | −35.40964448964132 | 2.59×10⁻⁵ |
| Central finite difference, h=1e-4 | −35.40967037451992 | **2.48×10⁻⁹** |
| Central finite difference, h=1e-6 | −35.40967043136334 | 5.44×10⁻⁸ (FD noise floor reached) |
| Closed-form analytic (`−170·sin(2π(year−13)/30)·(2π/30)`) | −35.409670376995855 | **7.11×10⁻¹⁵** (machine epsilon) |

The AD result matches the closed-form analytic derivative to machine
precision and the finite-difference approximation converges to it as `h`
shrinks (until FD's own truncation/round-off tradeoff floor around `h=1e-6`,
the classic FD noise pattern) — as strong a validation as a single-function
example can give. No source changes to `GREB.jl` were needed to get this
result.

## Summary: what's already true vs. what's blocked

| | Status |
|---|---|
| Differentiate a bare scalar argument a function reads directly (e.g. `forcing`'s `year`) | **Works today, zero changes**, validated above |
| Differentiate w.r.t. any `PhysicsConfig` field (CO2 level, hydrology params, orbital params) | **Blocked** — `MethodError` on `Float64` field assignment, confirmed empirically |
| Differentiate w.r.t. any value that must pass through a state struct (`CirculationWorkspace`, `ClimateFields`, etc.) | **Blocked** — same `Float64`-field wall, confirmed empirically |
| Differentiate through a `@turbo` loop, once the eltype issue above is fixed | Would **run correctly but silently lose SIMD** (confirmed — falls back to `@inbounds @fastmath`, not an error) |
| Ensemble mechanism: isolation via per-task `PhysicsConfig`/`ClimateFields` | **Works** — measured bit-identical results, sequential vs. 4-thread concurrent |
| Ensemble mechanism: linear speedup with thread count | **Does not hold as-is** — nested parallelism with `tendencies!`'s internal 2-way spawn caps realized speedup (measured 2.05x on 4 threads for 6 tasks); needs explicit handling |
| Ensemble mechanism: reusing `redirect_stdout(devnull)` inside spawned tasks | **Unsafe today** — measured output corruption; needs a different suppression strategy for concurrent use |

## What full-model differentiability would additionally require

1. Genericize the 6 structs / 120 array fields / 46 call-site annotations
   above to a type parameter `T` (additive, `T=Float64` default preserves
   today's behavior and performance for existing callers) — the single
   largest, most mechanical prerequisite.
2. Decide the `@turbo` story for the differentiated code path: accept the
   measured SIMD-loss fallback, or maintain a parallel non-`@turbo` kernel
   for non-`Float64` element types (extra code to keep in sync with the
   `@turbo` originals — a real maintenance cost, not just a one-time change).
3. Force `tendencies!`'s sequential branch (`ws_a === ws_q`) when
   differentiating — untested here, but consistent with "most AD frameworks
   don't safely trace across spawned tasks."
4. A story for `Vector{MonthlyRecord}`'s dynamic, heap-growing accumulation
   across the multi-year time loop under reverse-mode AD in particular
   (`push!`-based growth is a known Zygote pain point; less of an issue for
   `ForwardDiff`/`Enzyme` but still worth deciding up front).
5. A concrete choice between forward-mode (measured here: 2–5x memory cost
   per differentiated dimension, fine for a handful of scalar parameters,
   `1.46–3.64 GiB` per `ClimateFields` instance for 1–4 parameters) and
   reverse-mode (`Enzyme` — handles mutation well, cost independent of input
   dimension, but its own interaction with `@turbo` needs the same kind of
   direct probing done here before trusting it, and this session didn't have
   time to run that probe).
6. `Distributed.jl` as a documented, not-yet-implemented option for scaling
   the ensemble mechanism beyond one machine (per-worker process isolation
   trades away the low overhead measured above for true memory/fault
   isolation and multi-machine reach — worth revisiting if ensemble sizes
   grow past what a single machine's RAM can hold, given the ~745 MiB/member
   floor measured above).
7. A concurrency-safe replacement for the `redirect_stdout(devnull)`
   print-suppression idiom before any ensemble mechanism reuses it inside
   spawned tasks (per-task `IOBuffer` capture, or removing the `println`
   calls from `greb_model!`'s hot path in favor of `@debug`/`Logging`, which
   *is* task-safe unlike `redirect_stdout`).

## Small, low-risk follow-ups worth considering separately (not sized here beyond what's above)

- Stop `greb_model!` mutating the caller's `cfg.co2_scenario` in place
  (`src/model.jl:354,359`) — a real isolation hazard today, not just for a
  future ensemble mechanism (calling `greb_model!` twice with the same `cfg`
  for two different scenarios already silently misbehaves).
- Extract the sine-wave/linear CO2-forcing formulas already inside
  `forcing()` into small named scalar helper functions — zero behavior
  change, but turns today's one differentiable-by-accident code path into a
  deliberate, discoverable pattern for future examples.
