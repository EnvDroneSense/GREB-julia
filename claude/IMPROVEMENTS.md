# GREB.jl — Potential Improvements (Structure & Performance)

Observations and concrete suggestions for the `GREB` module. As of this pass
(branch `optimize/safe-mechanical-pass`) the model has been split from a
single 2245-line `src/GREB.jl` into topical files under `src/` — see §1.2.
Nothing here has been applied except where marked **✅ done** — this remains a
menu for future work, primarily the state-struct refactor in §1.1.

The single theme underneath most of what remains: **the model runs on ~40
mutable module-level globals** (grid arrays, climatology fields, work buffers,
flux corrections). That one design choice drives most of the remaining
structural *and* performance issues below — it's why the module eats ~750 MB
at load and why the code is hard to thread or test.

---

## 0. Known bugs — fixed this pass ✅

Three real, previously-shipped bugs were found and fixed while working through
this list (the first while implementing §0 from the prior version of this
doc; the other two while re-verifying §1.1's claims before splitting files —
finding which globals were live vs. dead required tracing every read/write,
which surfaced these). All three are fixed on `optimize/safe-mechanical-pass`
with regression tests added in `test/runtests.jl`.

### 0.1 `hydro!` referenced a field that doesn't exist on `CirculationWorkspace`
In `hydro!`'s `cfg.log_eva == 0` branch (now `src/physics/hydrology.jl`):
`ws.cE[i, j] = ...` — but `CirculationWorkspace` only declares `cE_buf`, not
`cE`. Threw `type CirculationWorkspace has no field cE` at runtime for any
`PhysicsConfig` with `log_eva = 0`. Invisible to the test suite because
`PhysicsConfig`'s default is `log_eva = -1`, which never reached that branch.
**Fixed**: renamed both occurrences to `ws.cE_buf`.

### 0.2 `HYDRO_PARAMS` was a `Tuple` of `Pair`s, not a `Dict`
`set_hydrology_parameters!` (now `src/config.jl`) did
`get(HYDRO_PARAMS, cfg.log_rain, default)` against a **tuple literal**
`(-1 => (...), 1 => (...), ...)`. `get` on a `Tuple` does *positional* integer
indexing, not key lookup, so:
- `log_rain = 0` (the default for every `PhysicsConfig`) → index 0 is out of
  tuple bounds → silently fell back to `(1.0, 0.0, 0.0, 0.0)` ("Original
  GREB") instead of the value tagged to key `0` ("Best GREB / ERA-Interim").
- `log_rain ∈ {1, 2, 3}` (valid tuple positions) → `get` returned the *`Pair`
  at that position* (wrong type, wrong association) → destructuring it into 4
  variables threw `BoundsError`. These modes have never worked.

**Fixed**: wrapped the literal in `Dict(...)` so `get` does real key lookup.

### 0.3 `set_hydrology_parameters!` never wrote back to `cfg`
Deeper than 0.2: the function computed the right values but assigned them to
**module-level globals** `c_q, c_rq, c_omega, c_omegastd` (via a function-local
`global` declaration) — never to `cfg.c_q` etc., the `PhysicsConfig` fields
that `hydro!` actually reads. The `@info` log line printed the *correct*,
freshly-computed values (misleadingly, since they were never delivered
anywhere else). **Net effect: `log_rain` had no effect on the model at all** —
every run used `hydro!` with `cfg.c_q=1.0, cfg.c_rq=0.0, cfg.c_omega=0.0,
cfg.c_omegastd=0.0` (the struct defaults, i.e. "Original GREB" parameters),
regardless of `log_rain`. Combined with 0.2, this is the one bug in this
batch that changes actual model output: after the fix, `log_rain=0` (default)
now genuinely uses "Best GREB (ERA-Interim)" parameters as documented, and
`log_rain ∈ {1,2,3}` finally do something.

**Fixed**: `set_hydrology_parameters!` now assigns `cfg.c_q = ...` etc.
directly; the disconnected module globals were removed.

**Audit performed**: grepped every other `global` reassignment inside a
function (all others — `TF_correct`, `dTrad`, `sw_solar_forcing_state` — are
genuine top-level module state, not this "computed but never delivered"
pattern), every `cfg.field = ...` assignment (all in `create_experiment_config`,
all correct), and every other `Dict`-vs-tuple usage passed to `get()` (the
other two `Dict(...)` literals in the file are real dicts). No sibling
instances of 0.2/0.3 found.

---

## 1. Structure

### 1.1 Replace module-global mutable state with state structs  ⭐ highest impact — deferred
Today the model's entire state lives as module globals: climatology fields
(`Tclim`, `uclim`, `qclim`, …), work buffers, masks (`co2_part`), and flux
corrections (`TF_correct`, …). Functions read/write them by name
(`load_greb_jld2!` fills them; `tendencies!`/`forcing` read them).

Problems this causes:
- **No reentrancy / concurrency** — you cannot run two experiments at once, or
  thread over ensembles; they'd clobber shared globals.
- **Hidden coupling** — the kind of thing that produced §0.1–0.3.
- **Type instability** (see §2.1) and **eager 750 MB allocation** (see §2.2).
- **Hard to test** — you can't construct a small, isolated state (this is also
  why §0.1's bug slipped through: there's no lightweight way to unit-test
  `hydro!` in isolation across all `log_eva` branches — though the new
  per-kernel `@benchmark` suite in `benchmark/run_benchmarks.jl` is a step in
  that direction for performance, if not yet for correctness).

Suggested shape:
```julia
struct ClimateFields          # the loaded climatology (immutable container, mutable arrays)
    Tclim::Array{Float64,3}; uclim::Array{Float64,3}; ...
end
mutable struct ModelState     # evolving state + reusable work buffers
    Ts::Matrix{Float64}; Ta::Matrix{Float64}; To::Matrix{Float64}; q::Matrix{Float64}
    ws::CirculationWorkspace; acc::MonthlyAccumulator; ...
end
greb_model(fields::ClimateFields, cfg::PhysicsConfig, run::RunSpec) -> Result
```
`CirculationWorkspace` and `MonthlyAccumulator` already do this well for their
slices — extend the same pattern to the rest (climatology fields, flux
corrections, annual/monthly accumulators, now consolidated in `src/state.jl`).
This is a large but mechanical refactor and unlocks almost everything else.
**Deferred to a follow-up pass**: it touches nearly every function signature,
and should be validated against real control-run output (now that
`greb_dataset_jld2/` is available locally — see the repo's `.gitignore`) rather
than only structurally.

### 1.2 Split the module into files  ✅ done
`src/GREB.jl` was a single 2245-line file; it's now a module shell with
`include()`s, in dependency order:
```
src/GREB.jl              # module, usings, includes, exports, PrecompileTools workload
src/constants.jl          # grid/time consts, physical constants, HYDRO_PARAMS, calendar lookup
src/config.jl             # PhysicsConfig, create_experiment_config, set_hydrology_parameters!
src/state.jl              # CirculationWorkspace, MonthlyAccumulator, TimeState, MonthlyRecord, mutable globals
src/io.jl                 # read_jld2, load_*_jld2!
src/physics/radiation.jl  # SWradiation!, LWradiation!
src/physics/hydrology.jl  # hydro!
src/physics/ocean.jl      # seaice!, deep_ocean!
src/circulation.jl        # diffusion!, advection!, convergence!, circulation!
src/tendencies.jl         # tendencies!, forcing
src/output.jl             # diagnostics!, output!, time_loop!
src/postprocess.jl        # build_monthly_climatology, apply_scenario_anomalies, compute_annual_ice_climatology
src/model.jl              # init_model!, qflux_correction!, greb_model!
```
Pure move — no function body, constant value, or struct definition changed.
Verified by `using GREB` loading cleanly and the full test suite passing
identically before/after.

### 1.3 Collapse giant positional argument lists into structs — deferred
Several kernels still thread ~10–20 positional arrays:
- `output!(it, irec, mon, Ts0, Ta0, To0, q0, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens, …)` (`src/output.jl`)
- `diagnostics!(it, year, CO2, Ts0, Ta0, To0, q0, albedo, sw, lw_surf, q_lat, q_sens, timestate)` (`src/output.jl`)
- `time_loop!(it, year, CO2, mon, irec, Ts, Ta, q, To, output_buf, ws, acc, timestate, cfg)` (`src/output.jl`)

These are error-prone (positional mixups — the exact failure mode that would
have hidden §0.1 even longer if it were a swapped argument instead of a
missing field) and hard to read. Passing a `ModelState`/`Diagnostics` struct
(from §1.1) removes most of them. Deferred alongside §1.1.

### 1.4 Turn run durations into a `RunSpec`, not positional ints — deferred
`greb_model!(time_flux, time_ctrl, time_scnr, cfg; jld2_dir="")` (`src/model.jl`)
takes three bare ints whose order is easy to swap. A `RunSpec(; flux=0, ctrl=1,
scnr=1)` keyword struct is self-documenting and pairs with §1.1.

### 1.5 ~~Reunite the two notebook-only helpers with the package~~ ✅ done
`examples/run_greb.jl` already provides a non-Pluto path — it builds a
`PhysicsConfig` via `create_experiment_config` and calls `greb_model!` directly.
No further action needed here.

### 1.6 Remove hidden/dead globals ✅ done
Removed, all confirmed dead by grep (no reads/writes anywhere else):
- `WZ_CACHE` — declared, never used.
- `circ_workspace` — a module-global `CirculationWorkspace()` instance that
  duplicated the one every `greb_model!` call builds for itself; never read.
- The entire **"Monthly-mean buffers"** block (`Tmm, Tamm, Tomm, qmm, apmm,
  icmm, prmm, evmm, qcrclmm, swmm, lwmm, qlatmm, qsensmm` — 13 arrays): every
  real use in the file goes through `acc.Tmm` etc. (the `MonthlyAccumulator`
  struct's identically-named fields) — these bare module globals were pure
  shadowed dead weight.
- The entire **"Control run monthly means"** block (`Tmn_ctrl, Tamn_ctrl,
  Tomn_ctrl, qmn_ctrl, icmn_ctrl, prmn_ctrl, evamn_ctrl, qcrclmn_ctrl` — 8
  arrays): `forcing`'s `icmn_ctrl` *parameter* shadowed the global name
  locally; the caller in `greb_model!` always passes a freshly-computed
  `ice_forcing` array, never this global.
- 4 dead entries (`icmn, prmn, evmn, qcrclmn`) inside the otherwise-live
  "annual-mean accumulators" block in `src/state.jl` (`Tsmn, Tamn, ...` are
  real and kept).

`sw_solar_forcing_state = Ref(1.0)` remains — it's a real, actively-mutated
global (set in `greb_model!`, read in `SWradiation!`). Belongs in run state
once §1.1 lands.

---

## 2. Performance

### 2.1 ~~`const`-ify the constants~~ ✅ done
Grid dimensions (`xdim`, `ydim`, `nstep_yr`, …) and the physical constants
(`σ`, `ρ_ocean`, `grav`, `cp_air`, …), now in `src/constants.jl`, are all
`const`. No further action needed.

### 2.2 Stop allocating ~750 MB at module load — deferred
There are still **28 top-level `zeros(Float64, xdim, ydim, nstep_yr)`
allocations** in `src/state.jl` — ~27 MB each ⇒ **~750 MB reserved even for
`using GREB` with no run**. Move these into the state structs from §1.1 so
they're allocated only when a model is actually set up.

### 2.3 Consider `Float32` for climatology fields — deferred
The JLD2 files store `Float32` data and are up-converted to `Float64` on load
(`src/io.jl`). Keeping the large 3D fields `Float32` (or making element type a
parameter) halves memory (~375 MB) and roughly doubles SIMD lane throughput
under `@turbo`. Now that `greb_dataset_jld2/` is available locally, this can
be benchmarked against real data (see `benchmark/run_benchmarks.jl`) rather
than estimated.

### 2.4 Thread the spatial operators — deferred
Time-stepping is inherently sequential, but per-step spatial kernels
(`diffusion!`, `advection!`, `SWradiation!`, `LWradiation!`, `hydro!`) are
independent across grid columns. Once state is passed explicitly (§1.1, no
shared globals) these become safe to `@threads` / `@batch` (Polyester) over
latitude — no threading is used anywhere in the module today.

### 2.5 ~~Cut per-call allocations in `forcing`~~ ✅ done
`forcing`'s `icmn_ctrl` is a required argument (`src/tendencies.jl`); the
caller in `greb_model!` passes a single `ice_forcing` array computed once via
`compute_annual_ice_climatology`, not fresh zeros every call.

**Newly measured** (via `benchmark/run_benchmarks.jl`, see `claude/BENCHMARKS.md`):
`diffusion!` and `advection!` each allocate several thousand times per call
(114–140 KiB) despite writing into pre-allocated workspace buffers, and
`circulation!`/`tendencies!` allocate several MiB per call (`circulation!`
calls `diffusion!`/`advection!` up to `ntime`=24 times per invocation, so
these compound). This wasn't previously measured or documented. Worth a
`@code_warntype`/allocation-source audit before the state-struct refactor —
the closures created by the `wz = if h_scl == z_air ... elseif ... end`
branch inside `diffusion!`/`advection!` are a reasonable first suspect.

### 2.6 Reduce time-to-first-run (TTFX) ✅ done
Added `PrecompileTools.@compile_workload` (`src/GREB.jl`, new `PrecompileTools`
dependency) running two tiny 1-year control-only runs (default config, and
`log_eva=0`) on unloaded fields, so the heavy kernels — `SWradiation!`,
`LWradiation!`, both `hydro!` branches, `circulation!`, `diffusion!`,
`advection!`, `convergence!`, `seaice!`, `deep_ocean!`, `output!`,
`diagnostics!`, `time_loop!`, `qflux_correction!` — compile at package build
time. Trade-off: precompile itself now takes ~2 minutes instead of seconds;
that cost is paid once (or in CI) rather than by every user's first call.

---

## 3. Testing, tooling & reproducibility

- **Reference/regression tests**: `test/runtests.jl` now includes
  branch-coverage loops over `log_eva ∈ {-1,0,1}` × `log_rain ∈ {-1,0,1,2,3}`
  (added specifically because §0.1 hid in an untested branch) and a direct
  unit test asserting `set_hydrology_parameters!` writes the right values
  into `cfg` for every `log_rain`/`log_clim` combination (added because of
  §0.2/§0.3). Still no golden/snapshot test against real numeric output —
  now that `greb_dataset_jld2/` is available locally (gitignored, not
  committed), a short real control run could be snapshotted and diffed with a
  tolerance. Not done this pass — deferred alongside §1.1, since that's when
  numeric-output validation becomes essential anyway.
- **Per-kernel unit tests**: not added as `Test.jl` unit tests, but
  `benchmark/run_benchmarks.jl` now exercises every hot kernel individually
  with synthetic inputs (see below) — a natural base to extend with
  correctness assertions (e.g. `diffusion!` conserves the field mean).
- **CI** ✅ done: `.github/workflows/ci.yml` runs `Pkg.test()` on a Julia
  `1.9`/`1` × `ubuntu-latest` matrix on push/PR to `main`.
- **`[compat]` bounds** ✅ done: `Project.toml` now bounds `julia`, `JLD2`,
  `LoopVectorization`, and `PrecompileTools` at their currently-resolved
  versions.
- **`NCDatasets` and `StaticArrays` dependencies** ✅ done: both were declared
  but unused in `src/GREB.jl` (`NCDatasets`: no `NCDataset(...)` call anywhere
  in the package; `StaticArrays`: no `SVector`/`SArray`/`SMatrix`/`MVector`
  anywhere — the "static longitude indices" comment was aspirational, they're
  plain `Vector{Int}` comprehensions). Both dropped from `Project.toml [deps]`
  and their `using` statements removed. (The notebook still uses `NCDatasets`
  in its own environment, untouched.)
- **Benchmark suite** ✅ done (new, beyond the original list): added
  `benchmark/run_benchmarks.jl` (own environment, mirroring `test/Project.toml`,
  so `BenchmarkTools` never becomes a package dependency) — `@benchmark`s
  every hot kernel against real data when `greb_dataset_jld2/` is present, and
  appends a dated, tagged, git-commit-stamped section to `claude/BENCHMARKS.md`
  every run, so before/after comparisons accumulate across optimization
  passes instead of being overwritten. A baseline entry (pre-this-pass) and a
  post-this-pass entry are both recorded there.
- **Data acquisition**: `scripts/convert_greb_to_jld2.jl` converts raw `.bin`
  files into the `.jld2` layout; the raw files themselves are still
  "available on request" per `DATA_README.md`. Unchanged this pass.
- **Docs**: a `Documenter.jl` site would suit the physics-heavy API; the
  notebook can become a rendered tutorial page. Unchanged this pass.

---

## 4. Suggested order (low-risk → high-value)

1. ~~Fix `§0` bugs~~ ✅ done — all three, with regression tests.
2. ~~Split into files (§1.2)~~ ✅ done — mechanical, improves everything after.
3. ~~Add `[compat]`, CI, dead-dependency cleanup, `log_eva`/`log_rain` branch
   coverage, benchmark suite (§3)~~ ✅ done.
4. **Introduce state structs** (§1.1, §1.3, §1.4) — the big one; enables
   §2.2, §2.4, and removes the global-coupling class of bugs for good. Next
   up — validate against real control-run output (`greb_dataset_jld2/` is now
   available locally) rather than structurally only.
5. Lazy allocation (§2.2), threading (§2.4), `Float32`/allocation audit
   (§2.3, §2.5's newly-measured `diffusion!`/`advection!`/`circulation!`
   allocations).

> ⚠️ Every performance change must be validated against a reference run —
> this is a numerical model, so "faster" only counts if the output is
> unchanged within tolerance. Use `benchmark/run_benchmarks.jl` +
> `claude/BENCHMARKS.md` to track this going forward.

---

## Changelog of this document

- **2026-08-06**: Safe/mechanical optimization pass (branch
  `optimize/safe-mechanical-pass`). Fixed and added regression tests for
  three real bugs (§0.1–0.3, the last two found while re-verifying this
  doc's own claims and newly discovered this pass). Split the module into
  files (§1.2). Removed all confirmed-dead globals (§1.6, expanded beyond
  what was previously documented). Dropped unused `NCDatasets` and
  `StaticArrays` dependencies, added `[compat]` bounds, added a
  `PrecompileTools` workload (§2.6), added GitHub Actions CI, and added a
  `BenchmarkTools`-based per-kernel benchmark suite with a running log
  (`claude/BENCHMARKS.md`) — which surfaced previously-unmeasured allocation
  hotspots in `diffusion!`/`advection!`/`circulation!` (§2.5). The
  state-struct refactor (§1.1/§1.3/§1.4) and the performance work that
  depends on it (§2.2–§2.4) are explicitly deferred to a follow-up pass, to
  be validated against real data now available at `greb_dataset_jld2/`
  (gitignored, not committed).
- **2026-08-05**: Refreshed against `src/GREB.jl` post-JLD2-conversion
  (commit `50a0bc2`). Marked §1.5, §2.1, §2.5 as done (already fixed in the
  codebase since this doc was first written). Confirmed §1.6/§2.2/§2.4 claims
  still hold via grep. Confirmed `NCDatasets` is unused. Added §0: a real
  `ws.cE` vs `ws.cE_buf` field-name bug in `hydro!`'s `log_eva == 0` branch,
  found while re-reading the file for this refresh.
