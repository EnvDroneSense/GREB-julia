# GREB.jl — Potential Improvements (Structure & Performance)

Observations and fixes for the `GREB` module, accumulated across seven audit
passes (2026-08-05 through 2026-08-12). The model was split from a single
2245-line `src/GREB.jl` into topical files under `src/` (§1.2) and its ~40
mutable module-level globals were replaced with explicit `ClimateFields`/
`ModelState` structs (§1.1) — that one refactor is why the module used to eat
~750 MB at load, why it was hard to thread or test, and why the same
"config flag silently disconnected from behavior" bug shape kept recurring.
Everything below is either **✅ done**, explicitly deferred with a reason, or
a documented finding that was investigated and intentionally not changed.
For the pass-by-pass discovery narrative and full benchmark methodology
behind any entry below, see [`CHANGELOG.md`](CHANGELOG.md).

---

## 0. Known bugs — fixed ✅

18 real, previously-shipped bugs found and fixed across five passes (a 19th
finding was investigated and reverted — see 0.14), each with a regression
test in `test/runtests.jl`. Passes 1–3 were found by grepping every
config-field write against its read sites; pass 4 was a direct line-by-line
comparison against the Fortran reference `greb.model.mscm.f90`, and every fix
there cites the exact Fortran line(s) confirmed by direct read (not recalled
from memory). Discovery narrative for each: `CHANGELOG.md`.

| # | Symptom | Fix | Fortran ref |
|---|---|---|---|
| 0.1 | `hydro!`'s `log_eva==0` branch referenced nonexistent `ws.cE`, threw at runtime (untested — default `log_eva=-1` never reached it) | Renamed to `ws.cE_buf` | — |
| 0.2 | `HYDRO_PARAMS` was a `Tuple` of `Pair`s; `get(...)` did positional indexing, silently wrong preset at default `log_rain=0` | Wrapped literal in `Dict(...)` | — |
| 0.3 | `set_hydrology_parameters!` wrote to disconnected globals, not `cfg` — `log_rain` had **zero effect** on any run | Assign `cfg.c_q` etc. directly | — |
| 0.4 | `circulation!`'s `do_conv` gate inverted vs. its four siblings — moisture-flux convergence off by default | All five flags now direct `Bool` truthiness | — |
| 0.5 | `co2_part` module-global leaked between runs (regional-CO₂ experiments) | `init_model!` resets `co2_part .= 1.0`; moot after §1.1 | — |
| 0.6 | Invalid `log_rain` silently fell back to a default instead of erroring | `error()` on unrecognized key | — |
| 0.7 | `hydro!`'s `log_eva` modes 1/2 didn't exist, silently duplicated mode `-1` | Implemented Fortran's real coefficients per mode + `error()` fallback | `greb.model.mscm.f90` |
| 0.8 | `init_model!`'s `_drsp` climatology overrides mismapped (`cldclim` missing, `qclim`/`mldclim` gated on wrong switch) | Added missing override, corrected switch names | `greb.model.mscm.f90` |
| 0.9 | `LWradiation!` zeroed both `LW_up`/`LW_down` on `log_atmos_dmc` off; Fortran only zeros `LW_down` | Removed extra `LW_up .= 0.0` | `greb.model.mscm.f90` |
| 0.10 | Paleo/orbital solar-forcing table swap loaded but never wired up | `greb_model!` swaps table at scenario start for the 4 affected experiments | `greb.model.mscm.f90` |
| 0.11 | `circulation!` never zeroed `ws.dX_conv` — every `Ta` sub-step in every default run picked up a stale humidity convergence term | Zero all 3 buffers once at top of `circulation!` | `greb.model.mscm.f90:856-858` |
| 0.12 | `deep_ocean!` disabled turbulent mixing under sea ice (wrong shared mask) | Split mask: turbulent mixing uses `z_topo<0` alone | `greb.model.mscm.f90:818-830` |
| 0.13 | `hydro!`'s rain limit used a single grid point (`wz_vapor[1,1]`) globally | Per-point `wz_vapor[i,j]` | `greb.model.mscm.f90:755` |
| 0.14 | `min_T_K=233.15` K vs. Fortran's raw `40` K floor — investigated, kept as-is | Tried matching Fortran, reverted: 40 K never physically binds, not worth copying literally | `greb.model.mscm.f90:470-477` |
| 0.15 | `seaice!` early `return` skipped the glacier→`cap_land` override when `log_ice` off | Removed early return | `greb.model.mscm.f90:786-792` |
| 0.16 | `hydro!` computed `rq` after `log_eva` dispatch, using an already-overwritten `ws.qs` | Moved `rq` computation before dispatch | `greb.model.mscm.f90:715-719` |
| 0.17 | `diffusion!`/`advection!` polar sub-stepping was Gauss-Seidel under `@turbo` (unsound), not Jacobi | Two-pass Jacobi via `ws.dTxh` scratch buffer | `greb.model.mscm.f90:983-1036,1227-1228` |
| 0.18 | Humidity update used `clamp` instead of Fortran's threshold-replace; `log_hydro_dmc==0` only zeroed eva/rain, not all of `dq` | Fixed both in `time_loop!`'s humidity update | `greb.model.mscm.f90:481-483,486` |
| 0.19 | 19 of 36 exported functions' docstrings silently detached from their bindings (notebook-cell comment between docstring and def) | Moved comment before docstring across 9 files; verified with `Docs.hasdoc` | — |

**Investigated, not changed**: `log_eva==1`'s land gust literal (`144.0`,
Fortran has `144.**2`) — kept unsquared, Fortran's own exponent looks like
the anomaly given every sibling mode's gust term is 2–12 m/s; confirmed
identical in a second Fortran copy, so not a transcription error.
`:earth_sun_distance`'s `rS0` rescale confirmed correctly scenario-only in
both Julia and Fortran. Full reasoning: `CHANGELOG.md`.

---

## 1. Structure

### 1.1 Module-global mutable state → state structs ✅ done
Replaced by `ClimateFields`/`ModelState` (`src/state.jl`), threaded
explicitly through every physics/circulation/tendencies/output/model
function. `greb_model!` defaults to a fresh `ClimateFields()` per call,
making cross-run leaks (§0.5) structurally impossible for that path.
Validated to ~7e-12 absolute agreement against pre-refactor output; the same
1-year run went from 49s/17.1 GiB to 7.0s/100 MiB (unplanned bonus that also
fixed §2.2 and unblocked §2.4).

### 1.2 Module split into files ✅ done
`src/GREB.jl` is now a shell with `include()`s in dependency order:
`constants.jl` → `config.jl` → `state.jl` → `io.jl` →
`physics/{radiation,hydrology,ocean}.jl` → `circulation.jl` →
`tendencies.jl` → `output.jl` → `postprocess.jl` → `model.jl`.

### 1.3 Positional argument lists → structs ✅ done
`SurfaceState` (`src/state.jl`) wraps `Ts0`/`Ta0`/`To0`/`q0` by reference;
the rest of `output!`/`diagnostics!`'s former 9–13 positional args already
lived on the existing `tend` `NamedTuple` or `ws::CirculationWorkspace`, so
those are passed through instead of a second struct. Zero new allocations.
Bonus: fixed a real ordering footgun between `time_loop!`'s and
`diagnostics!`/`output!`'s old positional argument order (`q`/`To` swapped).

### 1.4 Run durations → `RunSpec` ✅ done
`greb_model!(run::RunSpec, cfg; ...)`, `RunSpec(; flux=0, ctrl=1, scnr=1)`
(`src/config.jl`), replacing three easily-swapped positional ints. All call
sites updated.

### 1.5 Reunite notebook-only helpers with the package ✅ done
`examples/run_greb.jl` provides a non-Pluto path.

### 1.6 Remove hidden/dead globals/fields ✅ done
Removed: `WZ_CACHE`, `circ_workspace`, two shadowed monthly-mean buffer
blocks (21 arrays), 4 dead accumulator entries, `co2_part_scn`,
`sw_solar_forcing_data`, `CirculationWorkspace`'s `dX_crcl`/`ws.eva`/
`ws.rain`. All confirmed dead by grep before removal.

---

## 2. Performance

### 2.1 Constants ✅ done
Grid dimensions and physical constants in `src/constants.jl` are `const`.

### 2.2 Stop allocating ~750 MB at module load ✅ done
Direct consequence of §1.1 — climatology arrays are `ClimateFields` fields,
allocated only when constructed, not at `using GREB` time.

### 2.3 `Float32` for climatology/working fields — ✅ done (2026-08-13)
Retyping just the static (JLD2-native-`Float32`) climatology is safe, but
`init_model!`'s `copy(Tclim[:,:,end])`-style inits would silently downgrade
the *working* simulation state (`Ts`/`Ta`/`To`/`q`) too without an explicit
`Float64.(...)` guardrail. Measured the actual arithmetic payoff of a full
switch (standalone clone of the circulation kernels, both eltypes):

| | `circulation!(Ta)` | `circulation!(q)` | total |
|---|---|---|---|
| `Float64` | 1301–1385 µs | 1360–1455 µs | 2661–2840 µs |
| `Float32` | 1210–1218 µs | 1337–1370 µs | 2554–2580 µs |

Only **~5–10%** faster, not 2× — the zonal terms gather through permuted
neighbor-index arrays rather than reading contiguous memory, and the 96×48
grid already fits in cache at either precision. `CirculationWorkspace` is
also a concrete non-parametric `Float64` struct, so a real switch needs a
`CirculationWorkspace{T}` refactor first. **Not worth it for this kernel
shape** given §2.11 already gets 2.0× more cheaply; canceled by user
request. Could matter more for a GPU/ensemble-batched use case (§2.11),
untested here. Full methodology: `CHANGELOG.md`.

**Follow-up, benchmarked (`CircularArrays.jl` for the gather pattern itself) — tried, measured, rejected.**
An outside suggestion proposed `CircularArrays.jl` as a direct fix for the
gather pattern identified above: wrap `T1`/`wz` so periodic-longitude
neighbor reads (`T1[j-1,k]`/`T1[j+1,k]`) become ordinary offset accesses
instead of indexing through `lon_jm1..jp3`. Built a standalone clone of
`diffusion!`'s zonal stencil (`src/circulation.jl:73-94`) two ways —
today's `Matrix{Float64}` + index-array approach vs. a `CircularVector`
wrapping each row (`CircularArrays.jl` circularizes *every* dimension of a
`CircularArray`, so a whole-2D-array wrap would incorrectly wrap the
non-periodic latitude axis too — a per-row `CircularVector` is the only
structurally correct option) — and measured both, `@turbo`-annotated
identically:

| | time/call | vs. baseline |
|---|---|---|
| Baseline (`Matrix` + index arrays, `@turbo`) | 16–30 µs | 1× |
| `CircularVector`-per-row (`@turbo`) | 340–426 µs | **0.03–0.09× (10–20× slower)** |

Output was correct (diffs `1.6e-16`, float-reassociation only, 0 bytes/call
either way), but `@turbo` **silently fails to compile** over the
`CircularVector`-backed loop (`LoopVectorization.check_args` rejects it,
falling back to a plain `@inbounds @fastmath` loop with a runtime warning)
— confirmed directly in `code_native`: the baseline's inner loop is real
packed-SIMD code (`vgatherqpd`×12, `vmovupd`/`vmovapd`×26, `...pd` = packed
double), while the circular version's inner loop is **100% scalar**
(`vmovsd`/`vmulsd`/`vsubsd`/`vaddsd`, `...sd` = scalar double, zero packed
instructions of any kind). The fix removes the gather instructions
entirely, exactly as claimed — but takes the SIMD vectorization down with
it, which costs far more than the gather ever did. **Rejected**: this is
worse than doing nothing on every axis (correctness aside), confirming the
`@turbo`-compatibility risk flagged before running the experiment. Script:
`benchmark/circular_arrays_experiment.jl`.

**Re-review (2026-08-13) — different verdict for non-gather kernels, still
not implemented.** The original analysis above only ever tested
`circulation!` (gather-bound via the periodic-longitude index arrays). Redid
the Float32-vs-Float64 comparison fresh on this machine/Julia version,
*and* extended it to `hydro!`/`LWradiation!`-style kernels, which are pure
elementwise (no neighbor gather at all — every read is at `[i,j]`, nothing
offset):

| Kernel shape | `Float64` | `Float32` | Speedup |
|---|---|---|---|
| Gather-bound (`circulation!`-style zonal diffusion) | 210.5 ns | 220.9 ns | **0.95× (no gain, slightly slower)** |
| Pure elementwise, no gather (`hydro!`-style) | 8733 ns | 2889 ns | **3.02×** |
| Transcendental-heavy elementwise (`LWradiation!`-style, `log()` calls) | 29800 ns | 11300 ns | **2.64×** |

So the original "not worth it" conclusion still holds for `circulation!`
itself — on this machine Float32 doesn't even keep its earlier claimed
5-10%, confirming it's not the lever for the dominant 65-75%-of-runtime
kernel. But `hydro!`/`LWradiation!`/`ocean!`'s elementwise kernels (no
gather to fight) show a genuine **2.6-3× speedup**, never tested before —
plausibly double SIMD lane width plus faster `exp`/`log`/`sqrt` at single
precision. These three together are a meaningful chunk of the ~25-35% of
per-timestep cost that isn't circulation (§2.10), so a 2.6-3× speedup there
could plausibly be a double-digit-percent full-model win — bigger than
anything in §2.15 and in the same league as §2.11's threading win.
Precision check: accumulating 730 `Float32` additions of O(250) values (a
proxy for `diagnostics!`/`output!`'s annual-mean accumulators) against the
same in `Float64` gives a relative error of ~2e-7 — negligible if the
*accumulator* arrays stay `Float64` and only per-timestep working buffers
go to `Float32`.
**Not implemented, and not a quick follow-up**: this would need a real
design decision, not a drop-in change — `ClimateFields`/`CirculationWorkspace`
are concrete non-parametric `Float64` structs (same blocker as before), a
mixed-precision split (elementwise physics in `Float32`, circulation and
accumulators in `Float64`) means crossing precision boundaries at several
call sites, and any change here needs validation against the Fortran
reference over full multi-decade runs, not just a standalone kernel diff.
Worth a dedicated pass if the user wants to pursue it — flagging it here as
reopened, not reclosed.

**Implemented (2026-08-13) — full conversion, plus the fix that made the
gather kernel benefit too.** Investigating the "not worth it for
`circulation!`" result further (per user question) turned up the actual
cause: the gather-kernel comparison above changed the data to `Float32`
but left the periodic-longitude index arrays (`lon_jm1`/`lon_jp1`/`lon_jm2`/
`lon_jp2`/`lon_jm3`/`lon_jp3`) as `Int64` — a width mismatch.
`LoopVectorization` needs matching-width gather addresses; `Float32` data
gathered through `Int64` indices nearly *doubled* the gather-instruction
count (15→28) instead of halving it. Narrowing the index arrays to `Int32`
alongside the data dropped the count back to 15 and made the kernel
genuinely faster than the `Float64`/`Int64` baseline (~245ns → ~175ns,
confirmed reproducible across 3 runs) — reversing the "not worth it"
verdict for `circulation!` itself, not just the elementwise kernels.

Given that, converted the whole model to `Float32` end-to-end, per explicit
user decision on scope (not a mixed-precision split):
- `ClimateFields`/`CirculationWorkspace`/`SurfaceState`/`MonthlyAccumulator`/
  `ModelState`/`MonthlyRecord` (`state.jl`): every field retyped, including
  the public output (`MonthlyRecord` is now `Matrix{Float32}`, not promoted
  back to `Float64`).
- `PhysicsConfig`'s numeric fields (`co2_concentration`, `c_q`/`c_rq`/
  `c_omega`/`c_omegastd`, `earth_sun_distance_pct`, `co2_scenario::Dict{Int,
  Float32}`) — also converted, not just cast locally, so a stray `Float64`
  scalar can't re-promote an elementwise kernel back to double precision.
- `lon_jm1`/`lon_jp1`/`lon_jm2`/`lon_jp2`/`lon_jm3`/`lon_jp3` → `Int32` (the
  fix above), plus every other `constants.jl` physical/derived constant.
- Every bare `Float64` literal inside a `@turbo` loop across
  `circulation.jl`, `physics/hydrology.jl`, `physics/radiation.jl`,
  `physics/ocean.jl`, `output.jl`, `model.jl`, `tendencies.jl` got an `f0`
  suffix — the highest-risk, most mechanical part of the change, since one
  missed literal silently re-promotes that loop to `Float64` and reproduces
  the original false-negative. One simplification found along the way:
  plain integer literals (`3`, `4`, `20`) do **not** need suffixing —
  `Float32 op Int` promotes to `Float32` in Julia, not `Float64`; only
  literals with a decimal point/exponent do.
- `load_greb_jld2!` no longer upconverts at all — JLD2 data is already
  `Float32`, so loading is now a same-type copy instead of a promotion (a
  small free win on top).

**Verification**: full 438-test suite passes (6 tests needed tolerance
loosening — `atol=1e-9`/`rtol=1e-10` style assertions to `~1e-5`–`1e-4`,
an expected consequence of the precision switch, not a workaround; two
`==` comparisons against `Float64` literal tuples/dicts needed `isapprox`
since a round-tripped `Float32` literal doesn't compare exactly equal to
the original `Float64` literal). Numeric sanity check (1-year control run,
`Float32` vs. a saved `Float64` reference from before the change): max
absolute difference 0.0036 K on `Ts`, 0.00086 K on `Ta`, 0.003 K on `To`,
2e-7 kg/kg on `q` — physically negligible, and the printed annual
diagnostic line matched the `Float64` baseline to 2 decimal places.
`@code_native` on the real (not standalone) `diffusion!` confirmed no
`cvtss2sd`/`cvtsd2ss`-style conversion instructions anywhere in the
compiled loop — no literal was missed.

**Real full-model result**: `benchmark/run_benchmarks.jl` (real NCEP data,
`-t 3`) — **0.68-0.72s/year, mean 0.70s**, vs. the ~1.1-1.2s/year `Float64`
baseline measured earlier this session. **~1.6× faster**, matching the
weighted estimate from the kernel-level numbers above and landing in the
same league as §2.11's threading win — the single biggest lever found this
session.

### 2.4 Thread the spatial operators — tried, measured, reverted
`Threads.@threads :static` on `diffusion!`/`advection!`'s outer row loop,
validated bit-identical output, but **3–6× slower** end-to-end at `-t auto`/14
threads (`tendencies!`: 4.2ms→23.7ms) — the grid (96×48, ~1µs/row) is too
small for the scheduling overhead to pay off at any thread count tried.
Reverted to plain `for k in 1:ydim` loops. **Takeaway**: don't retry plain
`Threads.@threads` at this grid/kernel shape without a reason to expect a
different result; `Polyester.@batch` or a much larger grid might fare
better.

**Follow-up, benchmarked (`Polyester.jl` retry) — real but narrow, fragile win; not implemented.**
Retried the reverted idea with `Polyester.@batch` instead of
`Threads.@threads`, on a standalone clone of `diffusion!`'s *full* outer-`k`
loop (mid-latitude **and** polar branches — the polar branches reuse a
single shared `ws.T1h`/`ws.dTxh` scratch buffer across an inner sub-step
loop in the real code, `src/state.jl:10-11`, which would race under
concurrent `k`; the clone gives each lane a private `(xdim, nlanes)`
scratch matrix indexed by `Threads.threadid()` instead). Verified
bit-identical output vs. sequential across every thread count tested,
including 5 repeated trials per count to rule out an intermittent race —
none found. Timed at every thread count from 1 to the machine's full 14
logical cores:

| `-t` | time/call | speedup |
|---|---|---|
| 1 | 61.6 µs | 1.00× |
| 2 | 65.1 µs | 0.92× |
| 3 | 54.3 µs | 1.17× |
| **4** | **35.4 µs** | **1.79×** |
| 5 | 110.5 µs | 0.55× |
| 6 | 158.7 µs | 0.37× |
| 7 | 166.6 µs | 0.35× |
| 8 | 137.6 µs | 0.43× |
| 14 (all logical cores) | 16,900–18,000 µs | **0.003–0.004×** (270–340× *slower*) |

Real result, not a wash: **`-t 4` gives a genuine 1.79× speedup** with
correct output — a materially better outcome than `Threads.@threads`'s flat
3–6× slowdown at every thread count tried in the original attempt. But the
curve is sharply non-monotonic and **not safe to run at `-t auto`** — every
count above 4 degrades, and running at the machine's full logical-core
count (14) doesn't just lose the speedup, it collapses catastrophically
(reproduced twice, both ~300× slower than sequential), consistent with
Polyester's static per-lane scheduling overhead dominating once lane count
exceeds what a 48-row loop can actually use, compounded by likely
oversubscription past physical-core count on a hyperthreaded machine. **Not
implemented**: a real win exists, but only if the thread count is
explicitly pinned (e.g. `-t 4`, not `-t auto`) rather than left to the
caller's environment — that's a behavioral/deployment decision, not a pure
code change, and is left for a follow-up once a specific thread-count
policy is decided. Script: `benchmark/polyester_experiment.jl`.

### 2.5 Cut per-call allocations in `forcing` ✅ done
`icmn_ctrl` is a required argument computed once, not fresh zeros every call.

### 2.6 Reduce time-to-first-run (TTFX) ✅ done
`PrecompileTools.@compile_workload` compiles the heavy kernels at package
build time instead of on every user's first call.

### 2.7 Cut `SWradiation!`'s per-call allocation ✅ done
`@views`-wrapped the loop body — was materializing a fresh vector per
iteration (40704 bytes/48 allocs), the *entire* allocation footprint of
`tendencies!`. Every kernel now benchmarks at 0 bytes.

### 2.8 Hot-loop field/config localization + type stability — audited, clean
`fields.xyz`/`cfg.xyz` are already hoisted to locals before every hot loop;
no struct field is untyped/abstract. Nothing to change.

### 2.9 Precomputation/caching — audited, clean
Every candidate for hoisting is already hoisted (`regional_co2_*` masks
where possible — the other 2 genuinely can't be, they depend on
`icmn_ctrl`) or already a precomputed lookup table/`const` array. Nothing to
do.

**Flagged, not fixed, while auditing**: `output.jl`'s `global_mean` was an
unweighted grid mean — resolved, see §4.11.

### 2.10 Where per-timestep time goes vs. Fortran — profiled, documented
A 1-simulated-year `:full_model` control run takes ~2.7s wall-clock
(post-precompile). Profiling (`Profile.@profile` + manual `@elapsed`) over
5000 warm `time_loop!` calls attributes it as:

| Stage | Time | Share |
|---|---|---|
| `circulation!` (`Ta`, `z_air`, 24 substeps) | 1370 µs | ~34% |
| `circulation!` (`q`, `z_vapor`, 24 substeps) | 1498 µs | ~37% |
| everything else (radiation, hydrology, ocean, output) | ~1100–1800 µs | ~26–29% |

**~65–75% of every timestep is `circulation!`**, run twice per main step.
Ruled out as causes: allocations (96 bytes/step total) and the polar
sub-stepping (only the 2 outermost rows need >1 sub-step). The remaining gap
to native Fortran is likely structural — 24×2 separate `@turbo`-swept
function calls per step rarely fuse as tightly as one hand-written Fortran
subroutine.

**The Fortran-comparison assumptions were checked directly, not assumed:**
the reference build (`gfortran -Ofast -ffast-math -funroll-loops -fopenmp`,
found in `run.greb.*.csh` next to the Fortran source) isn't apples-to-apples
with a single-threaded Julia run. `greb.model.mscm.f90`'s `tendencies`
subroutine (lines 512–535) wraps `SW`/`LW`+sensible/`hydro`/`circulation`(Ta)/
`circulation`(q)/`deep_ocean` each in an `!$omp section` — the **same six
stages** §2.11 found and benchmarked independently — so the folklore "<1s/yr"
figure already assumes multi-threaded execution; the same folder's own
`readme.md` states unoptimized/unparallelized Fortran takes **~30s/simulated
year**, slower than this Julia port's current ~2.7s/year single-threaded.
Checked the other flags too: `julia -O3` vs. default `-O2` showed no
measurable difference (2.13–2.22s/year, within noise — `@turbo`'s codegen is
largely independent of the global `-O` level); `-march=native` is already
Julia's default (`Base.JLOptions().cpu_target == "native"`, confirmed on
this machine). `-ffast-math` has no Julia equivalent applied globally, but
see §2.12 for its per-kernel impact where it'd actually matter.

### 2.12 Missing `@turbo` spots — ✅ done (4 of 4 candidates)
Every hot kernel in `circulation.jl`/`hydrology.jl`/`ocean.jl` uses `@turbo`.
The 4 candidates below (all rewritten as explicit `@turbo` loops) were found
by checking whether `LWradiation!`'s fix generalizes, prioritizing by call
frequency (per-substep code compounds: `circulation!` runs 24 substeps × 2
fields = 48×/main-timestep). Verified against the full test suite (incl.
golden regression); real 1-year `greb_model!` run: 2.7s → 2.212s.

| Spot | Frequency | Baseline | `@turbo` | Speedup | ~Saved/year |
|---|---|---|---|---|---|
| `diffusion!` mid-lat branch (`circulation.jl:63-66`) | ~46 rows × 48/step | 9.2 µs/call | 2.3 µs/call | 75% | 0.22s |
| `circulation!` substep accumulate (`circulation.jl:326`) | 48/step | 208.8 µs/step | 45.6 µs/step | 78% | 0.12s |
| `LWradiation!` (`radiation.jl:82`, 3 `log()`/cell) | 1/step | 145.7 µs | 36.3 µs | 75% | 0.08s |
| `time_loop!` state update (`output.jl:137-169`) | 1/step | 63.3 µs | 12.3 µs | 80% | 0.04s |

All four verified numerically identical to baseline (diffs `1e-14`–`1e-18`,
pure float reassociation) — the first correctness pass hit a spurious `NaN`
traced to the *test harness* forgetting to call `init_model!`/`seaice!` first
(`cap_surf` still all-zero → division by zero), not a real bug. **Total
~0.46s/simulated year** across these four — composes with (doesn't
double-count) §2.11's ~1.0s/year, since §2.11 parallelizes *across* the two
`circulation!` calls while these speed up code *inside* each call.

**One negative result, useful as a boundary case**: `circulation!`'s final
subtraction (`circulation.jl:330`, a single elementwise op, only 2 calls/step
— not per-substep) was **17.6% slower** under `@turbo` — the macro's own
overhead isn't worth paying below some size/frequency threshold. Confirms
this isn't "add `@turbo` everywhere," it's "add it where the call is large
or frequent enough" — `LWradiation!`/mid-lat-diffusion/substep-accumulate all
clear that bar; this one line doesn't. `@fastmath` (tested only on
`LWradiation!`, 45% faster there) is the softer middle option when a spot is
too awkward to rewrite as an explicit `@turbo` loop.

**Follow-up sweep (2026-08-12), 4 more confirmed wins + 1 more negative
result** — a broader pass (2 parallel Explore agents across all of `src/`)
resolved the two candidates above and found 3 more of the same shape:

| Spot | Baseline | `@turbo` | Speedup |
|---|---|---|---|
| `diagnostics!` accumulate (`output.jl:11-21`, 11 ops, 1×/step) | 22.65 µs | 9.75 µs | 57% |
| `state.jl accumulate!` (13 ops, 1×/step) | 24.2 µs | 11.85 µs | 51% |
| `SWradiation!`'s `a_atmos`/albedo-combine/`sw`-flux-loop (1×/step) | 5.5 µs | 3.4 µs | 38% |
| `qflux_correction!`'s per-timestep update (12 ops, spin-up only) | 57.0 µs | 18.0 µs | 68% |

All four **implemented** (`output.jl`, `state.jl`, `radiation.jl`,
`model.jl`), verified numerically identical (diffs `1e-16`ish) and against
the full test suite incl. golden regression. Combined contribution to a
real `:full_model` `flux=0` run is small — the first 3 apply there (~20ms/
simulated year estimated from the per-call numbers × 730 steps); the 4th
only runs during the flux-correction spin-up (`run.flux > 0`, off by
default) so it contributes nothing to that particular run. This ~20ms/year
is genuine but small enough to be swamped by this machine's own run-to-run
wall-clock noise (~300ms swings observed across repeated 1-year timings) —
unlike §2.11/§2.12's fixes, which were large enough to show up cleanly
end-to-end, this batch's win is only visible in the isolated per-call
microbenchmarks above, not in an aggregate before/after timing. Documented
honestly rather than reporting a noisy or cherry-picked aggregate number.

**One more confirmed negative result**: `seaice!`'s trailing glacier-override
line (`ocean.jl:46`, single op) was flat-to-slower under `@turbo` (−2.5% to
−47% across repeated runs, i.e. noise-dominated but never a real win) —
**not implemented**, a second instance of the same "too small/infrequent"
pattern as `circulation.jl:330`. `tendencies.jl:39`'s `Q_sens = ct_sens *
(Ta - Ts)` (single op, 1×/step) was also consistently slower under `@turbo`
(−13% to −16%) — **not implemented**, third instance of the pattern. Both
below the profitability threshold: a single elementwise op over 96×48
doesn't carry enough work to amortize `@turbo`'s own dispatch/setup cost.

### 2.13 `forcing()`'s regional-CO₂ mask recomputed every scenario timestep — ✅ done
`:regional_co2_ocean`/`:regional_co2_land_ice` (`tendencies.jl`) recomputed
their `co2_part` mask from scratch on every scenario timestep, even though
it depends only on `z_topo` and `icmn_ctrl` — both fixed for the whole
scenario run. Not a `@turbo` candidate (branchy, not a clean elementwise
op) — fixed algorithmically by guarding the mask computation with `it ==
1`; `fields.co2_part` correctly holds the static mask for every later
timestep since nothing else touches it for these two experiments. Measured
in isolation (`forcing()` called directly, `:regional_co2_ocean`): 7.65 µs
baseline → 8.35 µs at `it==1` (pays the real cost once) → **0.1 µs at
`it>1`** (98.7% faster) — ~6ms/simulated year saved for these two
experiments (small in absolute terms, since the mask itself is cheap; the
point is skipping it 729/730 times per year is free). Verified via the full
test suite's `"greb_model! reaches every :experiment symbol..."` test,
which exercises both experiments directly.

### 2.15 Four small-scale fixes (2026-08-13) — ✅ done
All verified against the full 438-test suite and a full-model benchmark run
(`benchmark/run_benchmarks.jl`, real NCEP data, `-t 3`) before/after —
output unchanged, no regression. Each fix was also benchmarked **in
isolation** against its pre-fix version (old code pulled from git `HEAD`,
run side-by-side with the fixed version in the same process, real NCEP
data where applicable):

| Fix | Old | New | Speedup | Est. saving/year |
|---|---|---|---|---|
| `output!` copy elimination (per field) | 23.6 µs | 3.2 µs | **7.4×** | ~3.2 ms (×13 fields, ×12 months) |
| `LWradiation!` fused `LW_up` write | 36.6 µs | 34.3 µs | 1.07× | ~1.7 ms (×730 calls) |
| `hydro!` loop fusion | 16.8 µs | 11.5 µs | **1.46×** | ~3.9 ms (×730 calls) |
| polar sub-step constants (`diffusion!`) | 36.9 ns | 2.2 ns | **16.8×** | ~22 ms (×700,800 calls) |
| polar sub-step constants (`advection!`) | 33.3 ns | 2.4 ns | **17.7×** | ~22 ms (×700,800 calls) |

Combined estimated saving: **~53 ms/simulated year**, out of a ~1.1-1.2
s/year baseline (~4-5%) — real, but still inside the full-model benchmark's
own run-to-run noise band (0.94-1.22s/year measured on this machine), so it
doesn't show up as a clean before/after delta at the whole-run level; the
per-fix isolation above is what actually demonstrates each one's effect.
Detail per fix:
- **`output!`'s redundant double allocation** (`output.jl`): each of the 13
  monthly-record fields was built as `copy(acc.Tmm ./ ndm)` — `./` already
  returns a fresh, non-aliased array, so `copy(...)` allocated a second
  array and immediately discarded the first. Removed the redundant `copy`.
  Biggest *relative* win of the four (7.4×) but only runs at month
  boundaries (12×/year), so modest in absolute terms.
- **Extra full-grid pass in `LWradiation!`** (`radiation.jl`): `LW_up .=
  LW_down` copied the whole grid in a second pass right after the `@turbo`
  loop that computed `LW_down`. Merged the write into the same loop
  (`LW_up[i,j] = LW_down_val` alongside `LW_down[i,j] = LW_down_val`). Only
  a 1.07× win in isolation — copying a cache-resident 96×48 array is
  already fast (basically a 36 KB memcpy), so removing the second pass
  saves little; the value here is code cleanliness, not raw speed.
- **`hydro!`'s 6 separate `@turbo` loops over the same grid** (default
  config runs 5 of them; `log_rain==1` adds the 6th): saturation humidity,
  relative humidity, evaporation (one of 4 `log_eva` branches), precipitation,
  the optional rain-limit clamp, and water-vapor tendencies were each a
  separate full-grid pass, despite being purely elementwise (no
  cross-grid-point coupling) — nothing here needs neighbor values the way
  `diffusion!`/`advection!` do. Fused each `log_eva` branch into a single
  `@turbo` loop that walks the whole per-point pipeline (including the
  rain-limit clamp, now an inline `ifelse` gated on `cfg.log_rain==1` rather
  than a separate conditional pass) in one grid traversal. Output matches
  the pre-fusion version to floating-point roundoff (~1e-13, from
  reassociated `exp`/`sqrt` calls). `hydro!` runs once per main timestep
  (not per circulation sub-step), so this is ~730×/year.
- **Polar sub-step constants recomputed every call** (`circulation.jl`'s
  `diffusion!`/`advection!`): `dd`/`dtdff2`/`time2`/`ccx2` in each polar-row
  branch depend only on `k` (via `dxlat_grid[k]`) and fixed module
  constants — never on the per-call `T1` data — yet were recomputed from
  scratch on every call. There are **20 of 48 rows** flagged `IS_POLAR`
  (not a handful, as the raw formula count might suggest) — 20 rows × 2
  fields (Ta, q) × 24 sub-steps × 730 steps/year = **700,800 calls/year**
  for diffusion alone, and the same count again for advection (different
  formula, same waste). Precomputed as `POLAR_DIFF_TIME2`/`POLAR_DIFF_CCX2`/
  `POLAR_ADV_TIME2`/`POLAR_ADV_CCX2` (`constants.jl`, same pattern as the
  existing `ccx_diff`/`ccx_adv`/`dxlat_grid` lookup tables), indexed by `k`.
  Individually this was only a handful of scalar ops (~35 ns), but at
  700,800× the call count this turned out to be the **largest absolute
  saving of the four** (~44 ms/year combined) — the "small-scale fix" label
  undersold this one; the call-count multiplier is what makes it matter.

### 2.14 Other sweep findings — investigated, not implemented
A wider pass (2026-08-12) also looked beyond `@turbo` candidates, at
architecture/I-O patterns:
- **CMIP5/ENSO anomaly reload** (`model.jl`'s `greb_model!`,
  `load_cc_anomaly_jld2!`/`load_enso_anomaly_jld2!`): always re-reads from
  disk even when the same `fields` is reused across multiple `greb_model!`
  calls with the same experiment — real for a repeated-run/ensemble use
  case, but a safe fix needs a cache-key/invalidation design (what if
  `jld2_dir`/`cfg` differs between calls?) that's a genuine design decision,
  not a quick win — left for whoever actually hits this in an ensemble
  workflow.
- **`postprocess.jl`'s runtime-`Symbol` `getfield` loop**
  (`build_monthly_climatology`/`apply_scenario_anomalies`): a real dynamic
  dispatch per field per record, but dwarfed by the ~4608-element array op
  it wraps at current `MonthlyRecord` counts (12×`time_ctrl`/`time_scnr`) —
  not worth an `@generated`/`Val`-based rewrite at this scale.
- **`CirculationWorkspace()`'s ~40 separate allocations**: only happens 3×
  per `greb_model!` call (`ws`/`ws_a`/`ws_q`), not per-timestep — negligible
  (~1.44 MiB, low-µs) next to a multi-second run; consolidating into fewer
  larger arrays would save allocation count but not meaningfully change
  wall-clock time for the current usage pattern, and would hurt readability
  for no real gain absent a hypothetical high-frequency-construction use
  case that doesn't currently exist in this codebase.
- **Rust FFI for hot kernels** (2026-08-13): benchmarked a Rust `cdylib`
  (`opt-level=3`, `lto`, both default and `target-cpu=native`) against the
  existing `@turbo` mid-latitude zonal-diffusion kernel
  (`circulation.jl:72-94`), calling it via `ccall` with real NCEP climatology
  data and byte-identical output. Result: Rust was **slower**, not
  faster — 374-413 ns/call vs. 233 ns/call for the current `@turbo` version.
  The `ccall` round-trip itself is negligible (~2 ns, measured with a no-op
  control), so FFI overhead isn't the cause: this kernel is a gather-heavy
  stencil (5 pairs of non-contiguous periodic-neighbour lookups per
  element), and `LoopVectorization.@turbo` already emits SIMD-gather code
  tuned for exactly that pattern, which plain `rustc -O3` autovectorization
  doesn't match without hand-written SIMD intrinsics. At this grid size
  (96×48) and kernel granularity, adding a second-language build (Rust
  toolchain in CI, `ccall` bindings, binary-artifact distribution for a
  registered package) would add real complexity for a net slowdown — not
  pursued further.
- **Loop iteration order vs. column-major memory layout** (2026-08-13):
  `physics/ocean.jl:20,40,81` and `physics/radiation.jl:20,39` write
  `@turbo for i in 1:xdim, j in 1:ydim` (`j`, not the memory-contiguous
  `i`, is the fastest-varying/innermost index) while other kernels nest the
  other way (`i` innermost). Benchmarked both orderings standalone on two
  kernels: `seaice!`'s first loop body (cheap, few branches) and
  `deep_ocean!`'s loop body (heavier — divisions, more branches). For
  **both**, with `@turbo` present, source order made no reliable
  difference — medians swapped which order "won" between repeated runs
  (e.g. `deep_ocean!`-style: 9193 vs. 8100 ns one run, 4243 vs. 4329 ns the
  next), and `@code_native` showed packed SIMD loads (`vmovup*`), **no
  gather instructions**, for either order in both kernels.
  `LoopVectorization.@turbo` picks the true unit-stride axis for
  vectorization itself, independent of source loop nesting order — so this
  mismatch costs nothing today in the actual `@turbo`'d code.
  - However, the same `deep_ocean!`-style body **without** `@turbo` (plain
    `@inbounds` loop) is genuinely order-sensitive: correct
    (memory-contiguous) order measured 8400-9800 ns, mismatched order
    measured 18300-21200 ns — a reproducible **~2.2× slowdown** — and
    `@code_native` confirmed why: the mismatched plain loop compiles to
    `vgather`/scalar loads instead of packed `vmovup*`. So the general
    principle ("wrong order costs nothing") only holds *because* `@turbo`
    is present; it would be a real bug in any loop that lacks it.
  - Checked every double-nested spatial loop in `src/` for this: every
    site using the "wrong" `i`-outer/`j`-inner order already has `@turbo`
    (the two files above); every plain/`@inbounds`-only double loop
    (`model.jl:119` in one-time `init_model!` setup, `tendencies.jl:206-233`
    in rare `it==1` static-mask branches for `:regional_co2_ocean`/
    `:regional_co2_land_ice`) already uses the correct `j`-outer/`i`-inner
    order. No live instance of "no `@turbo` + wrong order" exists in the
    codebase today. Not changed; the `i`/`j` nesting inconsistency in
    `ocean.jl`/`radiation.jl` remains a readability nit only — but the
    ~2.2× plain-loop finding is worth keeping in mind for any *future*
    spatial loop that can't use `@turbo` (e.g. calls a function
    `LoopVectorization` doesn't support): get the nesting order right for
    those, since there's no compiler safety net without `@turbo`.

### 2.11 Parallelize `tendencies!`'s independent stages — ✅ done
All six stages (`SW`/`LW`/`hydro`/`deep_ocean`/`circulation!(Ta)`/
`circulation!(q)`) read only pre-timestep state and write disjoint buffers,
except the two `circulation!` calls which need separate workspace instances.

| Grouping | Time (2000-call avg) | Speedup |
|---|---|---|
| Sequential (all 6 stages) | 2747.5–2887.5 µs | 1× |
| 2-way: `circulation!(Ta)` \| `circulation!(q)` | 1540.0 µs | 1.88× |
| 6-way: every stage individually spawned | 1801.5 µs | 1.56× (worse) |
| **3-way**: `circulation!(Ta)` \| `circulation!(q)` \| rest bundled | **1376.0 µs** | **2.0×** |

Bare `Threads.@spawn`+`wait` overhead (~19.4 µs) exceeds several stages'
own cost, so spawning all 6 individually is a net loss. The 3-way split
captures the full ceiling at 2.0× (~1.0s/simulated year), verified
bit-identical vs. sequential.

**Implemented** in `tendencies!` (`src/tendencies.jl`): `ws_a`/`ws_q`
keyword args (default `=ws`, i.e. sequential — strictly opt-in), gated on
`Threads.nthreads() > 1 && ws_a !== ws_q`. `time_loop!`/`qflux_correction!`
forward the same kwargs; `greb_model!` is the only call site that allocates
distinct `ws_a`/`ws_q` and thus actually parallelizes. Verified: full suite
(319/319, incl. golden regression) passes unchanged at `-t 1`/`-t 2`/`-t 3`.

**Needs `-t 3`, not `-t 2`, to realize the full ceiling** — this is a true
3-lane split (`circulation!(Ta)` \| `circulation!(q)` \| rest), and with
only 2 OS threads the two spawned circulation tasks queue onto the single
non-main worker thread and run *sequentially* there, overlapping with "rest"
on the main thread but not with each other. Real 1-year `greb_model!`
timing (§2.12's turbo fixes already included in every number; `deepcopy`
of the large climatology arrays kept strictly outside every timed region —
an earlier pass of this benchmark hid this exact 2-vs-3-thread distinction
by leaving `deepcopy` inside the timed block, diluting the real
per-thread-count differences with a large constant):

| Threads | Time | vs. `-t 1` |
|---|---|---|
| `-t 1` | 1.912s | 1× |
| `-t 2` | 1.363s | 1.40× |
| **`-t 3`** | **1.170s** | **1.63×** |
| `-t 4` | 1.294s | no further gain (within noise) |

`-t 3` is the minimum thread count that gives the full benefit — `-t 4`+
adds nothing further since only 3 concurrent pieces of work exist per
timestep. **~2.31× total** vs. the original ~2.7s/year baseline at `-t 3`,
consistent with (slightly under) the ~2.8× the isolated integration test
below measured, the remaining gap being `output!`/`diagnostics!`/
bookkeeping outside `tendencies!` that this change doesn't touch (Amdahl's
law, given `tendencies!`'s ~70–75% share, §2.10).

Other identified, unimplemented opportunities (not benchmarked further):
ensemble/parameter-sweep parallelism via `Threads`/`Distributed`/GPU arrays
(likely the bigger practical payoff given how the model is actually used —
sensitivity sweeps, CO₂ scenarios); low-cost `Float32`/GPU-array swaps via
multiple dispatch (the physics kernels have no hardcoded `Float64`, only
`CirculationWorkspace`'s fields do, per §2.3); AD-based calibration via
`Enzyme.jl`/`ForwardDiff.jl` for gradient-based parameter tuning.

**Pre-implementation integration test** (kept for methodology reference):
combining the 3-way split with §2.12's turbo fixes in an isolated
`tendencies!`+`time_loop!` harness measured ~2.8× (2.83–2.88s →
1.00–1.04s) — higher than the ~1.78× (2.7s → 1.513s) the real
`greb_model!` shows post-implementation, since the isolated harness doesn't
carry `output!`/`diagnostics!`/bookkeeping overhead. One measurement gotcha
caught along the way, since fixed in the real implementation too: the
first, un-warmed-up timing showed the combo 2.5× *slower*, purely because
its new functions hadn't been JIT-compiled yet (unlike `greb_model!`'s
`PrecompileTools`-baked ones) — resolved with a warm-up call before timing.

---

## 3. Testing, tooling & reproducibility

- **Golden/snapshot regression test** ✅ done: real 1yr control + 1yr scenario
  `:full_model` run against `greb_dataset_jld2/`, checked against reference
  monthly global-mean `Ts`/`Ta`/`q` within `atol=1e-6`; `@test_skip`s when the
  dataset isn't present. `RUN_GOLDEN` env var (default on) additionally
  skips it locally even when the dataset is present.
- **Branch-coverage tests**: `log_eva`/`log_rain` sweep uses a zipped 5-run
  sweep instead of a 20-run cross product (independent branches, same
  coverage, half the runtime); a direct test confirms `log_eva` values
  produce genuinely different output.
- **CI** ✅, **`[compat]` bounds** ✅, **dead dependencies removed**
  (`NCDatasets`/`StaticArrays`/`Statistics`) ✅, **misplaced docstrings**
  fixed (§0.19) ✅, **benchmark suite** (`benchmark/run_benchmarks.jl` +
  `claude/BENCHMARKS.md`) ✅, **`Documenter.jl` site** ✅ (`docs/`, deployed
  to GitHub Pages on push to `main` via `.github/workflows/docs.yml`).
- **Testing gaps** ✅ closed: `qflux_correction!`'s asymmetric `Ta` handling
  confirmed intentional (matches Fortran, `greb.model.mscm.f90:540-597`),
  regression test added; the ~20 `:experiment` symbols unreachable via
  `create_experiment_config` now covered by a direct-construction sweep;
  `build_monthly_climatology`/`apply_scenario_anomalies` have direct unit
  tests; both JLD2 loaders' "file present" branches covered via synthetic
  `mktempdir`+`jldopen` datasets.

### 3.0 Parallelizing the test suite across processes — ✅ done
Splitting `runtests.jl`'s 24 testsets into 2 process-level shards (not
`Threads.@threads` — `Test.jl`'s testset tracking isn't thread-safe) along
an existing boundary:

| | Time |
|---|---|
| Sequential (one process) | 106.5s |
| Split A alone (17 cheap testsets) | 34.9s |
| Split B alone (7 heavy testsets, incl. golden regression) | 65.1s |
| **Parallel (2 concurrent `julia` processes)** | **68.1s** |

**1.56×** on 14 available cores — short of 2× because the split is
unbalanced (B takes ~2× A). A rebalanced/finer split would likely get
closer to the ceiling.

**Implemented**: `test/runtests.jl`'s 24 testsets are now split into
`run_light_tests()` (17) and `run_heavy_tests()` (7, incl. golden
regression) at this exact boundary, gated by `GREB_TEST_SHARD`
(`"all"`/`"light"`/`"heavy"`, default `"all"` — plain `Pkg.test()`/`]test`
is unaffected). `.github/workflows/ci.yml` gained a `shard: [light, heavy]`
matrix axis. Verified locally: 142 + 177 = 319, matching the unsharded
total exactly; both shards pass independently.

---

## 4. Findings needing domain/design input — documented, not fixed

### 4.1 Seven config switches wired to nothing ✅ resolved
`log_vapor_drsp`/`log_ice_drsp`/`log_ice_dmc`/`solar_multiplier` had zero
Fortran basis — **removed** from `PhysicsConfig`. `log_tsurf_ext`/
`log_hwind_ext`/`log_omega_ext` are declared in Fortran but genuinely dead
there too — **kept**, with a comment, for structural fidelity.

### 4.2 `log_clim` doesn't select the ERA/NCEP dataset its label implies ✅ resolved
`cfg.log_clim` (hydrology regression coefficients) and `load_greb_jld2!`'s
`dataset` kwarg (ERA vs. NCEP file selection) are fully disconnected. Per
user decision, kept orthogonal (a caller can combine `log_clim=1` with
`dataset=:era` for a sensitivity experiment) and documented the relationship
with cross-referencing comments on both (`src/config.jl`, `src/io.jl`).

### 4.6 More `const`-ify candidates — resolved by §1.1
The original target (module-level globals never rebound) no longer exists;
nothing to do.

### 4.8 `forcing()`'s dispatch style is inconsistent with the codebase — deferred by decision
Its ~180-line `if/elseif` chain over `cfg.experiment` is the odd one out
next to the Dict-based dispatch used elsewhere. **User decision: leave
as-is** — zero known bugs after two Fortran-audit passes, and refactoring
risks breaking the "missing branch = valid no-op" semantics §0.6 protects,
for a style-only win.

### 4.9 Testing gaps — ✅ closed, see §3

### 4.10 CI has no lint/format or doc-build check — deferred by decision
**User decision: skip for now** — `JuliaFormatter` isn't present anywhere in
the repo; adding one means picking a style and reformatting the whole
codebase, bigger scope than a mechanical addition.

### 4.11 `output.jl`'s `global_mean` was an unweighted grid mean ✅ fixed
Fortran's `gmean()` (`greb.model.mscm.f90:1497-1513`) area-weights by
`cos(lat)`; Julia's plain `sum/count` didn't. Fixed by reusing the existing
`dxlat_grid` (already `cos(lat)`-proportional) as the weight. Diagnostic-only
(`println`, never stored in `MonthlyRecord`) — golden regression test
unaffected.

---

## 5. Next up

Nothing blocking. Remaining open items:
- §4.8 (`forcing()` dispatch style) and §4.10 (CI lint/format) — both
  explicitly deferred by user decision, not forgotten.
- Ensemble/GPU parallelism and AD-based calibration (§2.11) remain
  documented ideas only, not benchmarked.
- §2.14's 3 investigated-not-implemented findings (CMIP5/ENSO anomaly
  reload caching, `postprocess.jl`'s Symbol dispatch, `CirculationWorkspace()`
  allocation count) — each needs a design decision or isn't worth it at
  current scale, not quick safe wins.
- §8's 5 newly found bugs remain documented, not fixed (explicit user
  decision to keep that pass report-only).

Landed across the last two passes: §2.12's original 4 `@turbo` fixes plus
its 4 follow-up fixes, §2.13's `forcing()` mask-recomputation fix, §3.0's
test-suite sharding, and §2.11's 3-way `tendencies!` thread split (needs
`-t 3`, not `-t 2`, to realize its full ceiling — see §2.11's correction
note) — see each section above and `CHANGELOG.md` for the real (not just
benchmarked) numbers.

> ⚠️ Every performance change must be validated against a reference run —
> "faster" only counts if output is unchanged within tolerance. Every
> physics change must be validated against the Fortran reference — grep the
> actual subroutine text, don't recall it from memory.

---

## 6. Data organization & JLD2 compression — investigated, findings documented

Full comparison doc: `claude/DATA_ORGANIZATION_OPTIONS.md`. Everything below
is measured (real `jldopen`/`tar` runs), not estimated.

### 6.1 `Data/`'s layout doesn't match its own docs
`DATA_README.md`/`convert_greb_to_jld2.jl` expect an `input/` subdirectory
that doesn't exist on disk (files sit directly under `Data/`). Two
root-level files (`erainterim.evaporation.clim`,
`erainterim.omega.vertmean.nomean.clim`, 25.7 MiB) were confirmed unread by
any Fortran variant (grepped both available source trees) — **not a porting
gap**. Both removed from `Data/`; `DATA_README.md` updated.

### 6.2 JLD2 compression: real numbers
Built the real `greb_dataset_jld2/` (605.4 MiB, 53 files) and measured
`compress=true` on 6 representative files:

| File | Compression ratio | Read slowdown |
|---|---|---|
| `ncep.tsurf...clim.jld2` | 80.1% | 9.5× |
| `cmip5.tsurf...forcing.new.jld2` | 86.6% | 7.9× |
| `erainterim.tsurf.elnino.forcing.jld2` | 93.2% | 20.4× |
| `global.topography.jld2` (static) | 29.5% | 4.0× |
| `solar_radiation.clim.jld2` | 40.4% | 2.0× |
| `solar_eccentricity.jld2` (grouped) | 40.8% | 11.4× |

Whole-tree `tar -I 'gzip -9'`: 530.6 MB from 634.8 MB (83.6%, 42s),
consistent with the per-file numbers. **Conclusion: compress only for
distribution** (a separate `tar`/`xz` archive), never the copy
`load_greb_jld2!` reads — the ~16–22% size win isn't worth an 8–20× read-time
risk.

### 6.3 JLD2 grouping: by load-pattern, not by source ✅ flux-correction merge implemented
Tested merging 5 CMIP5 fields into one file: combined-file single-field
reads cost 2.3× more (24.0ms vs 10.3ms, ~14ms absolute — irrelevant for a
once-at-startup load), but reading all 5 combined is *faster* (45.0ms vs
50.2ms). Rule: combine when a group is always read together, not by
upstream source — `ncep.*`/`erainterim.*` are read in varying subsets
depending on `dataset=:ncep`/`:era` and active experiment, so combining by
source would wrongly couple independently-varying load conditions.

**Implemented**: `load_flux_corrections_jld2!` always loads all 3 flux
files together, matching the "always read in full" case — merged into
`flux_corrections.jld2` (`convert_flux_corrections` in
`scripts/convert_greb_to_jld2.jl`), measured **~35% faster** (19.89ms vs
30.78ms combined vs. separate), no size penalty. `test/runtests.jl` updated
to match; full suite re-run clean (319/319 pass).

---

## 7. Fortran switch/experiment validation

Cross-referenced every `PhysicsConfig` field and `cfg.experiment` branch
against two `greb.shell.mscm.f90` variants and four `run.greb.*.csh`
scripts. 7.1–7.3 were initially report-only by explicit user decision, then
actually implemented in a follow-up pass (2026-08-12); 7.4–7.6 remain
documentation-only findings.

### 7.1 `:rcp26`/`:rcp45`/`:rcp60` ✅ implemented
Fortran `log_exp=96/97/98` load their CO2 files directly; Julia's
`forcing()` used to error unconditionally for these three. Fixed: they now
share the existing `:ssp*`/`:historical_co2` `cfg.co2_scenario[year]` lookup
dispatch in `forcing()` (`tendencies.jl`), loaded at scenario start via
`_CO2_SCENARIO_SYMBOLS`/`_CO2_SCENARIO_KEY` in `model.jl` (`:rcp60` maps to
the jld2 key `"rcp6"`, `:rcp26`/`:rcp45` map 1:1) from the already-working
`ipcc_scenarios.jld2`. `create_experiment_config` gained a matching preset.
Test: the former "not yet implemented" placeholder loop in
`test/runtests.jl` was replaced with a real synthetic-`ipcc_scenarios.jld2`
run for each of the three symbols.

### 7.2 `:custom_co2` (EXP=100) ✅ implemented
Fortran's `log_exp=100` lets a user point at an arbitrary `year CO2`
text file. Fixed: `PhysicsConfig` gained a `custom_co2_path::String` field;
`src/io.jl` gained `load_custom_co2_scenario(path)`, a plain-text parser
matching Fortran's format (skips blank/`#`-comment lines); `model.jl` loads
it into `cfg.co2_scenario` at scenario start (erroring clearly if the path
is unset), reusing the same `forcing()` dispatch as §7.1/`:ssp*`.
`create_experiment_config(:custom_co2; co2_path=...)` sets the path.

### 7.3 `_dmc`/`_drsp` switch families — combined presets added ✅ implemented
Every individual switch was already a real, wired `PhysicsConfig` field;
`create_experiment_config` had no preset that named a whole family at once
the way the Fortran run scripts do. Fixed (pure discoverability/convenience
addition, no dispatch-logic changes — every switch already worked):
`create_experiment_config(:decon_mean_climate; ...)` mirrors
`run.greb.decon_mean_climate.csh` (`log_exp=1`) — the 6 `_dmc` switches plus
the 5 circulation switches (`log_ice`/`log_hdif`/`log_hadv`/`log_vdif`/
`log_vadv`), each keyword-overridable, defaulting to `true` exactly as the
Fortran script's namelist defaults. `create_experiment_config(:decon_2xco2;
...)` mirrors `run.greb.decon2xco2.csh` (`log_exp=10`) — the 5 `_drsp`
switches (`log_crcl_drsp` excluded, no Fortran counterpart in this family)
plus the same 5 circulation switches, with `co2_concentration` fixed at
680.0.

### 7.4 `:sst_plus1` has no Fortran `log_exp` counterpart — confirmed intentional ✅
`model.jl:273,374-378` implements it correctly as a Julia-only addition
(logic lives entirely in `model.jl`, `forcing()` needs no special case).
Confirmed deliberate by explicit user decision, not a mistranslated
experiment number — no action taken.

### 7.5 `:a1b_scenario`/`:a1b_enhanced` — hardcoded formula, not a data file
No `ipcc.scenario.a1b*.txt` exists anywhere; Julia hardcodes a
piecewise-linear 1950→2000→2050→2100 (310→370→520 ppm) approximation
(`tendencies.jl:80-90,118-128`) — reasonable given no source table was
supplied, but an approximation, not a validated match.

### 7.6 Confirmed correct / already-resolved — no action
`log_clim`'s Fortran-convention mismatch (already resolved, §4.2);
obliquity/eccentricity range (Fortran comment is stale, glob-based file
discovery already handles the real file set correctly);
`earth_sun_distance_pct`'s `dradius` translation (confirmed correct);
regional CO2 experiments' dynamic-vs-static mask asymmetry between the two
halves of the family (both paths produce correct masks, architectural note
only).

---

## 8. Bugs found in the 2026-08-12 sweep — ✅ all fixed

A fresh Fortran-vs-Julia audit (run in parallel with the in-progress
performance work in §2.11) found 5 real bugs, initially left undocumented-
but-unfixed to avoid colliding with concurrent performance edits, then
implemented in a follow-up pass once that settled. Every finding was
verified directly against `greb.model.mscm.f90` before being fixed, and
each has a regression test in `test/runtests.jl`. Full discovery method:
`CHANGELOG.md`.

### 8.1 `forcing()`'s regional-CO₂ ice mask used only January, not the annual mean ✅ fixed
`tendencies.jl:200,216` (`:regional_co2_ocean`/`:regional_co2_land_ice`) took
`icmn_ctrl[:, :, 1]` — January's ice cover only — where Fortran
(`greb.model.mscm.f90:1279-1283`) averages `icmn_ctrl` across all 12 months
before thresholding. Fixed: `icmn_ctrl1` is now
`dropdims(sum(icmn_ctrl, dims=3), dims=3) ./ size(icmn_ctrl, 3)` at both
call sites. Test: `forcing() regional-CO2 ice mask uses the annual mean, not
January (§8.1)`.

### 8.2 `hydro!`'s `log_eva==1` branch dropped a carried-over base gust term ✅ fixed
`physics/hydrology.jl:101-113`. Fortran's `abswind` already carries `+2.0²`
(land)/`+3.0²` (ocean) into the `log_eva` dispatch
(`greb.model.mscm.f90:710-712`); the `log_eva==1` branch (`:727-728`) adds
its own `144.²`/`7.1²` *on top of* that already-modified variable, so the
true combined constant is `4+144²` land, `9+50.41` ocean — Julia recomputed
wind from scratch and added only `144.0`/`50.41`. Fixed: `gust_land_1`/
`gust_ocean_1` now add the existing `gust_land`/`gust_ocean` (`4.0`/`9.0`)
constants already used by the `log_eva==-1` branch. Test: `hydro!
log_eva==1 gust includes Fortran's carried-over +2.0²/+3.0² base term
(§8.2)`.

### 8.3 `hydro!` applied a spurious extra clamp to `dq_rain` alone ✅ fixed
`physics/hydrology.jl:150-153` clamped `dq_rain` to `-0.9*q/Δt` inside
`hydro!`, duplicating (and pre-empting, with a different formula) the real
Fortran clamp on the *combined* post-`dt` total, which Julia already
implements correctly and separately in `output.jl:163-164` (the §0.18 fix).
Fortran's `hydro` subroutine (`greb.model.mscm.f90:687-761`) has no such
clamp. Fixed: removed the `min_dq`/clamp lines from `hydro!`'s water-vapor-
tendency loop, leaving `output.jl`'s combined-`dq` clamp as the sole clamp.
Test: `hydro! doesn't apply an extra -0.9q clamp to dq_rain (§8.3)`.

### 8.4 Control-run climatology baseline averaged every control year instead of only the final year ✅ fixed
`model.jl:308-321` (`ctrl_output`) accumulates a `MonthlyRecord` for every
month of every control year; `build_monthly_climatology`/
`compute_annual_ice_climatology` (`postprocess.jl`) averaged across the
*entire* vector, where Fortran's `output` subroutine
(`greb.model.mscm.f90:1422-1446`) only ever assigns into
`Tmn_ctrl`/`icmn_ctrl`/etc. during the **final** control year — a straight
per-month value from the last year, not a multi-year average. Silently
correct only when `RunSpec.ctrl==1` (the Julia default); wrong for any
multi-year control run (the norm in every real Fortran run script). Fixed:
both functions now take only the final 12 records of the input vector
(`records[max(1, n-11):n]`), falling back to `records[1]`/zeros for the
`n<12` edge case exactly as before. Golden regression test unaffected
(`RunSpec()`'s default `ctrl=1` makes final-year-only and multi-year-average
identical). Tests: `build_monthly_climatology/apply_scenario_anomalies` and
`compute_annual_ice_climatology` (both updated to assert the final-year
value, not the old multi-year mean).

### 8.5 Flux-correction spin-up loaded precomputed files, then immediately overwrote them ✅ fixed
`model.jl:281-290`: when `!cfg.log_topo_drsp && cfg.log_qflux_dmc`, Julia
called `load_flux_corrections_jld2!` and then unconditionally called
`qflux_correction!` right after, which recomputed and overwrote the exact
same `fields.TF_correct`/`qF_correct`/`ToF_correct` arrays — the load had
zero effect on the run. Fortran's equivalent
(`greb.model.mscm.f90:357-390`) is a mutually exclusive `if/elseif/else`:
compute fresh **or** read precomputed files **or** zero — never both. Fixed:
the spin-up's `if`/`if` pair is now `if`/`elseif`/`else`, matching Fortran's
structure. (The separate, lower-confidence question of whether
`log_topo_drsp`/`log_qflux_dmc` are even the right switches for this gate —
Fortran's real gate is `log_exp` — is the already-documented §7.3 finding
about `_dmc`/`_drsp` families lacking a combined preset; not re-litigated
here.) Test: `greb_model! flux-correction spin-up: loaded files aren't
overwritten by qflux_correction! (§8.5)`.

---

## 9. Julia General Registry publication readiness

### 9.1 Data distribution (Artifacts.jl) — deferred by decision

The model needs ~1.16 GB of climatology/scenario data (`Data/` 581 MB +
`greb_dataset_jld2/` 580 MB) to actually run, but neither directory is
git-tracked — both are gitignored, so the git tree at any tag (what `Pkg`'s
package tarball is built from) already excludes them. That's correct: a
General-registry package must stay small, and nothing in the registry
process requires or expects large data to live inside the package tree
itself.

What's missing is the other half — an automated way for a `Pkg.add("GREB")`
user to actually *get* the data. Today it's "ask the maintainer, or
regenerate from raw `.bin` files via `scripts/convert_greb_to_jld2.jl`"
(README's Input Data section). There is no `Artifacts.toml` and no
`artifact"..."` usage anywhere in the repo — every loader in `src/io.jl`
(`load_greb_jld2!`, `load_cc_anomaly_jld2!`, `load_solar_forcing_jld2`,
`load_co2_scenario_jld2`, etc.) takes a plain required `jld2_dir::String`
positional argument, joined with hardcoded subdirectory names
(`climatology/`, `static/`, `solar/`, `scenario/`).

**Deferred plan, not implemented in this pass:**
1. Build the `greb_dataset_jld2/` tree, tar it, and attach it as an asset on
   a GitHub Release (the repo already publishes releases — see `v0.1.0`,
   `v1.0.0`).
2. Generate `Artifacts.toml` at the repo root via
   `Pkg.Artifacts.create_artifact`/`bind_artifact!`, pointing at that release
   asset URL with its SHA256 + git-tree-sha1.
3. Change `src/io.jl`'s loaders to resolve `artifact"greb_dataset"` instead
   of requiring a manually-supplied `jld2_dir`, falling back to an explicit
   path override for anyone who wants to point at a local/custom dataset
   (e.g. during development, or a different reanalysis product).
4. This also closes a real coverage gap noted during the registry audit:
   `.github/workflows/ci.yml` never fetches `greb_dataset_jld2/`, so the
   "heavy" test shard's data-gated tests — including the golden-regression
   check — currently always hit the `@test_skip "greb_dataset_jld2/ not
   present"` branch in `test/runtests.jl` and silently never run in CI.
   Once artifact-based loading exists, CI can resolve the artifact like any
   other user and actually exercise that coverage.

### 9.2 Repo rename (`GREB-julia` → `GREB.jl`) — deferred by decision

Registry convention expects the GitHub repo backing a package named `GREB`
to be called `GREB.jl`; this repo is `GREB-julia`. Not an AutoMerge
blocker, but commonly expected by reviewers. Rename is deferred for now —
revisit before actually submitting the registration PR
(`docs/make.jl`'s `deploydocs(repo = "github.com/EnvDroneSense/GREB-julia.git", ...)`
and any other hardcoded repo URLs would need updating alongside the
GitHub-side rename).

### 9.3 Registry audit summary

A full AutoMerge-requirements pass (2026-08-12) otherwise found the package
in good shape: `Project.toml` has compat bounds on every dependency
including `julia` itself; `LICENSE` (MIT) exists; `test/runtests.jl` runs
cleanly with zero external data required (data-gated tests skip gracefully,
per §9.1 above); and a live check against `JuliaRegistries/General`
confirmed the name `GREB` is not yet taken. Fixed as part of this pass:
version reset to `1.0.0` for a clean first registration (the package had
briefly been bumped to `1.1.0`), a new `CITATION.cff`, a
`Pkg.add("GREB")` forward-reference in the README, and a git-history
rewrite (`git filter-repo`) to strip ~26 MB of now-deleted `Data/input/*.bin`
files that were still permanently baked into old commits.
