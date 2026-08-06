# GREB.jl — Potential Improvements (Structure & Performance)

Observations and concrete suggestions for the `GREB` module. The model has
been split from a single 2245-line `src/GREB.jl` into topical files under
`src/` — see §1.2. Three bug-hunting passes have gone through this list so far
(2026-08-05/06, 2026-08-06 §4, and 2026-08-06 third pass in §0); nothing else
here has been applied except where marked **✅ done**.

The single theme underneath most of what used to remain: **the model ran on
~40 mutable module-level globals** (grid arrays, climatology fields, work
buffers, flux corrections). That one design choice drove most of the
structural *and* performance issues below — it's why the module used to eat
~750 MB at load, why the code was hard to thread or test, and why the same
"config flag silently disconnected from behavior" bug shape kept recurring
across all three bug-hunting passes. **§1.1 (the state-struct refactor) is
done** — see that section for what changed and how it was validated. §1.1
also unblocked §2.4 (threading); §2.4 itself was attempted, benchmarked, and
reverted — see that section for the numbers.

---

## 0. Known bugs — fixed ✅

Six real, previously-shipped bugs have been found and fixed across two
bug-hunting passes over this codebase. All are fixed with regression tests
added in `test/runtests.jl`.

### First pass (2026-08-05/06)

Three bugs, found while implementing §0 from the prior version of this doc
and while re-verifying §1.1's claims before splitting files (tracing every
global's reads/writes to know which were live vs. dead surfaced the other two).

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

### Second pass (2026-08-06)

A fresh, independent read through every file in the (by now split) `src/`
tree turned up three more real bugs — all of the same "a config flag is
silently disconnected from the behavior it's supposed to gate" shape as
0.2/0.3.

#### 0.4 `circulation!`'s `do_conv` gate was inverted
`src/circulation.jl` precomputes five flags before the per-substep loop; four
read `cfg.log_X == 1`, but `do_conv` alone read `cfg.log_conv == 0`:
```julia
do_diff_v = cfg.log_vdif == 1 && h_scl == z_vapor
do_diff_h = cfg.log_hdif == 1 && h_scl == z_air
do_adv_v  = cfg.log_vadv == 1 && h_scl == z_vapor
do_adv_h  = cfg.log_hadv == 1 && h_scl == z_air
do_conv   = cfg.log_conv == 0 && h_scl == z_vapor   # inverted
```
`log_conv::Bool` defaults to `true`, so moisture-flux-convergence
(Eq. 18, Stassen et al. 2019) was silently **off** for every
default-configured run — backwards from the flag's own name and default.
Confirmed identical in `notebooks/GREB_julia.jl`, so this predates the
package extraction; it just had never been read closely before. **Fixed**:
converted all five lines to direct boolean truthiness (`cfg.log_vdif`, …,
`cfg.log_conv`) — this fixes the inversion and a style inconsistency (mixing
numeric `==1`/`==0` comparisons with truthy tests elsewhere in the same file)
in one change. Regression test constructs a nonzero synthetic `omegaclim` and
asserts `log_conv=true` and `log_conv=false` now produce different
`circulation!` output (previously identical).

#### 0.5 `co2_part` module-global leaked between runs
`co2_part` (`src/state.jl`) is read every timestep in `src/physics/radiation.jl`
but only ever written by `forcing()`'s `regional_co2_*` branches
(`src/tendencies.jl`), and nothing reset it to `1.0` at the start of a run.
Running a `:regional_co2_nh`-style experiment and then any other experiment
in the same Julia session silently carried the first run's halved CO₂ mask
into the second run's flux-correction, control, *and* scenario phases —
exactly the kind of hazard repeated re-runs in the interactive Pluto notebook
would trigger. **Fixed**: `init_model!` (`src/model.jl`) now resets
`co2_part .= 1.0` at the start of every `greb_model!` call. Regression test
runs a regional experiment, then a plain one, and asserts the second run's
mask is back to all-`1.0`.

#### 0.6 Silent fallback on an invalid `log_rain`
`set_hydrology_parameters!` (`src/config.jl`) did
`get(HYDRO_PARAMS, cfg.log_rain, (1.0,0.0,0.0,0.0))` — a typo'd `log_rain`
silently behaved like "Original GREB" with no warning, unlike
`create_experiment_config`'s own `error()` on an unknown `experiment` symbol
in the same file. **Fixed**: `error()`s now if `cfg.log_rain` isn't a key of
`HYDRO_PARAMS`, instead of silently defaulting.

*(Considered but **not** fixed: adding a similar `error()` for an
unrecognized `cfg.experiment` symbol in `forcing()`/`init_model!`'s
`if/elseif` chains. Ruled out — both chains are missing branches for many
*valid* experiments, including the default `:full_model`, which intentionally
fall through and do nothing; a blanket final `else error(...)` would break
those, not just catch typos. A real fix needs a canonical list of valid
`experiment` symbols to check against up front, which is a design decision,
not a mechanical one — tracked in §4.)*

Second-pass audit performed the same way as the first: grepped every
`cfg.field` write against every read site, and every function-local `global`.
Two more dead globals and one more unused dependency turned up as a side
effect — see §1.6 and §3.

### Third pass (2026-08-06): direct audit against the Fortran reference

The user asked for a line-by-line comparison against the original Fortran
model (`greb.model.mscm.f90`) to confirm the Julia port hadn't drifted
conceptually — either from the notebook→package extraction or from the first
two passes' own fixes. That audit confirmed §0.4's `do_conv`/`log_conv` fix
was correct (Fortran's own switch-family convention and the original
notebook's preset groupings both independently confirm `log_conv::Bool=true`
means "on") and turned up four more real bugs, all fixed with regression
tests and checked against the actual Fortran subroutine text (not recalled
from memory):

#### 0.7 `hydro!`'s `log_eva` modes 1 and 2 didn't exist
`src/physics/hydrology.jl`'s evaporation section had branches for
`log_eva ∈ {-1, 0}` (both correct) but silently reused mode `-1`'s formula
for `1`/`2`/anything else, instead of Fortran's two genuinely distinct
parameterizations (different wind-gust terms and land/ocean coefficients).
**Fixed**: implemented the real mode-1 (`0.04`/`0.73` land/ocean coefficients,
`144.0`/`50.41` gust terms, wind from `uclim`/`vclim`) and mode-2 (`0.56`/
`0.79` coefficients, `81.0`/`16.0` gust terms, wind from `wsclim`) formulas,
plus an `error()` fallback for any other value. Regression test extended the
existing `log_eva`×`log_rain` branch-coverage sweep to include mode `2`, and
added a direct assertion that all four `log_eva` values produce different
`hydro!` output for the same synthetic input (the branch-coverage sweep alone
only checks output *shape*, which is exactly why this went undetected).

#### 0.8 `init_model!`'s `_drsp` climatology overrides were mismapped
Fortran's decon-2xCO2 sensitivity overrides (`if(log_cloud_drsp==0)
cldclim=0.7`, `log_humid_drsp` → `qclim=0.0052`, `log_ocean_drsp` →
`mldclim=d_ocean`) were, in `src/model.jl`'s `init_model!`: the `cldclim=0.7`
override missing entirely, the `qclim` override gated on `log_vapor_dmc`
instead of `log_humid_drsp`, and the `mldclim` override gated on
`log_ocean_dmc` instead of `log_ocean_drsp`. A user running the documented
"deconstruct 2×CO2" experiment and toggling the real `_drsp` switches got
silently wrong or no-op behavior. **Fixed**: added the missing override,
corrected the other two switch names. Regression test constructs a config
with all three `_drsp` switches off and asserts each climatology array
matches its documented constant.

#### 0.9 `LWradiation!`'s `log_atmos_dmc` asymmetry
Fortran only zeros `LWair_down` when `log_atmos_dmc==0`; `LWair_up` is
snapshotted *before* that zeroing and keeps its full computed value (the
switch decouples the surface from atmospheric downwelling feedback without
touching the atmosphere's own emission term). Julia's `LWradiation!` zeroed
both. **Fixed**: removed the extra `LW_up .= 0.0`. Only affects experiments
that explicitly disable `atmos_dmc`; no change to `:full_model`. Regression
test confirms `LW_down==0` but `LW_up` unchanged when the switch is off.

#### 0.10 Paleo/orbital solar forcing was never wired up
Fortran swaps in an entirely different `(ydim, nstep_yr)` solar table for
paleo/orbital experiments (`sw_solar = sw_solar_scnr`, log_exp 30/31/35/36).
Julia had the loader (`load_solar_forcing_jld2`) but nothing ever called it
from `greb_model!` — `:paleo_231kyr`/`:paleo_solar_modern_co2`/`:obliquity`/
`:eccentricity` silently ran on the ordinary modern seasonal solar cycle.
**Fixed**: `greb_model!` now swaps in the right table at the start of the
scenario run for those four experiments (via a new `orbital_index`
`PhysicsConfig` field for selecting an eccentricity/obliquity table row), and
`:earth_sun_distance` now applies Fortran's scalar rescale (`rS0 =
(1/(1+0.01·dradius))²`) via a new `earth_sun_distance_pct` field, the same
mechanism already used for `:solar_plus27`. Regression test writes a
synthetic solar-table JLD2 fixture and confirms the swap branch runs and
`sw_solar` is restored afterward (see §1.1 — this hazard is structurally
gone for the default fresh-`fields`-per-call path now anyway).

Audited the same way as the first two passes (grep every switch write vs.
read site) plus a fresh, independent read of the Fortran subroutine bodies
for every function these four bugs touch. No further deviations found in
`hydro!`, `LWradiation!`, `SWradiation!`, `seaice!`, `deep_ocean!`,
`circulation!`, or `init_model!`.

---

## 1. Structure

### 1.1 Replace module-global mutable state with state structs  ✅ done
The model's ~40 mutable module globals (climatology fields, derived grid
fields, flux corrections, the `co2_part` mask, `sw_solar`) are gone. Two
structs now hold everything, `src/state.jl`:
```julia
mutable struct ClimateFields   # loaded climatology, grid fields, flux corrections, masks — 37 fields
    z_topo::Matrix{Float64}; ...; Tclim::Array{Float64,3}; uclim::Array{Float64,3}; ...
end
mutable struct ModelState      # runtime solar-forcing scalar + annual-mean diagnostic accumulators
    sw_solar_forcing::Float64; Tsmn::Matrix{Float64}; ...
end
```
Every physics/circulation/tendencies/output/model function gained an explicit
`fields::ClimateFields` parameter (and `state::ModelState` where it touches
the solar-forcing scalar or the annual-mean accumulators) instead of closing
over module globals. `load_greb_jld2!`/`load_flux_corrections_jld2!` now
*return*/*fill* a `ClimateFields` rather than mutating globals in place.
`greb_model!(time_flux, time_ctrl, time_scnr, cfg; jld2_dir="", fields=ClimateFields())`
defaults to a fresh instance per call — the old cross-run leak hazard
(§0.5's `co2_part` bug, and a similar risk for a paleo run's swapped
`sw_solar` table) is now structurally impossible for that default path;
`init_model!`/`greb_model!` still defensively reset `co2_part`/`sw_solar` for
callers who explicitly reuse one `fields` instance across multiple runs (e.g.
to avoid reloading real climatology).

Migrated one file at a time, bottom-up (`physics/*.jl` → `circulation.jl` →
`tendencies.jl` → `output.jl` → `model.jl` → `io.jl`), running the full test
suite after each. `CirculationWorkspace`/`MonthlyAccumulator` already followed
this pattern for their own slices — the same shape was just extended to
everything else.

**Validated against real data** (`greb_dataset_jld2/`), not just structurally:
a 1-year control+scenario `:full_model` run's `MonthlyRecord` output before
and after the refactor agrees to within **~7e-12 absolute** on fields with
O(1–300) magnitude (Ts, precip, SW, …) — floating-point reassociation noise
from a different compiled instruction order (typed struct fields vs. globals
changes what LLVM/`@turbo` can vectorize), not a behavior change. Exactly
bit-identical turned out to be the wrong bar for a change that also alters
codegen; agreement to 13+ significant figures over 1460 accumulated timesteps
is the right one.

**Unplanned bonus**: this was pursued for correctness/testability, not speed,
but the same real-data run went from 49s/17.1 GiB allocated to **7.0s/100 MiB
allocated** — a ~7x runtime, ~170x allocation reduction. Global mutable
bindings in Julia are type-unstable at every access; struct fields with
concrete types are not. This also incidentally fixes §2.2 (see below) and
unlocks §2.4 (threading) — see those sections.

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

### 1.3 Collapse giant positional argument lists into structs — still deferred
Several kernels still thread ~10–20 positional arrays — `output!`,
`diagnostics!`, `time_loop!` (`src/output.jl`) all still take Ts0/Ta0/To0/q0/
albedo/ice/precip/... as bare positional args, now *alongside* the new
`fields`/`state` parameters from §1.1 rather than instead of them. §1.1
deliberately only moved what was actually **global** (climatology, grid
fields, flux corrections, the annual-mean accumulators) — `Ts`/`Ta`/`To`/`q`
and the per-month working arrays were already ordinary function-local
variables, not part of the global-state bug class, so leaving them positional
was in scope but not required by §1.1. Collapsing them into a
`SurfaceState`-style struct is still a real readability win and still
deferred — same error-prone-positional-mixup rationale as before.

### 1.4 Turn run durations into a `RunSpec`, not positional ints — still deferred
`greb_model!(time_flux, time_ctrl, time_scnr, cfg; jld2_dir="", fields=ClimateFields())`
(`src/model.jl`) still takes three bare ints whose order is easy to swap — §1.1
added a `fields` keyword alongside them but didn't touch this. A
`RunSpec(; flux=0, ctrl=1, scnr=1)` keyword struct remains a good follow-up.

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

`sw_solar_forcing_state = Ref(1.0)` — now `ModelState.sw_solar_forcing`, moved
as part of §1.1.

**Second-pass addendum**: two more, found the same way — `co2_part_scn` and
`sw_solar_forcing_data` (`src/state.jl`), both declared and never read or
written anywhere else. Removed.

---

## 2. Performance

### 2.1 ~~`const`-ify the constants~~ ✅ done
Grid dimensions (`xdim`, `ydim`, `nstep_yr`, …) and the physical constants
(`σ`, `ρ_ocean`, `grav`, `cp_air`, …), now in `src/constants.jl`, are all
`const`. No further action needed.

### 2.2 Stop allocating ~750 MB at module load ✅ done
Fixed as a direct consequence of §1.1: the ~28 large arrays that used to be
top-level `zeros(Float64, xdim, ydim, nstep_yr)` calls in `src/state.jl` (run
once at module load, ~27 MB each ⇒ ~750 MB just from `using GREB`) are now
fields inside `ClimateFields`, allocated only when `ClimateFields()`/
`load_greb_jld2!`/`greb_model!` is actually called.

### 2.3 Consider `Float32` for climatology fields — deferred
The JLD2 files store `Float32` data and are up-converted to `Float64` on load
(`src/io.jl`). Keeping the large 3D fields `Float32` (or making element type a
parameter) halves memory (~375 MB) and roughly doubles SIMD lane throughput
under `@turbo`. Now that `greb_dataset_jld2/` is available locally, this can
be benchmarked against real data (see `benchmark/run_benchmarks.jl`) rather
than estimated.

### 2.4 Thread the spatial operators — tried, measured, reverted
Implemented `Threads.@threads :static for k in 1:ydim` on `diffusion!`'s and
`advection!'s` outer latitude-row loop (`src/circulation.jl`) — each row
writes a disjoint output column from read-only input, genuinely
embarrassingly parallel in principle. One real race had to be fixed first:
both functions' polar sub-stepping branches reused a single shared scratch
vector (`ws.T1h`), read-modify-written across several `k` rows near each
pole; fixed with a per-thread `T1h_threads::Matrix{Float64}` column, sized by
`Threads.maxthreadid()` (not `Threads.nthreads()` — Julia 1.9+'s
default+interactive threadpool split means `threadid()` values inside a
`:static` `@threads` loop can exceed `nthreads()`; sizing by `nthreads()`
produced a real `BoundsError` on the first `-t 4` run, caught before it went
further). Validated bit-identical (not just within tolerance) across thread
counts via a dedicated subprocess-based test — confirmed no race remained.

**Benchmarked, and reverted** — the grid is too small (`xdim=96, ydim=48`,
~1µs of actual work per row) for `Threads.@threads`'s per-call scheduling
overhead to pay off at any thread count tried:

| kernel | `-t 1` (before → after) | `-t auto`/14 threads (before → after) |
|---|---|---|
| `diffusion!` | 45–117µs, 0 allocs → 54µs, **7 allocs** | → **220µs, 72 allocs (4x slower)** |
| `advection!` | 23–74µs, 0 allocs → 68µs, **7 allocs** | → **206µs, 72 allocs (3x slower)** |
| `circulation!` | 1.9ms, 0 allocs → 1.9ms, **336 allocs** | → **10.9ms, 3456 allocs (5.7x slower)** |
| `tendencies!` | 4.2ms, 48 allocs → 4.2ms, **720 allocs** | → **23.7ms, 6960 allocs (5.6x slower)** |

Full numbers in `claude/BENCHMARKS.md`'s "-t 1"/"-t auto" entries for this
pass. Even at `-t 1`, `Threads.@threads` itself allocates (task/partition
bookkeeping) even with nothing to parallelize across — a real, if small,
regression. At `-t auto`, `circulation!`'s `ntime`-iteration sub-step loop
calls `diffusion!`/`advection!` repeatedly, so the per-call threading
overhead compounds every sub-step: 5–6x slower end-to-end, not a marginal
loss. The plan going into this explicitly flagged threading overhead as
possibly not paying off at this grid size and included a benchmark gate for
exactly this reason — the numbers came back unambiguously negative, and
unlike the plan's anticipated fallback (gate `@threads` behind
`nthreads() > 1`), that wouldn't actually help here: `nthreads() > 1` is
precisely the regime that's *worse*, not better. Reverted `diffusion!`/
`advection!` to their original plain `for k in 1:ydim` loops; `ws.T1h`
reverted to a single `Vector{Float64}` (no per-thread buffer needed without
threading). Kept the unrelated find made along the way: `CirculationWorkspace`'s
`dTxh` field was dead code (zero reads/writes anywhere in `src/`) and stayed
removed.

**Takeaway for any future threading attempt**: this grid (96×48) and these
kernels' per-row cost (~1µs) are the wrong shape for `Threads.@threads`.
A lower-overhead approach (e.g. `Polyester.@batch`, which avoids the
task-scheduling machinery `Threads.@threads` uses) might fare better, or
threading might only pay off on a substantially larger grid. Don't re-attempt
with plain `Threads.@threads` at this resolution without a reason to expect
a different result.

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
  branch-coverage loops over `log_eva`/`log_rain` (added specifically because
  §0.1 hid in an untested branch) and a direct unit test asserting
  `set_hydrology_parameters!` writes the right values into `cfg` for every
  `log_rain`/`log_clim` combination (added because of §0.2/§0.3).
- **Golden/snapshot regression test** ✅ done: a real 1yr control + 1yr
  scenario `:full_model` run against `greb_dataset_jld2/`, snapshotting
  monthly global-mean `Ts`/`Ta`/`q` and asserting they match reference values
  within `atol=1e-6` — well above the ~1e-12 float-reassociation noise §1.1
  already measured, tight enough to catch a real regression. Skips (via
  `@test_skip`) when the dataset isn't present. Added as an end-to-end
  tripwire before attempting §2.4's threading change; kept afterward as a
  general-purpose regression test for any future kernel change.
- **Test-suite runtime**: the `log_eva`/`log_rain` branch-coverage sweep used
  to run a full `greb_model!` year for every one of the 4×5=20 combinations,
  ~65% of all simulated model-years in the suite, for no extra bug-catching
  power — `log_eva` and `log_rain` gate independent branches of `hydro!`, so
  the cross product exercises no interaction the two axes don't already cover
  separately. Changed to a zipped sweep (cycling the shorter list): 5 runs
  instead of 20, still hitting every value of both axes at least once — same
  regression-catching guarantee, 75% less work. Also dropped a redundant full
  `greb_model!` call from the thread-equivalence test's subprocess worker
  (added, then found unnecessary, in the same pass) — the `diffusion!`/
  `advection!` bit-identity check alone already covers the actual race risk.
- **Per-kernel unit tests**: not added as `Test.jl` unit tests, but
  `benchmark/run_benchmarks.jl` now exercises every hot kernel individually
  with synthetic inputs (see below) — a natural base to extend with
  correctness assertions (e.g. `diffusion!` conserves the field mean).
- **CI** ✅ done: `.github/workflows/ci.yml` runs `Pkg.test()` on a Julia
  `1.9`/`1` × `ubuntu-latest` matrix on push/PR to `main`.
- **`[compat]` bounds** ✅ done: `Project.toml` now bounds `julia`, `JLD2`,
  `LoopVectorization`, and `PrecompileTools` at their currently-resolved
  versions.
- **`NCDatasets`, `StaticArrays`, and `Statistics` dependencies** ✅ done: all
  three were declared but unused (`NCDatasets`: no `NCDataset(...)` call
  anywhere in the package; `StaticArrays`: no `SVector`/`SArray`/`SMatrix`/
  `MVector` anywhere — the "static longitude indices" comment was
  aspirational, they're plain `Vector{Int}` comprehensions; `Statistics`,
  found in the second pass: no `mean`/`std`/`var`/`median`/`cor`/`cov` call
  anywhere in `src/`). All three dropped from `Project.toml [deps]` and their
  `using` statements removed. (The notebook still uses `NCDatasets` in its own
  environment, untouched.)
- **Misplaced docstrings** ✅ done (second pass): three exported functions —
  `load_flux_corrections_jld2!`, `load_greb_jld2!` (`src/io.jl`) and
  `convergence!` (`src/circulation.jl`) — had their docstring written as a
  bare string literal *inside* the function body (after `function foo(...)`)
  instead of immediately before the `function` keyword. A string in that
  position is just a discarded statement, not a real docstring — `?convergence!`
  showed nothing. Moved all three above their `function` line.
- **Benchmark suite** ✅ done (new, beyond the original list): added
  `benchmark/run_benchmarks.jl` (own environment, mirroring `test/Project.toml`,
  so `BenchmarkTools` never becomes a package dependency) — `@benchmark`s
  every hot kernel against real data when `greb_dataset_jld2/` is present, and
  appends a dated, tagged, git-commit-stamped section to `claude/BENCHMARKS.md`
  every run, so before/after comparisons accumulate across optimization
  passes instead of being overwritten. One baseline entry exists so far
  (2026-08-05, pre-first-pass); a post-fixes comparison entry is added once
  this second pass's fixes are committed.
- **Data acquisition**: `scripts/convert_greb_to_jld2.jl` converts raw `.bin`
  files into the `.jld2` layout; the raw files themselves are still
  "available on request" per `DATA_README.md`. Unchanged this pass.
- **Docs**: a `Documenter.jl` site would suit the physics-heavy API; the
  notebook can become a rendered tutorial page. Unchanged this pass.

---

## 4. Second-pass findings — documented, not fixed

Found during the same second-pass read-through as §0.4–0.6. These either need
domain physics expertise to fix correctly, touch public API/behavior more
broadly than a mechanical change should, or naturally belong with the
already-deferred state-struct refactor (§1.1) — so they're recorded here
rather than fixed inline.

### 4.1 Seven config switches are wired to nothing
`PhysicsConfig` fields that are set (including by `create_experiment_config`)
but never read by any physics function, so toggling them has zero effect on
model output: `solar_multiplier` (the actual `+27 W/m²` solar-forcing effect
is computed independently inside `forcing()`, `src/tendencies.jl`, not via
this field), `log_vapor_drsp`/`log_ice_drsp` ("CO₂ Response Switches" —
`log_clouds_drsp`/`log_humid_drsp` are now wired, fixed as §0.8),
`log_tsurf_ext`/`log_hwind_ext`/`log_omega_ext` ("External Forcing"
switches), and `log_ice_dmc` (distinct from `log_ice`, which *is* wired in
`radiation.jl`/`ocean.jl`). All are documented, notebook-exposed controls —
this is the same "config field disconnected from behavior" bug class as
§0.2–0.6/§0.8, but wiring 7 features up (or removing them) is a scope/design
decision, not a mechanical fix.

### 4.2 `log_clim` doesn't select the ERA/NCEP dataset its label implies
The notebook documents a Hydrology control "Climatology (ERA/NCEP)" bound to
`log_clim`. In the module, `cfg.log_clim` only ever affects one thing
(`src/config.jl`: swaps four regression coefficients when
`log_rain==0 && log_clim==1`) — it does **not** select which dataset gets
loaded. The actual ERA-vs-NCEP choice is the separate, disconnected `dataset`
keyword to `load_greb_jld2!` (`src/io.jl`), set before `greb_model!` runs and
never connected to `cfg` at all. Fixing this properly means deciding whether
`cfg` should drive data loading (a real API question) — deferred alongside §1.1.

### 4.3 ~~`hydro!`'s 4 documented evaporation modes are really 2~~ ✅ done
Fixed as §0.7 — modes `1`/`2` now implement Fortran's actual distinct
formulas instead of duplicating mode `-1`.

### 4.4 `forcing()` per-timestep overhead — partially fixed ✅
`forcing()` (`src/tendencies.jl`) now short-circuits with an early
`if cfg.experiment == :full_model; return (...); end` before the `elseif`
chain, avoiding the allocating `startswith(string(cfg.experiment),
"regional_co2_")` check on every timestep of the default, most common
experiment. Pure short-circuit — returns exactly what the old fallthrough
already computed for `:full_model`, verified no branch sets anything for it.

Still deferred: the regional-CO2 mask (`co2_part`) is still fully recomputed
every timestep inside the `regional_co2_*` branches even though it depends
only on static topography — hoisting that needs a new init-time hook and
only benefits those niche experiments, not the default path this fix
targeted.

### 4.5 `hydro!` recomputes per-call "constants" as locals ✅ done
`const_factor1/2/3`, `gust_land/ocean`, `cE_land/ocean`, `const_latent`
(`src/physics/hydrology.jl`) are now module-level `const`s
(`_HYDRO_CONST_FACTOR1`, etc.), computed once instead of as local variables
recomputed every timestep. Zero behavior change — same values, same source
constants (`ce`, `cq_latent`, `ρ_air`).

### 4.6 More `const`-ify candidates
Beyond the constants already made `const` (§2.1), several module globals are
verifiably never *rebound* (only mutated via `.+=`/`fill!`/broadcasting) —
the 11 `diagnostics!` accumulators (`Tsmn, Tamn, Tomn, qmn, amn, swmn, lwmn,
qlatmn, qsensmn, ftmn, fqmn`, `src/state.jl`) plus climatology fields like
`z_topo`/`Tclim`/`cldclim`. Marking these `const` is individually low-risk
(Julia errors immediately if a hidden rebind exists) but spans many locations
— deserves its own dedicated, careful pass rather than a piecemeal one.

### 4.7 Missing docstrings on the public API ✅ done
All ~22 exported functions/types that had no docstring now do: `SWradiation!`,
`LWradiation!`, `hydro!`, `seaice!`, `deep_ocean!`, `diffusion!`, `advection!`,
`circulation!`, `tendencies!`, `forcing`, `diagnostics!`, `output!`,
`time_loop!`, `init_model!`, `qflux_correction!`, `build_monthly_climatology`,
`apply_scenario_anomalies`, `compute_annual_ice_climatology`, `PhysicsConfig`,
`TimeState`, `MonthlyRecord`, `load_solar_forcing_jld2`. `greb_model!` itself
already had one from §1.1. Mechanical — brief description + argument/return
notes, following the existing `convergence!`-style format, no behavior change.

### 4.8 `forcing()`'s dispatch style is inconsistent with the rest of the codebase
Its ~183-line/~25-branch `if/elseif` chain over `cfg.experiment`
(`src/tendencies.jl`) is the odd one out next to the Dict-based dispatch used
elsewhere (`HYDRO_PARAMS` in `constants.jl`, `file_map` in `io.jl`). Natural
to address alongside the §1.1 state-struct refactor, since `forcing`'s
signature changes then anyway — not a standalone mechanical fix, and (per
§0.6's note) the chain's "missing branch = valid no-op" semantics need to be
preserved carefully by whatever replaces it.

### 4.9 Testing gaps
- **Scenario-run coverage**: was zero except §0.5's regression test; §0.10's
  new regression test also now runs a `time_scnr=1` scenario for
  `:paleo_231kyr`. Still mostly uncovered: `build_monthly_climatology`,
  `apply_scenario_anomalies`, and the `is_forced_boundary`/`is_sst_plus1`
  branches in `greb_model!` (only `is_orbital_exp`'s path is exercised so far).
- **`qflux_correction!`'s loop body has zero coverage**: every test uses
  `time_flux = 0`, so `for it in 1:(time_flux*ndt_days*ndays_yr)` is always an
  empty range. README's own "Known Issues" section flags this exact module as
  suspect ("not sure that it is fully working"), yet nothing exercises it.
- **Most `experiment` symbols are unreachable in tests**: `create_experiment_config`
  only covers 9 of ~25+ symbols `forcing()`/`init_model!`/`greb_model!`
  special-case (e.g. `:a1b_scenario`, `:co2_10x`, `:solar_cycle_11yr`,
  `:rcp26/45/60`) — only reachable via `PhysicsConfig(experiment=...)`
  directly. `:regional_co2_nh`, `:paleo_231kyr` are now exercised (§0.5, §0.10
  regression tests); the rest still aren't.
- **Branch-coverage tests assert shape, not distinctness**: the `log_eva`/
  `log_rain` sweep still only checks `length(result.ctrl) == 12` for its
  main pass, but a separate direct test now added alongside §0.7's fix
  confirms `hydro!`'s four `log_eva` values produce genuinely different
  output, not just the same shape.
- **Synthetic JLD2 fixtures**: §0.10's regression test now writes one (a
  `solar_paleo.jld2` in a tempdir) to exercise `load_solar_forcing_jld2`
  without the real, gitignored dataset. `load_greb_jld2!`/
  `load_flux_corrections_jld2!`'s "file exists" branches remain untested this way.

### 4.10 CI has no lint/format or doc-build check
`.github/workflows/ci.yml` only runs `julia-runtest`. No `JuliaFormatter`
check and no docs build step, so style drift or documentation regressions
(like §3's misplaced-docstring bug) aren't caught by CI. Adding one means
picking a formatting style first — a decision, not a mechanical addition.

---

## 5. Suggested order (low-risk → high-value)

1. ~~Fix `§0` bugs~~ ✅ done — all ten, across three passes, with regression tests.
2. ~~Split into files (§1.2)~~ ✅ done — mechanical, improves everything after.
3. ~~Add `[compat]`, CI, dead-dependency cleanup, `log_eva`/`log_rain` branch
   coverage, benchmark suite (§3)~~ ✅ done.
4. ~~Introduce state structs (§1.1)~~ ✅ done — validated against real
   control-run output (`greb_dataset_jld2/`): agrees with the pre-refactor
   baseline to ~7e-12 absolute (floating-point reassociation, not a behavior
   change), and incidentally fixed §2.2 (lazy allocation) and unlocked §2.4
   (threading) as a side effect. §1.3/§1.4 (argument-list/`RunSpec` collapsing)
   remain deferred — real but smaller wins than §1.1 was.
5. Threading the spatial operators (§2.4) — **attempted, benchmarked, reverted**.
   `Threads.@threads` on `diffusion!`/`advection!` was 3-6x *slower* at
   `-t auto` (grid too small, per-row work too cheap for the scheduling
   overhead) and even added allocations at `-t 1`. Validated the
   implementation was race-free (bit-identical across thread counts) before
   the benchmarks made clear it wasn't worth keeping. Picked up several other
   low-risk §4 items in the same pass instead: `forcing()`'s `:full_model`
   fast path (§4.4), `hydro!`'s constant hoist (§4.5), missing docstrings
   (§4.7), and a golden/snapshot regression test (§3).
6. **Next up**: `Float32`/allocation audit (§2.3, §2.5's `diffusion!`/
   `advection!`/`circulation!` allocations — worth re-measuring now that §1.1
   changed codegen), §1.3/§1.4's remaining argument-collapsing. If threading
   is revisited, §2.4 now has a documented reason to try `Polyester.@batch`
   instead of `Threads.@threads`, not just retry the same approach.
7. The §4 findings that need domain/design input before they can even be
   scoped as mechanical work (§4.1's remaining unwired switches, §4.2
   `log_clim`/dataset split, §4.8 `forcing()` dispatch) — bring in whoever
   knows the intended physics/API shape.

> ⚠️ Every performance change must be validated against a reference run —
> this is a numerical model, so "faster" only counts if the output is
> unchanged within tolerance. Use `benchmark/run_benchmarks.jl` +
> `claude/BENCHMARKS.md` to track this going forward.

---

## Changelog of this document

- **2026-08-06 (threading attempt + low-risk cleanup)**: The user asked to
  plan and implement threading (§2.4, unblocked by §1.1) plus some of this
  doc's other low-risk items, keeping this doc current as work landed.
  Added a golden/snapshot regression test against real data (§3) as a
  tripwire, then threaded `diffusion!`/`advection!` over latitude rows
  (`Threads.@threads :static`); found and fixed a real race in the process
  (`ws.T1h`'s shared polar scratch buffer needed a per-thread copy) plus a
  sizing gotcha (`Threads.threadid()` can exceed `Threads.nthreads()` — must
  size by `maxthreadid()`), and validated bit-identical output across thread
  counts. The benchmark gate the plan called for then showed the threading
  itself was a net loss — 3-6x *slower* at `-t auto` (14 threads), plus new
  allocations even at `-t 1` — so it was reverted; §2.4 now documents the
  actual numbers and why a plain `Threads.@threads` retry isn't worth it at
  this grid size. Also implemented §4.4's `:full_model` fast path, §4.5's
  `hydro!` constant hoist, and §4.7's missing docstrings (~22 exports).
  Separately, investigated and cut the test suite's own runtime: the
  `log_eva`/`log_rain` branch-coverage sweep ran a full model-year for all 20
  cross-product combinations for no extra coverage (the two switches gate
  independent branches) — reduced to 5 zipped runs covering every value of
  both axes. Updated §2.4/§4.4/§4.5/§4.7/§3/§5 to reflect all of the above.
- **2026-08-06 (third pass + state-struct refactor)**: The user asked for a
  direct line-by-line comparison against the original Fortran reference
  (`greb.model.mscm.f90`) to confirm no conceptual drift, then to plan and
  execute both the fixes it found and the previously-deferred §1.1 refactor.
  Found and fixed four more real bugs against the Fortran text (§0.7–0.10):
  `hydro!`'s missing `log_eva` modes 1/2, `init_model!`'s mismapped `_drsp`
  climatology overrides, `LWradiation!`'s `log_atmos_dmc` asymmetry, and
  paleo/orbital solar forcing never being wired up — all with regression
  tests, all committed separately before starting the refactor so it had a
  known-good baseline to validate against. Then replaced the ~40 mutable
  module globals with `ClimateFields`/`ModelState` (§1.1), threaded explicitly
  through every physics/circulation/tendencies/output/model function,
  migrated one file at a time with the full test suite run after each.
  Validated against a real `greb_dataset_jld2/` control+scenario run: matches
  the pre-refactor baseline to ~7e-12 absolute (floating-point reassociation,
  not a behavior change) and — unplanned — cut runtime ~7x (49s→7.0s) and
  allocations ~170x (17.1 GiB→100 MiB) on that same run, plus fixed §2.2
  (lazy allocation) and unlocked §2.4 (threading) as direct side effects.
  Updated §4.1/§4.3/§4.9 to reflect what §0.8/§0.7 fixed, and §1.3/§1.4/§2.4
  to reflect what §1.1 changed but didn't fully subsume.
- **2026-08-06 (second pass)**: A fresh, independent read through the (now
  split) `src/` tree found three more real bugs of the same "disconnected
  config flag" shape as §0.2/§0.3 — fixed with regression tests: `circulation!`'s
  inverted `do_conv` gate (§0.4), a `co2_part` module-global leak between runs
  (§0.5), and a silent fallback on invalid `log_rain` (§0.6). Also fixed as
  mechanical cleanup: two more dead globals (§1.6 addendum), the unused
  `Statistics` dependency (§3), and three exported functions whose docstrings
  were silently discarded due to placement after `function` instead of before
  it (§3). Everything from the first pass that had been sitting uncommitted
  (the split, the first three bug fixes, dependency/compat/precompile
  changes, CI) was committed. Added §4, a new set of documented-but-not-fixed
  findings — 9 unwired config switches, the `log_clim`/dataset disconnect,
  `hydro!`'s incomplete evaporation modes, `forcing()` per-timestep overhead,
  more `const`-ify candidates, missing docstrings on ~20 more exports,
  `forcing()`'s dispatch style, testing gaps, and CI lint/doc-build gaps —
  each deferred for a stated reason (domain expertise, API design, or scope)
  rather than fixed blind. Re-ran the benchmark suite for a genuine
  post-fixes comparison entry in `claude/BENCHMARKS.md`.
- **2026-08-06 (first pass)**: Safe/mechanical optimization pass. Fixed and
  added regression tests for
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
