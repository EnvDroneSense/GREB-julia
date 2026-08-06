# GREB.jl — Potential Improvements (Structure & Performance)

Observations and concrete suggestions for the `GREB` module. The model has
been split from a single 2245-line `src/GREB.jl` into topical files under
`src/` — see §1.2. Two bug-hunting passes have gone through this list so far
(2026-08-05/06 and 2026-08-06 §4); nothing else here has been applied except
where marked **✅ done** — this remains a menu for future work, primarily the
state-struct refactor in §1.1.

The single theme underneath most of what remains: **the model runs on ~40
mutable module-level globals** (grid arrays, climatology fields, work buffers,
flux corrections). That one design choice drives most of the remaining
structural *and* performance issues below — it's why the module eats ~750 MB
at load and why the code is hard to thread or test.

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

**Second-pass addendum**: two more, found the same way — `co2_part_scn` and
`sw_solar_forcing_data` (`src/state.jl`), both declared and never read or
written anywhere else. Removed.

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

### 4.1 Nine config switches are wired to nothing
`PhysicsConfig` fields that are set (including by `create_experiment_config`)
but never read by any physics function, so toggling them has zero effect on
model output: `solar_multiplier` (the actual `+27 W/m²` solar-forcing effect
is computed independently inside `forcing()`, `src/tendencies.jl`, not via
this field), `log_clouds_drsp`/`log_vapor_drsp`/`log_ice_drsp`/`log_humid_drsp`
("CO₂ Response Switches"), `log_tsurf_ext`/`log_hwind_ext`/`log_omega_ext`
("External Forcing" switches), and `log_ice_dmc` (distinct from `log_ice`,
which *is* wired in `radiation.jl`/`ocean.jl`). All are documented,
notebook-exposed controls — this is the same "config field disconnected from
behavior" bug class as §0.2–0.6, but wiring 9 features up (or removing them)
is a scope/design decision, not a mechanical fix.

### 4.2 `log_clim` doesn't select the ERA/NCEP dataset its label implies
The notebook documents a Hydrology control "Climatology (ERA/NCEP)" bound to
`log_clim`. In the module, `cfg.log_clim` only ever affects one thing
(`src/config.jl`: swaps four regression coefficients when
`log_rain==0 && log_clim==1`) — it does **not** select which dataset gets
loaded. The actual ERA-vs-NCEP choice is the separate, disconnected `dataset`
keyword to `load_greb_jld2!` (`src/io.jl`), set before `greb_model!` runs and
never connected to `cfg` at all. Fixing this properly means deciding whether
`cfg` should drive data loading (a real API question) — deferred alongside §1.1.

### 4.3 `hydro!`'s 4 documented evaporation modes are really 2
`log_eva` branches for `1`/`2`/anything-else are byte-identical to the
`-1` branch (`src/physics/hydrology.jl`) — only 2 of the 4 documented
evaporation formulas actually exist. Confirmed the same gap exists in the
original notebook, so this is a pre-existing physics-completeness question,
not a coding bug — needs someone who knows the intended formulas for modes
`1`/`2`, not a mechanical fix.

### 4.4 `forcing()` per-timestep overhead
`forcing()` (`src/tendencies.jl`) has no `:full_model` branch (the default,
most common experiment), so every timestep in the scenario loop falls through
the entire ~25-branch chain, including `startswith(string(cfg.experiment),
"regional_co2_")`, which allocates a new `String` on every call regardless of
match. Additionally, the regional-CO2 mask (`co2_part`) is fully recomputed
every timestep even though it depends only on static topography and never
changes within a run — easy to hoist to a one-time setup once §1.1 gives
`init_model!` a natural place to precompute it.

### 4.5 `hydro!` recomputes per-call "constants" as locals
`const_factor1`, `const_factor2/3`, `gust_land/ocean`, `cE_land/ocean`,
`const_latent` (`src/physics/hydrology.jl`) are plain local variables
recomputed every timestep, inconsistent with `src/constants.jl`'s pattern of
precomputing everything else once at module load. Cheap individually, adds up
across tens of thousands of timesteps; low priority.

### 4.6 More `const`-ify candidates
Beyond the constants already made `const` (§2.1), several module globals are
verifiably never *rebound* (only mutated via `.+=`/`fill!`/broadcasting) —
the 11 `diagnostics!` accumulators (`Tsmn, Tamn, Tomn, qmn, amn, swmn, lwmn,
qlatmn, qsensmn, ftmn, fqmn`, `src/state.jl`) plus climatology fields like
`z_topo`/`Tclim`/`cldclim`. Marking these `const` is individually low-risk
(Julia errors immediately if a hidden rebind exists) but spans many locations
— deserves its own dedicated, careful pass rather than a piecemeal one.

### 4.7 Missing docstrings on the public API
Beyond the 3 fixed in §3, ~20 more exported functions/types have no
docstring at all — including `greb_model!` itself, the package's main entry
point (its `(time_flux, time_ctrl, time_scnr, cfg)` argument order is only
explained by one README usage line). Also: `SWradiation!`, `LWradiation!`,
`hydro!`, `seaice!`, `deep_ocean!`, `diffusion!`, `advection!`, `circulation!`,
`tendencies!`, `forcing`, `diagnostics!`, `output!`, `time_loop!`, `init_model!`,
`qflux_correction!`, `build_monthly_climatology`, `apply_scenario_anomalies`,
`compute_annual_ice_climatology`, `PhysicsConfig`, `TimeState`, `MonthlyRecord`.

### 4.8 `forcing()`'s dispatch style is inconsistent with the rest of the codebase
Its ~183-line/~25-branch `if/elseif` chain over `cfg.experiment`
(`src/tendencies.jl`) is the odd one out next to the Dict-based dispatch used
elsewhere (`HYDRO_PARAMS` in `constants.jl`, `file_map` in `io.jl`). Natural
to address alongside the §1.1 state-struct refactor, since `forcing`'s
signature changes then anyway — not a standalone mechanical fix, and (per
§0.6's note) the chain's "missing branch = valid no-op" semantics need to be
preserved carefully by whatever replaces it.

### 4.9 Testing gaps
- **Zero scenario-run coverage**: every `greb_model!` call in
  `test/runtests.jl` (aside from the new §0.5 regression test) uses
  `time_scnr = 0`, so the scenario loop's body — all of `forcing()`,
  `build_monthly_climatology`, `apply_scenario_anomalies`, and the
  `is_forced_boundary`/`is_sst_plus1`/`is_orbital_exp` branches in
  `greb_model!` — is otherwise untested.
- **`qflux_correction!`'s loop body has zero coverage**: every test uses
  `time_flux = 0`, so `for it in 1:(time_flux*ndt_days*ndays_yr)` is always an
  empty range. README's own "Known Issues" section flags this exact module as
  suspect ("not sure that it is fully working"), yet nothing exercises it.
- **Most `experiment` symbols are unreachable in tests**: `create_experiment_config`
  only covers 9 of ~25+ symbols `forcing()`/`init_model!`/`greb_model!`
  special-case (e.g. `:a1b_scenario`, `:co2_10x`, `:solar_cycle_11yr`,
  `:obliquity`, `:rcp26/45/60`) — only reachable via `PhysicsConfig(experiment=...)`
  directly, and none are exercised in the pre-existing test suite (this is
  exactly how §0.5's `co2_part` leak went undetected; the new regression test
  for it is the first test to reach a `regional_co2_*` branch at all).
- **Branch-coverage tests assert shape, not distinctness**: the `log_eva`/
  `log_rain` sweep only checks `length(result.ctrl) == 12` — it wouldn't have
  caught §4.3 (branches producing identical output) because it never compares
  outputs across branch values.
- **No synthetic JLD2 fixture**: `read_jld2`'s happy path and
  `load_greb_jld2!`/`load_flux_corrections_jld2!`'s "file exists" branches are
  untested; a minimal file written via `JLD2` to a `tempname()` path would
  exercise these without needing the real, gitignored dataset.

### 4.10 CI has no lint/format or doc-build check
`.github/workflows/ci.yml` only runs `julia-runtest`. No `JuliaFormatter`
check and no docs build step, so style drift or documentation regressions
(like §3's misplaced-docstring bug) aren't caught by CI. Adding one means
picking a formatting style first — a decision, not a mechanical addition.

---

## 5. Suggested order (low-risk → high-value)

1. ~~Fix `§0` bugs~~ ✅ done — all six, across two passes, with regression tests.
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
6. The §4 findings that need domain/design input before they can even be
   scoped as mechanical work (§4.1 unwired switches, §4.2 `log_clim`/dataset
   split, §4.3 `hydro!`'s missing evaporation formulas, §4.8 `forcing()`
   dispatch) — bring in whoever knows the intended physics/API shape.

> ⚠️ Every performance change must be validated against a reference run —
> this is a numerical model, so "faster" only counts if the output is
> unchanged within tolerance. Use `benchmark/run_benchmarks.jl` +
> `claude/BENCHMARKS.md` to track this going forward.

---

## Changelog of this document

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
