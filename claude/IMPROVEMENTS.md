# GREB.jl — Potential Improvements (Structure & Performance)

Observations and fixes for the `GREB` module, accumulated across seven audit
passes (2026-08-05 through 2026-08-06). The model was split from a single
2245-line `src/GREB.jl` into topical files under `src/` (§1.2) and its ~40
mutable module-level globals were replaced with explicit `ClimateFields`/
`ModelState` structs (§1.1) — that one refactor is why the module used to eat
~750 MB at load, why it was hard to thread or test, and why the same
"config flag silently disconnected from behavior" bug shape kept recurring.
Everything below is either **✅ done**, explicitly deferred with a reason, or
a documented finding that was investigated and intentionally not changed.

---

## 0. Known bugs — fixed ✅

18 real, previously-shipped bugs found and fixed across five passes (a 19th
finding was investigated and reverted — see 0.14), each with a regression
test in `test/runtests.jl`. Passes 1–3 were found by grepping every
config-field write against its read sites; pass 4 (fourth pass, below) was a
direct line-by-line comparison against the Fortran
reference `greb.model.mscm.f90`, and every fix there cites the exact Fortran
line(s) confirmed by direct read (not recalled from memory).

### First pass (2026-08-05/06)
- **0.1 `ws.cE` typo** (`hydro!`, `log_eva==0`): referenced a nonexistent
  `CirculationWorkspace` field, threw at runtime; invisible to tests because
  the default `log_eva=-1` never reached that branch. Fixed: renamed to
  `ws.cE_buf`.
- **0.2 `HYDRO_PARAMS` was a `Tuple` of `Pair`s**, not a `Dict`:
  `set_hydrology_parameters!`'s `get(HYDRO_PARAMS, cfg.log_rain, default)` did
  positional tuple indexing — `log_rain=0` (the default) silently fell back to
  the wrong preset, `log_rain∈{1,2,3}` threw `BoundsError`. Fixed: wrapped the
  literal in `Dict(...)`.
- **0.3 `set_hydrology_parameters!` never wrote to `cfg`**: assigned to
  disconnected module globals `c_q/c_rq/c_omega/c_omegastd` instead of
  `cfg.c_q` etc. — **`log_rain` had zero effect on any run**, always used
  struct defaults. Fixed: assign `cfg.c_q = ...` directly; removed the
  disconnected globals.

### Second pass (2026-08-06)
- **0.4 `circulation!`'s `do_conv` gate was inverted**: read
  `cfg.log_conv == 0` while all four sibling flags read `== 1` — moisture-flux
  convergence silently off for every default run. Fixed: all five flags now
  direct `Bool` truthiness. *(Correction, found in the fourth pass: Fortran's
  own `log_conv` default is `0`, and `log_conv==0` is what triggers
  convergence — it's the one switch in the family with an inverted
  convention, not an instance of a shared one as originally stated here. The
  Julia fix is still correct; only the reasoning was wrong.)*
- **0.5 `co2_part` module-global leaked between runs**: written only by
  `forcing()`'s regional-CO₂ branches, never reset — a regional experiment run
  once could contaminate any later, unrelated run in the same session. Fixed:
  `init_model!` resets `co2_part .= 1.0` per run (moot by construction once
  §1.1 made `ClimateFields` per-call rather than global).
- **0.6 Silent fallback on invalid `log_rain`**: `get(HYDRO_PARAMS, ..., default)`
  swallowed typos silently. Fixed: `error()`s on an unrecognized key instead.

### Third pass (2026-08-06): direct audit against the Fortran reference
- **0.7 `hydro!`'s `log_eva` modes 1/2 didn't exist**: silently duplicated
  mode `-1`'s formula. Fixed: implemented Fortran's real coefficients (mode 1:
  `0.04`/`0.73` land/ocean, gust terms `144.0`/`50.41`, wind from
  `uclim`/`vclim`; mode 2: `0.56`/`0.79`, gust `81.0`/`16.0`, wind from
  `wsclim`) plus an `error()` fallback for other values.
- **0.8 `init_model!`'s `_drsp` climatology overrides mismapped**: the
  `cldclim=0.7` override was missing entirely; `qclim`/`mldclim` overrides
  were gated on the wrong (`_dmc`) switches instead of `_drsp`. Fixed: added
  the missing override, corrected the two switch names.
- **0.9 `LWradiation!`'s `log_atmos_dmc` asymmetry**: zeroed both `LW_up` and
  `LW_down` when off; Fortran only zeros `LW_down` (`LW_up` is snapshotted
  first). Fixed: removed the extra `LW_up .= 0.0`.
- **0.10 Paleo/orbital solar-forcing table swap never wired up**: Fortran
  swaps in an alternate `(ydim, nstep_yr)` solar table for paleo/orbital
  experiments (`sw_solar = sw_solar_scnr`); Julia had the loader but never
  called it. Fixed: `greb_model!` swaps the table for the four affected
  experiments at scenario start, restores it after.

### Fourth pass (2026-08-06): second Fortran audit + threading attempt
- **0.11 `circulation!` never zeroed `ws.dX_conv` itself**: Fortran zeroes
  `dx_diffuse`/`dx_advec`/`dx_conv` once per call before its sub-step loop
  (`greb.model.mscm.f90:856-858`); Julia relied on each kernel to zero its own
  output, but `convergence!` only writes when called — structurally never for
  `h_scl==z_air`. **Every temperature (`Ta`) sub-step in every default run**
  picked up a stale moisture-convergence term left over from the *previous
  timestep's* humidity call. Fixed: zero all three buffers once at the top of
  `circulation!`.
- **0.12 `deep_ocean!` disabled turbulent mixing under sea ice**: gated
  turbulent mixing on the same `Ts>=To_ice2` mask as entrainment/detrainment;
  Fortran gates turbulent mixing on ocean alone, using `Tx=max(To_ice2,Ts)`
  specifically so it stays well-defined under ice
  (`greb.model.mscm.f90:818-830`). Cut off ocean-atmosphere heat exchange
  under sea ice and high-latitude winters. Fixed: split the mask so turbulent
  mixing uses `z_topo<0` alone.
- **0.13 `hydro!`'s rain limit used a single grid point**: `wz_vapor[1, 1]`
  (the Antarctic coast) applied globally instead of the full spatially-varying
  field (`greb.model.mscm.f90:755`). Fixed: per-point `wz_vapor[i, j]`.
- **0.14 `min_T_K` investigated, kept as-is**: Fortran's `Tmin_limit=40` is a
  raw-Kelvin floor (`greb.model.mscm.f90:470-477`), which looked at first like
  a units mismatch against Julia's `233.15` K (-40°C). Tried matching it
  (`min_T_K=40.0`), then reverted: 40 K is colder than anywhere on Earth ever
  gets, so it's a floor that can *never* physically bind — not a real
  numerical safety net worth copying literally. Kept `233.15` K, a real
  cold-extreme floor, on the judgment that Fortran isn't always the right
  reference to match exactly.
- **0.15 `seaice!` skipped the glacier override when `log_ice` was off**: an
  early `return` inside the `!log_ice` branch skipped the
  glacier→`cap_land` override Fortran applies unconditionally afterward
  (`greb.model.mscm.f90:786-792`). Fixed: removed the early return.
- **0.16 `hydro!` computed `rq` after the `log_eva` dispatch**: Fortran
  computes `rq = q/qs` (feeding `dq_rain`) *before* dispatching on `log_eva`
  (`greb.model.mscm.f90:715-719`); Julia computed it afterward, by which point
  `log_eva==0`'s branch had already overwritten `ws.qs` with a
  Tskin-based value, corrupting `rq` for that mode. Fixed: moved the `rq`
  computation before the `log_eva` dispatch.
- **0.17 `diffusion!`/`advection!`'s polar sub-stepping used Gauss-Seidel,
  not Jacobi**: wrote `ws.T1h[j] += dq` inside the same `@turbo` loop that
  computed `dq`, so later `j` read values already updated by earlier `j` in
  the same sweep. Fortran computes the whole row's increment from the OLD
  row first, then applies it all at once (`greb.model.mscm.f90:983-1036`,
  `:1227-1228`) — Julia's version was also unsound under `@turbo`, which
  assumes no cross-iteration dependency. Fixed: two-pass Jacobi (compute into
  a re-added `ws.dTxh` scratch buffer, then apply) for both kernels' polar
  branches.
- **0.18 Humidity update used `clamp`, not Fortran's threshold-replace**:
  `clamp(dq, -0.9q, 0.02)` pulls *any* `dq` below `-0.9q` up to it; Fortran
  only replaces `dq` with `-0.9q1` when `dq<=-q1`, leaving the range between
  `-q1` and `-0.9q1` untouched (`greb.model.mscm.f90:481-483`). Separately,
  `log_hydro_dmc==0` only zeroed the eva/rain terms; Fortran zeroes the
  *entire* `dq` including circulation and flux correction
  (`greb.model.mscm.f90:486`). Fixed both in `time_loop!`'s humidity update.

### Fifth pass (2026-08-06): allocation fix + a real doc-tooling bug
- **0.19 19 exported functions' docstrings were silently detached from their
  bindings**: the codebase's convention of a `# ── notebook cell ... ──`
  comment line placed *between* a docstring and the `function`/`struct` it
  documents breaks Julia's `@doc` binding — a bare comment in that position
  is enough to disconnect the string literal from the definition, even
  though it reads as attached in the source. Confirmed with a minimal
  repro and with `Docs.hasdoc` before/after across every exported symbol:
  `init_model!`, `qflux_correction!`, `diffusion!`/`advection!`/`circulation!`,
  `diagnostics!`/`output!`/`time_loop!`, `build_monthly_climatology`/
  `apply_scenario_anomalies`/`compute_annual_ice_climatology`,
  `SWradiation!`/`LWradiation!`, `seaice!`/`deep_ocean!`, `tendencies!`/
  `forcing`, `load_solar_forcing_jld2`, `hydro!` — 19 of the module's 36
  exports had a docstring in the source that `@doc`/Documenter's `@autodocs`
  would never actually see. Found while wiring up §3's Documenter.jl site
  (an `@autodocs` page would have silently rendered blank for most of the
  API). Fixed by moving each notebook-cell comment to *before* its
  docstring instead of after (a mechanical, mass find-and-fix across the 9
  affected files, verified by re-running `Docs.hasdoc` against all 36
  exports — 0 missing after). Also added docstrings to `xdim`/`ydim`/
  `nstep_yr`, the 3 exports that genuinely never had one.

**Investigated, not changed**: `log_eva==1`'s land gust literal (`144.0`) —
Fortran's text has `144.**2` (=20736), but every sibling mode's gust term is
2–12 m/s while 144 m/s is physically absurd for a turbulent gust; 12 m/s
(`sqrt(144)`) fits the pattern, suggesting Fortran's own `**2` is the
anomaly, not Julia's value. Kept Julia's `144.0` unsquared per user
direction (use physical judgment when the reference itself looks wrong);
confirmed identical in a second, independently-maintained copy of the
Fortran source, so this isn't a transcription error. Also verified
`:earth_sun_distance`'s `rS0` rescale is correctly scenario-only in both
Julia and Fortran (`greb.model.mscm.f90:406-416` applies it in the same
scenario-start block as the paleo/obliquity/eccentricity swaps) — no fix
needed, a prior finding suggesting an init-time fix was mistaken. `min_T_K`
(§0.14) is a third case of the same pattern: tried matching Fortran, then
reverted after judging the Fortran value itself wasn't worth copying.

---

## 1. Structure

### 1.1 Replace module-global mutable state with state structs ✅ done
The ~40 mutable module globals (climatology, grid fields, flux corrections,
`co2_part`, `sw_solar`) are gone, replaced by `ClimateFields`/`ModelState`
(`src/state.jl`), threaded explicitly through every physics/circulation/
tendencies/output/model function. `greb_model!` defaults to a fresh
`ClimateFields()` per call, making the old cross-run leak hazards (§0.5)
structurally impossible for that path; reused-`fields` callers still get a
defensive reset.

Validated against real data: a 1-year control+scenario run's output before
and after the refactor agreed to ~7e-12 absolute (floating-point
reassociation from different codegen, not a behavior change). Unplanned
bonus: the same run went from 49s/17.1 GiB to 7.0s/100 MiB (global mutable
bindings are type-unstable at every access; struct fields aren't) — this also
fixed §2.2 and unblocked §2.4.

### 1.2 Split the module into files ✅ done
`src/GREB.jl` (2245 lines) is now a module shell with `include()`s in
dependency order: `constants.jl` → `config.jl` → `state.jl` → `io.jl` →
`physics/{radiation,hydrology,ocean}.jl` → `circulation.jl` →
`tendencies.jl` → `output.jl` → `postprocess.jl` → `model.jl`. Pure move, no
logic changed.

### 1.3 Collapse giant positional argument lists into structs ✅ done
`output!`/`diagnostics!` took 9–13 positional `(xdim,ydim)` matrices each.
Checked which of those args were genuinely new state to wrap vs. things
already sitting on an existing object: `Ts0/Ta0/To0/q0` are one cohesive,
persistently-allocated "current surface state" (new `SurfaceState` struct,
`src/state.jl`), but `albedo`/`sw`/`lw_surf`/`q_lat`/`q_sens`/`ice`/`precip`/
`evap`/`qcrcl` already live on the `tend` `NamedTuple` `tendencies!` returns
every call, or on `ws::CirculationWorkspace` (`ws.precip_out`/`evap_out`/
`qcrcl_out`) — so the fix passes `surf`/`tend`/`ws` through instead of
manually unpacking their fields at each call site, rather than inventing a
second struct that would just duplicate what already exists. Zero new
allocations: both `tend` and `ws` already exist at every call site,
`SurfaceState` just wraps existing arrays by reference. Bonus: fixed a real
footgun in the process — `time_loop!`'s old positional order was
`Ts,Ta,q,To` while `diagnostics!`/`output!`'s was `Ts0,Ta0,To0,q0` (`q`/`To`
swapped) — named struct fields make that ordering irrelevant.

### 1.4 Turn run durations into a `RunSpec` ✅ done
`greb_model!(time_flux, time_ctrl, time_scnr, cfg; ...)`'s three
easily-swapped positional ints are now `greb_model!(run::RunSpec, cfg; ...)`,
`RunSpec(; flux=0, ctrl=1, scnr=1)` (`src/config.jl`). All call sites
(`examples/run_greb.jl`, `test/runtests.jl`, `src/GREB.jl`'s precompile
workload, `README.md`) updated.

### 1.5 ~~Reunite the two notebook-only helpers with the package~~ ✅ done
`examples/run_greb.jl` already provides a non-Pluto path.

### 1.6 Remove hidden/dead globals/fields ✅ done
Removed across three passes, all confirmed dead by grep: `WZ_CACHE`,
`circ_workspace`, the shadowed "Monthly-mean buffers" and "Control run
monthly means" blocks (21 arrays), 4 dead accumulator entries, `co2_part_scn`,
`sw_solar_forcing_data`. Fourth-pass addendum: `CirculationWorkspace`'s
`dX_crcl`, `ws.eva`, `ws.rain` fields (zero reads/writes after §0.18's
humidity-update rewrite made the `ws.eva`/`ws.rain` zero-stand-ins
unnecessary); `ws.temp_buf` — previously also dead — is now genuinely used as
that rewrite's scratch buffer, so it stayed.

---

## 2. Performance

### 2.1 ~~`const`-ify the constants~~ ✅ done
Grid dimensions and physical constants in `src/constants.jl` are all `const`.

### 2.2 Stop allocating ~750 MB at module load ✅ done
Direct consequence of §1.1: climatology arrays are now `ClimateFields` fields,
allocated only when actually constructed, not at `using GREB` time.

### 2.3 Consider `Float32` for climatology fields — analyzed, not implemented
JLD2 stores `Float32`; up-converting to `Float64` on load gains no precision
that wasn't already lost on disk, so retyping the static climatology loses
nothing in isolation. **But** `init_model!` does `Ts_ini = copy(Tclim[:,:,end])`
etc. — `copy` preserves `eltype`, so naively retyping `ClimateFields` to
`Float32` would silently downgrade the entire simulation's *working* state
(`Ts`/`Ta`/`To`/`q`, accumulated over the whole run), not just the read-only
climatology. Safe implementation would need an explicit `Float64.(...)`
guardrail at those four `init_model!` lines. Analyzed this pass but not
implemented — feature canceled by user request before landing.

### 2.4 Thread the spatial operators — tried, measured, reverted
Implemented `Threads.@threads :static` on `diffusion!`/`advection!`'s outer
latitude-row loop; fixed a real race in the process (`ws.T1h`'s shared polar
scratch buffer needed a per-thread copy, sized by `Threads.maxthreadid()` not
`Threads.nthreads()` — Julia 1.9+'s threadpool split means `threadid()` can
exceed `nthreads()`). Validated bit-identical across thread counts.

**Benchmarked, and reverted**: the grid (96×48, ~1µs of work per row) is too
small for `Threads.@threads`'s scheduling overhead to pay off at any thread
count tried — `-t auto`/14 threads was 3–6x *slower* end-to-end
(`tendencies!`: 4.2ms→23.7ms), and even `-t 1` picked up new allocations from
the threading machinery itself. Full numbers in `claude/BENCHMARKS.md`'s
"-t 1"/"-t auto" entries. Reverted to plain `for k in 1:ydim` loops; kept the
unrelated dead-field removal (`dTxh`) found while touching the struct — later
un-reverted for a different reason, see §0.17. **Takeaway**: this grid/kernel
shape is wrong for `Threads.@threads`; a lower-overhead approach
(`Polyester.@batch`) or a much larger grid might fare better — don't retry
plain `Threads.@threads` here without a reason to expect a different result.

### 2.5 ~~Cut per-call allocations in `forcing`~~ ✅ done
`icmn_ctrl` is a required argument computed once via
`compute_annual_ice_climatology`, not fresh zeros every call.

### 2.6 Reduce time-to-first-run (TTFX) ✅ done
`PrecompileTools.@compile_workload` compiles the heavy kernels at package
build time (~2 min once, e.g. in CI) instead of on every user's first call.

### 2.7 Cut `SWradiation!`'s per-call allocation ✅ done
`for j in 1:ydim; sf = ...; @. sw[:, j] = sf * (1.0 - albedo[:, j]); end` —
`albedo[:, j]` on the RHS of `@.` is plain `getindex`, not a `dotview` like
the LHS, so it materialized a fresh `Vector{Float64}` every iteration: 40704
bytes / 48 allocs per call, confirmed with `@allocated`, and the *entire*
allocation footprint of `tendencies!` (which calls it once per timestep) —
every other kernel already benchmarked at 0 bytes. Fixed by wrapping the
loop body in `@views`. Found via the existing per-kernel benchmark harness
(`benchmark/run_benchmarks.jl`) rather than a fresh `@profview` session,
since it already isolates every hot kernel individually and pointed
straight at the one non-zero entry.

### 2.8 Hot-loop field/config localization + type stability — audited, clean
Checked whether `fields.xyz`/`cfg.xyz` get re-dereferenced inside
`@turbo`/`@inbounds` loops (they don't — already hoisted to locals before
every hot loop across `circulation.jl`/`hydrology.jl`/`ocean.jl`) and
whether any struct field is untyped/abstract (none are — `PhysicsConfig`,
`ClimateFields`, `CirculationWorkspace`, `ModelState`, `RunSpec` are all
concretely typed, so `@code_warntype` has nothing to flag). No changes
needed; recorded so this doesn't get re-audited from scratch next pass.

### 2.9 Precomputation/caching audit — clean, nothing left to hoist
Re-checked every candidate for a static/precomputable value: the
`regional_co2_*` masks are 4-of-6 hoisted to `init_model!` already (§4.4);
the other 2 (`_ocean`/`_land_ice`) genuinely can't be, since they depend on
`icmn_ctrl` (the control run's own output), which doesn't exist at init
time. Solar forcing is already a precomputed lookup table
(`fields.sw_solar[lat, day-of-year]`) — no per-timestep trig anywhere.
Latitude/area weights (`lat_grid`, `dxlat_grid`, `ccx_diff`, `ccx_adv`,
`IS_POLAR`) are already `const` arrays in `constants.jl`, referenced
directly by the hot loops that use them, not recomputed. Nothing to do.

**Flagged, not fixed** (incidental finding while auditing latitude
weights): `output.jl:40`'s `global_mean = sum(state.Tsmn) / (xdim*ydim)` is
an *unweighted* grid mean, not `cos(lat)`-area-weighted — every grid cell
counts equally regardless of its actual area, which shrinks toward the
poles. This is a scientific-accuracy question (would change diagnostic
output numbers), not a performance one — needs the same
Fortran-comparison-plus-judgment treatment as `min_T_K`/the gust literal,
not a silent change made in passing during a perf-focused pass.

---

## 3. Testing, tooling & reproducibility

- **Golden/snapshot regression test** ✅ done: a real 1yr control + 1yr
  scenario `:full_model` run against `greb_dataset_jld2/`, checked against
  reference monthly global-mean `Ts`/`Ta`/`q` within `atol=1e-6`. Skips via
  `@test_skip` when the dataset isn't present. Reference values recomputed
  after this pass's §0.11–0.18 bug fixes (real behavior changes).
- **Branch-coverage tests**: `log_eva`/`log_rain` sweep (added because §0.1
  hid in an untested branch) now uses a zipped 5-run sweep instead of a
  4×5=20-run cross product — the two switches gate independent branches, so
  the cross product bought no extra coverage; cut steady-state suite runtime
  roughly in half. A separate direct test confirms `hydro!`'s `log_eva` values
  produce genuinely different output, not just the same shape (coverage
  gap that let §0.7 hide).
- **CI** ✅ done, **`[compat]` bounds** ✅ done, **dead dependencies**
  (`NCDatasets`/`StaticArrays`/`Statistics`) removed ✅ done, **misplaced
  docstrings** (3 functions, string literal after `function` instead of
  before, plus 19 more with the *comment*-after-docstring variant of the
  same shape — §0.19) fixed ✅ done, **benchmark suite**
  (`benchmark/run_benchmarks.jl` + `claude/BENCHMARKS.md`) added ✅ done,
  **`Documenter.jl` site** ✅ done (`docs/`, deployed to GitHub Pages on
  push to `main` via `.github/workflows/docs.yml`).
- **Test suite cost audit** ✅ done: a live timed `Pkg.test()` run (164
  passed / 2 failed, 2m36.9s) found and fixed a real, currently-failing
  assertion — `@test GREB.min_T_K == 233.15` fails because
  `273.15 - 40.0` rounds to `233.14999999999998` in Float64, not the
  literal; fixed by comparing with `≈` instead of `==`. Re-read the whole
  suite and confirmed no testset is a genuine duplicate of another (each
  traces to a distinct bug), but found 4 non-`DATA_DIR`-gated testsets that
  each run one or more full 12-month `greb_model!` simulations on synthetic
  data and together ate 119.1s of the 156.9s total (76%) for
  crash/shape-only signal; one had a safe cut (the `co2_part` mask-reset
  test's first run only needs `init_model!`'s mask, not a completed
  scenario year — trimmed `RunSpec()` to `RunSpec(scnr=0)`), the rest had
  no safe cut without losing real coverage. Added a `RUN_GOLDEN` env var
  (default on) to additionally skip the golden snapshot test locally even
  when `greb_dataset_jld2/` is present, alongside the existing `DATA_DIR`
  gate.
- **Testing gaps** ✅ closed:
  - `qflux_correction!`'s loop body had zero coverage. Checked the Fortran
    reference directly for the one thing that looked like a real bug — Julia
    gives `Ts`/`To`/`q` climatology-forced correction fields but never `Ta`.
    **Confirmed intentional, not a bug**: Fortran's own `qflux_correction`
    (`greb.model.mscm.f90:540-597`) does the identical thing — no
    `TaF_correct` array exists anywhere in the Fortran source either. Added a
    test asserting `Ts`/`To`/`q` land exactly on their climatology targets
    (the correction's algebraic fixed point) plus a code comment citing the
    Fortran lines, so nobody "fixes" the asymmetry by mistake later.
  - Added a sweep test covering the ~20 `:experiment` symbols `forcing()`/
    `init_model!` handle but `create_experiment_config` doesn't build presets
    for — reachable by constructing `PhysicsConfig(experiment=:foo)`
    directly. Caught a real gotcha while writing it: `forcing()` only runs
    during the *scenario* loop, never control, so the sweep needs
    `RunSpec()`'s default `scnr=1` — `RunSpec(scnr=0)` would have made the
    whole sweep silently vacuous.
  - Added direct unit tests for `build_monthly_climatology`/
    `apply_scenario_anomalies` (multi-year averaging, the non-12-multiple
    `records[1]` fallback branch, empty/mismatched inputs) — previously only
    exercised indirectly through 2 full-model-run tests.
  - Added synthetic-dataset tests (`mktempdir`+`jldopen`, same pattern as the
    paleo-swap test) for `load_greb_jld2!`/`load_flux_corrections_jld2!`'s
    "file present" branches, previously only covered when the real
    (gitignored) dataset happened to be present locally.

---

## 4. Findings needing domain/design input — documented, not fixed

### 4.1 Seven config switches wired to nothing ✅ resolved
Checked each of `solar_multiplier`/`log_vapor_drsp`/`log_ice_drsp`/
`log_tsurf_ext`/`log_hwind_ext`/`log_omega_ext`/`log_ice_dmc` against
`greb.model.mscm.f90` directly, then asked the user what to do with each
category found:
- `log_vapor_drsp`/`log_ice_drsp`/`log_ice_dmc`/`solar_multiplier` have zero
  Fortran basis (no equivalent variable exists there at all, or — for
  `solar_multiplier` — the real effect is already computed correctly
  elsewhere, in `forcing()`'s `:solar_plus27` branch). **Removed** from
  `PhysicsConfig`, including the one dead assignment to `solar_multiplier`
  in `create_experiment_config(:solar_plus27)`.
- `log_tsurf_ext`/`log_hwind_ext`/`log_omega_ext` ARE declared and
  namelist-read in Fortran (`greb.model.mscm.f90:152-154,242`) but never
  used there either — genuinely dead upstream, not a missing-wiring bug.
  **Kept**, with a comment explaining why, for structural fidelity with the
  Fortran switch family.

### 4.2 `log_clim` doesn't select the ERA/NCEP dataset its label implies ✅ resolved
Confirmed `cfg.log_clim` only swaps 4 hydrology regression coefficients
(`set_hydrology_parameters!`) while `load_greb_jld2!`'s separate `dataset`
kwarg picks the actual ERA-vs-NCEP files — the two are fully disconnected.
Per the user's decision, kept them orthogonal (a caller can legitimately
combine `log_clim=1` with `dataset=:era` for a sensitivity experiment) and
just documented the relationship with cross-referencing comments on both
(`src/config.jl`'s `log_clim` field, `src/io.jl`'s `load_greb_jld2!`
docstring) — no behavior change.

### 4.6 More `const`-ify candidates — resolved by §1.1
The original target (module-level globals never rebound) no longer exists:
§1.1 moved everything into `ClimateFields`/`ModelState` struct fields, and
Julia has no `const`-field mechanism for mutable structs. Nothing to do here.

### 4.8 `forcing()`'s dispatch style is inconsistent with the rest of the codebase — flagged future improvement, kept as-is
Its ~180-line `if/elseif` chain over `cfg.experiment` is the odd one out next
to the Dict-based dispatch used elsewhere (`HYDRO_PARAMS`, `io.jl`'s
`file_map`, the solar-swap table). Asked the user: explicit decision was to
**leave it as-is for now** — zero known bugs after two Fortran-audit passes,
and refactoring risks silently breaking the "missing branch = valid no-op"
semantics a past bug (§0.6) taught this codebase to protect, for a
style-only win. Recorded here as a real, not-forgotten future improvement,
not silently closed out.

### 4.9 Testing gaps — ✅ closed, see §3

### 4.10 CI has no lint/format or doc-build check
Asked the user: explicit decision was to skip this for now (checked —
`JuliaFormatter` isn't present anywhere in the repo; adding one means
picking a style and reformatting the whole existing codebase, bigger scope
than a mechanical addition).

### 4.11 `output.jl`'s `global_mean` is an unweighted grid mean ✅ fixed
Unlike `min_T_K`/the gust literal, this one didn't stay a judgment call —
checked Fortran's own diagnostic directly and it settled the question.
`greb.model.mscm.f90`'s `gmean()` (`:1497-1513`) does `sum(data*w)/sum(w)`
with `w = cos(lat)`; Julia's `sum(state.Tsmn)/(xdim*ydim)` was a plain
unweighted mean, diverging from Fortran, not matching it. Fixed by reusing
the existing `dxlat_grid` (already `cos(lat)`-proportional, `constants.jl`)
as the weight — the proportionality constant cancels in the ratio, so this
gives an identical result to Fortran's bare `cos(lat)` weight. Diagnostic-
only (`global_mean` is a `println` value, never stored in `MonthlyRecord`),
so this couldn't change any returned/tested model output — confirmed by the
golden regression test passing unchanged.

---

## 5. Suggested order (low-risk → high-value)

1. ~~§0 bug fixes~~ ✅ done — 18 across five passes, all with regression tests.
2. ~~Split into files (§1.2)~~, ~~state structs (§1.1)~~, ~~CI/compat/dead-deps/
   benchmark suite (§3)~~, ~~`RunSpec` (§1.4)~~ ✅ all done.
3. ~~Thread the spatial operators (§2.4)~~ — attempted, benchmarked, reverted
   (3–6x slower at this grid size). ~~`Float32` (§2.3)~~ — analyzed, canceled
   before implementation.
4. ~~Cut `SWradiation!`'s allocation (§2.7)~~, ~~hot-loop/type-stability audit
   (§2.8)~~, ~~precomputation audit (§2.9)~~, ~~test-suite cost audit (§3)~~,
   ~~`Documenter.jl` site (§3)~~ ✅ all done.
5. ~~`SurfaceState` struct (§1.3)~~, ~~testing gaps (§3)~~, ~~§4.1/4.2/4.11
   resolved~~ ✅ all done. §4.8/§4.10 got an explicit user decision to defer,
   not silently dropped.
6. **Next up**: nothing blocking — remaining open items are §4.8 (forcing()
   dispatch style) and §4.10 (CI lint/format), both explicitly deferred by
   user decision rather than left unaddressed.

> ⚠️ Every performance change must be validated against a reference run —
> "faster" only counts if output is unchanged within tolerance. Every
> physics change must be validated against the Fortran reference — grep the
> actual subroutine text, don't recall it from memory.

---

## 6. Data organization & JLD2 compression — investigated, findings documented

Full comparison doc: `claude/DATA_ORGANIZATION_OPTIONS.md`. Triggered by the
CMIP5/ENSO forcing wiring work (§4's `:rcp85`/`:elnino`/`:lanina` stubs) making
it worth asking how `Data/` and its converted `greb_dataset_jld2/` output should
be organized and shared. Everything below is measured, not estimated — the
conversion script was actually run (`Data` → a scratch `greb_dataset_jld2/`) and
real `jldopen`/`tar` calls were timed, rather than assuming behavior.

### 6.1 `Data/`'s layout doesn't match its own docs
`DATA_README.md` and `convert_greb_to_jld2.jl`'s default (`input_path =
Data/input`) both expect an `input/` subdirectory; on disk, all ~90 files sit
directly under `Data/` (confirmed — `Data/input/` doesn't exist, had to point
the conversion explicitly at `Data`). Recommended fix: reorganize by category
(mirroring `convert_all`'s own static/solar/climatology routing) rather than
just used-vs-unused, since with the CMIP5/ENSO wiring in progress there's almost
nothing left to quarantine — only 2 of 49 root-level files
(`erainterim.evaporation.clim`, `erainterim.omega.vertmean.nomean.clim`, 25.7
MiB) are genuinely unreferenced by any code path.

**✅ Confirmed not a porting gap, files removed.** Grepped every Fortran
source variant available (`greb.model.mscm.f90`/`greb.shell.mscm.f90`, both
the NCEP-only Downloads copy and the ERA-Interim/CMIP5/ENSO-capable
`greb-official-official` copy) for `evap`/`evaporation`/`nomean` — zero
`open()` statements and zero filename occurrences for either file in any
variant checked. Neither was ever consumed by the original model; the
Julia port's non-use correctly mirrors it, not an oversight. Both raw files
removed from `Data/`, and `DATA_README.md` updated to stop listing
`erainterim.evaporation.clim.bin` as required.

### 6.2 JLD2 compression: real numbers, not assumed ones
Built the actual `greb_dataset_jld2/` output once (605.4 MiB, 53 files — JLD2's
container format adds ~0.04% overhead over raw `.bin`) and measured
`jldopen(path, "w"; compress=true)` against six representative files:

| File | Compression ratio | Read slowdown vs. uncompressed |
|---|---|---|
| `ncep.tsurf...clim.jld2` | 80.1% | **9.5×** |
| `cmip5.tsurf...forcing.new.jld2` | 86.6% | **7.9×** |
| `erainterim.tsurf.elnino.forcing.jld2` | 93.2% | **20.4×** |
| `global.topography.jld2` (static) | 29.5% | 4.0× |
| `solar_radiation.clim.jld2` | 40.4% | 2.0× |
| `solar_eccentricity.jld2` (grouped) | 40.8% | 11.4× |

The 3D climate fields (619 of 635 MB) barely compress — 80–93% of original —
because float32 reanalysis/ensemble-mean data has limited redundancy, while
costing 8–20× longer to read every time. A real whole-tree `tar -I 'gzip -9'`
run got the same story at dataset scale: only 530.6 MB from 634.8 MB (83.6%,
measured, 42s), cross-checked against the per-category ratios above to within
0.6%. **Conclusion: compress only for distribution** (a plain `tar`/`xz`
archive of `greb_dataset_jld2/`, decoupled from the JLD2 read path entirely),
never the copy `load_greb_jld2!` actually reads — the size win (~16–22%
dataset-wide) doesn't justify an 8–20× read-time risk if anyone points the
loader at the compressed copy by mistake.

### 6.3 JLD2 grouping: combine by load-pattern, not by source
Initial instinct was "combine per-field files by their upstream source"
(one `cmip5.jld2`, one `erainterim.jld2`, one `ncep.jld2`). Tested this
directly instead of assuming: merged 5 real CMIP5 fields (67.3 MB) into one
multi-key file and timed reads.

| | Combined (5-key) file | Separate files |
|---|---|---|
| Read only 1 field | 24.0 ms | 10.3 ms — combined costs **2.3×** more |
| Read all 5 fields | 45.0 ms | 50.2 ms — combined is *faster* |

So combining helps exactly when a group is always read in full, and hurts when
only a subset is needed — but the 2.3× penalty is ~14 ms in absolute terms,
irrelevant for a `load_greb_jld2!` call that runs once at startup, not per
timestep. The real reason to avoid combining by source isn't performance:
- `ncep.*`/`erainterim.*` are each read in different subsets depending on
  `dataset=:ncep`/`:era` and which experiment is active — combining by source
  would couple independently-varying load conditions into one file.
- `cmip5.*` is 16 raw files but probably only ~5 distinct fields will feed
  `*_anom_cc`, with an undecided `.bin`-vs-`.new.bin` duplicate for most of
  them (see the earlier session's byte-level diff — same shape, genuinely
  different values, no code or docs pick a winner). Combining now bakes both
  variants into a shared blob that would need unpacking and rewriting later
  just to drop the loser.
- The one case where combining *is* recommended — merging the three
  flux-correction files into `flux_corrections.jld2` — works precisely because
  `load_flux_corrections_jld2!` (`src/io.jl:82-100`) already always loads all
  three together, matching the "always read in full" condition the benchmark
  above says combining rewards. `solar_eccentricity.jld2`/`ipcc_scenarios.jld2`
  (already-combined, pre-existing) fit the same rule: a scenario index or CO2
  lookup table is inherently read as one shape, never partially.

  **✅ Implemented.** Measured directly on the real 3 files (38.51 MiB
  combined) before implementing: reading all 3 from one combined file took
  19.89 ms vs. 30.78 ms from 3 separate files (**~35% faster**, matching the
  "always read in full" case above), at 40,367,196 vs. 40,382,127 bytes
  (no size penalty). Implemented in `scripts/convert_greb_to_jld2.jl`
  (`convert_flux_corrections`, excluded from `convert_all`'s per-field loop
  via `FLUX_CORRECTION_NAMES`) and `src/io.jl`'s `load_flux_corrections_jld2!`
  (now reads 3 keys from one file instead of 3 files); `test/runtests.jl`'s
  synthetic-dataset test updated to match. Full suite re-run: 319/319 pass.

---

## 7. Fortran switch/experiment validation — report-only, no code changes

Cross-referenced every `PhysicsConfig` field and `cfg.experiment` branch
against the actual Fortran reference: two `greb.shell.mscm.f90` variants (an
older one without ERA-Interim/CMIP5/ENSO support, and the newer
`greb-official-official` one with `log_clim` dataset-switching plus CMIP5
`log_exp=230` and ENSO `log_exp=240/241` forcing) and four `run.greb.*.csh`
scripts (`decon_mean_climate` = `log_exp=1`, `decon2xco2` = `log_exp=10`,
`scenarios` = the full `log_exp=20-105` sensitivity table, `hydro` = the
CMIP5/ENSO forced-boundary script). Per explicit user decision, this pass is
documentation-only — the gaps below are real but deliberately left unfixed.

### 7.1 `:rcp26`/`:rcp45`/`:rcp60` — real gap, data is ready but code isn't ✅ documented
Fortran `log_exp=96/97/98` load `ipcc.scenario.rcp26/45/6.forcing.txt`
directly. Julia's `forcing()` explicitly errors instead:
`tendencies.jl:159-160` (`:rcp26`), `:162-163` (`:rcp45`), `:165-166`
(`:rcp60`) — all three say "requires external CO₂ data file. Not yet
implemented." But `ipcc_scenarios.jld2` already contains working
`"rcp26"`/`"rcp45"`/`"rcp6"` keys (confirmed by actually running
`convert_greb_to_jld2.jl` earlier this session — the CO2 data converts
cleanly), and the exact same `(:ssp119,...,:historical_co2) =>
cfg.co2_scenario[yr]` dispatch already at `tendencies.jl:175-180` would work
unchanged for these three. `test/runtests.jl:441-444` confirms the gap is
intentional-and-tested (`@test_throws ErrorException`), not an oversight.
Closing this is a ~10-line change (add the 3 symbols to
`_CO2_SCENARIO_SYMBOLS`/`create_experiment_config`, delete the 3 `error()`
branches) — left undone per explicit user decision to keep this pass
report-only.

### 7.2 `:custom_co2` (EXP=100) — no file-loading path
Fortran's `log_exp=100` lets a user point at an arbitrary
`../input/<name>.txt` CO2 file. Julia's `:custom_co2` also just errors
(`tendencies.jl:171-172`) — there's no `PhysicsConfig` field for a
user-supplied file path at all, so this isn't a small dispatch-table gap
like 7.1; it needs a new field plus a loader. Lower priority than 7.1.

### 7.3 `_dmc`/`_drsp` switch families are fully wired, but have no combined preset
Every individual switch from `run.greb.decon_mean_climate.csh`
(`log_cloud_dmc`/`log_ocean_dmc`/`log_atmos_dmc`/`log_co2_dmc`/
`log_hydro_dmc`/`log_qflux_dmc`) and `run.greb.decon2xco2.csh`
(`log_topo_drsp`/`log_cloud_drsp`/`log_humid_drsp`/`log_ocean_drsp`/
`log_hydro_drsp`) exists as a real `PhysicsConfig` field and is read at a
real gate — `model.jl:50,55,60,67,71,75,79,120,141,284-285`,
`radiation.jl:116`, `circulation.jl:301`, `hydrology.jl:33`, `ocean.jl:15,63`,
`output.jl:154,158` — so none of these switches are dead. But
`create_experiment_config` (`config.jl:104-157`) has no preset that flips a
whole family at once the way each Fortran run script does; only
`:constant_topo` sets a single `_drsp` field. Reproducing either Fortran run
script today means using the bare `PhysicsConfig(log_cloud_dmc=false, ...)`
keyword constructor directly instead of a named experiment. A
discoverability/convenience gap, not a correctness bug.

### 7.4 `:sst_plus1` has no Fortran `log_exp` counterpart — confirm intentional
`model.jl:273,374-378` implements a `:sst_plus1` experiment (forces
ocean-point `Ts = Tclim+1.0`) matching no `log_exp` value in any of the four
run scripts or their comment tables. `forcing()` has no case for it either —
it silently falls through to the generic default at `tendencies.jl:67-68`
(`CO2 = cfg.co2_concentration`), which only produces correct behavior
because `:sst_plus1`'s real logic lives entirely in `model.jl`. Not a bug as
implemented, but worth an explicit confirmation that this is a deliberate
Julia-only addition rather than a mistranslated experiment number from a
Fortran source file not included in what was reviewed here.

### 7.5 `:a1b_scenario`/`:a1b_enhanced` — hardcoded formula, not a data file
Fortran's `log_exp=95` ("IPCC A1B scenario") is data-driven like the other
IPCC scenarios, but no `ipcc.scenario.a1b*.txt` exists anywhere in `Data/`
(confirmed absent, both in this repo and in the reference `input/` folder).
Julia's `:a1b_scenario` (`tendencies.jl:80-90`) and `:a1b_enhanced`
(`tendencies.jl:118-128`, byte-identical duplicated code) instead hardcode a
piecewise-linear 1950→2000→2050→2100 (310→370→520 ppm) approximation. A
reasonable stand-in given the source table was never supplied, but it's an
approximation, not a reproduction of Fortran's real A1B forcing — worth
flagging so it isn't mistaken for a validated match later.

### 7.6 Confirmed correct / already-resolved — no action
- **`log_clim==1` in `set_hydrology_parameters!`** (`config.jl:175-177`)
  doesn't match Fortran's `log_clim ∈ {0 (ERA), -1 (NCEP)}` convention
  (`greb.shell.mscm.f90:25,40` in the newer reference) — but this exact
  discrepancy was **already investigated and resolved** in §4.2: `log_clim`
  and the `dataset=:era/:ncep` file-selection kwarg were deliberately kept
  orthogonal by user decision. Not re-opened here, just cross-referenced so
  this validation pass isn't mistaken for having missed it.
- **Obliquity/eccentricity range**: the Fortran run-script comments describe
  `OBL∈[-250,900]`/`ECC∈[-30,30]`, but the actual `solar_forcing_scenarios/`
  files use a step-5 convention (0..230) — already correctly handled by
  `convert_greb_to_jld2.jl`'s glob-based discovery (its own docstring
  documents fixing a prior hardcoded-`0:25:230`-stride bug), and Julia's
  generic `orbital_index::Int` reads whatever `coords` values the files
  actually contain. Confirmed already correct; the Fortran comment's stated
  range is simply stale relative to the real file set.
- **`dradius`/`earth_sun_distance_pct`**: naming differs, but
  `(1/(1+0.01*pct))^2` (`tendencies.jl:156`) is a faithful translation of
  Fortran's `log_exp=37` radius-change physics. Confirmed correct.
- **Regional CO2 experiments (`log_exp=40-47`)**: all 8 symbols exist and
  dispatch correctly. Noted only: `:regional_co2_ocean`/`_land_ice`/
  `_winter`/`_summer`'s masks are set dynamically inside `forcing()`
  (`tendencies.jl:187-232`) while the other four are reset in `init_model!`
  (`model.jl:21-38`) — an architectural asymmetry between the two halves of
  the same experiment family, not a bug (both paths produce correct masks),
  worth a note so it doesn't read as an oversight later.

---

## Changelog of this document

- **2026-08-10 (flux-correction combine implemented + evaporation/nomean
  resolved, §6.1/§6.3)**: two follow-ups on §6's data-organization findings.
  First, confirmed via direct grep of every available Fortran source variant
  (Downloads' NCEP-only `greb.shell.mscm.f90`/`greb.model.mscm.f90`, and
  `greb-official-official`'s ERA-Interim/CMIP5/ENSO-capable copies) that
  `erainterim.evaporation.clim.bin`/`erainterim.omega.vertmean.nomean.clim.bin`
  are unread by any variant — not a Julia-porting gap; both raw files have
  since been removed from `Data/`, `DATA_README.md` updated to match. Second,
  actually implemented §6.3's flux-correction merge (previously recommended
  by analogy only): measured the real 3 files directly first (19.89ms vs.
  30.78ms, ~35% faster, no size penalty), then implemented
  `convert_flux_corrections`/`FLUX_CORRECTION_NAMES` in
  `scripts/convert_greb_to_jld2.jl` and rewrote `load_flux_corrections_jld2!`
  (`src/io.jl`) to read the combined `flux_corrections.jld2`; updated
  `test/runtests.jl`'s synthetic-dataset test and both `README.md`/
  `DATA_README.md` to match. Full suite re-run clean: 319/319 pass.
- **2026-08-10 (Fortran switch/experiment validation, §7)**: cross-referenced
  every `PhysicsConfig` field and `cfg.experiment` branch against 8 supplied
  Fortran reference files (two `greb.shell.mscm.f90` variants, plus
  `run.greb.decon_mean_climate.csh`/`decon2xco2.csh`/`scenarios.csh`/
  `hydro.csh`). Found and documented (report-only, by explicit user
  decision): `:rcp26`/`:rcp45`/`:rcp60` error out despite their CO2 data
  already converting cleanly (§7.1); `:custom_co2` has no file-loading path
  (§7.2); the Fortran `_dmc`/`_drsp` switch families are fully wired but have
  no combined experiment preset, unlike their Fortran run scripts (§7.3);
  `:sst_plus1` has no Fortran `log_exp` counterpart (§7.4); `:a1b_scenario`/
  `:a1b_enhanced` hardcode a formula since no source data file exists
  (§7.5). Confirmed as already-correct or already-resolved: the `log_clim`
  discrepancy (already closed in §4.2), the obliquity/eccentricity numeric
  range (already fixed by the glob-based discovery in
  `convert_greb_to_jld2.jl`), `earth_sun_distance_pct`, and the regional CO2
  experiment family's dynamic-vs-static mask asymmetry (§7.6).
- **2026-08-10 (data organization & JLD2 compression investigation, §6)**:
  added `claude/DATA_ORGANIZATION_OPTIONS.md`, prompted by the CMIP5/ENSO
  forcing wiring work. Actually ran `convert_greb_to_jld2.jl` against `Data`
  into a scratch output (605.4 MiB, 53 files) to get measured numbers instead
  of estimates; measured `compress=true` size/read-time tradeoffs on 6
  representative files (8–20× slower reads for 80–93% compression ratios on
  the 3D climate fields); measured real `tar+gzip`/`tar+xz` whole-tree
  archives (only 16–22% smaller, cross-checked against per-category ratios to
  within 0.6%); and — in response to the user directly challenging the initial
  recommendation — tested combining 5 CMIP5 fields into one multi-key file
  head-to-head against separate files, finding a real but small (2.3×, ~14 ms)
  single-key-read penalty that's irrelevant for a once-at-startup load. Revised
  conclusion: group by load-pattern (things always read together, like the
  flux-correction trio) not by upstream source (`cmip5`/`erainterim`/`ncep`),
  since those sources are read in different subsets depending on
  `dataset=:ncep`/`:era` and active experiment, and CMIP5 additionally has an
  undecided `.bin`-vs-`.new.bin` duplicate that shouldn't be baked into a
  shared blob yet.
- **2026-08-06 (seventh pass: `SurfaceState` struct, testing gaps, §4
  decisions, README)**: implemented §1.3's argument-struct collapse — a
  `SurfaceState` struct for `Ts`/`Ta`/`To`/`q`, reusing the already-existing
  `tend` `NamedTuple`/`ws::CirculationWorkspace` for everything else instead
  of inventing a second struct (zero new allocations). Closed all four §3
  testing gaps: `qflux_correction!` (verified the `Ta`-correction asymmetry
  matches Fortran exactly, not a bug), the ~20 unreachable `:experiment`
  symbols (direct-construction sweep — caught a real footgun while writing
  it: `forcing()` only runs during the scenario loop, so `RunSpec(scnr=0)`
  would have made the sweep vacuous), `build_monthly_climatology`/
  `apply_scenario_anomalies` (direct unit tests), and both JLD2 loaders'
  untested "file present" branches (synthetic `mktempdir`+`jldopen`
  datasets). Got the user's decisions on §4.1/4.2/4.8/4.10 and resolved
  4.1/4.2 accordingly (removed 4 config switches with zero Fortran basis,
  documented 3 that mirror genuinely-dead Fortran knobs; documented the
  `log_clim`/`dataset` relationship without linking them); kept 4.8/4.10
  deferred by explicit decision, not silently dropped. Resolved 4.11 myself
  via direct Fortran comparison (no longer a judgment call once checked) —
  `global_mean` now area-weights by `cos(lat)` like Fortran's own `gmean()`,
  diagnostic-only so it can't affect any returned/tested output. Refreshed
  `README.md`: fixed stale TOC anchors (added real `Prerequisites`/
  `Installation` sections, renamed the second, duplicate "Quick Start" to
  "Running the Model", added a real `Project Structure` tree), modernized
  the "Running the Model" section to the actual plain-Julia API instead of
  Pluto-only instructions, updated "Known Issues" with the confirmed
  `qflux_correction!` finding, and dropped a stale `Statistics` dependency
  row (removed as dead in an earlier pass, still listed here).
- **2026-08-06 (sixth pass: test-suite cost audit + allocation fix + doc
  tooling)**: live-timed the full test suite (2m36.9s) and found a real,
  currently-failing assertion (`min_T_K`'s float-literal equality, fixed
  with `≈`); confirmed no testset duplicates another, trimmed the one safe
  cut found (`co2_part` mask-reset test), and added a `RUN_GOLDEN` env var
  to skip the golden snapshot test locally. Used the existing benchmark
  harness plus `@allocated` (no `@profview` needed) to find and fix
  `SWradiation!`'s 40704-byte/48-alloc hotspot — the entire allocation
  footprint of `tendencies!` — with `@views`; every kernel now benchmarks
  at 0 bytes. Audited hot-loop field/config localization, type stability,
  and remaining caching/precomputation opportunities — both audits came
  back clean, nothing left to change, documented as such (§2.8, §2.9).
  Flagged (not fixed) `global_mean`'s missing `cos(lat)` area weighting as
  a scientific-judgment question, same treatment as `min_T_K`/the gust
  literal (§4.11). While wiring up a `Documenter.jl` site, found a real,
  previously-unknown bug (§0.19): a `# ── notebook cell ── ` comment
  between a docstring and its function/struct silently detaches the
  docstring from Julia's doc system — affected 19 of 36 exports. Fixed by
  moving the comment before the docstring everywhere it occurred (verified
  with `Docs.hasdoc` across all exports, 0 missing after) and added the 3
  docstrings that were genuinely absent (`xdim`/`ydim`/`nstep_yr`). Shipped
  `docs/` (Documenter site: home page, tutorial adapted from
  `examples/run_greb.jl`/README's Quick Start, `@autodocs` API reference)
  and `.github/workflows/docs.yml` (deploys to GitHub Pages on push to
  `main`).
- **2026-08-06 (fifth pass: bug sweep + small fixes + compaction)**: fresh
  Fortran-audit sweep found 8 more real findings (§0.11–0.18, including one
  active in every default run — the `circulation!` `dX_conv` leak, §0.11), 7
  of which were fixed with regression tests and verified against
  `greb.model.mscm.f90` directly (5 of 8 independently re-verified, not just
  trusted from the sub-agent that found them). The 8th (§0.14, `min_T_K`) was
  tried, then reverted on reflection that Fortran's own value doesn't
  actually make physical sense to copy literally — kept alongside the
  `log_eva==1` gust literal (§0's other "Investigated, not changed" entry) as
  a second example of the same judgment call. Analyzed
  §2.3's `Float32` proposal, found a real correctness hazard (§2.3), then
  canceled implementation per user request. Implemented §1.4 (`RunSpec`) and
  partially implemented §4.4's `co2_part` hoist (4 of 6 regional experiments;
  the other 2 depend on data that doesn't exist yet at init time). Resolved
  §4.6 as obsolete (its premise no longer applies post-§1.1). Recomputed the
  golden regression test's reference values (real output changed). Compacted
  this entire document — §0 and the changelog were the heaviest offenders,
  each bug/entry cut from ~15-20 lines to ~3-6 while preserving every exact
  field name, formula, and trigger condition needed to avoid re-introducing
  a fixed bug.
- **2026-08-06 (fourth pass: threading attempt + low-risk cleanup)**: added
  a golden regression test, then implemented and benchmarked threading for
  `diffusion!`/`advection!` — reverted after benchmarks showed a 3-6x
  slowdown at `-t auto` (§2.4). Implemented §4.4's `:full_model` fast path,
  §4.5's `hydro!` constant hoist, §4.7's missing docstrings (~22 exports).
  Cut the test suite's own runtime (20-run cross product → 5-run zipped
  sweep, same coverage).
- **2026-08-06 (third pass: Fortran audit + state-struct refactor)**: direct
  line-by-line comparison against `greb.model.mscm.f90` found and fixed 4
  bugs (§0.7–0.10). Then replaced ~40 mutable module globals with
  `ClimateFields`/`ModelState` (§1.1), validated against real data to
  ~7e-12 absolute, with an unplanned ~7x runtime / ~170x allocation cut.
- **2026-08-06 (second pass)**: independent re-read of the split `src/` tree
  found 3 more bugs (§0.4–0.6). Mechanical cleanup: dead globals, unused
  `Statistics` dependency, 3 misplaced docstrings. Added §4 (documented,
  not fixed).
- **2026-08-06 (first pass)**: fixed 3 bugs (§0.1–0.3). Split the module into
  files (§1.2). Removed dead globals (§1.6). Dropped unused `NCDatasets`/
  `StaticArrays`, added `[compat]` bounds, `PrecompileTools` workload (§2.6),
  CI, benchmark suite.
- **2026-08-05**: initial refresh against post-JLD2-conversion `src/GREB.jl`;
  found the `ws.cE`/`ws.cE_buf` bug (§0.1) while re-reading for this pass.
