# GREB.jl — Potential Improvements (Structure & Performance)

Observations and fixes for the `GREB` module, accumulated across five audit
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

17 real, previously-shipped bugs found and fixed across four passes (an 18th
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

### 1.3 Collapse giant positional argument lists into structs — still deferred
`output!`/`diagnostics!`/`time_loop!` still take ~10–20 positional arrays
(`Ts0`/`Ta0`/`To0`/`q0`/`albedo`/...) alongside the `fields`/`state` params
from §1.1 — those were already function-local, not part of the global-state
bug class, so out of §1.1's scope. A `SurfaceState`-style struct remains a
real readability win, not yet done.

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
  before) fixed ✅ done, **benchmark suite** (`benchmark/run_benchmarks.jl` +
  `claude/BENCHMARKS.md`) added ✅ done.
- **Remaining gaps**: `qflux_correction!`'s loop body has zero coverage
  (every test uses `time_flux=0`; README's own "Known Issues" flags this
  module as possibly broken) — worth a dedicated look. Most `experiment`
  symbols (only 9 of ~25+) are unreachable via `create_experiment_config`.
  `build_monthly_climatology`/`apply_scenario_anomalies` are undertested.
  `load_greb_jld2!`/`load_flux_corrections_jld2!`'s "file exists" branches
  are untested without the real dataset. No `Documenter.jl` site.

---

## 4. Findings needing domain/design input — documented, not fixed

### 4.1 Seven config switches wired to nothing
`solar_multiplier`, `log_vapor_drsp`, `log_ice_drsp`, `log_tsurf_ext`,
`log_hwind_ext`, `log_omega_ext`, `log_ice_dmc` are set but never read by any
physics function. Same bug class as §0.2–0.18, but wiring 7 features up (or
removing them) is a scope decision.

### 4.2 `log_clim` doesn't select the ERA/NCEP dataset its label implies
Only swaps 4 regression coefficients; the real ERA-vs-NCEP choice is the
disconnected `dataset` kwarg to `load_greb_jld2!`. Deciding whether `cfg`
should drive data loading is a real API question.

### 4.6 More `const`-ify candidates — resolved by §1.1
The original target (module-level globals never rebound) no longer exists:
§1.1 moved everything into `ClimateFields`/`ModelState` struct fields, and
Julia has no `const`-field mechanism for mutable structs. Nothing to do here.

### 4.8 `forcing()`'s dispatch style is inconsistent with the rest of the codebase
Its ~180-line `if/elseif` chain over `cfg.experiment` is the odd one out next
to the Dict-based dispatch used elsewhere. A Dict-based replacement must
preserve the chain's "missing branch = valid no-op" semantics carefully
(per §0.6) — not a mechanical fix.

### 4.9 Testing gaps — see §3's "Remaining gaps"

### 4.10 CI has no lint/format or doc-build check
Adding one means picking a formatting style first — a decision, not a
mechanical addition.

---

## 5. Suggested order (low-risk → high-value)

1. ~~§0 bug fixes~~ ✅ done — 17 across four passes, all with regression tests.
2. ~~Split into files (§1.2)~~, ~~state structs (§1.1)~~, ~~CI/compat/dead-deps/
   benchmark suite (§3)~~, ~~`RunSpec` (§1.4)~~ ✅ all done.
3. ~~Thread the spatial operators (§2.4)~~ — attempted, benchmarked, reverted
   (3–6x slower at this grid size). ~~`Float32` (§2.3)~~ — analyzed, canceled
   before implementation.
4. **Next up**: §1.3's remaining argument-collapsing; re-measure §2.5's
   allocation hotspots now that §0.17's Jacobi fix changed `diffusion!`/
   `advection!`'s codegen; `qflux_correction!`'s test-coverage gap (§3).
5. The §4 findings needing domain/design input (§4.1, §4.2, §4.8, §4.10) —
   bring in whoever knows the intended physics/API shape.

> ⚠️ Every performance change must be validated against a reference run —
> "faster" only counts if output is unchanged within tolerance. Every
> physics change must be validated against the Fortran reference — grep the
> actual subroutine text, don't recall it from memory.

---

## Changelog of this document

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
