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

### 2.3 `Float32` for climatology/working fields — analyzed and benchmarked, not implemented
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

### 2.4 Thread the spatial operators — tried, measured, reverted
`Threads.@threads :static` on `diffusion!`/`advection!`'s outer row loop,
validated bit-identical output, but **3–6× slower** end-to-end at `-t auto`/14
threads (`tendencies!`: 4.2ms→23.7ms) — the grid (96×48, ~1µs/row) is too
small for the scheduling overhead to pay off at any thread count tried.
Reverted to plain `for k in 1:ydim` loops. **Takeaway**: don't retry plain
`Threads.@threads` at this grid/kernel shape without a reason to expect a
different result; `Polyester.@batch` or a much larger grid might fare
better.

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

### 2.12 Missing `@turbo` spots — measured across 4 candidates, not implemented
Every hot kernel in `circulation.jl`/`hydrology.jl`/`ocean.jl` uses `@turbo`
except a handful of plain-`@.`-broadcast spots. Checked whether the same fix
generalizes, prioritizing by call frequency (per-substep code compounds:
`circulation!` runs 24 substeps × 2 fields = 48×/main-timestep):

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

Not benchmarked but plausible same-scale candidates for later: `state.jl:
206-218`'s per-timestep monthly accumulation (9 elementwise adds, runs every
step) and `SWradiation!`'s remaining combine/flux lines
(`radiation.jl:59,66`). None of these are implemented — real, safe wins
whenever someone next touches this code, alongside §2.11's larger threading
opportunity.

### 2.11 Parallelize `tendencies!`'s independent stages — benchmarked, not yet wired in
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
bit-identical vs. sequential. **Not yet implemented in `tendencies.jl`** —
requires `-t 2`+, pending a decision on the `-t 1`-default compatibility
tradeoff.

Other identified, unimplemented opportunities (not benchmarked further):
ensemble/parameter-sweep parallelism via `Threads`/`Distributed`/GPU arrays
(likely the bigger practical payoff given how the model is actually used —
sensitivity sweeps, CO₂ scenarios); low-cost `Float32`/GPU-array swaps via
multiple dispatch (the physics kernels have no hardcoded `Float64`, only
`CirculationWorkspace`'s fields do, per §2.3); AD-based calibration via
`Enzyme.jl`/`ForwardDiff.jl` for gradient-based parameter tuning.

**Combined with §2.12's turbo fixes and integration-tested over a full
1-year run: ~2.8× (2.83–2.88s → 1.00–1.04s), better than either alone —
they compose rather than interfere. Output verified bit-identical (mod.
float noise) against both the golden-regression reference and a fresh
baseline run. One measurement gotcha caught along the way: the first,
un-warmed-up timing showed the combo 2.5× *slower*, purely because its new
functions hadn't been JIT-compiled yet (unlike `greb_model!`'s
`PrecompileTools`-baked ones) — resolved with a warm-up call; a real
implementation would need this folded into `@compile_workload`.**

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

### 3.0 Parallelizing the test suite across processes — benchmarked, not implemented
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
closer to the ceiling. Not implemented as an actual CI change (would need a
load-balanced shard split plus a CI matrix or `ReTestItems.jl`-style
runner) — a benchmarked, real opportunity for whoever next touches CI
turnaround time.

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
- Landing §2.11's benchmarked 3-way `tendencies!` parallelization, pending a
  decision on the `-t 1`-default compatibility tradeoff.
- Landing §2.12's 4 benchmarked `@turbo` fixes (~0.46s/year combined, low
  risk, independent of §2.11 and composes with it).
- Ensemble/GPU parallelism and AD-based calibration (§2.11) remain
  documented ideas only, not benchmarked.

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

## 7. Fortran switch/experiment validation — report-only, no code changes

Cross-referenced every `PhysicsConfig` field and `cfg.experiment` branch
against two `greb.shell.mscm.f90` variants and four `run.greb.*.csh`
scripts. Per explicit user decision, documentation-only — gaps below are
real but deliberately left unfixed.

### 7.1 `:rcp26`/`:rcp45`/`:rcp60` — data is ready, code isn't
Fortran `log_exp=96/97/98` load their CO2 files directly; Julia's
`forcing()` explicitly errors (`tendencies.jl:159-166`), even though
`ipcc_scenarios.jld2` already has working `"rcp26"`/`"rcp45"`/`"rcp6"` keys
and the existing scenario dispatch (`tendencies.jl:175-180`) would work
unchanged. `test/runtests.jl:441-444` confirms the gap is
intentional-and-tested. Closing this is a ~10-line change — left undone per
explicit decision to keep this pass report-only.

### 7.2 `:custom_co2` (EXP=100) — no file-loading path
Fortran's `log_exp=100` lets a user point at an arbitrary CO2 file; Julia's
`:custom_co2` just errors (`tendencies.jl:171-172`) with no
`PhysicsConfig` field for a user-supplied path — needs a new field plus a
loader, not just a dispatch-table entry. Lower priority than 7.1.

### 7.3 `_dmc`/`_drsp` switch families have no combined preset
Every individual switch is a real, wired `PhysicsConfig` field, but
`create_experiment_config` has no preset that flips a whole family at once
the way each Fortran run script does — a discoverability/convenience gap,
not a correctness bug.

### 7.4 `:sst_plus1` has no Fortran `log_exp` counterpart — confirm intentional
`model.jl:273,374-378` implements it correctly as a Julia-only addition
(logic lives entirely in `model.jl`, `forcing()` needs no special case) —
worth an explicit confirmation this is deliberate, not a mistranslated
experiment number.

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
