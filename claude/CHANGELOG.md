# GREB.jl — Audit Log

Forensic, pass-by-pass companion to [`IMPROVEMENTS.md`](IMPROVEMENTS.md). That
document tracks the *current* state of the code; this one is an append-only
diary of how it got there — every fix's discovery method, the reasoning
behind corrections and reverted attempts, and benchmark methodology. Entries
are never edited retroactively; a later pass that revises an earlier finding
adds a new entry and cross-references the old one rather than rewriting it.

---

## Bug-discovery narrative (passes 1–5, 2026-08-05/06)

18 real, previously-shipped bugs found and fixed across five passes (a 19th
finding was investigated and reverted — see the third-pass-correction note
under 0.4 below, and 0.14), each with a regression test in
`test/runtests.jl`. Passes 1–3 were found by grepping every config-field
write against its read sites; pass 4 was a direct line-by-line comparison
against the Fortran reference `greb.model.mscm.f90`, and every fix there
cites the exact Fortran line(s) confirmed by direct read (not recalled from
memory). Final, current-state summary of each fix lives in
`IMPROVEMENTS.md` §0's table; this section is the "how it was found and why"
for each one.

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
(0.14 above) is a third case of the same pattern: tried matching Fortran, then
reverted after judging the Fortran value itself wasn't worth copying.

---

## Changelog of this document

- **2026-08-12 (implemented §7.1/§7.2/§7.3)**: §7's Fortran switch/experiment
  validation was originally report-only by explicit decision; asked to
  actually implement the 3 actionable gaps (§7.4/§7.5/§7.6 stay
  documentation-only — §7.4 separately reconfirmed intentional).
  Cross-checked the real Fortran run scripts directly
  (`run.greb.decon_mean_climate.csh`, `run.greb.decon2xco2.csh`,
  `run.greb.hydro.csh`, found in the sibling `greb-official-official/`
  folder) rather than relying on the earlier audit's summary alone — this
  confirmed the exact family membership for §7.3 (`decon_mean_climate.csh`'s
  namelist sets 6 `_dmc` fields + 5 circulation switches;
  `decon2xco2.csh`'s sets 5 `_drsp` fields, notably *not* `log_crcl_drsp` +
  the same 5 circulation switches) and the Fortran text-file mechanism
  behind §7.1/§7.2 (both `log_exp=96-98` and `log_exp=100` read the same
  `year CO2` unit-26 text file format, differing only in which file gets
  symlinked to it).
  - §7.1: `forcing()` (`tendencies.jl`) now dispatches `:rcp26`/`:rcp45`/
    `:rcp60` through the same `cfg.co2_scenario[year]` lookup already used
    by `:ssp*`/`:historical_co2`, replacing the 3 unconditional `error()`
    branches. `model.jl`'s `_CO2_SCENARIO_SYMBOLS`/`_CO2_SCENARIO_KEY` gained
    the 3 symbols (`:rcp60 => :rcp6`, since the jld2 key doesn't match the
    symbol name the way `:rcp26`/`:rcp45` do). `create_experiment_config`
    gained a matching preset.
  - §7.2: `PhysicsConfig` gained `custom_co2_path::String`; `src/io.jl`
    gained `load_custom_co2_scenario`, a plain-text `year CO2` parser
    (skips blank/`#` lines) matching the Fortran file format exactly;
    `model.jl` loads it into `cfg.co2_scenario` at scenario start (clear
    error if the path is unset), reusing §7.1's `forcing()` dispatch
    unchanged. `create_experiment_config(:custom_co2; co2_path=...)` sets it.
  - §7.3: added `create_experiment_config(:decon_mean_climate; ...)` and
    `create_experiment_config(:decon_2xco2; ...)`, each keyword-argument per
    family member (default `true`, matching the Fortran namelist defaults),
    mirroring the two Fortran run scripts. Pure discoverability/convenience
    addition — no `forcing()`/dispatch changes needed, since every
    individual switch already worked and both experiments' CO2 handling
    falls through to existing, already-correct code paths (`:decon_2xco2`
    is architecturally identical to the existing `:co2_double` pattern).
  Replaced the old "not yet implemented" placeholder test
  (`test/runtests.jl`, previously asserting `@test_throws` for all 4
  symbols) with real synthetic-`ipcc_scenarios.jld2`/temp-CO2-file runs for
  `:rcp26`/`:rcp45`/`:rcp60`/`:custom_co2`, plus new tests for the two decon
  presets (default-kwarg values, one-override propagation, and a
  `greb_model!` smoke run each). Full suite re-run clean: 249 light + 189
  heavy = 438 pass (up from 319 pre-pass, reflecting the new tests added
  across this and the intervening §8 pass).
- **2026-08-12 (implemented the §8 bug-sweep findings)**: §8's 5 findings
  were originally left documented-but-unfixed to avoid colliding with the
  in-progress performance work; asked to implement them once that settled.
  Re-checked each finding's cited lines against the current source first
  (concurrent performance-pass edits had touched `hydrology.jl`/`model.jl`/
  `tendencies.jl`/`state.jl` in the meantime) — all 5 were confirmed still
  present, unchanged in substance, before fixing:
  - §8.1: `tendencies.jl`'s two `icmn_ctrl1 = @view icmn_ctrl[:,:,1]` sites
    now average over the full 12-month dimension.
  - §8.2: `hydro!`'s `log_eva==1` gust constants now add the existing
    `gust_land`/`gust_ocean` (`4.0`/`9.0`) constants, matching Fortran's
    carried-over `abswind` base term.
  - §8.3: removed `hydro!`'s spurious `min_dq`/`-0.9q/Δt` clamp on
    `dq_rain` alone; `output.jl`'s combined-`dq` clamp (§0.18) is now the
    only clamp, matching Fortran.
  - §8.4: `build_monthly_climatology`/`compute_annual_ice_climatology`
    (`postprocess.jl`) now return only the final 12 records of their input,
    not a multi-year average — required updating 2 existing tests
    (`build_monthly_climatology/apply_scenario_anomalies`,
    `compute_annual_ice_climatology`) whose 2-year-average assertions
    (`m+6`) encoded the old buggy behavior; both now assert the
    Fortran-matching final-year value (`m+12`). Confirmed the golden
    regression test is unaffected (`RunSpec()`'s default `ctrl=1` makes
    final-year-only and multi-year-average numerically identical).
  - §8.5: `model.jl`'s flux-correction spin-up `if`/`if` pair (load, then
    unconditionally also compute-fresh) is now `if`/`elseif`/`else`,
    matching Fortran's mutually-exclusive dispatch.
  Added one new regression test per fix (`§8.1`/`§8.2`/`§8.3`/`§8.5` —
  `§8.4` reused/updated the 2 existing postprocess tests instead of adding
  new ones, since they already covered the exact functions being changed).
  One test-authoring mistake caught and fixed along the way: the first
  `§8.3` test draft used the wrong sign for `c_q` (produced a small
  positive `dq_rain` instead of a large negative one, so it never actually
  reached clamp territory) and separately failed to account for
  `init_model!` calling `set_hydrology_parameters!`, which clobbers a
  manually-set `c_q`/`c_rq`/`c_omega`/`c_omegastd` back to the
  `log_rain`-indexed preset — fixed by setting those fields *after*
  `init_model!` instead of before. Full suite re-run clean: 210/210 light +
  182/182 heavy.
- **2026-08-12 (Polyester.jl + CircularArrays.jl trials, §2.3/§2.4)**: an
  outside reviewer's package suggestion list named two "five-minute
  experiments" for open performance items — `Polyester.jl` to retry §2.4's
  reverted `Threads.@threads` row-threading with a lower-overhead scheduler,
  and `CircularArrays.jl` to remove the `vgather`-vs-`vmovup` gather pattern
  §2.3/§2.12 diagnosed as the reason `Float32`/further `@turbo` tuning
  underdeliver. First pass at this request answered structurally (read the
  code, reasoned about compatibility) without running anything — user
  correctly rejected that as not what was asked, and asked for real
  benchmark code instead. Built a new `benchmark/` directory (`Project.toml`
  with `CircularArrays`/`Polyester`/`LoopVectorization` as benchmark-only
  deps, kept separate from GREB's own `Project.toml` per the user's explicit
  choice) with two standalone scratch scripts, same "verbatim-logic clone,
  correctness-then-timing" discipline as every prior benchmark in this
  project:
  - `benchmark/circular_arrays_experiment.jl` — cloned `diffusion!`'s zonal
    stencil (`circulation.jl:73-94`) two ways. Found `CircularArrays.jl`
    circularizes every dimension of a `CircularArray` (no per-dimension
    option), so the only structurally valid use here is a `CircularVector`
    per row, not a whole-array wrap (which would incorrectly wrap the
    non-periodic latitude axis). Result: correct output, but `@turbo`
    silently fails `LoopVectorization.check_args` on the `CircularVector`
    version and falls back to a scalar (non-SIMD) loop — confirmed directly
    in `code_native` (baseline: `vgatherqpd`×12 + real packed-`pd` SIMD ops;
    circular: zero packed ops, 100% scalar `sd` instructions) — making it
    10–20× *slower* (16–30µs → 340–426µs), not faster. Rejected.
  - `benchmark/polyester_experiment.jl` — cloned `diffusion!`'s full
    outer-`k` loop (mid-latitude + polar branches), giving each `@batch`
    lane a private `(xdim, nlanes)` scratch matrix instead of the real
    code's single shared `ws.T1h`/`ws.dTxh` (which would race under
    concurrent `k`). Verified bit-identical output at every thread count
    tried, with 5 repeated trials per count to rule out an intermittent
    race (none found). Swept `-t 1` through `-t 14` (the machine's full
    logical-core count): found a real 1.79× speedup at `-t 4`, but a sharply
    non-monotonic curve — every count above 4 degrades, and `-t 14`
    collapses catastrophically (270–340× *slower*, reproduced twice). A real
    win exists but only at a specific, pinned thread count, not at `-t
    auto` — left undecided (not implemented) pending a thread-count policy
    decision, and to avoid landing a change in `circulation.jl` while the
    user's own performance work is in flight there (matching §2.11's
    "benchmarked, not yet wired in" precedent). Documented in `IMPROVEMENTS.md`
    §2.3/§2.4 as follow-ups; `benchmark/` and its scripts committed so this
    is the first real reusable benchmark harness in the repo (every prior
    benchmark was an uncommitted scratch script).
- **2026-08-12 (performance sweep: 4 more `@turbo` wins, 1 algorithmic fix,
  2 more negative results, 3 investigated-not-implemented)**: user asked for
  a wider sweep beyond §2.12's 4 already-landed spots, explicitly including
  its two never-benchmarked candidates. Ran 2 parallel Explore agents across
  all of `src/` (one for remaining per-timestep broadcast/allocation/type-
  instability patterns, one for architecture/I-O/control-flow patterns
  outside the inner kernels), each briefed not to re-propose anything
  already covered as done/reverted/deferred. Combined findings → benchmarked
  6 `@turbo` candidates + 1 algorithmic candidate with a single scratch
  script (`bench_sweep2.jl`), same methodology as §2.12 (baseline vs.
  explicit `@turbo` loop, correctness diff, only land real wins):
  - **4 wins, implemented**: `diagnostics!`'s accumulate (`output.jl`,
    57% faster, 22.65→9.75 µs), `state.jl`'s `accumulate!` (51%, 24.2→11.85
    µs), `SWradiation!`'s remaining `a_atmos`/albedo-combine/`sw`-flux-loop
    (38%, 5.5→3.4 µs), `qflux_correction!`'s per-timestep update (68%,
    57.0→18.0 µs). The `qflux_correction!` correctness check first showed a
    `NaN` — same root cause as an earlier session's gotcha: the benchmark's
    synthetic `cap_surf` was all-zero (real `cap_surf` only gets populated
    by `init_model!`), division by zero. Fixed by adding a nonzero offset in
    the benchmark; re-checked clean (`1e-16`-level diffs).
  - **2 negative results, not implemented**: `seaice!`'s trailing glacier-
    override line (`ocean.jl:46`, single op) was flat-to-slower across
    repeated runs (−2.5% to −47%, noise-dominated, never a real win);
    `tendencies.jl:39`'s `Q_sens` (single op) was consistently slower
    (−13% to −16%). Both are a 2nd/3rd instance of the already-documented
    "single op, too small for `@turbo`'s overhead to pay off" pattern
    (`circulation.jl:330`).
  - **1 algorithmic fix, implemented**: `forcing()`'s `:regional_co2_ocean`/
    `:regional_co2_land_ice` branches recomputed a *static* mask (depends
    only on `z_topo`/`icmn_ctrl`, both fixed for the whole scenario run)
    every scenario timestep instead of once. Not a `@turbo` candidate
    (branchy) — fixed by guarding with `it == 1`. Measured in isolation:
    7.65 µs baseline → 0.1 µs at `it>1` (98.7% faster), ~6ms/simulated year
    saved for these two experiments (the mask itself is cheap; the win is
    skipping it 729/730 times).
  - **3 investigated, not implemented** (documented as findings, §2.14):
    CMIP5/ENSO anomaly reload always re-reads from disk even when `fields`
    is reused across `greb_model!` calls (real for a repeated-run/ensemble
    use case, needs a cache-invalidation design decision — out of scope);
    `postprocess.jl`'s runtime-`Symbol` `getfield` loop (real but dwarfed by
    the array op it wraps at current record counts); `CirculationWorkspace()`'s
    ~40-allocation constructor (only called 3×/run, not per-timestep —
    negligible).
  - Full suite (319/319, incl. golden regression and the `:regional_co2_*`
    sweep test) verified passing after landing all 5 implemented changes.
    Attempted an aggregate before/after 1-year timing at `-t 3` but this
    batch's total contribution to a `:full_model`/`flux=0` run (~20ms/year,
    estimated from the per-call numbers × 730 steps — `qflux_correction!`'s
    and `forcing()`'s fixes don't even apply to that run at all) is smaller
    than this machine's own run-to-run wall-clock noise (four repeated
    timings: 1.15s/1.37s/1.41s/1.46s, a ~300ms spread) — unlike §2.11/§2.12's
    wins, which were large enough to show up cleanly end-to-end. Reported
    the isolated per-call microbenchmarks as the real signal instead of a
    noisy or cherry-picked aggregate number.
- **2026-08-12 (correction: the "-t 2, no gain past -t 4" timing below was a
  benchmark artifact)**: the user asked directly — "did you do the 3-way
  split and not the 2-way? would you expect it to need 3 threads?" —
  which was the right question. Re-examined the timing script
  (`time_1yr_after_impl.jl`) and found `deepcopy(fields)` (copying ~24 large
  3D climatology arrays, needed because `greb_model!` mutates `fields` in
  place) was sitting *inside* the `@elapsed` block as an argument
  expression, adding the same large thread-count-*independent* constant to
  every timing — which is exactly what hid the real 2-vs-3-thread
  distinction the architecture predicts. Fixed by moving both `deepcopy`
  calls (warm-up and timed) strictly outside their respective timed
  regions, then re-ran `-t 1`/`-t 2`/`-t 3`/`-t 4`:

  | Threads | Time | vs. `-t 1` |
  |---|---|---|
  | `-t 1` | 1.912s | 1× |
  | `-t 2` | 1.363s | 1.40× |
  | `-t 3` | **1.170s** | **1.63×** |
  | `-t 4` | 1.294s | no further gain (within noise) |

  `-t 3` *is* meaningfully faster than `-t 2` (confirming the two spawned
  `circulation!` tasks were queuing sequentially onto the single extra
  worker thread at `-t 2`, only truly running as 3 concurrent lanes at
  `-t 3`+), and `-t 4` adds nothing further — the correct signature of a
  genuine 3-way split, and the correct operational recommendation is
  `-t 3`, not `-t 2`. Total vs. the original ~2.7s/year baseline is
  **~2.31×** at `-t 3`, much closer to (and consistent with, modulo
  Amdahl's-law overhead outside `tendencies!`) the ~2.8× the isolated
  integration test below measured — the previous entry's "~1.78× total,
  lower than the isolated test's ~2.8×" framing is superseded by this one;
  the earlier entry is left as-is below per this document's append-only
  convention, not edited. Full suite (319/319, incl. golden regression)
  re-verified passing at `-t 3`.
- **2026-08-12 (§2.11 landed for real: 3-way thread split of `tendencies!`)**:
  implemented the benchmarked-not-implemented split from the entry below.
  `tendencies!` (`src/tendencies.jl`) gained `ws_a`/`ws_q` keyword args
  (default `=ws`, so every existing call site keeps its exact sequential
  behavior); when a caller passes distinct workspaces *and*
  `Threads.nthreads() > 1`, the two `circulation!` calls run via
  `Threads.@spawn` while `SW`/`LW`/`Q_sens`/`hydro!`/`deep_ocean!` run on the
  main `ws`. `time_loop!` (`src/output.jl`) and `qflux_correction!`
  (`src/model.jl`) forward the same kwargs; `greb_model!` is the only call
  site that actually opts in, allocating `ws_a`/`ws_q` once alongside `ws`.
  Full suite (319/319, including golden regression) passes unchanged at
  `-t 1` and at `-t 2` (parallel path actually exercised — confirmed via a
  direct `include("test/runtests.jl")` at `-t 2`, 1m42s vs. 2m49s at `-t 1`,
  itself a real-world echo of the win below). Real 1-year `greb_model!`
  timing (turbo fixes from the entry below already included in both):
  `-t 1` 2.212s → `-t 2` 1.513s (**1.46×** from threading alone on top of
  turbo's already-landed ~1.22×; **~1.78×** vs. the original ~2.7s/year
  baseline). `-t 4` gives no further gain (1.559s, within noise of `-t 2`) —
  expected, only 2 concurrent `circulation!` tasks exist. This real,
  full-pipeline number is lower than the ~2.8× from the isolated integration
  test below: that test measured `tendencies!`+`time_loop!`'s hot loop in
  isolation, while `greb_model!`'s ~2.7s/year baseline also includes
  `output!`/`diagnostics!`/JLD2-adjacent bookkeeping that this change doesn't
  touch — consistent with Amdahl's law given tendencies!'s ~70-75% share
  (§2.10). Measured, not re-estimated, so this number supersedes the
  integration test's for planning purposes.
- **2026-08-12 (§3.0 landed: test suite sharding)**: split
  `test/runtests.jl`'s 24 testsets into `run_light_tests()` (17) and
  `run_heavy_tests()` (7, incl. golden regression), matching the boundary
  already benchmarked in §3.0 below. Gated by `GREB_TEST_SHARD`
  (`"all"`/`"light"`/`"heavy"`, default `"all"`) so plain `Pkg.test()`/`]test`
  is unaffected. `.github/workflows/ci.yml` gained a `shard: [light, heavy]`
  matrix axis crossed with the existing `julia-version`/`os` axes, setting
  `GREB_TEST_SHARD` as a step env var. Verified locally: light shard 142
  tests, heavy shard 177 tests, 142+177=319 matches the unsharded total
  exactly; both shards pass independently.
- **2026-08-12 (§2.12 landed: 4 `@turbo` fixes)**: implemented all 4 rows
  from the table below (`diffusion!`'s mid-lat branch and `circulation!`'s
  substep accumulate in `src/circulation.jl`; `LWradiation!` in
  `src/physics/radiation.jl`; `time_loop!`'s state update in
  `src/output.jl`) as explicit `@turbo for j in 1:ydim; for i in 1:xdim`
  loops, matching each one's already-benchmarked prototype. Left
  `circulation!`'s final difference (line 330, the documented negative
  result) as plain `@.`. `time_loop!`'s rewrite additionally had to fold in
  two details the original prototype didn't cover: the `dq_crcl_use`
  buffer-selection branch (kept outside the loop, unchanged) and
  `log_hydro_dmc`'s force-zero of the humidity increment (replaced with a
  precomputed `0.0`/`1.0` multiplier inside the loop — same result, since
  the multiplied value is always finite). Full suite (319/319, incl. golden
  regression against the real NCEP dataset) passes unchanged. Real 1-year
  `greb_model!` timing: 2.7s → 2.212s, in line with the ~0.46s/year
  predicted below.
- **2026-08-12 (sixth Fortran audit — 3 parallel Explore agents, §8)**: a
  fresh bug sweep requested while performance work (§2.11/§2.12) was in
  progress elsewhere, run independently and left undone by explicit user
  decision (documentation only, no concurrent behavioral edits). Split the
  codebase across 3 parallel read-only Explore agents — (1)
  `circulation.jl`/`tendencies.jl`/`state.jl`, (2)
  `physics/{hydrology,ocean,radiation}.jl`, (3) `model.jl`/`output.jl`/
  `postprocess.jl`/`io.jl`/`config.jl`/`constants.jl` — each briefed on all
  18 previously-fixed bugs (this document's §0/changelog) so they wouldn't
  re-report known issues, and each instructed to verify every claim by
  reading both the Julia and Fortran text directly rather than from memory.
  Every finding the agents returned was then independently re-verified a
  second time by directly reading the cited Fortran lines myself before
  being kept — consistent with this codebase's established practice (see
  the fifth-pass entry below, "5 of 8 independently re-verified"). 5 new
  findings confirmed and documented in §8, none fixed:
  - `forcing()`'s regional-CO₂ ice mask reads only January's `icmn_ctrl`
    slice instead of the Fortran-mandated 12-month average
    (`greb.model.mscm.f90:1279-1283`).
  - `hydro!`'s `log_eva==1` branch is missing a `+2.0²`/`+3.0²` gust
    carryover that Fortran's shared `abswind` variable picks up before the
    `log_eva` dispatch even starts (`greb.model.mscm.f90:710-712,727-728`)
    — distinct from the already-investigated-and-kept `144.0` vs `144.**2`
    magnitude question in §0's "Investigated, not changed" note.
  - `hydro!` has a second, spurious `-0.9*q/Δt` clamp on `dq_rain` alone
    with no Fortran counterpart, duplicating/interfering with the already-
    correct combined-`dq` clamp in `output.jl` (the §0.18 fix) at an
    earlier, wrong point in the pipeline.
  - The control-run climatology baseline (`ctrl_output` → `build_monthly_
    climatology`/`compute_annual_ice_climatology`) averages every control
    year, where Fortran's `output` subroutine only ever retains the final
    control year (`it/ndt_days > 365*(time_ctrl-1)`,
    `greb.model.mscm.f90:1422-1446`) — silently correct only when
    `RunSpec.ctrl==1`, silently wrong for any multi-year control run (the
    norm in every real Fortran run script).
  - The flux-correction spin-up's "load precomputed files" branch
    (`model.jl:281-290`) is dead code: the unconditional `qflux_correction!`
    call right after it always overwrites whatever was just loaded, where
    Fortran's equivalent dispatch is a mutually exclusive
    load-or-compute-or-zero (`greb.model.mscm.f90:357-390`).
- **2026-08-12 (§2.11+§2.12 combined, integration-tested)**: built a full
  1-year control-run integration test combining the 3-way thread split with
  3 of the 4 turbo fixes, to check they don't interfere. Result: ~2.8×
  (2.83–2.88s → 1.00–1.04s), better than either alone; output bit-identical
  to both the golden reference and a fresh baseline run. First attempt
  showed a spurious 2.5× *slowdown* — traced to the new functions not being
  JIT-warmed the way `greb_model!`'s precompiled path is; fixed with an
  explicit warm-up call before timing. Not implemented for real yet.
- **2026-08-12 (does the `@turbo` fix generalize? — 3 more wins, 1 loss,
  §2.12)**: follow-up question after `LWradiation!`'s `@turbo` benchmark —
  does the same fix apply elsewhere? Checked every remaining plain-`@.`
  broadcast in the hot path, prioritizing by call frequency since
  `circulation!` sub-steps (24×2 = 48/main-timestep) make even small
  per-call costs compound fast. Found three more real wins:
  `diffusion!`'s mid-latitude branch (`circulation.jl:63-66`, ~46 rows ×
  48/step, 75% faster, ~0.22s/year), `circulation!`'s per-substep
  accumulate (`circulation.jl:326`, 48/step, 78% faster, ~0.12s/year), and
  `time_loop!`'s per-step state update (`output.jl:137-169`, 80% faster,
  ~0.04s/year) — ~0.46s/simulated year total alongside the already-found
  `LWradiation!` fix, verified numerically identical to baseline in all
  cases (diffs `1e-14`–`1e-18`). Hit one real snag along the way: the first
  correctness check for the state-update rewrite returned `NaN`, traced to
  the *test harness* — forgot to call `init_model!`/`seaice!` before
  benchmarking, so `fields.cap_surf` was still all-zero, dividing by zero;
  not a bug in the `@turbo` rewrite itself, confirmed by re-running with
  `cap_surf` properly initialized (clean `0.0`–`5.7e-14` diffs). Also found
  a useful negative result: `circulation!`'s final subtraction
  (`circulation.jl:330`), called only twice per step (not per-substep), was
  **17.6% slower** under `@turbo` — the macro's overhead isn't free below
  some size/frequency threshold, so this isn't a blanket "add `@turbo`
  everywhere" finding, it's a frequency/size-gated one. Documented as an
  expanded §2.12 (folding in the earlier `LWradiation!`-only version);
  flagged two more same-scale candidates (`state.jl:206-218`'s monthly
  accumulation, `SWradiation!`'s remaining combine/flux lines) as
  not-yet-benchmarked for a future pass. None of the four implemented yet.
- **2026-08-12 (Fortran build-flag comparison + `LWradiation!` fast-math gap,
  §2.10/§2.12)**: the Julia-vs-Fortran speed comparison (§2.10) had never
  actually checked whether the Fortran reference's build/run assumptions
  matched the single-threaded Julia baseline it was compared against — the
  user surfaced the real Fortran compile line
  (`gfortran -Ofast -ffast-math -funroll-loops -fopenmp greb.model.mscm.f90
  greb.shell.mscm.f90 -o greb.x`, found via a sub-agent search across the
  Fortran source tree) and asked whether we'd applied the same assumptions.
  Found: (1) `-fopenmp` isn't decorative — `greb.model.mscm.f90`'s
  `tendencies` subroutine (lines 512–535) wraps `SW`/`LW`+sensible/`hydro`/
  `circulation`(Ta)/`circulation`(q)/`deep_ocean` each in its own
  `!$omp section`, the *same six stages* §2.11 found and benchmarked
  independently in Julia — meaning the "Fortran runs a year in <1s" folklore
  already assumes multi-threaded execution, not evidence of an inherent
  single-thread speed advantage. (2) The same folder's `readme.md` states
  unoptimized/unparallelized Fortran takes ~30s/simulated year — *slower*
  than this Julia port's current ~2.7s/year single-threaded — flipping the
  original framing: naive Fortran is much slower than this port; Fortran's
  speed comes from optimization+parallelization this comparison hadn't
  matched. (3) Tested `julia -O3` vs. default `-O2` on the full 1-year
  timing: no measurable difference (2.13–2.22s, within noise) — `@turbo`'s
  LoopVectorization-generated code is largely independent of Julia's global
  `-O` level. (4) Checked the `-march=native` assumption directly rather
  than asserting it: `Base.JLOptions().cpu_target == "native"` and LLVM's
  JIT reports targeting this exact CPU (`alderlake`) by default — already
  matched, no flag needed. (5) `-ffast-math` has no direct Julia equivalent
  applied globally, but its *per-kernel* impact is real where a kernel isn't
  already `@turbo`-annotated — checked by benchmarking three variants of
  `LWradiation!` (the one hot kernel using plain `@.` broadcasting instead of
  `@turbo`, unlike every sibling kernel): baseline 145.7µs, `@fastmath`-
  wrapped 79.7µs (45% faster), rewritten as an explicit `@turbo` loop 36.3µs
  (75% faster), all numerically identical to baseline within float
  reassociation (`1e-13`–`1e-16`). Real and safe, but only ~3% of total step
  cost (`LWradiation!`'s share), so ~80ms/simulated-year saved — a good small
  follow-up alongside §2.11's ~1s/year 3-way split, not a substitute for it.
  Documented as §2.10's addendum and new §2.12; neither implemented yet.
- **2026-08-12 (assembly-level confirmation of the Float32 gather finding,
  §2.3)**: follow-up question on the same day's Float32 benchmark — asked
  whether `@code_llvm`/`@code_native` would confirm the gather-vs-contiguous
  hypothesis directly instead of just inferring it from timing. Isolated the
  gather-pattern zonal-diffusion loop (indexed through
  `lon_jm1`/`jp1`/`jm2`/`jp2`/`jm3`/`jp3`) vs. a plain contiguous elementwise
  loop as a control, and counted instruction mnemonics in `code_native`
  output (not `code_llvm` — LLVM IR is pre-instruction-selection, so the
  `vgather`/`vmovup` distinction only appears after instruction selection,
  in native asm). The contiguous control compiled to the *identical*
  instruction count in both types (34 `vmovup`/2 `vmovap`) — confirming the
  expected win where it should occur (same instructions, double the data per
  instruction). The real gather-pattern loop instead needed *more*
  instructions under `Float32` (24 `vgather`/44 total memory ops vs. 12
  `vgather`/26 for `Float64`) — gather throughput doesn't scale favorably
  with narrower elements on this hardware, directly explaining the measured
  ~5–10% (not 2×) real-world speedup rather than leaving it as an inference
  from timing alone. Added to §2.3 as confirming evidence.
- **2026-08-12 (Float32 kernel benchmark + test-suite process parallelism,
  §2.3/§3.0)**: two follow-ups requested after the `tendencies!`-stage
  parallelization benchmark. First, measured (not just analyzed) the actual
  compute payoff of `Float32` for the circulation kernels — built a
  standalone verbatim clone of `diffusion!`/`advection!`/`convergence!`/
  `circulation!` with every constant/climatology array converted to the
  tested eltype (avoids the mixed-precision trap where a stray `Float64`
  constant silently promotes everything back), and found only a 5–10%
  speedup, not the 2× a naive "half the width, double the throughput" SIMD
  argument would suggest — attributed to the permuted-neighbor-index gather
  pattern in the zonal terms and the grid already fitting in cache at either
  precision, plus the discovery that `CirculationWorkspace` is a concrete
  (non-parametric) `Float64` struct, so a real switch needs a
  `CirculationWorkspace{T}` refactor on top of §2.3's original
  `init_model!`-eltype hazard — concluded not worth it for this kernel shape
  given §2.11's 3-way split already gets 2.0× more cheaply. Confirmed
  directly in generated assembly, not just inferred from timing: isolated
  the gather-pattern zonal-diffusion loop vs. a plain contiguous control
  loop and counted `code_native` mnemonics. The contiguous control compiles
  to the identical instruction count in both types (34 `vmovup`/2 `vmovap`)
  — the expected win, same instructions moving 2× the data. The real
  gather-pattern loop instead needed more instructions under `Float32` (24
  `vgather`/44 total memory ops vs. 12 `vgather`/26 for `Float64`) — gather
  throughput doesn't scale favorably with narrower elements on this
  hardware. Second, tested process-level (not thread-level — `Test.jl`'s
  testset tracking isn't thread-safe) parallelization of the test suite:
  split `runtests.jl`'s 24 testsets into 2 standalone scripts along an
  existing boundary and ran them as concurrent `julia` processes, measuring
  1.56× (106.5s → 68.1s) despite an unbalanced split (one half takes ~2× the
  other) — a real, benchmarked opportunity for CI turnaround, not
  implemented as an actual pipeline change.
- **2026-08-12 (benchmarked parallelizing `tendencies!`'s stages, §2.11)**:
  followed up on the same day's profiling pass by actually testing
  `Threads.@spawn` on `tendencies!`'s independent stages. Confirmed via
  buffer audit (`state.jl:8-68`) that all six stages
  (`SW`/`LW`/`hydro`/`deep_ocean`/`circulation!(Ta)`/`circulation!(q)`) read
  only pre-timestep state and write disjoint buffers except the two
  `circulation!` calls, which need separate workspace instances. Naively
  spawning all 6 individually measured *worse* (1.56×) than just splitting
  the two `circulation!` calls (1.88×) — traced to bare `Threads.@spawn`+
  `wait` overhead (~19.4µs) exceeding `deep_ocean!`/`SWradiation!`/`hydro!`'s
  own cost (5.5–18.5µs each), so spawning them individually is a net loss.
  The right grouping — `circulation!(Ta)` and `circulation!(q)` each on their
  own spawned task, `SW`+`LW`+`hydro`+`deep_ocean` bundled sequentially into
  a third — hits 2.0×, verified bit-identical against sequential output for
  every stage. Documented as benchmarked-but-not-yet-implemented, same as
  the rest of §2.11.
- **2026-08-12 (per-timestep profiling vs. Fortran, §2.10/§2.11)**: timed a
  1-simulated-year `:full_model` control run post-precompile (~2.7s
  wall-clock) to compare against the Fortran reference's speed, then
  profiled where that time actually goes (`Profile.@profile` over 5000 warm
  `time_loop!` calls, plus manual `@elapsed` timing of each kernel).
  `circulation!` (`circulation.jl:299`) — run once for temperature, once for
  humidity, each cycling 24 diffusion+advection substeps — accounts for
  ~65–75% of every timestep; allocations (96 bytes/step) and the polar
  Jacobi sub-stepping (only 2 of 48 rows actually need >1 sub-step) were both
  ruled out as the cause, contrary to initial suspicion. Documented (§2.11,
  not implemented) four Julia-specific opportunities found in the process:
  parallelizing the two independent `circulation!` calls (coarser-grained
  than §2.4's already-reverted row-level threading attempt, so not assumed
  to fail the same way); ensemble/parameter-sweep parallelism via
  `Threads`/`Distributed`/GPU arrays reusing the same broadcast kernels
  unchanged (likely bigger practical payoff than single-run speed, given how
  this model actually gets used); low-cost `Float32`/GPU-array swaps via
  multiple dispatch; and AD-based parameter calibration via
  `Enzyme.jl`/`ForwardDiff.jl`.
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
- **2026-08-12 (this reorganization)**: split the single 1080-line
  `IMPROVEMENTS.md` into this file (pass-by-pass forensic narrative) and a
  rewritten, present-tense `IMPROVEMENTS.md` (current state of the code,
  organized by topic, §0 converted to a bug table). No information was
  dropped — every bug number, Fortran citation, and benchmark result from
  the original still appears in one of the two files.
