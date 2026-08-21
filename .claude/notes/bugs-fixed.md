# Bugs found and fixed

**Status:** Done — all entries below are fixed and covered by tests.

Two sweeps of bug discovery: the original pass (§0) and the 2026-08-12 sweep (§8).
The pass-by-pass forensic narrative lives in [`audit-history.md`](audit-history.md).

<!-- Split out of the former monolithic claude/IMPROVEMENTS.md on 2026-08-21.
     Section numbers in cross-references (§N.M) refer to that original file;
     see INDEX.md for where each section now lives. -->

---

## 0. Known bugs — fixed ✅

18 real, previously-shipped bugs found and fixed across five passes (a 19th
finding was investigated and reverted — see 0.14), each with a regression
test in `test/runtests.jl`. Passes 1–3 were found by grepping every
config-field write against its read sites; pass 4 was a direct line-by-line
comparison against the Fortran reference `greb.model.mscm.f90`, and every fix
there cites the exact Fortran line(s) confirmed by direct read (not recalled
from memory). Discovery narrative for each: `AUDIT_LOG.md`.

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
both Julia and Fortran. Full reasoning: `AUDIT_LOG.md`.

---

## 8. Bugs found in the 2026-08-12 sweep — ✅ all fixed

A fresh Fortran-vs-Julia audit (run in parallel with the in-progress
performance work in §2.11) found 5 real bugs, initially left undocumented-
but-unfixed to avoid colliding with concurrent performance edits, then
implemented in a follow-up pass once that settled. Every finding was
verified directly against `greb.model.mscm.f90` before being fixed, and
each has a regression test in `test/runtests.jl`. Full discovery method:
`AUDIT_LOG.md`.

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
