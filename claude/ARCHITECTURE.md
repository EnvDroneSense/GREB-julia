# GREB.jl — How the code works

**Purpose of this document:** a shared reference (human + Claude) for understanding
and improving the codebase. It describes `src/GREB.jl` as it *is* — including its
warts — not as it should be. Companion docs: [PACKAGING_WORKPLAN.md](PACKAGING_WORKPLAN.md)
(what to change), [IMPROVEMENTS.md](IMPROVEMENTS.md), [WORKPLAN.md](WORKPLAN.md)
(history of the notebook→package extraction).

*Last verified against `src/GREB.jl` on the `package-reorg` branch, 2026-07-17.*

---

## 1. TL;DR

GREB is a coarse (96×48, 3.75°) global energy-balance climate model with four
prognostic 2-D fields — **surface temperature `Ts`**, **atmosphere temperature
`Ta`**, **deep-ocean temperature `To`**, and **specific humidity `q`** — stepped
forward with explicit Euler at a 12-hour time step (730 steps/year). Everything
else (winds, clouds, soil moisture, mixed-layer depth, solar radiation) is
**prescribed climatology** read from binary "JDAL2" files into module-level
global arrays. A model run has three phases: an optional **flux-correction
spin-up** (computes nudging fields so the control climate matches observations),
a **control run**, and a **scenario run** (perturbed CO₂/solar/boundary
conditions), returning monthly-mean fields; scenario output is post-processed
into **anomalies relative to the control climatology**.

The code is a single 2 250-line module, mechanically extracted from a Pluto
notebook. Its two defining architectural traits:

1. **Global mutable state.** Climatology, static fields, flux corrections, and
   some parameters live as module-level globals (mostly non-`const`). Functions
   read them implicitly; `load_greb_jdal2!` and `init_model!` write them.
2. **Zero-allocation hot loop.** All per-step temporaries live in a
   pre-allocated `CirculationWorkspace`; kernels are `@turbo` (LoopVectorization)
   SIMD loops over the 96×48 grid.

---

## 2. Grid, calendar, and units

Defined at the top of the module (currently **non-`const`** globals — Phase 2 target):

| Constant | Value | Meaning |
|:--|:--|:--|
| `xdim × ydim` | 96 × 48 | lon × lat grid, 3.75° spacing; `j=1` is the North-Pole row, `j=48` the South-Pole row (lat from `lat_grid[k] = dlat*k − dlat/2 − 90` → −88.125…+88.125, so index 1 = −88° S… **note:** k=1 is *southernmost*) |
| `Δt` | 43 200 s (12 h) | main (physics) time step |
| `Δt_crcl` | 1 800 s (30 min) | transport sub-step; `ntime = 24` sub-steps per main step |
| `nstep_yr` | 730 | main steps per 365-day year (no leap years) |
| `jday_mon` | 31,28,… | month lengths; monthly output triggers on cumulative day boundaries |

All temperatures in **K**, humidity in **kg/kg**, fluxes in **W/m²**.
Precip/evap are converted for output to **mm/day** via `conv_factor = r_qviwv * 86400`.
Longitude is periodic; neighbour indices are precomputed (`lon_jm1 … lon_jp3`,
the only `const`s in the transport section). Latitude rows where the zonal grid
spacing `dxlat_grid[k] ≤ 250 km` are flagged `IS_POLAR` and get CFL
sub-sub-stepping in transport.

---

## 3. Data model: who owns what state

### 3.1 Prognostic state (owned by the caller)

`Ts, Ta, To, q :: Matrix{Float64}(96×48)` are created in `greb_model!` from
climatology initial conditions and threaded explicitly through
`time_loop! → tendencies! → kernels`. They are the *only* honestly-functional
part of the data flow.

### 3.2 Module-level globals (the implicit environment)

Everything below is a module global that kernels read directly (not passed as
arguments). Grouped by writer:

**Filled by `load_greb_jdal2!(dir; dataset=:ncep|:era)` from JDAL2 files:**

| Global | Shape | Content |
|:--|:--|:--|
| `z_topo`, `glacier` | 96×48 | topography [m] (<0 = ocean), glacier mask (>0.5) |
| `Tclim, uclim, vclim, qclim, swetclim` | 96×48×730 | surface temp, 850 hPa winds, humidity, soil-moisture climatologies |
| `cldclim, mldclim, Toclim, omegaclim, omegastdclim, wsclim` | 96×48×730 | clouds, mixed-layer depth, deep-ocean temp, vertical velocity (mean & std), wind speed |
| `sw_solar` | 48×730 | daily-mean incoming solar per latitude |
| `uclim_p/m, vclim_p/m` | 96×48×730 | wind sign-splits precomputed for upwind advection |
| `TF_correct, qF_correct, ToF_correct` | 96×48×730 | flux-correction fields (also *written* by `qflux_correction!`) |

**Derived by `init_model!(cfg)`:** `dTrad` (radiation-temperature offset
`−0.16·Tclim − 5`), `z_ocean` (3× max mixed-layer depth), `wz_air`/`wz_vapor`
(topographic pressure weights `exp(−z_topo/H)`), `cap_surf` (per-cell heat
capacity, land vs ocean·MLD), plus experiment-dependent *in-place edits of the
climatologies themselves* (see hazard H2).

**Anomaly fields** `*_anom_enso`, `*_anom_cc` (96×48×730): allocated as zeros,
**never loaded from any file in the package** — see bug B3.

**Runtime knobs:** `sw_solar_forcing_state::Ref{Float64}` (solar multiplier read
by `SWradiation!`, set per-step in the scenario loop), `co2_part` /
`co2_part_scn` (regional CO₂ masks, mutated by `forcing` for `regional_co2_*`
experiments).

**Dead/vestigial globals:** `WZ_CACHE`, `circ_workspace`, the standalone
monthly/annual accumulator blocks `Tmm…qsensmm`, `Tmn_ctrl…qcrclmn_ctrl`,
`co2_part_scn` — defined but unused by the model driver (the annual set
`Tsmn…fqmn` *is* used, by `diagnostics!`). Removal candidates.

### 3.3 Per-run scratch

- `CirculationWorkspace` (~40 pre-allocated 96×48 buffers + 2 vectors): every
  kernel writes its outputs into dedicated `*_buf` fields and returns views of
  them in NamedTuples. One instance per `greb_model!` call.
- `MonthlyAccumulator`: 13 running-sum fields + count; `reset!` after each month.

---

## 4. Configuration: `PhysicsConfig` and experiments

`PhysicsConfig` (mutable, `@kwdef`) is a bag of ~35 switches in four families,
following the MSCM (Monash Simple Climate Model) deconstruction philosophy:

- **`log_*_dmc`** — "deconstruct mean climate": turn a process off in the *base*
  climate (clouds, vapor, ice, circulation, hydrology, atmosphere, CO₂, ocean, q-flux).
- **`log_*_drsp`** — "deconstruct response": keep the process in the mean climate
  but freeze its *response* to forcing (topo, humidity, ocean, …).
- **Circulation flags** — `log_ice, log_hdif, log_hadv, log_vdif, log_vadv, log_conv`
  (ints/bools selecting which transport terms run).
- **Hydrology modes** — `log_rain ∈ {-1,0,1,2,3}` selects the precipitation
  parameterization (`HYDRO_PARAMS` lookup: original GREB, +RH, +ω, +RH&ω, "best
  GREB"), `log_eva ∈ {-1,0,…}` the evaporation formulation, `log_clim` the
  tuning dataset (0 = ERA-Interim, 1 = NCEP).

`create_experiment_config(sym)` returns presets (`:full_model`, `:co2_double` →
680 ppm, `:co2_quadruple`, `:solar_plus27`, `:elnino`, `:lanina`,
`:paleo_231kyr`, `:rcp85`, `:constant_topo`). Many more experiment symbols are
recognized downstream by `forcing()` (CO₂ ramps, sine/step, regional masks,
orbital) — the preset list and the `forcing()` list are **not in sync**.

`set_hydrology_parameters!(cfg)` maps `(log_rain, log_clim)` to precipitation
coefficients `(c_q, c_rq, c_omega, c_omegastd)` — but see **bug B1**: it writes
them to module globals while `hydro!` reads `cfg.c_q…`, so the lookup currently
has **no effect** on the physics.

---

## 5. Run orchestration

```
greb_model!(time_flux, time_ctrl, time_scnr, cfg; jdal2_dir="")
│
├─ 1. init_model!(cfg)              # derive fields, apply experiment edits,
│                                   # → (Ts_ini, Ta_ini, To_ini, q_ini, CO2_ctrl)
├─ 2. qflux_correction!(…time_flux) # spin-up: run the model nudged to climatology,
│                                   # storing per-step corrections TF/qF/ToF_correct
├─ 3. control run                   # time_ctrl years at fixed CO2_ctrl
│      for it in 1:time_ctrl*730:  time_loop!(…, ctrl_output, …)
├─ 4. scenario run                  # reset state to initial conditions, then
│      for it in 1:time_scnr*730:
│        forcing(it, year, cfg, ice_forcing)   # → per-step CO2, solar multiplier
│        (forced-boundary exps overwrite Ts from Tclim)
│        time_loop!(…, scnr_output, …)
└─ 5. post-process: scnr_output .-= monthly climatology of ctrl_output
       (skipped for orbital experiments)
→ returns (ctrl = Vector{MonthlyRecord}, scnr = Vector{MonthlyRecord})
```

Notes:

- **Flux correction** (`qflux_correction!`) is the GREB trademark: at every step
  it computes what the *uncorrected* model would predict, sets the correction to
  exactly the gap to climatology (`TF = (Tclim − Ts_pred)·cap_surf/Δt`, and
  analogous additive corrections for `q`, `To`), then applies it. After
  `time_flux` years the arrays `TF/qF/ToF_correct(·,·,730)` hold a
  seasonally-resolved nudging field used in phases 3–4. The spin-up *mutates*
  `Ts_ini…q_ini` in place, so control/scenario start from the spun-up state.
  (README flags this module as "being debugged" — treat with suspicion.)
- **Control run** output feeds `compute_annual_ice_climatology` → passed to
  `forcing` (used only by `regional_co2_ocean/land_ice` masks).
- **Scenario output is anomalies**, not absolute fields, for all non-orbital
  experiments — a frequent gotcha when plotting.
- `time_loop!` per step: calendar lookup (precomputed table for 200 yr) →
  `tendencies!` → Euler updates → clamps (`Ts,Ta ≥ 233.15 K`; `Δq` clamped to
  `[−0.9q, +0.02]`) → `seaice!` (updates `cap_surf` for next step) → unit
  conversions → `output!` (monthly means) + `diagnostics!` (annual means,
  printed to stdout).

---

## 6. The physics step: `tendencies!`

Called once per 12-h step; each kernel writes into workspace buffers and the
results are combined in `time_loop!`'s Euler update:

```
Ts += dT_ocean + Δt·(SW + LW_surf − LW_down + Q_lat + Q_sens + TF_corr)/cap_surf
Ta += dTa_crcl + Δt·(LW_up + LW_down − em·LW_surf + Q_lat_air − Q_sens)/cap_air
To += dTo + ToF_corr
q  += clamp(Δt·(dq_eva + dq_rain) + dq_crcl + qF_corr, …)
```

| Kernel | What it does |
|:--|:--|
| `SWradiation!` | Ice-cover fraction from `Ts` (linear ramp between land/ocean freeze thresholds, glacier override); surface albedo (0.1 ice-free, +0.25 ice); atmospheric albedo `= 0.35·cloud`; combined albedo; `SW = S·multiplier·(1−albedo)` with `S = sw_solar[lat, step]`. Returns `(SW, albedo, ice_cover)`. |
| `LWradiation!` | Greenhouse via a 10-parameter (`p_emi`) log-regression emissivity in effective CO₂ & vapor columns (topography-scaled by `wz_air`, regional mask `co2_part`), cloud-adjusted. `LW_surf = −σTs⁴`; `LW_down = LW_up = −em·σ(Ta + dTrad)⁴`. |
| sensible heat | `Q_sens = 22.5·(Ta − Ts)` (inline in `tendencies!`). |
| `hydro!` | Saturation humidity via Magnus-type formula on `Ts` (or a skin temperature for `log_eva==0`); evaporation `Q_lat ∝ (q − qs)·wind·swet` with land/ocean gustiness; precipitation `dq_rain = (c_q + c_rq·RH + c_ω·ω + c_ωstd·ωstd)·cq_rain·q`; converts to vapor tendencies and atmospheric latent heating. Early-outs to zeros if hydrology/atmosphere switched off. |
| `circulation!` ×2 | Transport of `Ta` (scale height `z_air`) and `q` (`z_vapor`): 24 sub-steps of `diffusion!` + `advection!` (+ `convergence!` = `−q·ω·const` for vapor, Stassen eq. 18), accumulating `X_work`; returns total change `dX`. |
| `deep_ocean!` | Mixed-layer ↔ deep-ocean exchange: entrainment/detrainment from seasonal MLD change plus turbulent mixing (`co_turb = 5 W/K/m²`), active only for ice-free ocean cells. Returns `(dT_ocean, dTo)`. |
| `seaice!` | Not a tendency: rewrites `cap_surf` in place from current ice fraction (frozen ocean cells get land-like heat capacity). |

**Transport numerics** (`diffusion!`/`advection!`): 2-D stencils with
topographic weights `wz` at every access; zonal terms use a 10-4-1-weighted
wide stencil over ±3 periodic neighbours; meridional terms use asymmetric
upwind stencils with special cases for the 4 rows nearest each pole; advection
is upwind via the precomputed sign-split wind fields. Polar rows re-integrate
with a locally reduced CFL-stable sub-step and a positivity clamp
(`dq ≥ −0.9·X`). Everything is `@turbo`.

---

## 7. Output pipeline

Every step, `output!` adds the current fields into `MonthlyAccumulator`; at each
month boundary it pushes a `MonthlyRecord` — a `NamedTuple` of 13 **copied**
96×48 matrices (`Ts, Ta, To, q, albedo, ice, precip, evap, qcrcl, sw, lw, qlat,
qsens`) — onto the output vector. So a run returns `12 × years` records per
phase; memory ≈ 0.48 MB/month. `diagnostics!` separately accumulates annual
means and **prints** a `year / global-mean / tropical-Pacific / North-Europe`
line each model year (stdout is the only progress reporting).

---

## 8. Performance architecture

- **No allocations in the hot loop**: workspace buffers everywhere; kernels
  return NamedTuples of buffer references. Consequence: kernel outputs are
  *aliased* — never call a kernel twice and keep the first result.
- **`@turbo`** (LoopVectorization) on all grid loops, with branch-free `ifelse`
  chains instead of `if`.
- **Precomputation**: periodic neighbour index tables, wind sign-splits,
  diffusion/advection coefficients per latitude, `calendar_lookup`
  (146 000-entry table built at module load), inverse ice-threshold ranges,
  `dTrad`.
- **Known perf hazard**: nearly all of the globals the kernels read
  (`z_topo`, `wz_air`, climatologies, grid constants…) are **non-`const`**,
  so every access is a type-unstable global load. `@turbo` loops mostly touch
  them through local views/bindings which limits the damage, but plain
  broadcasts and scalar reads pay full price. Fixing this is Phase 2 of the
  packaging workplan and should be benchmarked with `benchmark/benchmarks.jl`.

---

## 9. Known bugs, hazards, and smells (improvement backlog)

Verified against the source on 2026-07-17. **B = likely bug, H = hazard/footgun, S = smell.**

- **B1 — Hydrology parameters never reach the physics.**
  `set_hydrology_parameters!(cfg)` assigns `c_q, c_rq, c_omega, c_omegastd` to
  *module globals* (`global c_q, …`, [src/GREB.jl:273](../src/GREB.jl#L273)), but `hydro!` reads
  `cfg.c_q…` ([src/GREB.jl:1095](../src/GREB.jl#L1095)), which keep their `PhysicsConfig` defaults
  `(1.0, 0, 0, 0)` = the plain "original GREB" scheme. Net effect: `log_rain`
  and `log_clim` selections are silently ignored. Fix: write into `cfg` fields
  (and drop the globals). May be an extraction artifact or inherited from the
  notebook — check against the notebook and the Fortran reference before fixing;
  changes results.
- **B2 — `log_eva == 0` path crashes.** `hydro!` writes `ws.cE[i,j]`
  ([src/GREB.jl:1159](../src/GREB.jl#L1159)) but the workspace field is named `cE_buf`; any run with
  `cfg.log_eva = 0` throws. Untested code path (default is −1).
- **B3 — ENSO / climate-change anomaly fields are never loaded.**
  `Tclim_anom_enso`, `*_anom_cc`, etc. are allocated as zeros and no loader
  fills them, so `:elnino`, `:lanina`, and `:rcp85`'s boundary forcing add zero
  anomalies — these experiments currently degenerate to (roughly) control runs.
  The original data source for these fields needs to be identified and a loader
  added.
- **H1 — `qflux_correction!` correctness is doubted** by the original author
  (README "Known Issues"). Anything downstream of `TF/qF/ToF_correct` inherits
  that doubt.
- **H2 — `init_model!` destructively edits the loaded climatology.**
  `Tclim .+= anomalies` (rcp85/ENSO), `z_topo = min(z_topo, 1)`
  (`constant_topo`), `qclim .= 0`, `cldclim .= 0`, `mldclim .= d_ocean`, and
  flux-correction zeroing all mutate the shared globals **in place, with no
  restore**. Running a second experiment in the same session without re-calling
  `load_greb_jdal2!` operates on corrupted inputs (and `:rcp85` twice would
  double-apply anomalies). This is the single biggest obstacle to a clean API.
- **H3 — Scenario output is anomalies, control output is absolute.** Easy to
  misread downstream; consider returning both or tagging the record.
- **H4 — Kernel outputs alias workspace buffers** (§8); storing them without
  `copy` gives silently time-varying "results". `output!` copies; user code must too.
- **S1 — Non-`const` module globals** throughout (perf + safety; Phase 2).
- **S2 — Dead code**: `WZ_CACHE`, `circ_workspace` global instance, the
  standalone `Tmm…`/`Tmn_ctrl…` accumulator blocks, `co2_part_scn`; also
  `forcing`'s `icmn_ctrl` zeros-default allocates when omitted.
- **S3 — `println`/`@info` sprinkled through the model core** (`diagnostics!`,
  loaders, `greb_model!`); no logging control, awkward for library use.
- **S4 — `forcing()`'s experiment list ≫ `create_experiment_config`'s presets**,
  several experiments error("not yet implemented"), and `:sst_plus1`/`:a1b_*`
  exist only in the driver; the experiment registry wants unifying.
- **S5 — Misleading names**: `min_humidity_change = 0.9` is actually "max
  removable fraction of q per step"; `uclim_m` holds the *positive* wind part
  (comment admits the swap); `lat_grid` runs south→north while the README calls
  row 1 "north".
- **S6 — Latitude convention is implicit** (see §2) and only documented in
  scattered comments; sample-point indices in `diagnostics!` are magic numbers.

---

## 10. File map

| Path | Role |
|:--|:--|
| `src/GREB.jl` | everything (module, 2 253 lines) — split planned in Phase 2 |
| `test/runtests.jl` | smoke tests, no data needed (31 asserts) |
| `examples/run_greb.jl` | plain-Julia driver: load data → run → PNG maps |
| `benchmark/benchmarks.jl` | per-kernel + full-year BenchmarkTools suite (synthetic fields) |
| `notebooks/GREB_julia.jl` | original Pluto notebook (source of truth for the extraction) |
| `DATA_README.md` | JDAL2 format spec & dataset layout |

The JDAL2 binary format (`read_jdal2`): magic `"JDAL2"` + version byte, Int32
dimension-name table, Int32 shape, a type byte (Float32 only), then raw data —
column-major, reshaped to the header shape.
