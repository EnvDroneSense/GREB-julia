### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 1cf6d892-c827-48fd-8e77-bb1117d03c20
begin
	import Pkg
	Pkg.develop(path=joinpath(@__DIR__, ".."))
	using GREB
	using Plots
	using ProgressLogging
end

# ╔═╡ 1165eeb4-10d8-4fb5-859d-7fe410189608
begin
	using Statistics
	using PlutoUI
	using NCDatasets
	using StaticArrays  # For optimization: static longitude indices
	using LoopVectorization  # For SIMD vectorization with @turbo
	using JLD2
end

# ╔═╡ 3f63ce64-effd-49c7-9906-92cda40d59f0
md"""
# The Globally Resolved Energy Balance (GREB) Model 

Authors: Dietmar Dommenget, Janine Flöter, Tobias Bayr and Christian Stassen

Translated + optimized for Julia: Thomas Struys

**References:**

- Dommenget, D., and Flöter, J. (2011). Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model. *Journal of Climate Dynamics*, 37: 2143. doi:[10.1007/s00382-011-1026-0](https://doi.org/10.1007/s00382-011-1026-0)

- Stassen, C., Dommenget, D., and Loveday, N. (2019). A hydrological cycle model for the Globally Resolved Energy Balance (GREB) model v1.0. *Geoscientific Model Development*, 12, 425-440. doi:[10.5194/gmd-12-425-2019](https://doi.org/10.5194/gmd-12-425-2019)

- Dommenget, D., Nice, K., Bayr, T., Kasang, D., Stassen, C., and Rezny, M. The Monash Simple Climate Model Experiments: An interactive database of the mean climate, climate change and scenarios simulations. *Geoscientific Model Development* (submitted)


input fields: The GREB model needs the following fields to be specified before

| Variable | Dimensions | Description |
|:---------|:-----------|:------------|
| `z_topo(xdim,ydim)` | 2D | topography (<0 are ocean points) [m] |
| `glacier(xdim,ydim)` | 2D | glacier mask (>0.5 are glacier points) |
| `Tclim(xdim,ydim,nstep_yr)` | 3D | mean Tsurf [K] |
| `uclim(xdim,ydim,nstep_yr)` | 3D | mean zonal wind speed [m/s] |
| `vclim(xdim,ydim,nstep_yr)` | 3D | mean meridional wind speed [m/s] |
| `qclim(xdim,ydim,nstep_yr)` | 3D | mean atmospheric humidity [kg/kg] |
| `mldclim(xdim,ydim,nstep_yr)` | 3D | mean ocean mixed layer depth [m] |
| `Toclim(xdim,ydim,nstep_yr)` | 3D | mean deep ocean temperature [K] |
| `swetclim(xdim,ydim,nstep_yr)` | 3D | soil wetness, fraction of total [0-1] |
| `sw_solar(ydim,nstep_yr)` | 2D | 24hrs mean solar radiation [W/m²] |
"""

# ╔═╡ 0780e492-cdcd-43fb-af85-f9a4bab37450
PlutoUI.TableOfContents(title="📋 Table of Contents", depth=2)

# ╔═╡ d0ebd34f-a83a-4dc1-9c66-24a3878b681d
md"""
## 💾 Input/Output

### Data Loading Module

Functions to read the `.jld2` input files produced by `scripts/convert_greb_to_jld2.jl`
(standard JLD2 container: each file stores plain Julia values under the keys
`"data"`, `"dim_names"`, and optionally `"coords"`/`"ctl"`).
"""

# ╔═╡ ebca5877-1305-40f1-ac52-66356dc17661
md"""
### Solar Forcing Storage

Global variable for orbital/paleoclimate solar forcing. Loaded on-demand using grouped JLD2 files.

**Usage:**
```julia
# Load eccentricity forcing (index 0-60)
global sw_solar_forcing_data = load_solar_forcing_jld2(jld2_dir, :eccentricity, 36)

# Load obliquity forcing
global sw_solar_forcing_data = load_solar_forcing_jld2(jld2_dir, :obliquity, 230)

# Load paleoclimate forcing
global sw_solar_forcing_data = load_solar_forcing_jld2(jld2_dir, :paleo)
```
"""

# ╔═╡ 8165ed3f-ede0-45c9-83ac-e4e49262457c
md"""
## 🕹️ Advanced Physics Control (MSCM)

**Process Isolation Capabilities**

The MSCM physics switches enable decomposition experiments by selectively enabling/disabling physical processes:

### Mean Climate Switches (`*_dmc`)
- `log_clouds_dmc`: Cloud feedback in mean climate
- `log_vapor_dmc`: Water vapor feedback in mean climate  
- `log_ice_dmc`: Ice-albedo feedback in mean climate
- `log_crcl_dmc`: Atmospheric circulation in mean climate
- `log_hydro_dmc`: Hydrological cycle in mean climate

### 2×CO₂ Response Switches (`*_drsp`) 
- `log_clouds_drsp`: Cloud feedback in CO₂ response
- `log_vapor_drsp`: Water vapor feedback in CO₂ response
- `log_crcl_drsp`: Atmospheric circulation in CO₂ response
- `log_hydro_drsp`: Hydrological cycle in CO₂ response  

**Implementation Patterns**:
1. **Early returns**: Functions exit early when switches are disabled
2. **Climatology control**: Mean climate switches modify initialization
3. **Tendency zeroing**: Computed changes are zeroed based on switches

This enables precise process attribution studies and decomposition experiments.
"""

# ╔═╡ 07a88c94-93bb-4564-8737-980144c4af43
md"""
## ⚙️ Configuration

### 📐 Grid Dimensions & Time-Stepping

The model uses a **3.75° × 3.75°** latitude-longitude grid with **96** longitude points × **48** latitude points.

Time-stepping uses a **12-hour** main time step (`Δt`) and a **30-minute** sub-time step for circulation (`Δt_crcl`), giving **730** time steps per year.

Calendar arrays track days per month for output functions.
"""

# ╔═╡ 813d59f7-eeea-402a-96e8-e5cbe2fd4582
md"""
## 𝚿 `mo_physics` - Physical Parameters & Climate Fields

This module defines all physical constants, model parameters, emissivity coefficients, and declares the climate field arrays.
"""

# ╔═╡ 54e63724-44ca-42b9-9246-7afb271b8154
md"""
### Hydrology Parameters
Optimized parameter system with lookup table for different parameterization modes and climatology datasets.
"""

# ╔═╡ caec8065-9874-4d8c-8e7e-598a97506f06
md"""
## 🔧 Workspace Structs

Pre-allocated structures to eliminate runtime allocations.
"""

# ╔═╡ 12a89062-bc27-4e36-8b5e-1670b4bd0b28
md"""
### Spatial CO₂ Fields
Spatial masking arrays for regional/partial CO₂ experiments.
"""

# ╔═╡ 1879fd40-0117-478b-acdd-18112738b81f
md"""
### Natural Constants
Fundamental physical constants used throughout the model.
"""

# ╔═╡ 0308940f-16bc-4b03-8321-9ea8c91a21c7
md"""
### Emissivity Regression Parameters
10-element vector used in the LW radiation subroutine to compute atmospheric emissivity from CO₂, water vapor, and cloud cover.
"""

# ╔═╡ 78605b5f-b03d-4851-b339-aea63bad3688
md"""
### Climate Field Arrays (Placeholders)

These arrays will be populated with real data once available. For now we initialize them with zeros so the notebook is runnable.

**2D fields** `(xdim, ydim)` = (96, 48):
- `z_topo` — topography [m] (< 0 = ocean)
- `glacier` — glacier mask (> 0.5 = glaciated)
- `z_ocean` — derived ocean depth (3× max MLD)
- `cap_surf` — surface heat capacity [J/K/m²]

**3D fields** `(xdim, ydim, nstep_yr)` = (96, 48, 730):
- `Tclim` — climatological surface temperature [K]
- `uclim`, `vclim` — zonal/meridional wind [m/s]
- `qclim` — atmospheric humidity [kg/kg]
- `mldclim` — mixed-layer depth [m]
- `Toclim` — deep ocean temperature [K]
- `cldclim` — cloud cover fraction
- `swetclim` — soil wetness [0–1]
- `omegaclim` — vertical velocity [Pa/s] (MSCM)
- `omegastdclim` — omega standard deviation [Pa/s] (MSCM)
- `wsclim` — wind speed climatology [m/s] (MSCM)
- `uclim_m/p`, `vclim_m/p` — negative/positive wind components (MSCM)
- `TF_correct`, `qF_correct`, `ToF_correct` — flux corrections
- `dTrad` — Tatmos–radiation offset

**2D solar** `(ydim, nstep_yr)` = (48, 730):
- `sw_solar` — 24-hr mean solar radiation [W/m²]
"""

# ╔═╡ 22d5e751-f095-42d7-9c24-78e546b3ce37
md"""
### Time-Step Indices
`jday` and `ityr` track the current calendar day and time-step within the year during integration. In Fortran these were module-level integers; here they start at 1.
"""

# ╔═╡ d403b01e-7d9f-4778-813d-bbbb0bfdbfb9
md"""
## ⚙️ `mo_diagnostics` - Diagnostic & Output Fields

This module declares the accumulator arrays for annual-mean diagnostics and monthly-mean output.
"""

# ╔═╡ c94558a5-0ae2-4563-93d7-0609ea3cce41
md"""
### Annual-Mean Accumulators
These fields accumulate values over a year and are divided by `nstep_yr` at year-end in the `diagnostics` subroutine, then reset to zero.

| Variable | Description |
|:---------|:------------|
| `Tsmn` | Surface temperature |
| `Tamn` | Air temperature |
| `Tomn` | Deep ocean temperature |
| `qmn` | Atmospheric humidity |
| `amn` | Albedo |
| `icmn` | Ice-cover fraction (MSCM) |
| `prmn` | Precipitation tendency (MSCM) |
| `evmn` | Evaporation tendency (MSCM) |
| `qcrclmn` | Humidity circulation tendency (MSCM) |
| `swmn` | Shortwave radiation |
| `lwmn` | Longwave radiation |
| `qlatmn` | Latent heat flux |
| `qsensmn` | Sensible heat flux |
| `ftmn` | Temperature flux correction |
| `fqmn` | Humidity flux correction |
"""

# ╔═╡ 2d02d36f-efee-411e-befb-52d61b3fb16e
md"""
### Monthly-Mean Output Buffers
Accumulated within each calendar month, divided by month length at month-end.

| Variable | Description |
|:---------|:------------|
| `Tmm` | Surface temperature |
| `Tamm` | Air temperature |
| `Tomm` | Deep ocean temperature |
| `qmm` | Atmospheric humidity |
| `apmm` | Albedo |
| `icmm` | Ice-cover fraction (MSCM) |
| `prmm` | Precipitation tendency (MSCM) |
| `evmm` | Evaporation tendency (MSCM) |
| `qcrclmm` | Humidity circulation tendency (MSCM) |
"""

# ╔═╡ 1b38b24d-bb82-4d52-a26d-850d142e25d5
md"""
## 🎮 `greb_model` - Main Driver

The driver orchestrates four phases:

1. **Initialization** - Derived fields and sensitivity overrides
2. **Flux correction** - Spin-up to compute correction fields
3. **Control run** - Integration at fixed CO₂
4. **Scenario run** - Integration with time-varying CO₂

Output is stored as in-memory `Vector` of `MonthlyRecord` named tuples.
"""

# ╔═╡ 0ff4be7c-580c-4db1-a874-175badc8c11a
md"""
### Monthly Output Record Type
Each monthly snapshot stores nine 2-D fields as a `NamedTuple`: 
`Ts`, `Ta`, `To`, `q`, `albedo`, `ice`, `precip`, `evap`, `qcrcl`.
This replaces the Fortran direct-access binary writes to units 21/22.
"""

# ╔═╡ 3beff4df-89da-4f25-a0e4-be03048c5f2c
md"""
## `init_model!`

Initializes the model state from climatology and experiment configuration:

- Sets hydrology parameters (`c_q`, `c_rq`, `c_omega`, `c_omegastd`) via `set_hydrology_parameters!`
- Computes derived fields: `dTrad` (radiation offset), `z_ocean` (3× max mixed-layer depth), `wz_air`/`wz_vapor` (pressure weights)
- Applies experiment overrides based on `cfg` switches: zeroes clouds, fixes topography, adjusts humidity/ocean climatology
- Handles special experiments (El Niño, La Niña, RCP8.5) by adding/subtracting anomaly fields
- Sets `cap_surf` per grid cell (land vs. ocean heat capacity, with ice-aware ocean capacity)
- Extracts initial conditions from the last time step of climatology: `Ts`, `Ta` (= `Ts`), `To`, `q`
- Determines control CO₂ level based on experiment type (default 340 ppm, IPCC scenarios 280 ppm, or 0 for deconstruction)

Returns `(Ts_ini, Ta_ini, To_ini, q_ini, CO2_ctrl)`.
"""

# ╔═╡ fd193b13-336f-42f8-8f93-3ca3bea2e620
md"""
## 🔄 `time_loop` - Single Time-Step Integration

Advances the model by one main time step (12 hours):

1. Looks up calendar day and time-step index from iteration counter
2. Computes all physics tendencies via `tendencies!`
3. Updates `Ts` (surface) and `Ta` (air) from energy budgets, with clamps at 233.15 K
4. Updates `To` (deep ocean) with mixing tendencies and flux correction
5. Updates `q` (humidity) with evaporation, rain, circulation, and flux correction (clamped)
6. Adjusts ocean heat capacity via `seaice!` based on current ice cover
7. Converts moisture tendencies to mm/day for output
8. Calls `output!` (monthly accumulation) and `diagnostics!` (annual accumulation)

Returns updated `(mon, irec)`.
"""

# ╔═╡ 241496f6-580b-4e81-acbb-3295696fe40e
md"""
## ⚔️ `tendencies` - Physics Aggregation

Calls all physical process modules in sequence for one time step:

| Step | Process | Key Output |
|:-----|:--------|:-----------|
| 1 | Shortwave radiation | `albedo`, `SW`, `ice_cover` |
| 2 | Longwave radiation | `LW_surf`, `LW_down`, `LW_up`, `em` |
| 3 | Sensible heat flux | `Q_sens` (zero if atmosphere disabled) |
| 4 | Hydrological cycle | `Q_lat`, `dq_eva`, `dq_rain` |
| 5 | Circulation (heat) | `dTa_crcl` |
| 6 | Circulation (moisture) | `dq_crcl` |
| 7 | Deep ocean coupling | `dT_ocean`, `dTo` |

Processes respect MSCM isolation switches - disabled processes return zero tendencies.

Returns a NamedTuple with all tendencies and diagnostics.
"""

# ╔═╡ a4ab662c-5876-4861-aecf-2f695d15e15d
md"""
## ☀️ `SWradiation` - Shortwave Radiation

Computes ice cover, surface albedo, and net shortwave flux.

- **Ice fraction**: linear ramp between ice-free and fully frozen thresholds, separate for land (`Tl_ice1` → `Tl_ice2`) and ocean (`To_ice1` → `To_ice2`); glaciers always fully iced
- **Atmospheric albedo**: `a_cloud × cloud_fraction`
- **Surface albedo**: ramps from `a_no_ice` to `a_no_ice + da_ice` based on ice fraction; fixed to `a_no_ice` if ice feedback disabled (`log_ice = false`)
- **Combined albedo**: `α = α_surf + α_atmos − α_surf × α_atmos`
- **Net SW**: `SW = solar_flux × (1 − α)` with optional solar multiplier

Returns `(SW, albedo, ice_cover)`.
"""

# ╔═╡ 2a859d43-5768-4b98-88fa-978b12203513
md"""
## 🌙 `LWradiation` - Long-Wave Radiation

Computes atmospheric emissivity and longwave fluxes.

- **Effective columns**: CO₂ and water vapor scaled by pressure weight `wz_air`; CO₂ also scaled by spatial mask `co2_part`
- **Emissivity**: log-regression with 10 coefficients combining CO₂, vapor, and cloud effects
- **Fluxes**:  
  `LW_surf = −σ·Ts⁴`  
  `LW_down = −ε·σ·(Ta + dTrad)⁴`  
  `LW_up = LW_down`
- Zeroed if atmosphere disabled (`log_atmos_dmc = false`)

Returns `(LW_surf, LW_up, LW_down, em)`.
"""

# ╔═╡ 0febe534-237e-4921-b39c-3828dbae9d19
md"""
## 💧 `hydro` - Hydrological Cycle

Computes latent heat flux and water vapor tendencies (evaporation, precipitation).

- **Saturation humidity**: Clausius-Clapeyron approximation with pressure scaling `wz_air`
- **Evaporation**: bulk aerodynamic formula with gustiness; three modes selectable via `log_eva`:
  - `-1` - original GREB (wind from `uclim`, `vclim`)
  - `0`  - skin temperature method with wind speed climatology
  - `1`  - enhanced
- **Precipitation**: `dq_rain = (c_q + c_rq·RH + c_omega·ω + c_omegastd·σ_ω) · cq_rain · q`
- **Vapor tendencies**: `dq_eva = −Q_lat / cq_latent / r_qviwv`, rain limited to `−0.9·q/Δt`
- Returns zero if hydrology or atmosphere disabled

Returns `(Q_lat, Q_lat_air, dq_eva, dq_rain)`.
"""

# ╔═╡ a1c04c52-f7c9-430a-8791-7fea15650b2c
md"""
## 🌫️ `convergence` - Moisture Flux Convergence

Computes moisture tendency from vertical velocity (omega) convergence:

`dX_conv = −T1 · ω · (Δt_crcl / z_vapor · 2.5 / (ρ_air·g))`

Implements Eq. 18 from Stassen et al. (2019). Called during circulation sub-stepping for water vapor only.
"""

# ╔═╡ 34fb479a-9834-4354-8f36-2680bedec798
md"""
## ❄️ `seaice` - Sea-Ice Heat Capacity

Updates `cap_surf` in-place for ocean points based on ice cover fraction.

- Ocean heat capacity blends linearly between `cap_land` (fully frozen, `Ts ≤ To_ice1`) and `cap_ocean × mld` (ice-free, `Ts ≥ To_ice2`)
- Land and glacier points keep `cap_land`
- Skipped entirely if ocean disabled (`log_ocean_dmc = false`)
- If ice feedback disabled (`log_ice = false`): land → `cap_land`, ocean → `cap_ocean × mld`
"""

# ╔═╡ ae38f814-fe9c-443f-ae3b-42fa7a7d199a
md"""
## 🌊 `deep_ocean` - Deep Ocean Coupling

Computes heat exchange between mixed layer and deep ocean.

- **Entrainment** (`mld` shallowing): deep ocean warms, surface cools — scaled by `c_effmix`
- **Detrainment** (`mld` deepening): surface cools, deep ocean warms
- **Turbulent mixing**: proportional to `co_turb × (max(Ts, To_ice2) − To)`
- Active only for ocean points with `Ts ≥ To_ice2`; returns zeros if ocean disabled

Returns `(dT_ocean, dTo)`.
"""

# ╔═╡ a4cd8d40-eedd-4a4b-ac17-51e68334328c
md"""
## 🌀 `circulation` - Atmospheric Circulation

Runs sub‑stepped diffusion + advection at the finer time step `Δt_crcl`.

- Loops `ntime` sub‑steps, each calling selected processes: horizontal/vertical diffusion, horizontal/vertical advection, moisture convergence
- Process selection controlled by switches: `log_hdif`, `log_vdif`, `log_hadv`, `log_vadv`, `log_conv`
- Early return with zeros if atmosphere or circulation disabled
- Returns `dX_crcl = X_final − X_in`
"""

# ╔═╡ cc7db228-0a06-4812-8b32-6541988b2115
md"""
## 🌀 `diffusion` - 3rd-Order Diffusion

Topography‑weighted horizontal diffusion on a periodic lat‑lon grid.

- **Meridional**: 2nd‑order centred differences; zero‑flux at poles
- **Zonal (mid‑latitudes)**: 3rd‑order stencil with ±1, ±2, ±3 neighbours (weights 10/4/1, normalised by 20)
- **Zonal (polar, `dxlat ≤ 250 km`)**: sub‑time‑stepping for stability, with clamp `dX ≥ −0.9·X`
- Scaled by pressure weight `exp(−z_topo / h_scl)`

Writes directly to `ws.dX_diff`.
"""

# ╔═╡ 99f0e9ad-438f-4fa5-be9d-36d0fa78d89c
md"""
## ➡️ `advection` - Upwind Advection

Upwind transport driven by climatological wind fields.

- **Meridional**: 2nd‑order upwind with ±1/±2 neighbours; one‑sided at pole boundaries
- **Zonal (mid‑latitudes)**: 2nd‑order upwind (±1/±2)
- **Zonal (polar)**: sub‑time‑stepping with 3rd‑order stencil (weights 10/4/1, normalised by 20), clamp `dX ≥ −0.9·X`
- Returns zeros if advection disabled for the transported quantity (`log_hadv` for heat, `log_vadv` for vapor)

Writes directly to `ws.dX_adv`.
"""

# ╔═╡ 8fdcfefd-0490-434a-a2bb-d171557b6ae7
md"""
## 🔄 `qflux_correction` - Flux Correction Spin-Up

Runs a `time_flux`-year simulation to compute correction fields that keep the model close to observed climatology.

Each time step:
1. Computes all physics tendencies
2. Calculates uncorrected state for `Ts`, `Ta`, `To`, `q`
3. Derives corrections: `TF_correct = (Tclim − Ts_uncorrected) × cap_surf / Δt` (and similarly for ocean and humidity)
4. Applies corrections to get corrected state
5. Updates sea ice and accumulates diagnostics

Corrections are stored in `TF_correct`, `ToF_correct`, `qF_correct` (3D arrays: lon × lat × time‑step‑in‑year) and applied during the main simulation to nudge the model toward climatology.
"""

# ╔═╡ 7e514484-9d81-4c5c-83ab-191d68c13043
md"""
---
## 🧪 `forcing` — Experiment Forcing

Returns `(CO2, sw_solar_forcing)` for the current year/time step based on experiment type.

| Category | Experiments | Description |
|:---------|:------------|:------------|
| **CO₂ scaling** | `:co2_double`, `:co2_quadruple`, `:co2_10x`, `:co2_half`, `:co2_zero` | Fixed CO₂ multipliers |
| **Solar** | `:solar_plus27`, `:solar_cycle_11yr` | Solar constant +27 W/m² or 11-year cycle |
| **Time‑varying CO₂** | `:co2_sine_wave`, `:co2_step`, `:a1b_scenario`, `:a1b_enhanced` | CO₂ as function of year |
| **Paleoclimate** | `:paleo_231kyr`, `:paleo_solar_modern_co2`, `:modern_solar_paleo_co2` | Paleo solar + CO₂ combinations |
| **Orbital** | `:obliquity`, `:eccentricity`, `:earth_sun_distance` | Solar forcing loaded externally |
| **Regional CO₂** | `:regional_co2_nh/sh/tropics/extratropics/ocean/land_ice/winter/summer` | Spatial masks via `co2_part` |
| **Forced boundary** | `:elnino`, `:lanina`, `:rcp85` | CO₂ = 340 ppm; anomalies applied in `init_model!` |
"""

# ╔═╡ 7b8e4130-7d8e-492a-a949-bd7fb6808dd6
md"""
---
## 📊 `diagnostics` - Annual Diagnostics

Accumulates annual‑mean fields and prints a summary at year end.

- Accumulates 12 fields every time step: `Ts`, `Ta`, `To`, `q`, `albedo`, `SW`, `LW`, `Q_lat`, `Q_sens`, and flux corrections `TF_correct`, `qF_correct`
- At `ityr == nstep_yr`: divides by number of steps, prints global mean temperature and two sample points (°C), then resets all accumulators to zero
"""

# ╔═╡ 7281dc60-e0e2-4a34-b4b6-c3f63a97f60f
md"""
---
## 💾 `output` - Monthly Output

Accumulates monthly‑mean fields and pushes a `MonthlyRecord` to the output buffer at each month boundary.

- Accumulates 13 fields every time step via `MonthlyAccumulator`
- At the last time step of each calendar month: averages by `days_in_month × steps_per_day`, pushes a copy as a NamedTuple, resets accumulator, advances month counter
- Returns updated `(irec, mon)` for bookkeeping
"""

# ╔═╡ db95cf44-4d37-463c-8fa7-8f084c2a20e8
md"""
---
## 🔄 Data Loading

Before running the model, load the climate input data from the `input/` directory:
"""

# ╔═╡ 4995d3d8-1f95-41d8-be6c-50663edbce10
begin
	jld2_dir = joinpath(@__DIR__, "..", "greb_dataset_jld2")  # dataset lives at repo root, not notebooks/
	fields = load_greb_jld2!(jld2_dir; dataset=:ncep)
end

# ╔═╡ 83e54812-2291-44e6-9b6e-0c0be57865d3
md"""
---
## `greb_model` - Main Integration Driver

Orchestrates the three simulation phases:

1. **Initialisation** - calls `init_model!` to set up derived fields, experiment overrides, initial conditions, and control CO₂
2. **Flux correction spin‑up** (`time_flux` years) - computes `TF_correct`, `ToF_correct`, `qF_correct` by nudging the model toward climatology; skipped if topography and q‑flux corrections are both disabled
3. **Control run** (`time_ctrl` years) - integrates at constant CO₂ (or constant forcing for orbital experiments); stores monthly output; builds ice climatology from control results
4. **Scenario run** (`time_scnr` years) - integrates with time‑varying CO₂ from `forcing()`; for forced‑boundary experiments (`:elnino`, `:lanina`, `:rcp85`), `Ts` is overwritten with climatology each step; for orbital experiments, starts at year 1 instead of 1950

Returns `(ctrl, scnr)` - vectors of `MonthlyRecord`. Scenario output is converted to anomalies relative to control climatology (except for orbital experiments).
"""

# ╔═╡ d17c570f-ad56-4d43-98cb-ccd03a8fcb4a
md"""
---
## 🎛️ Interactive Control Panel

Configure and run GREB experiments interactively via Pluto widgets.

**Controls:**
- **Experiment** - preset experiment type (full model, 2×CO₂, El Niño, etc.)
- **Configuration preset** - process isolation presets (full physics, no feedbacks, MSCM, sensitivity test, custom)
- **Mean climate / CO₂ response / circulation panels** - individual process toggles for custom configurations
- **Hydrology** - rain parameterization mode, evaporation mode, climatology dataset
- **External forcing** - toggles for forced surface temperature, wind, and vertical velocity

**Run controls:**
- Flux correction years, control run years, scenario run years (sliders)
- Run toggle checkbox to execute the model
- Results stored in `last_run` (control and scenario `MonthlyRecord` vectors)
"""

# ╔═╡ 71676720-4685-4216-a55e-5918829d258f
md"""
### Experiment Selection
"""

# ╔═╡ b7cef546-f0f4-4f42-a0ed-00ffb92a422e
begin
@bind experiment_type Select([
    :full_model => "Full Model (default)",
    :constant_topo => "Constant Topography",
    :co2_double => "2×CO₂",
    :co2_quadruple => "4×CO₂",
    :solar_plus27 => "Solar +27W/m²",
    :elnino => "El Niño",
    :lanina => "La Niña",
    :paleo_231kyr => "Paleo 231kyr",
    :rcp85 => "RCP8.5 Climate Change"
], default=:full_model)
end

# ╔═╡ 7969e897-121e-4fe0-9b75-36d21931357f
md"""
#### 🎛️ Configuration Preset
"""

# ╔═╡ dec9de54-eede-45ae-94cd-5bcc477298bf
@bind config_preset Select([
	"full" => "Full Physics (All processes active)",
	"no_feedbacks" => "No Feedbacks (Fixed albedo, clouds, vapor)",
	"mscm" => "MSCM Original (Mean State Climate Model)",
	"sensitivity" => "Sensitivity Test (Minimal processes)",
	"custom" => "Custom Configuration"
], default="full")

# ╔═╡ 2d2aa153-1b51-4972-89a0-de0f5d803838
md"""
#### 🌡️ Mean Climate Processes
"""

# ╔═╡ 66527a77-1926-41f7-b86f-84e0d65d9f30
@bind mean_climate_panel PlutoUI.combine() do Child
	md"""
	$(Child("clouds", CheckBox(default=true))) Clouds  
	$(Child("vapor", CheckBox(default=true))) Water Vapor  
	$(Child("ice", CheckBox(default=true))) Ice-Albedo Feedback  
	$(Child("crcl", CheckBox(default=true))) Circulation  
	$(Child("hydro", CheckBox(default=true))) Hydrology  
	$(Child("atmos", CheckBox(default=true))) Atmosphere  
	$(Child("co2", CheckBox(default=true))) CO₂  
	$(Child("ocean", CheckBox(default=true))) Ocean  
	$(Child("qflux", CheckBox(default=true))) Q-Flux Corrections
	"""
end

# ╔═╡ 393a8bae-c5ee-4a46-85f3-2c6a7cead357
md"""
#### 🔥 CO₂ Response Processes
"""

# ╔═╡ 506e51a4-d054-4e5d-8741-c4d36a9f2714
@bind co2_response_panel PlutoUI.combine() do Child
	md"""
	$(Child("clouds", CheckBox(default=true))) Clouds  
	$(Child("vapor", CheckBox(default=true))) Water Vapor  
	$(Child("crcl", CheckBox(default=true))) Circulation  
	$(Child("hydro", CheckBox(default=true))) Hydrology  
	$(Child("topo", CheckBox(default=true))) Topography  
	$(Child("humid", CheckBox(default=true))) Humidity Climatology
	"""
end

# ╔═╡ 3ac049a9-50e2-4db6-acb9-aacc20586ca4
md"""
#### ⚡ Circulation Components
"""

# ╔═╡ 26b4d172-72ec-4dd6-83c0-0779904634ae
@bind circulation_panel PlutoUI.combine() do Child
	md"""
	$(Child("ice", CheckBox(default=true))) Ice-Albedo Feedback  
	$(Child("hdif", CheckBox(default=true))) Horizontal Diffusion  
	$(Child("hadv", CheckBox(default=true))) Horizontal Advection  
	$(Child("vdif", CheckBox(default=true))) Vertical Diffusion  
	$(Child("vadv", CheckBox(default=true))) Vertical Advection  
	$(Child("conv", CheckBox(default=true))) Moisture Convergence
	"""
end

# ╔═╡ 2d4f533e-a240-459d-9575-f9982c03184b
md"""
#### 💧 Hydrology Parameterizations
"""

# ╔═╡ 8ef38e96-42f4-4bb2-9063-88db9842fb47
@bind hydro_rain_mode Select([
	-1 => "Original GREB",
	0 => "Best Performance",
	1 => "+Relative Humidity",
	2 => "+Omega Convergence",
	3 => "+Both RH & Omega"
], default=0)

# ╔═╡ d5331962-6bb1-485d-b538-ff35221e14e6
@bind hydro_eva_mode Select([
	-1 => "Original GREB",
	0 => "Skin Temperature",
	1 => "Enhanced",
	2 => "Wind Speed Climatology"
], default=-1)

# ╔═╡ 54e19b14-4bc1-4086-9015-ff2cd8f4afbf
@bind hydro_clim_dataset Select([
	0 => "ERA-Interim",
	1 => "NCEP"
], default=0)

# ╔═╡ efb23f98-fb75-42e6-8681-c5147a82a09c
md"""
#### 📡 External Forcing
"""

# ╔═╡ d74882fa-3a38-4cc7-b476-e52cb5a7db31
@bind external_forcing_panel PlutoUI.combine() do Child
	md"""
	$(Child("tsurf", CheckBox(default=false))) Surface Temperature  
	$(Child("hwind", CheckBox(default=false))) Horizontal Wind  
	$(Child("omega", CheckBox(default=false))) Vertical Velocity
	"""
end

# ╔═╡ dc0f1a5e-1e7b-4a9d-816e-731d2747ac4e
begin
	# Apply preset configurations or use custom settings
	if config_preset == "full"
		# Full Physics: All processes active
		log_clouds_dmc = log_vapor_dmc = log_ice_dmc = true
        log_crcl_dmc = log_hydro_dmc = true
        log_atmos_dmc = log_co2_dmc = log_ocean_dmc = log_qflux_dmc = true
        
        log_clouds_drsp = log_vapor_drsp = true
        log_crcl_drsp = log_hydro_drsp = true
        log_topo_drsp = log_humid_drsp = true
        
        log_ice = log_hdif = log_hadv = true
        log_vdif = log_vadv = log_conv = true
        
        log_rain, log_eva, log_clim = 0, -1, 0
        log_tsurf_ext = log_hwind_ext = log_omega_ext = false
		
	elseif config_preset == "no_feedbacks"
		# No Feedbacks: Fixed albedo, clouds, vapor
		log_clouds_dmc = log_vapor_dmc = log_ice_dmc = false
        log_crcl_dmc = log_hydro_dmc = true
        log_atmos_dmc = log_co2_dmc = log_ocean_dmc = log_qflux_dmc = true
        
        log_clouds_drsp = log_vapor_drsp = false
        log_crcl_drsp = log_hydro_drsp = true
        log_topo_drsp = log_humid_drsp = true
        
        log_ice = false
        log_hdif = log_hadv = true
        log_vdif = log_vadv = log_conv = true
        
        log_rain, log_eva, log_clim = 0, -1, 0
        log_tsurf_ext = log_hwind_ext = log_omega_ext = false
		
	elseif config_preset == "mscm"
		# MSCM Original: Mean State Climate Model configuration
		log_clouds_dmc = log_vapor_dmc = log_ice_dmc = true
        log_crcl_dmc = log_hydro_dmc = true
        log_atmos_dmc = log_co2_dmc = log_ocean_dmc = log_qflux_dmc = true
        
        log_clouds_drsp = log_vapor_drsp = true
        log_crcl_drsp = log_hydro_drsp = true
        log_topo_drsp = log_humid_drsp = true
        
        log_ice = log_hdif = log_hadv = true
        log_vdif = log_vadv = log_conv = true
        
        log_rain, log_eva, log_clim = 0, -1, 0
        log_tsurf_ext = log_hwind_ext = log_omega_ext = false
		
	elseif config_preset == "sensitivity"
		# Sensitivity Test: Minimal processes for testing
		log_clouds_dmc = log_vapor_dmc = log_ice_dmc = false
        log_crcl_dmc = log_hydro_dmc = false
        log_atmos_dmc = log_co2_dmc = log_ocean_dmc = true
        log_qflux_dmc = false
        
        log_clouds_drsp = log_vapor_drsp = false
        log_crcl_drsp = log_hydro_drsp = false
        log_topo_drsp = log_humid_drsp = false
        
        log_ice = log_hdif = log_hadv = false
        log_vdif = log_vadv = log_conv = false
        
        log_rain, log_eva, log_clim = -1, -1, 0
        log_tsurf_ext = log_hwind_ext = log_omega_ext = false
		
	else # custom
		log_clouds_dmc = mean_climate_panel.clouds
        log_vapor_dmc = mean_climate_panel.vapor
        log_ice_dmc = mean_climate_panel.ice
        log_crcl_dmc = mean_climate_panel.crcl
        log_hydro_dmc = mean_climate_panel.hydro
        log_atmos_dmc = mean_climate_panel.atmos
        log_co2_dmc = mean_climate_panel.co2
        log_ocean_dmc = mean_climate_panel.ocean
        log_qflux_dmc = mean_climate_panel.qflux
        
        log_clouds_drsp = co2_response_panel.clouds
        log_vapor_drsp = co2_response_panel.vapor
        log_crcl_drsp = co2_response_panel.crcl
        log_hydro_drsp = co2_response_panel.hydro
        log_topo_drsp = co2_response_panel.topo
        log_humid_drsp = co2_response_panel.humid
        
        log_ice = circulation_panel.ice
        log_hdif = circulation_panel.hdif
        log_hadv = circulation_panel.hadv
        log_vdif = circulation_panel.vdif
        log_vadv = circulation_panel.vadv
        log_conv = circulation_panel.conv
        
        log_rain = hydro_rain_mode
        log_eva = hydro_eva_mode
        log_clim = hydro_clim_dataset
        
        log_tsurf_ext = external_forcing_panel.tsurf
        log_hwind_ext = external_forcing_panel.hwind
        log_omega_ext = external_forcing_panel.omega
	end
end;

# ╔═╡ fee2639d-d8cf-4ee3-b824-73a0b02fef4a
md"""
### Run Duration Configuration (years)
"""

# ╔═╡ 952d5e5d-119a-4ed9-9e22-0b0fe66ae04d
@bind time_flux Slider(0:10, default=0, show_value=true)

# ╔═╡ 82bc0ee6-d0ca-4f47-a43f-67ae6702bb6a
md"""
**Flux correction years:** $(time_flux)
"""

# ╔═╡ 13647b0a-561c-4f65-a64c-b9a4ea76c9f0
@bind time_ctrl Slider(0:100, default=10, show_value=true)

# ╔═╡ 28371449-d04e-4d64-8b4a-cd37ceb25ef5
md"""
**Control run years:** $(time_ctrl)
"""

# ╔═╡ 4d3bedd5-1d4a-4150-828f-c5352c70874b
@bind time_scnr Slider(0:100, default=0, show_value=true)

# ╔═╡ 0ab27970-986e-4359-bed1-8908c5fd109b
md"""
**Scenario run years:** $(time_scnr)
"""

# ╔═╡ 10e20296-cc4f-4c4a-80bf-c39ae6450e85
md"""
### Execute Model

**Current Configuration:**
- Experiment: $(experiment_type)
- Flux correction: $time_flux years
- Control run: $time_ctrl years  
- Scenario run: $time_scnr years
"""

# ╔═╡ 9b7a847e-a0d2-4421-9fb9-ff01b73c4194
@bind run_toggle CheckBox(default=false)

# ╔═╡ f7fb4846-5ee8-4ca5-a76f-1edcd58a0559
md"""
### Accessing Results

After running the model, results are available in the `last_run` variable:

```julia
# Control run monthly records
last_run.ctrl  # Vector{MonthlyRecord}

# Scenario run monthly records  
last_run.scnr  # Vector{MonthlyRecord}

# Each MonthlyRecord contains:
# .Ts      - Surface temperature (96×48)
# .Ta      - Atmospheric temperature (96×48)
# .To      - Ocean temperature (96×48)
# .q       - Atmospheric humidity (96×48)
# .albedo  - Surface albedo (96×48)
# .ice     - Ice coverage (96×48)
# .precip  - Precipitation (96×48)
# .evap    - Evaporation (96×48)
# .qcrcl   - Circulation heat flux (96×48)
# .sw      - Shortwave radiation (96×48)
# .lw      - Longwave radiation (96×48)
# .qlat    - Latent heat flux (96×48)
# .qsens   - Sensible heat flux (96×48)
```

**Example: Plot global mean temperature evolution**
```julia
Ts_global_mean = [mean(rec.Ts) for rec in last_run.ctrl]
plot(Ts_global_mean, xlabel="Month", ylabel="Global Mean Ts [K]")
```
"""

# ╔═╡ 39afec69-f90d-4ff2-99d1-5cfec84d09a1
function current_physics_config()
	return PhysicsConfig(
		# Mean Climate Switches
		log_clouds_dmc = log_clouds_dmc,
		log_vapor_dmc = log_vapor_dmc,
		log_crcl_dmc = log_crcl_dmc,
		log_hydro_dmc = log_hydro_dmc,
		log_atmos_dmc = log_atmos_dmc,
		log_co2_dmc = log_co2_dmc,
		log_ocean_dmc = log_ocean_dmc,
		log_qflux_dmc = log_qflux_dmc,
		# CO₂ Response Switches
		# NOTE: log_ice_dmc / log_vapor_drsp are set below in the config-preset
		# cell but no longer exist on GREB.PhysicsConfig (removed upstream in
		# the package); the corresponding "ice"/"vapor" checkboxes above are
		# kept for UI parity with the mean-climate/CO₂-response panels but are
		# currently no-ops here.
		log_clouds_drsp = log_clouds_drsp,
		log_crcl_drsp = log_crcl_drsp,
		log_hydro_drsp = log_hydro_drsp,
		log_topo_drsp = log_topo_drsp,
		log_humid_drsp = log_humid_drsp,
		# Circulation Components
		log_ice = log_ice,
		log_hdif = log_hdif,
		log_hadv = log_hadv,
		log_vdif = log_vdif,
		log_vadv = log_vadv,
		log_conv = log_conv,
		# Hydrology Parameters
		log_rain = log_rain,
		log_eva = log_eva,
		log_clim = log_clim,
		# External Forcing
		log_tsurf_ext = log_tsurf_ext,
		log_hwind_ext = log_hwind_ext,
		log_omega_ext = log_omega_ext,
		experiment = experiment_type,
	)
end

# ╔═╡ cf2a51eb-a661-4c66-9e0a-d1730110e4bc
begin
	if !@isdefined(last_run)
		global last_run = nothing
	end

	if run_toggle
		cfg = current_physics_config()
		run = RunSpec(flux=time_flux, ctrl=time_ctrl, scnr=time_scnr)
		local result
		@progress name="GREB run" for _ in 1:1
			result = greb_model!(run, cfg; jld2_dir=jld2_dir, fields=fields)
		end
		global last_run = result
		md"""**✅ Run complete!**
		Control months: $(length(result.ctrl))
		Scenario months: $(length(result.scnr))"""
	else
		md"""**⏸️ Toggle ON to run model**"""
	end
end

# ╔═╡ b0d6a9e3-5c3a-4408-ab87-1067c2dfaeb7
begin
    if @isdefined(last_run) && !isnothing(last_run) && !isempty(last_run.ctrl)
        println("="^50)
        println("GREB MODEL OUTPUT")
        println("="^50)
        
        rec = last_run.ctrl[1]
        
        println("\n🌡️ TEMPERATURE (K):")
        println("   Surface (Ts): min=$(round(minimum(rec.Ts), digits=1)), mean=$(round(mean(rec.Ts), digits=1)), max=$(round(maximum(rec.Ts), digits=1))")
        println("   Air (Ta):     min=$(round(minimum(rec.Ta), digits=1)), mean=$(round(mean(rec.Ta), digits=1)), max=$(round(maximum(rec.Ta), digits=1))")
        println("   Ocean (To):   min=$(round(minimum(rec.To), digits=1)), mean=$(round(mean(rec.To), digits=1)), max=$(round(maximum(rec.To), digits=1))")
        
        println("\n💧 HYDROLOGY (mm/day):")
        println("   Precip: min=$(round(minimum(rec.precip), digits=2)), mean=$(round(mean(rec.precip), digits=2)), max=$(round(maximum(rec.precip), digits=2))")
        println("   Evap:   min=$(round(minimum(rec.evap), digits=2)), mean=$(round(mean(rec.evap), digits=2)), max=$(round(maximum(rec.evap), digits=2))")
        
        println("\n☀️ RADIATION (W/m²):")
        println("   SW: min=$(round(minimum(rec.sw), digits=1)), mean=$(round(mean(rec.sw), digits=1)), max=$(round(maximum(rec.sw), digits=1))")
        println("   LW: min=$(round(minimum(rec.lw), digits=1)), mean=$(round(mean(rec.lw), digits=1)), max=$(round(maximum(rec.lw), digits=1))")
        
        println("\n🪞 ALBEDO:")
        println("   min=$(round(minimum(rec.albedo), digits=2)), mean=$(round(mean(rec.albedo), digits=2)), max=$(round(maximum(rec.albedo), digits=2))")
        
        println("\n❄️ ICE FRACTION:")
        println("   min=$(round(minimum(rec.ice), digits=2)), mean=$(round(mean(rec.ice), digits=2)), max=$(round(maximum(rec.ice), digits=2))")
        
        println("\n✅ All finite: $(all(isfinite, rec.Ts))")
        
    else
        println("❌ No model output found. Run the model first!")
        println("   Set time_ctrl >= 1 and toggle run_toggle ON")
    end
end

# ╔═╡ 9bf4c37a-461c-4b03-abfc-e599f583244f
md"### 📈 Visualization"

# ╔═╡ e64b2bf2-9bf4-4bd5-a5ab-89dc085220fb
@bind plot_field Select([:Ts=>"Surface Temp (Ts)", :Ta=>"Air Temp (Ta)", :precip=>"Precipitation", :evap=>"Evaporation", :ice=>"Ice fraction", :sw=>"Shortwave", :lw=>"Longwave", :albedo=>"Albedo"], default=:Ts)

# ╔═╡ 5c9612a4-ead3-42a3-a73b-74a4e2d96575
@bind plot_run Select([:ctrl=>"Control run", :scnr=>"Scenario run"], default=:ctrl)

# ╔═╡ a9089bf1-3241-4038-81d6-b953a00c3cef
md"""
---
## ⏱️ Benchmarking

Benchmarking components to find optimizations.
"""

# ╔═╡ 71667e2a-79ef-47dc-ab74-8606e5e5699c
@bind run_benchmarks CheckBox(default=false)

# ╔═╡ 698a990b-cd78-43a1-a747-e79102767d62
# ╠═╡ skip_as_script = true
begin
	using BenchmarkTools, Profile, ProfileSVG

	# Create test data for benchmarking
	function setup_benchmark()
	    xdim_local = 96
	    ydim_local = 48

	    # Create workspace
	    ws = CirculationWorkspace()

	    # Ensure dX_crcl is properly sized
	    if size(ws.dX_crcl) != (xdim_local, ydim_local)
	        ws.dX_crcl = zeros(Float64, xdim_local, ydim_local)
	    end

	    # Time state
	    timestate = TimeState(1, 1)

	    # Configuration
	    cfg = PhysicsConfig()

	    # ClimateFields/ModelState — every current physics function takes
	    # `fields::ClimateFields` (and SWradiation!/tendencies! also take
	    # `state::ModelState`); a fresh (all-zero) instance is fine here since
	    # benchmarks only measure timing, not correctness.
	    fields = ClimateFields()
	    fields.omegaclim .= rand(Float64, xdim_local, ydim_local, nstep_yr) .* 0.1
	    state = ModelState()

	    # Climate fields
	    Ts = rand(Float64, xdim_local, ydim_local) .* 80 .+ 233.15
	    Ta = copy(Ts)
	    To = rand(Float64, xdim_local, ydim_local) .* 30 .+ 270.0
	    q = rand(Float64, xdim_local, ydim_local) .* 0.02
	    CO2 = 340.0

	    # For advection and convergence tests
	    T_test = rand(Float64, xdim_local, ydim_local)

	    return (ws=ws, timestate=timestate, cfg=cfg, fields=fields, state=state,
	            Ts=Ts, Ta=Ta, To=To, q=q, CO2=CO2,
	            z_air=8400.0, z_vapor=5000.0,
	            T_test=T_test)
	end

	if run_benchmarks
		# Setup once
		bench_data = setup_benchmark()
	else
		md"*Benchmarks skipped — enable the checkbox above to run.*"
	end
end

# ╔═╡ 8de366ff-1753-4ac5-84e7-b4cdaef18a39
md"☐ Enable to run BenchmarkTools/Profile cells below — skipped by default, several take minutes."

# ╔═╡ 9fa7c129-2289-4979-aa1f-8f8112e7875f
md"""### SWradiation:"""

# ╔═╡ caf443c0-5fbc-4fd4-8ef9-3565417d40ab
# ╠═╡ skip_as_script = true
#=╠═╡
if run_benchmarks
	# Benchmark SWradiation
	@benchmark SWradiation!($bench_data.Ts, $bench_data.fields, $bench_data.state, $bench_data.timestate, $bench_data.cfg, $bench_data.ws)
else
	md"*Benchmarks skipped — enable the checkbox above to run.*"
end
  ╠═╡ =#

# ╔═╡ 1d46e4a8-41fc-4b76-bc9d-a115b6b6d9c4
md"""### LWradiation:"""

# ╔═╡ e20657c6-cf03-4f84-a782-2a3e5d6d303c
# ╠═╡ skip_as_script = true
#=╠═╡
if run_benchmarks
	# Benchmark LWradiation
	@benchmark LWradiation!($bench_data.Ts, $bench_data.Ta, $bench_data.q, $bench_data.CO2, $bench_data.fields, $bench_data.timestate, $bench_data.cfg, $bench_data.ws)
else
	md"*Benchmarks skipped — enable the checkbox above to run.*"
end
  ╠═╡ =#

# ╔═╡ e2227570-2fcb-4090-8bbb-ec85af3a519c
md"""### hydro:"""

# ╔═╡ 1b94212a-a5b9-46dc-81ea-896e4c839755
# ╠═╡ skip_as_script = true
#=╠═╡
if run_benchmarks
	# Benchmark hydrological cycle
	@benchmark hydro!($bench_data.Ts, $bench_data.q, $bench_data.fields, $bench_data.timestate,
	                  $bench_data.cfg, $bench_data.ws)
else
	md"*Benchmarks skipped — enable the checkbox above to run.*"
end
  ╠═╡ =#

# ╔═╡ 40717ac9-ff3f-4835-ba7e-3cde10b2c3f0
md"""### circulation (temperature):"""

# ╔═╡ fb0719e8-b9fe-4ca3-8c5a-de3586d193a9
# ╠═╡ skip_as_script = true
#=╠═╡
if run_benchmarks
	@benchmark circulation!($bench_data.Ta, 8400.0, $bench_data.ws.dX_crcl, $bench_data.fields, $bench_data.ws, $bench_data.timestate, $bench_data.cfg)
else
	md"*Benchmarks skipped — enable the checkbox above to run.*"
end
  ╠═╡ =#

# ╔═╡ 0c201c53-0db2-4805-beec-e7705082509d
md"""### circulation (vapor):"""

# ╔═╡ eb10715c-247d-40e4-b11b-aceabb60628c
# ╠═╡ skip_as_script = true
#=╠═╡
if run_benchmarks
	@benchmark circulation!($bench_data.q, 5000.0, $bench_data.ws.dX_crcl, $bench_data.fields, $bench_data.ws, $bench_data.timestate, $bench_data.cfg)
else
	md"*Benchmarks skipped — enable the checkbox above to run.*"
end
  ╠═╡ =#

# ╔═╡ f9da0fed-a0a3-4bfb-8bac-6a839177bc73
md"""### seaice:"""

# ╔═╡ 2e05f3f6-687b-4096-a869-d64564db60da
# ╠═╡ skip_as_script = true
#=╠═╡
if run_benchmarks
	@benchmark seaice!($bench_data.Ts, $bench_data.fields, $bench_data.timestate, $bench_data.cfg)
else
	md"*Benchmarks skipped — enable the checkbox above to run.*"
end
  ╠═╡ =#

# ╔═╡ 09cbb86e-64e1-45e5-b0f6-1e8b868506e8
md"""### deep_ocean:"""

# ╔═╡ c186464e-e705-49be-ad5d-4352ce719c1f
# ╠═╡ skip_as_script = true
#=╠═╡
if run_benchmarks
	@benchmark deep_ocean!($bench_data.Ts, $bench_data.To, $bench_data.fields, $bench_data.timestate, $bench_data.cfg, $bench_data.ws)
else
	md"*Benchmarks skipped — enable the checkbox above to run.*"
end
  ╠═╡ =#

# ╔═╡ 395242dc-cb73-481f-ad25-71eba93d341f
md"""### diffusion:"""

# ╔═╡ 718e4f4b-af27-4271-a37f-83ea6e30bdbe
# ╠═╡ skip_as_script = true
#=╠═╡
if run_benchmarks
	begin
		T1 = rand(Float64, xdim, ydim)
		h_scl = 8400.0  # z_air
		@benchmark diffusion!($T1, $h_scl, $bench_data.fields, $bench_data.ws, $bench_data.timestate)
	end
else
	md"*Benchmarks skipped — enable the checkbox above to run.*"
end
  ╠═╡ =#

# ╔═╡ 59c97f2b-de11-4f42-9fcf-aaf2924d8677
md"""### advection:"""

# ╔═╡ b9b3564e-c59b-4f3a-80d7-4431fd8ce2db
# ╠═╡ skip_as_script = true
#=╠═╡
if run_benchmarks
	@benchmark advection!($bench_data.T_test, 8400.0, $bench_data.fields, $bench_data.ws, $bench_data.timestate, $bench_data.cfg)
else
	md"*Benchmarks skipped — enable the checkbox above to run.*"
end
  ╠═╡ =#

# ╔═╡ 435bae7c-e0d3-46af-b185-5221f687d6bf
md"""### convergence:"""

# ╔═╡ 1229bf7d-204c-41a2-97dd-9d75e79b18ce
# ╠═╡ skip_as_script = true
#=╠═╡
if run_benchmarks
	@benchmark convergence!($bench_data.T_test, $bench_data.fields, $bench_data.timestate, $bench_data.ws)
else
	md"*Benchmarks skipped — enable the checkbox above to run.*"
end
  ╠═╡ =#

# ╔═╡ 426c1511-f912-4e10-9b8b-f858471019bc
md"""### tendencies:"""

# ╔═╡ 50777c8a-ee6f-449b-b704-bba7e9b36490
# ╠═╡ skip_as_script = true
#=╠═╡
if run_benchmarks
	begin
		@benchmark tendencies!($bench_data.CO2, $bench_data.Ts, $bench_data.Ta, $bench_data.To, $bench_data.q, $bench_data.fields, $bench_data.state, $bench_data.ws, $bench_data.timestate, $bench_data.cfg)
	end
else
	md"*Benchmarks skipped — enable the checkbox above to run.*"
end
  ╠═╡ =#

# ╔═╡ d571c448-cea6-49c0-81a4-2376e4fc3f67
# ╠═╡ skip_as_script = true
#=╠═╡
if run_benchmarks
	begin
	### 📊 Comprehensive Performance Profiling Report (full width, high samples)
	# Run this after the benchmark setup (bench_data must exist).
	# This will take ~2‑3 minutes – be patient.

	if !@isdefined(bench_data)
	    error("bench_data not defined. Please run the benchmark setup cell first.")
	end

	println("Warming up...")
	for _ in 1:10
	    tendencies!(bench_data.CO2, bench_data.Ts, bench_data.Ta, bench_data.To,
	                bench_data.q, bench_data.fields, bench_data.state, bench_data.ws, bench_data.timestate, bench_data.cfg)
	end

	# -------- CPU Profile (time) – 500 iterations --------
	println("Running CPU profile (500 iterations)...")
	Profile.clear()
	@profile for _ in 1:500
	    tendencies!(bench_data.CO2, bench_data.Ts, bench_data.Ta, bench_data.To,
	                bench_data.q, bench_data.fields, bench_data.state, bench_data.ws, bench_data.timestate, bench_data.cfg)
	end

	cpu_io = IOBuffer()
	cpu_ctx = IOContext(cpu_io, :maxwidth => 0)          # full width
	Profile.print(cpu_ctx, format=:tree, mincount=1, C=false)
	cpu_profile = String(take!(cpu_io))

	# Flat profile (sorted by count) – also full width
	cpu_flat_io = IOBuffer()
	cpu_flat_ctx = IOContext(cpu_flat_io, :maxwidth => 0)
	Profile.print(cpu_flat_ctx, format=:flat, mincount=5, C=false)
	cpu_flat = String(take!(cpu_flat_io))

	# -------- Allocation Profile (memory) – 100 iterations --------
	println("Running allocation profile (100 iterations)...")
	Profile.Allocs.@profile sample_rate=0.01 begin
	    for _ in 1:100
	        tendencies!(bench_data.CO2, bench_data.Ts, bench_data.Ta, bench_data.To,
	                    bench_data.q, bench_data.fields, bench_data.state, bench_data.ws, bench_data.timestate, bench_data.cfg)
	    end
	end

	allocs_io = IOBuffer()
	allocs_ctx = IOContext(allocs_io, :maxwidth => 0)    # full width
	# No maxwidth keyword – width is controlled by IOContext
	Profile.Allocs.print(allocs_ctx, format=:tree, mincount=100)
	allocs_profile = String(take!(allocs_io))

	# -------- Benchmark summary --------
	println("Running benchmark...")
	trial = @benchmark tendencies!($bench_data.CO2, $bench_data.Ts, $bench_data.Ta, $bench_data.To,
	                               $bench_data.q, $bench_data.fields, $bench_data.state, $bench_data.ws, $bench_data.timestate, $bench_data.cfg)
	bench_summary = sprint(show, trial)

	# -------- Assemble Markdown report --------
	Markdown.MD([
	    Markdown.Header("Performance Profiling Report (Full Width)", 1),
	    Markdown.Header("Benchmark Summary (per `tendencies!` call)", 2),
	    Markdown.Code(bench_summary, "text"),
	    Markdown.Header("CPU Profile – Flat (sorted by sample count, mincount=5)", 2),
	    Markdown.Code(cpu_flat, "text"),
	    Markdown.Header("CPU Profile – Tree (call hierarchy, mincount=1)", 2),
	    Markdown.Code(cpu_profile, "text"),
	    Markdown.Header("Allocation Profile – Tree (mincount=100)", 2),
	    Markdown.Code(allocs_profile, "text"),
	    Markdown.Header("What to look for", 2),
	    Markdown.MD("""
	    - **Flat profile**: shows which functions have the most samples → primary bottlenecks.
	    - **Tree profile**: shows the call stack, so you see *why* those functions are called.
	    - **Allocation profile**: reveals lines that allocate memory → causes GC overhead.
	    - **Benchmark summary**: average time & memory per call.

	    Focus on:
	    - `getindex` / `setindex!` inside loops → non‑contiguous memory access.
	    - `broadcast` / `materialize!` → often creates temporary arrays.
	    - `log`, `exp`, `pow` → expensive math; consider pre‑computation or approximations.
	    - `invokelatest` or `wait` → dynamic dispatch or threading issues.
	    """)
	])

	# Save to file (full width, no truncation)
	open("profile_report.txt", "w") do io
	    println(io, "=== BENCHMARK SUMMARY ===\n", bench_summary)
	    println(io, "\n=== CPU FLAT PROFILE ===\n", cpu_flat)
	    println(io, "\n=== CPU TREE PROFILE ===\n", cpu_profile)
	    println(io, "\n=== ALLOCATION PROFILE ===\n", allocs_profile)
	end
	println("\nReport also saved to profile_report.txt")
	end
else
	md"*Benchmarks skipped — enable the checkbox above to run.*"
end
  ╠═╡ =#

# ╔═╡ 3dd14a9f-130f-4431-a2f3-51df38483a57
# ╠═╡ skip_as_script = true
#=╠═╡
function export_to_netcdf(records::Vector{MonthlyRecord}, filename::String;
                          lon_shift::Bool = true)  # true → shift to [-180,180]
    nt = length(records)
    xdim, ydim = size(records[1].Ts)  # 96, 48

    # Build coordinate arrays
    dlon = 360.0 / xdim
    dlat = 180.0 / ydim
    lon = (0:xdim-1) * dlon .+ dlon/2
    lat = (1:ydim) .* dlat .- dlat/2 .- 90.0

    if lon_shift
        lon = mod.(lon .+ 180, 360) .- 180  # wrap to [-180,180]
    end

    # Create NetCDF file
    ds = Dataset(filename, "c")
    defDim(ds, "lon", xdim)
    defDim(ds, "lat", ydim)
    defDim(ds, "time", nt)

    # Write coordinates
    v = defVar(ds, "lon", Float64, ("lon",))
    v[:] = lon
    v = defVar(ds, "lat", Float64, ("lat",))
    v[:] = lat
    v = defVar(ds, "time", Float64, ("time",))
    v[:] = 1:nt  # or use actual dates if available

    # Define all output variables
    field_names = propertynames(records[1])
    vars = Dict()
    for f in field_names
        vars[f] = defVar(ds, String(f), Float32, ("lon", "lat", "time"),
                         deflatelevel=1)
    end

    # Fill data
    for (t, rec) in enumerate(records)
        for f in field_names
            vars[f][:, :, t] = getfield(rec, f)
        end
    end

    close(ds)
    println("Exported to $filename")
end
  ╠═╡ =#

# ╔═╡ 5713d5d6-c95f-4a34-b450-c240d4a2148b
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	# Assume last_run is the output from greb_model!
ctrl_records = last_run.ctrl
scnr_records = last_run.scnr

export_to_netcdf(ctrl_records, "greb_control.nc")
export_to_netcdf(scnr_records, "greb_scenario.nc")
end
  ╠═╡ =#

# ╔═╡ ff9e5cd0-8424-44ac-81fb-dfc195598677
begin
	if @isdefined(last_run) && !isnothing(last_run) && !isempty(last_run.ctrl)
		records = plot_run == :ctrl ? last_run.ctrl : last_run.scnr
		if isempty(records)
			md"*No months in the $(plot_run) run — nothing to plot.*"
		else
			field_data = getfield(records[end], plot_field)
			heatmap(field_data'; xlabel="lon index", ylabel="lat index",
				title="$(plot_field) — last month of $(plot_run) run", c=:viridis)
		end
	else
		md"*Run the model first (toggle above) to see plots.*"
	end
end

# ╔═╡ 3e6c76e7-22d7-4de6-a833-90c8b7567781
begin
	if @isdefined(last_run) && !isnothing(last_run) && !isempty(last_run.ctrl)
		records = plot_run == :ctrl ? last_run.ctrl : last_run.scnr
		if isempty(records)
			md"*No months in the $(plot_run) run — nothing to plot.*"
		else
			series = [mean(getfield(rec, plot_field)) for rec in records]
			plot(series; xlabel="Month", ylabel=String(plot_field), legend=false,
				title="Global mean $(plot_field) — $(plot_run) run")
		end
	else
		md"*Run the model first (toggle above) to see plots.*"
	end
end

# ╔═╡ Cell order:
# ╠═1165eeb4-10d8-4fb5-859d-7fe410189608
# ╟─1cf6d892-c827-48fd-8e77-bb1117d03c20
# ╟─3f63ce64-effd-49c7-9906-92cda40d59f0
# ╟─0780e492-cdcd-43fb-af85-f9a4bab37450
# ╟─d0ebd34f-a83a-4dc1-9c66-24a3878b681d
# ╟─ebca5877-1305-40f1-ac52-66356dc17661
# ╟─8165ed3f-ede0-45c9-83ac-e4e49262457c
# ╟─07a88c94-93bb-4564-8737-980144c4af43
# ╟─813d59f7-eeea-402a-96e8-e5cbe2fd4582
# ╟─54e63724-44ca-42b9-9246-7afb271b8154
# ╟─caec8065-9874-4d8c-8e7e-598a97506f06
# ╟─12a89062-bc27-4e36-8b5e-1670b4bd0b28
# ╟─1879fd40-0117-478b-acdd-18112738b81f
# ╟─0308940f-16bc-4b03-8321-9ea8c91a21c7
# ╟─78605b5f-b03d-4851-b339-aea63bad3688
# ╟─22d5e751-f095-42d7-9c24-78e546b3ce37
# ╟─698a990b-cd78-43a1-a747-e79102767d62
# ╟─d403b01e-7d9f-4778-813d-bbbb0bfdbfb9
# ╟─c94558a5-0ae2-4563-93d7-0609ea3cce41
# ╟─2d02d36f-efee-411e-befb-52d61b3fb16e
# ╟─1b38b24d-bb82-4d52-a26d-850d142e25d5
# ╟─0ff4be7c-580c-4db1-a874-175badc8c11a
# ╟─3beff4df-89da-4f25-a0e4-be03048c5f2c
# ╟─fd193b13-336f-42f8-8f93-3ca3bea2e620
# ╟─241496f6-580b-4e81-acbb-3295696fe40e
# ╟─a4ab662c-5876-4861-aecf-2f695d15e15d
# ╟─2a859d43-5768-4b98-88fa-978b12203513
# ╟─0febe534-237e-4921-b39c-3828dbae9d19
# ╟─a1c04c52-f7c9-430a-8791-7fea15650b2c
# ╟─34fb479a-9834-4354-8f36-2680bedec798
# ╟─ae38f814-fe9c-443f-ae3b-42fa7a7d199a
# ╟─a4cd8d40-eedd-4a4b-ac17-51e68334328c
# ╟─cc7db228-0a06-4812-8b32-6541988b2115
# ╟─99f0e9ad-438f-4fa5-be9d-36d0fa78d89c
# ╟─8fdcfefd-0490-434a-a2bb-d171557b6ae7
# ╟─7e514484-9d81-4c5c-83ab-191d68c13043
# ╟─7b8e4130-7d8e-492a-a949-bd7fb6808dd6
# ╟─7281dc60-e0e2-4a34-b4b6-c3f63a97f60f
# ╟─db95cf44-4d37-463c-8fa7-8f084c2a20e8
# ╟─4995d3d8-1f95-41d8-be6c-50663edbce10
# ╟─83e54812-2291-44e6-9b6e-0c0be57865d3
# ╟─d17c570f-ad56-4d43-98cb-ccd03a8fcb4a
# ╟─71676720-4685-4216-a55e-5918829d258f
# ╟─b7cef546-f0f4-4f42-a0ed-00ffb92a422e
# ╟─7969e897-121e-4fe0-9b75-36d21931357f
# ╟─dec9de54-eede-45ae-94cd-5bcc477298bf
# ╟─2d2aa153-1b51-4972-89a0-de0f5d803838
# ╟─66527a77-1926-41f7-b86f-84e0d65d9f30
# ╟─393a8bae-c5ee-4a46-85f3-2c6a7cead357
# ╟─506e51a4-d054-4e5d-8741-c4d36a9f2714
# ╟─3ac049a9-50e2-4db6-acb9-aacc20586ca4
# ╟─26b4d172-72ec-4dd6-83c0-0779904634ae
# ╟─2d4f533e-a240-459d-9575-f9982c03184b
# ╟─8ef38e96-42f4-4bb2-9063-88db9842fb47
# ╟─d5331962-6bb1-485d-b538-ff35221e14e6
# ╟─54e19b14-4bc1-4086-9015-ff2cd8f4afbf
# ╟─efb23f98-fb75-42e6-8681-c5147a82a09c
# ╟─d74882fa-3a38-4cc7-b476-e52cb5a7db31
# ╟─dc0f1a5e-1e7b-4a9d-816e-731d2747ac4e
# ╟─fee2639d-d8cf-4ee3-b824-73a0b02fef4a
# ╟─952d5e5d-119a-4ed9-9e22-0b0fe66ae04d
# ╟─82bc0ee6-d0ca-4f47-a43f-67ae6702bb6a
# ╟─13647b0a-561c-4f65-a64c-b9a4ea76c9f0
# ╟─28371449-d04e-4d64-8b4a-cd37ceb25ef5
# ╟─4d3bedd5-1d4a-4150-828f-c5352c70874b
# ╟─0ab27970-986e-4359-bed1-8908c5fd109b
# ╟─10e20296-cc4f-4c4a-80bf-c39ae6450e85
# ╟─9b7a847e-a0d2-4421-9fb9-ff01b73c4194
# ╟─f7fb4846-5ee8-4ca5-a76f-1edcd58a0559
# ╟─39afec69-f90d-4ff2-99d1-5cfec84d09a1
# ╟─cf2a51eb-a661-4c66-9e0a-d1730110e4bc
# ╟─b0d6a9e3-5c3a-4408-ab87-1067c2dfaeb7
# ╟─9bf4c37a-461c-4b03-abfc-e599f583244f
# ╟─e64b2bf2-9bf4-4bd5-a5ab-89dc085220fb
# ╟─5c9612a4-ead3-42a3-a73b-74a4e2d96575
# ╟─3e6c76e7-22d7-4de6-a833-90c8b7567781
# ╟─ff9e5cd0-8424-44ac-81fb-dfc195598677
# ╟─a9089bf1-3241-4038-81d6-b953a00c3cef
# ╟─71667e2a-79ef-47dc-ab74-8606e5e5699c
# ╟─8de366ff-1753-4ac5-84e7-b4cdaef18a39
# ╟─9fa7c129-2289-4979-aa1f-8f8112e7875f
# ╟─caf443c0-5fbc-4fd4-8ef9-3565417d40ab
# ╟─1d46e4a8-41fc-4b76-bc9d-a115b6b6d9c4
# ╟─e20657c6-cf03-4f84-a782-2a3e5d6d303c
# ╟─e2227570-2fcb-4090-8bbb-ec85af3a519c
# ╟─1b94212a-a5b9-46dc-81ea-896e4c839755
# ╟─40717ac9-ff3f-4835-ba7e-3cde10b2c3f0
# ╟─fb0719e8-b9fe-4ca3-8c5a-de3586d193a9
# ╟─0c201c53-0db2-4805-beec-e7705082509d
# ╟─eb10715c-247d-40e4-b11b-aceabb60628c
# ╟─f9da0fed-a0a3-4bfb-8bac-6a839177bc73
# ╟─2e05f3f6-687b-4096-a869-d64564db60da
# ╟─09cbb86e-64e1-45e5-b0f6-1e8b868506e8
# ╟─c186464e-e705-49be-ad5d-4352ce719c1f
# ╟─395242dc-cb73-481f-ad25-71eba93d341f
# ╟─718e4f4b-af27-4271-a37f-83ea6e30bdbe
# ╟─59c97f2b-de11-4f42-9fcf-aaf2924d8677
# ╟─b9b3564e-c59b-4f3a-80d7-4431fd8ce2db
# ╟─435bae7c-e0d3-46af-b185-5221f687d6bf
# ╟─1229bf7d-204c-41a2-97dd-9d75e79b18ce
# ╟─426c1511-f912-4e10-9b8b-f858471019bc
# ╟─50777c8a-ee6f-449b-b704-bba7e9b36490
# ╟─d571c448-cea6-49c0-81a4-2376e4fc3f67
# ╟─3dd14a9f-130f-4431-a2f3-51df38483a57
# ╟─5713d5d6-c95f-4a34-b450-c240d4a2148b
