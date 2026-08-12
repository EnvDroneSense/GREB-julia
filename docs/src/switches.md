# Physics Switches

A "switch" here means a [`PhysicsConfig`](@ref) field that turns a physical
process on/off or selects between parameterizations. They're grouped below
by the same families used in `src/config.jl`, since a single physical
process (e.g. clouds) is often controlled by *two different* switches with
different mechanisms - one for the mean-climate state, one for the response
to CO₂ forcing - so grouping by code family avoids duplicating that
explanation. Build a config with [`create_experiment_config`](@ref) and
override individual fields afterwards, e.g. `cfg.log_rain = 1`.

Two naming suffixes recur throughout: `_dmc` ("deconstruct mean climate")
and `_drsp` ("deconstruct response"). Setting a `_dmc` switch to `false`
removes that process from the model entirely and the model is normally
re-run to a new equilibrium - this quantifies how much the process shapes
the *mean* climate. Setting a `_drsp` switch to `false` keeps the control
climate fixed (flux corrections are re-estimated so the baseline is
unaffected) and only disables that process's contribution to the *change*
under CO₂ forcing - this isolates the process's role in climate
*sensitivity* rather than the mean state. This is the same deconstruction
approach used in the MSCM experiments (see [References](#references)
below, MSCM §2.1–2.2).

## Overview

| Switch | Type | Default | Family | Effect |
|:-------|:-----|:--------|:-------|:-------|
| `log_clouds_dmc` | Bool | `true` | Mean climate | Zeroes cloud climatology at init when `false` |
| `log_vapor_dmc`[^1] | Bool | `true` | Mean climate | No effect - see [Known Limitations](#known-limitations) |
| `log_crcl_dmc` | Bool | `true` | Mean climate | Disables circulation (advection+diffusion+convergence) when `false` |
| `log_hydro_dmc` | Bool | `true` | Mean climate | Disables hydrology (evap/rain/latent heat) when `false` |
| `log_atmos_dmc` | Bool | `true` | Mean climate | Master atmosphere switch - decouples atmosphere from surface when `false` |
| `log_co2_dmc` | Bool | `true` | Mean climate | Forces control-run CO₂ to 0 ppm when `false` |
| `log_ocean_dmc` | Bool | `true` | Mean climate | Master ocean-coupling switch - no mixed-layer/deep-ocean exchange when `false` |
| `log_qflux_dmc` | Bool | `true` | Mean climate | Governs whether flux corrections are computed, loaded, or zeroed |
| `log_clouds_drsp` | Bool | `true` | CO₂ response | Freezes cloud cover at a constant 0.7 when `false` |
| `log_crcl_drsp` | Bool | `true` | CO₂ response | Disables circulation's contribution to the CO₂ response when `false` |
| `log_hydro_drsp` | Bool | `true` | CO₂ response | Disables hydrology's contribution to the CO₂ response when `false` |
| `log_topo_drsp` | Bool | `true` | CO₂ response | Caps topography at 1.0 m ("constant topography") when `false` |
| `log_humid_drsp` | Bool | `true` | CO₂ response | Freezes humidity at a constant 0.0052 kg/kg when `false` |
| `log_ocean_drsp` | Bool | `true` | CO₂ response | Constant mixed-layer depth, no deep-ocean response, when `false` |
| `log_ice` | Bool | `true` | Circulation | Ice-albedo feedback - constant albedo/heat capacity when `false` |
| `log_hdif` | Bool | `true` | Circulation | Horizontal (heat) diffusion |
| `log_hadv` | Bool | `true` | Circulation | Horizontal (heat) advection |
| `log_vdif` | Bool | `true` | Circulation | Vertical/meridional diffusion (water vapor) |
| `log_vadv` | Bool | `true` | Circulation | Advection (water vapor) |
| `log_conv`[^2] | Bool | `true` | Circulation | Moisture convergence sub-step |
| `log_rain` | Int | `0` | Hydrology | Selects the precipitation parameterization |
| `log_eva` | Int | `-1` | Hydrology | Selects the evaporation/exchange-coefficient parameterization |
| `log_clim` | Int | `0` | Hydrology | NCEP-tuned hydrology coefficient override |
| `log_tsurf_ext` | Bool | `false` | External forcing | Loads surface-temperature anomaly forcing |
| `log_hwind_ext` | Bool | `false` | External forcing | Loads horizontal-wind anomaly forcing |
| `log_omega_ext` | Bool | `false` | External forcing | Loads vertical-velocity anomaly forcing |

[^1]: Declared in `src/config.jl` but never read anywhere in `src/` - has no effect. See [Known Limitations](#known-limitations).
[^2]: Behaves like the other circulation switches but is not exposed as a keyword by `create_experiment_config` - set `cfg.log_conv` directly. See [Known Limitations](#known-limitations).

## Mean-Climate Deconstruction Switches (`_dmc`)

These switches, forwarded as keywords by `create_experiment_config(:decon_mean_climate; ...)`,
remove a process from the model entirely so it can be re-equilibrated
without it - quantifying that process's contribution to the mean climate
(DF11 §3.8; MSCM §2.1).

| Switch | Default | Effect when `false` |
|:-------|:--------|:---------------------|
| `log_clouds_dmc` | `true` | Cloud climatology zeroed at init - no cloud shortwave reflection or longwave-emissivity contribution (DF11 §3.1–3.2). Cloud cover is a *prescribed* seasonal field in GREB, not a simulated/responding one, so this removes the process outright rather than disabling a dynamic feedback. |
| `log_vapor_dmc` | `true` | No effect (dead switch - see [Known Limitations](#known-limitations)). |
| `log_crcl_dmc` | `true` | `circulation!` returns zero tendencies - no advection, diffusion, or convergence. GREB's circulation is a fixed climatological wind field, not a responding dynamical core (DF11 §3.5), so this isolates how much of the spatial temperature/humidity pattern depends on that fixed transport. |
| `log_hydro_dmc` | `true` | Humidity climatology zeroed at init; `hydro!` (latent heat, evaporation, rain) disabled entirely - removes the hydrological cycle described in Hydro19 from the mean climate. |
| `log_atmos_dmc` | `true` | Master switch: disables sensible heat flux, freezes the downward-longwave feedback, and also gates circulation and hydrology - effectively decouples the atmosphere from the surface. |
| `log_co2_dmc` | `true` | Control-run CO₂ forced to 0 ppm - removes the greenhouse-gas baseline from the mean climate (DF11 §3.2, Eq. 4). |
| `log_ocean_dmc` | `true` | Master ocean-coupling switch: surface heat capacity uses the land value everywhere (no mixed layer), sea-ice heat-capacity blending is skipped, and the deep ocean contributes zero heat exchange (DF11 §3.7–3.8) - removes ocean heat storage/transport from the mean climate. |
| `log_qflux_dmc` | `true` | Governs whether the empirical flux-correction terms (the residual added to force the model's climatology to match observations, DF11 §3.8) are computed via spin-up, loaded from a precomputed file, or zeroed, depending on its combination with `log_topo_drsp`. |

## CO₂-Response Deconstruction Switches (`_drsp`)

These switches, forwarded as keywords by `create_experiment_config(:decon_2xco2; ...)`,
keep the control/mean climate fixed and instead disable a process's
contribution to the *change* under a CO₂ doubling - isolating each
process's role in climate sensitivity, since feedback strength is
mean-state dependent and a fair "no-feedback" test must not also change the
climate the feedback acts on (MSCM §2.1–2.2). `:decon_2xco2` always applies
the CO₂ doubling itself as the forcing (`co2_concentration = 680.0`); there
is no `log_co2_drsp` because the response experiments' entire purpose is to
measure feedbacks *to* that forcing, not to disable it.

| Switch | Default | Effect when `false` |
|:-------|:--------|:---------------------|
| `log_clouds_drsp` | `true` | Cloud cover frozen at a constant 0.7 rather than the seasonal climatology - removes the cloud contribution to the CO₂ response while keeping *some* cloud cover present (MSCM p. 2161). |
| `log_crcl_drsp` | `true` | Removes circulation's contribution to the forced response (combined with `log_crcl_dmc` to fully gate `circulation!`). |
| `log_hydro_drsp` | `true` | Removes hydrology's contribution to the forced response (combined with `log_hydro_dmc` to fully gate `hydro!`). |
| `log_topo_drsp` | `true` | Topography capped at 1.0 m ("constant topography") - removes topographic modulation of the response; also interacts with `log_qflux_dmc` in the flux-correction branching. |
| `log_humid_drsp` | `true` | Humidity frozen at a constant 0.0052 kg/kg - removes the water-vapor feedback to CO₂ (DF11 §3.2, Eq. 5) specifically, while `log_hydro_drsp` covers the rest of the hydrological cycle's response. |
| `log_ocean_drsp` | `true` | Mixed-layer depth held constant (no seasonal variation) and the deep ocean contributes zero response tendency - removes the ocean heat-uptake response to CO₂ forcing (DF11 §3.7). |

## Circulation Components

GREB's atmospheric transport uses a fixed climatological wind field - it
does not itself respond to forcing (DF11 §3.5). Diffusion approximates the
smearing effect of weather disturbances with a roughly one-week lifetime;
advection moves heat/moisture with the mean wind; both are scaled by
topography. The vertical terms additionally let the model represent
convergence-driven precipitation (e.g. the ITCZ) from the prescribed
vertical-velocity field under continuity/hydrostatic assumptions
(Hydro19 §3.1, §3.3, Eqs. 11, 17–18).

| Switch | Default | Effect when `false` |
|:-------|:--------|:---------------------|
| `log_ice` | `true` | Disables the ice-albedo feedback: albedo becomes a constant value regardless of surface temperature, and heat capacity ignores ice fraction (DF11 §3.6, Fig. 3). GREB has no explicit sea-ice model - albedo and heat capacity near freezing are simple functions of surface temperature, which is what this switch removes. |
| `log_hdif` | `true` | Disables horizontal (heat) diffusion. |
| `log_hadv` | `true` | Disables horizontal (heat) advection. |
| `log_vdif` | `true` | Disables vertical/meridional diffusion of water vapor. |
| `log_vadv` | `true` | Disables advection of water vapor. |
| `log_conv`[^2] | `true` | Disables the moisture-convergence sub-step (water vapor only). Not exposed via `create_experiment_config` - set `cfg.log_conv` directly. |

## Hydrology Parameterization

Evaporation in GREB is a bulk formula on the saturation deficit, wind speed,
and a *prescribed* (climatological, not dynamically simulated) soil-wetness
fraction. Precipitation scales with relative humidity and vertical velocity
(both mean and variability), which is what lets the model reproduce the
ITCZ and midlatitude storm tracks rather than only a diffuse pattern
(Hydro19, Abstract & §3.1–3.3, Eqs. 11–18).

**`log_rain`** - selects the precipitation regression coefficients:

| Value | Meaning |
|:------|:--------|
| `-1` | Original GREB |
| `0` (default) | Best GREB fit (ERA-Interim) |
| `1` | + relative humidity |
| `2` | + omega (vertical velocity) convergence |
| `3` | + relative humidity & omega |

When `log_rain == 1`, a rain-rate limiter is additionally applied.

**`log_eva`** - selects the evaporation/exchange-coefficient parameterization:

| Value | Meaning |
|:------|:--------|
| `-1` (default) | Original GREB (climatological wind + fixed gust term) |
| `0` | `Ts`-derived skin temperature + climatological wind speed + land/ocean exchange coefficients |
| `1` | Modified gust/coefficient variant of `-1` |
| `2` | Modified variant of `0` with different gust/coefficients |

**`log_clim`** - NCEP-vs-ERA coefficient override:

| Value | Meaning |
|:------|:--------|
| `0` (default) | No override |
| `1` | If also `log_rain == 0`, overrides the precipitation coefficients to a hardcoded NCEP-tuned set - independent of which climatology *files* were loaded via `load_greb_jld2!` |

> **Note:** `log_rain == 0 && log_clim == 1` is the only combination where
> `log_clim` has any effect - it's easy to set `log_clim = 1` and see no
> change if `log_rain` isn't also `0`.

## External-Forcing Gates

**These three switches only matter for `:elnino`, `:lanina`, and `:rcp85`.**
They gate whether the corresponding anomaly field is *loaded from disk*;
the anomaly is then applied unconditionally in `init_model!` based on
`cfg.experiment` alone, not on these switches. `create_experiment_config`
already sets all three to `true` for those three experiments - you should
not normally need to touch them directly.

| Switch | Default | Effect when `true` |
|:-------|:--------|:---------------------|
| `log_tsurf_ext` | `false` | Loads the surface-temperature anomaly forcing field (RCP8.5/ENSO). |
| `log_hwind_ext` | `false` | Loads the horizontal wind (u/v) anomaly field. |
| `log_omega_ext` | `false` | Loads the vertical-velocity (omega) anomaly field. |

## Experiment-Level Fields

Scenario parameters rather than physics switches - set via
`create_experiment_config` keywords or directly on the returned `cfg`:

| Field | Type | Default | Purpose |
|:------|:-----|:--------|:--------|
| `experiment` | `Symbol` | `:full_model` | Selects the experiment branch - see [Experiment Presets](#experiment-presets) below |
| `co2_concentration` | `Float64` | `340.0` | Static CO₂ (ppm) for experiments without a time-varying scenario table |
| `orbital_index` | `Int` | `0` | Row index into the `solar_scenarios` table for `:obliquity`/`:eccentricity` |
| `earth_sun_distance_pct` | `Float64` | `0.0` | Percent change in orbital radius for `:earth_sun_distance` |
| `co2_scenario` | `Dict{Int,Float64}` | `Dict()` | Year→ppm lookup, auto-populated for IPCC RCP/SSP/historical/custom-CO2 experiments |
| `custom_co2_path` | `String` | `""` | User-supplied "year CO2" text file path for `:custom_co2` |

`c_q`, `c_rq`, `c_omega`, and `c_omegastd` also live on `PhysicsConfig` but
are not meant to be set manually - `set_hydrology_parameters!` derives them
from `log_rain`/`log_clim` at the start of every run.

## Experiment Presets

What each `create_experiment_config` preset changes relative to `:full_model`:

| Preset | vs. `:full_model` |
|:-------|:-------------------|
| `:full_model` | Nothing - the baseline |
| `:constant_topo` | `log_topo_drsp = false` |
| `:co2_double` | `co2_concentration = 680.0` |
| `:co2_quadruple` | `co2_concentration = 1360.0` |
| `:solar_plus27` | Solar constant +27 W/m² |
| `:elnino` / `:lanina` | The three `log_*_ext` switches set `true`; adds/subtracts ERA-Interim ENSO anomalies |
| `:paleo_231kyr` | `co2_concentration = 200.0`; paleo solar-forcing table |
| `:rcp85` | `log_*_ext` switches `true`; loads CMIP5 RCP8.5 anomaly fields |
| `:rcp26` / `:rcp45` / `:rcp60` / `:ssp119` / `:ssp126` / `:ssp245` / `:ssp460` / `:ssp585` | No switches change; loads a year→CO2 lookup table |
| `:historical_co2` | Year counter starts at 1850; loads the observed CO₂ record |
| `:custom_co2` | `custom_co2_path` set from the `co2_path` keyword |
| `:decon_mean_climate` | Exposes `log_clouds_dmc`, `log_ocean_dmc`, `log_atmos_dmc`, `log_co2_dmc`, `log_hydro_dmc`, `log_qflux_dmc`, `log_ice`, `log_hdif`, `log_hadv`, `log_vdif`, `log_vadv` as keywords (all default `true`, i.e. behaves like `:full_model` unless overridden) |
| `:decon_2xco2` | `co2_concentration = 680.0`; exposes `log_topo_drsp`, `log_clouds_drsp`, `log_humid_drsp`, `log_ocean_drsp`, `log_hydro_drsp`, `log_ice`, `log_hdif`, `log_hadv`, `log_vdif`, `log_vadv` as keywords |

### Advanced / manually-constructed experiments

The following experiment symbols are dispatched in `src/model.jl`/`src/tendencies.jl`
but are not wrapped by `create_experiment_config` - construct them with
`PhysicsConfig(experiment=:symbol, ...)` directly and read the corresponding
branch in source for exact behavior:

| Category | Symbols |
|:---------|:--------|
| CO₂ variants | `:a1b_scenario`, `:a1b_enhanced`, `:co2_10x`, `:co2_half`, `:co2_zero`, `:co2_sine_wave`, `:co2_step` |
| Orbital/paleo | `:solar_cycle_11yr`, `:paleo_solar_modern_co2`, `:modern_solar_paleo_co2`, `:obliquity`, `:eccentricity`, `:earth_sun_distance` |
| Other | `:sst_plus1` |
| Regional CO₂ (static mask) | `:regional_co2_nh`, `:regional_co2_sh`, `:regional_co2_tropics`, `:regional_co2_extratropics` |
| Regional CO₂ (dynamic/seasonal mask) | `:regional_co2_ocean`, `:regional_co2_land_ice`, `:regional_co2_winter`, `:regional_co2_summer` |

## Known Limitations

- **`log_vapor_dmc`** is declared on `PhysicsConfig` but is never read
  anywhere in `src/` - setting it currently has no effect on the model.
- **`log_conv`** works like the other circulation-component switches
  (gates the moisture-convergence sub-step) but, unlike them, is not
  exposed as a keyword by `create_experiment_config` - it must be set with
  `cfg.log_conv = ...` after construction.

## References

Abbreviations used above: **DF11** = Dommenget & Flöter (2011); **Hydro19**
= Stassen, Dommenget & Loveday (2019); **MSCM** = Dommenget et al. (2019).

1. **DF11** - Dommenget, D., and Flöter, J. (2011). Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model. *Journal of Climate Dynamics*, 37: 2143. [doi:10.1007/s00382-011-1026-0](https://doi.org/10.1007/s00382-011-1026-0)
2. **Hydro19** - Stassen, C., Dommenget, D., and Loveday, N. (2019). A hydrological cycle model for the Globally Resolved Energy Balance (GREB) model v1.0. *Geoscientific Model Development*, 12, 425-440. [doi:10.5194/gmd-12-425-2019](https://doi.org/10.5194/gmd-12-425-2019)
3. **MSCM** - Dommenget, D., Nice, K., Bayr, T., Kasang, D., Stassen, C., and Rezny, M. The Monash Simple Climate Model Experiments: An interactive database of the mean climate, climate change and scenarios simulations. *Geoscientific Model Development*, 12, 2155-2179. [doi:10.5194/gmd-12-2155-2019](https://doi.org/10.5194/gmd-12-2155-2019)

See the repository [README](https://github.com/EnvDroneSense/GREB-julia#references)'s References section for the original GREB model homepage link.
