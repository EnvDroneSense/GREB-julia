begin
    # - Grid dimensions ──────────────────────────────────
    "Number of longitude grid points."
    const xdim = 96
    "Number of latitude grid points."
    const ydim = 48
    const dlon = 360.0f0 / xdim                     # longitude spacing [degrees]
    const dlat = 180.0f0 / ydim                     # latitude spacing  [degrees]

    # - Time stepping ────────────────────────────────────
    const ndays_yr = 365                            # days per year (no leap years)
    const Δt = 12.0f0 * 3600.0f0                    # main time step [s] (12 hours)
    const Δt_crcl = 1800.0f0                        # circulation sub-time step [s] (30 min)
    const ndt_days = Int(round(24 * 3600 / Δt))     # time steps per day
    "Time steps per year (`ndays_yr * ndt_days` = 730)."
    const nstep_yr = Int(ndays_yr * ndt_days)
    const ntime = max(1, Int(round(Δt / Δt_crcl)))  # Number of sub-steps within one main time step

    # - Calendar constants ───────────────────────────────
    const cjday_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    const jday_mon_cumsum = cumsum(cjday_mon)

    # - Physical Constants & Numerical Limits ────────────
    const min_T_K = 273.15f0 - 40.0f0   # -40°C, minimum allowed surface/air temperature [K]
    const max_humidity_change = 0.020f0 # Maximum humidity increment [kg/kg]
    const min_humidity_change = 0.9f0   # Fraction of humidity that can be removed
end;

begin
    # - Optimized Hydrology Parameter Lookup Table ────────────────────────
    const HYDRO_PARAMS = Dict(
        -1 => (1.0f0, 0.0f0, 0.0f0, 0.0f0),                         # Original GREB
        1 => (-1.391649f0, 3.018774f0, 0.0f0, 0.0f0),               # +Relative humidity
        2 => (0.862162f0, 0.0f0, -29.02096f0, 0.0f0),               # +Omega convergence
        3 => (-0.2685845f0, 1.4591853f0, -26.9858807f0, 0.0f0),     # +RH & Omega
        0 => (-1.88f0, 2.25f0, -17.69f0, 59.07f0)                   # Best GREB (ERA-Interim)
    )
end;

begin
    # ── Natural constants ────────────────────────────────────────────
    const const_pi = Float32(pi)   # π (model precision)
    const σ = 5.6704f-8            # Stefan-Boltzmann constant [W/m²/K⁴]
    const ρ_ocean = 999.1f0        # density of water at T=15°C [kg/m³]
    const ρ_land = 2600.0f0        # density of solid rock [kg/m³]
    const ρ_air = 1.2f0            # density of air at 20°C at sea
    const grav = 9.80665f0         # gravitational acceleration [m/s²]
    const cp_ocean = 4186.0f0      # specific heat of water at T=15°C [J/kg/K]
    const cp_land = cp_ocean / 4.5f0 # specific heat of dry land [J/kg/K]
    const cp_air = 1005.0f0        # specific heat of air [J/kg/K]
    const ε = 1.0f0                # emissivity for IR
end;

begin
    # ── Column depths [m] ───────────────────────────────────────────
    const d_ocean = 50.0f0                      # ocean column
    const d_land = 2.0f0                        # land column
    const d_air = 5000.0f0                      # air column

    # ── Heat capacities [J/K/m²] ────────────────────────────────────
    const cap_ocean = cp_ocean * ρ_ocean        # 1m ocean
    const cap_land = cp_land * ρ_land * d_land  # land column
    const cap_air = cp_air * ρ_air * d_air      # air column

    # ── Sensible heat [W/K/m²] ──────────────────────────────────────
    const ct_sens = 22.5f0                       # sensible heat coupling

    # ── Albedo parameters ───────────────────────────────────────────
    const da_ice = 0.25f0                        # albedo increase for ice-cover
    const a_no_ice = 0.1f0                       # albedo ice-free
    const a_cloud = 0.35f0                       # cloud albedo

    # ── Ice/snow temperature thresholds [K] ─────────────────────────
    const Tl_ice1 = 273.15f0 - 10.0f0            # land: full ice albedo
    const Tl_ice2 = 273.15f0                     # land: no ice albedo
    const To_ice1 = 273.15f0 - 7.0f0             # ocean: full ice
    const To_ice2 = 273.15f0 - 1.7f0             # ocean: no ice albedo

    # Precomputed inverse ranges (avoids division in hot loops)
    const inv_To_ice_range = 1.0f0 / (To_ice2 - To_ice1)
    const inv_Tl_ice_range = 1.0f0 / (Tl_ice2 - Tl_ice1)

    # ── Deep ocean ──────────────────────────────────────────────────
    const co_turb = 5.0f0                        # turbulent mixing coefficient [W/K/m²]
    const c_effmix = 0.5f0                       # mixing efficiency
    const turb_coeff = Δt * co_turb / cap_ocean  # precomputed mixing coefficient

    # ── Atmospheric transport ───────────────────────────────────────
    const κ = 8f5                          # diffusion coefficient [m²/s]

    # ── Latent heat / hydrology ─────────────────────────────────────
    const ce = 2f-3                        # latent heat transfer coefficient
    const cq_latent = 2.257f6              # latent heat of evaporation [J/kg]
    const cq_rain = -0.1f0 / 24.0f0 / 3600.0f0   # rain-related vapor decrease [1/s]

    # ── Scaling heights [m] ─────────────────────────────────────────
    const z_air = 8400.0f0                 # heat & CO₂ scaling height
    const z_vapor = 5000.0f0               # water vapor scaling height
    const const_factor = Δt_crcl / z_vapor * 2.5f0 / (ρ_air * grav)

    # ── Regression factor [kg/m³] ───────────────────────────────────
    const r_qviwv = 2.6736f3               # VIWV ↔ q_air regression factor
    const conv_factor = r_qviwv * 86400.0f0  # kg/kg → mm/day conversion

    # ── solar factor [%] ────────────────────────────────────────────
    const S0_var = 100.0f0        # default 100%

    # ── transport geometry & coefficients ───────────────────────────
    const deg_grid = 2.0f0 * const_pi * 6.371f6 / 360.0f0
    const dyy_grid = dlat * deg_grid
    const lat_grid = Float32[dlat * k - dlat / 2.0f0 - 90.0f0 for k in 1:ydim]
    const dxlat_grid = Float32[dlon * deg_grid * cos(2.0f0 * const_pi / 360.0f0 * lat_grid[k]) for k in 1:ydim]

    # ── Diffusion coefficients ──────────────────────────────────────
    const ccy_diff = κ * Δt_crcl / dyy_grid^2
    const ccx_diff = Float32[κ * Δt_crcl / dxlat_grid[k]^2 for k in 1:ydim]

    # ── Advection coefficients ──────────────────────────────────────
    const ccy_adv = Δt_crcl / dyy_grid / 2.0f0
    const ccx_adv = Float32[Δt_crcl / dxlat_grid[k] / 2.0f0 for k in 1:ydim]

    # ── periodic longitude neighbour indices ───────────────────────
    const lon_jm1 = Int32[mod1(i-1, xdim) for i in 1:xdim]
    const lon_jp1 = Int32[mod1(i+1, xdim) for i in 1:xdim]
    const lon_jm2 = Int32[mod1(i-2, xdim) for i in 1:xdim]
    const lon_jp2 = Int32[mod1(i+2, xdim) for i in 1:xdim]
    const lon_jm3 = Int32[mod1(i-3, xdim) for i in 1:xdim]
    const lon_jp3 = Int32[mod1(i+3, xdim) for i in 1:xdim]
end

const ΔT_AIR_FACTOR = Δt / cap_air

const p_emi = Float32[9.0721, 106.7252, 61.5562, 0.0179, 0.0028,
                      0.0570, 0.3462, 2.3406, 0.7032, 1.0662]

begin
    # Precomputed Calendar Lookup
    const max_timesteps = 200 * nstep_yr  # Support up to 200-year runs
    const calendar_lookup = [(
        day=mod((it - 1) ÷ ndt_days, ndays_yr) + 1,
        step=mod(it - 1, nstep_yr) + 1
    ) for it in 1:max_timesteps]

    const polar_treshold = 2.5f5  # 250 km in meters
    const IS_POLAR = [dxlat_grid[k] <= polar_treshold for k in 1:ydim]
end;

begin
    # ── Polar sub-stepping constants (diffusion!/advection!) ──────────────
    function _polar_diff_step(k)
        dd = max(1, round(Int, Δt_crcl / (dxlat_grid[k]^2 / κ)))
        dtdff2 = Δt_crcl / dd
        time2 = max(1, round(Int, Δt_crcl / dtdff2))
        return (time2=time2, ccx2=κ * dtdff2 / dxlat_grid[k]^2)
    end
    function _polar_adv_step(k)
        dd = max(1, round(Int, Δt_crcl / (dxlat_grid[k] / 10.0f0)))
        dtdff2 = Δt_crcl / dd
        time2 = max(1, round(Int, Δt_crcl / dtdff2))
        return (time2=time2, ccx2=dtdff2 / dxlat_grid[k] / 2.0f0)
    end

    const POLAR_DIFF_TIME2 = [_polar_diff_step(k).time2 for k in 1:ydim]
    const POLAR_DIFF_CCX2 = Float32[_polar_diff_step(k).ccx2 for k in 1:ydim]
    const POLAR_ADV_TIME2 = [_polar_adv_step(k).time2 for k in 1:ydim]
    const POLAR_ADV_CCX2 = Float32[_polar_adv_step(k).ccx2 for k in 1:ydim]
end;
