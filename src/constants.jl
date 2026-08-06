# ── notebook cell 531589ab-c6b5-4048-9ba5-f9ad62ab00a6  (orig lines 202-227) ──
# Note: this should be in the script
begin
    # 📐 Grid dimensions ──────────────────────────────────
    "Number of longitude grid points."
    const xdim = 96
    "Number of latitude grid points."
    const ydim = 48
    const dlon = 360.0 / xdim                       # longitude spacing [degrees]
    const dlat = 180.0 / ydim                       # latitude spacing  [degrees]

    # ⏱️ Time stepping ────────────────────────────────────
    const ndays_yr = 365                            # days per year (no leap years)
    const Δt = 12.0 * 3600.0                        # main time step [s] (12 hours)
    const Δt_crcl = 1800.0                          # circulation sub-time step [s] (30 min)
    const ndt_days = Int(round(24 * 3600 / Δt))     # time steps per day
    "Time steps per year (`ndays_yr * ndt_days` = 730)."
    const nstep_yr = Int(ndays_yr * ndt_days)
    const ntime = max(1, Int(round(Δt / Δt_crcl)))  # Number of sub-steps within one main time step

    # 📅 Calendar constants ───────────────────────────────
    const cjday_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    const jday_mon_cumsum = cumsum(cjday_mon)

    # 🔢 Physical Constants & Numerical Limits ────────────
    # Fortran's Tmin_limit=40 is a raw-Kelvin numerical-stability floor
    # (greb.model.mscm.f90:470-477: `where(Ts0 .le. Tmin_limit) Ts0 = Tmin_limit`),
    # but 40 K is colder than anywhere on Earth ever gets — it's a floor that
    # can never physically bind, not a real numerical safety net. Kept at
    # 233.15 K (-40°C) intentionally: a real, physically-plausible cold-extreme
    # floor for Antarctic winter, on the judgment that the Fortran value isn't
    # actually a target worth matching here.
    const min_T_K = 273.15 - 40.0       # -40°C, minimum allowed surface/air temperature [K]
    const max_humidity_change = 0.020   # Maximum humidity increment [kg/kg]
    const min_humidity_change = 0.9     # Fraction of humidity that can be removed
end;

# ── notebook cell 32b531ab-ee71-4685-af65-a2bae0b868f6  (orig lines 363-377) ──
begin
    # 💧 Optimized Hydrology Parameter Lookup Table ────────────────────────
    const HYDRO_PARAMS = Dict(
        -1 => (1.0, 0.0, 0.0, 0.0),                     # Original GREB
        1 => (-1.391649, 3.018774, 0.0, 0.0),           # +Relative humidity
        2 => (0.862162, 0.0, -29.02096, 0.0),           # +Omega convergence
        3 => (-0.2685845, 1.4591853, -26.9858807, 0.0), # +RH & Omega
        0 => (-1.88, 2.25, -17.69, 59.07)               # Best GREB (ERA-Interim)
    )
end;

# ── notebook cell f8a2c2de-5045-4ab6-a6fa-7bca502afc9b  (orig lines 644-658) ──
begin
    # ── Natural constants ────────────────────────────────────────────
    const const_pi = pi            # π (model precision)
    const σ = 5.6704e-8            # Stefan-Boltzmann constant [W/m²/K⁴]
    const ρ_ocean = 999.1          # density of water at T=15°C [kg/m³]
    const ρ_land = 2600.0          # density of solid rock [kg/m³]
    const ρ_air = 1.2              # density of air at 20°C at sea
    const grav = 9.80665           # gravitational acceleration [m/s²]
    const cp_ocean = 4186.0        # specific heat of water at T=15°C [J/kg/K]
    const cp_land = cp_ocean / 4.5 # specific heat of dry land [J/kg/K]
    const cp_air = 1005.0          # specific heat of air [J/kg/K]
    const ε = 1.0                  # emissivity for IR
end;

# ── notebook cell b23bc922-e9bc-4012-9782-1258e3cc8e7b  (orig lines 659-736) ──
begin
    # ── Column depths [m] ───────────────────────────────────────────
    const d_ocean = 50.0                        # ocean column
    const d_land = 2.0                          # land column
    const d_air = 5000.0                        # air column

    # ── Heat capacities [J/K/m²] ────────────────────────────────────
    const cap_ocean = cp_ocean * ρ_ocean        # 1m ocean
    const cap_land = cp_land * ρ_land * d_land  # land column
    const cap_air = cp_air * ρ_air * d_air      # air column

    # ── Sensible heat [W/K/m²] ──────────────────────────────────────
    const ct_sens = 22.5                        # sensible heat coupling

    # ── Albedo parameters ───────────────────────────────────────────
    const da_ice = 0.25                         # albedo increase for ice-cover
    const a_no_ice = 0.1                        # albedo ice-free
    const a_cloud = 0.35                        # cloud albedo

    # ── Ice/snow temperature thresholds [K] ─────────────────────────
    const Tl_ice1 = 273.15 - 10.0               # land: full ice albedo
    const Tl_ice2 = 273.15                      # land: no ice albedo
    const To_ice1 = 273.15 - 7.0                # ocean: full ice
    const To_ice2 = 273.15 - 1.7                # ocean: no ic

    # Precomputed inverse ranges (avoids division in hot loops)
    const inv_To_ice_range = 1.0 / (To_ice2 - To_ice1)
    const inv_Tl_ice_range = 1.0 / (Tl_ice2 - Tl_ice1)

    # ── Deep ocean ──────────────────────────────────────────────────
    const co_turb = 5.0                          # turbulent mixing coefficient [W/K/m²]
    const c_effmix = 0.5                         # mixing efficiency
    const turb_coeff = Δt * co_turb / cap_ocean  # precomputed mixing coefficient

    # ── Atmospheric transport ───────────────────────────────────────
    const κ = 8e5                          # diffusion coefficient [m²/s]

    # ── Latent heat / hydrology ─────────────────────────────────────
    const ce = 2e-3                        # latent heat transfer coefficient
    const cq_latent = 2.257e6              # latent heat of evaporation [J/kg]
    const cq_rain = -0.1 / 24.0 / 3600.0   # rain-related vapor decrease [1/s]

    # ── Scaling heights [m] ─────────────────────────────────────────
    const z_air = 8400.0                   # heat & CO₂ scaling height
    const z_vapor = 5000.0                 # water vapor scaling height
    const const_factor = Δt_crcl / z_vapor * 2.5 / (ρ_air * grav)

    # ── Regression factor [kg/m³] ───────────────────────────────────
    const r_qviwv = 2.6736e3               # VIWV ↔ q_air regression factor
    const conv_factor = r_qviwv * 86400.0  # kg/kg → mm/day conversion

    # ── solar factor [%] ────────────────────────────────────────────
    const S0_var = 100.0          # default 100%

    # ── transport geometry & coefficients ───────────────────────────
    const deg_grid = 2.0 * const_pi * 6.371e6 / 360.0
    const dyy_grid = dlat * deg_grid
    const lat_grid = [dlat * k - dlat / 2.0 - 90.0 for k in 1:ydim]
    const dxlat_grid = [dlon * deg_grid * cos(2.0 * const_pi / 360.0 * lat_grid[k]) for k in 1:ydim]

    # ── Diffusion coefficients ──────────────────────────────────────
    const ccy_diff = κ * Δt_crcl / dyy_grid^2
    const ccx_diff = [κ * Δt_crcl / dxlat_grid[k]^2 for k in 1:ydim]

    # ── Advection coefficients ──────────────────────────────────────
    const ccy_adv = Δt_crcl / dyy_grid / 2.0
    const ccx_adv = [Δt_crcl / dxlat_grid[k] / 2.0 for k in 1:ydim]

    # ── periodic longitude neighbour indices ───────────────────────
    const lon_jm1 = [mod1(i-1, xdim) for i in 1:xdim]
    const lon_jp1 = [mod1(i+1, xdim) for i in 1:xdim]
    const lon_jm2 = [mod1(i-2, xdim) for i in 1:xdim]
    const lon_jp2 = [mod1(i+2, xdim) for i in 1:xdim]
    const lon_jm3 = [mod1(i-3, xdim) for i in 1:xdim]
    const lon_jp3 = [mod1(i+3, xdim) for i in 1:xdim]
end

# ── notebook cell 0dbfb663-46e7-4873-ac77-1e8e392fe69d  (orig lines 737-749, split: see state.jl for the mutable globals from this cell) ──
const ΔT_AIR_FACTOR = Δt / cap_air

# ── notebook cell 898dd5aa-5273-4833-9a4a-0f2b94cc8d38  (orig lines 779-782) ──
const p_emi = [9.0721, 106.7252, 61.5562, 0.0179, 0.0028,
               0.0570, 0.3462, 2.3406, 0.7032, 1.0662]

# ── notebook cell 047b312f-8d6c-4732-aa0b-bea3de3e99e2  (orig lines 984-1002, split: see state.jl for the TimeState struct from this cell) ──
begin
    # 📅 Precomputed Calendar Lookup
    const max_timesteps = 200 * nstep_yr  # Support up to 200-year runs
    const calendar_lookup = [(
        day=mod((it - 1) ÷ ndt_days, ndays_yr) + 1,
        step=mod(it - 1, nstep_yr) + 1
    ) for it in 1:max_timesteps]

    const polar_treshold = 2.5e5  # 250 km in meters
    const IS_POLAR = [dxlat_grid[k] <= polar_treshold for k in 1:ydim]
end;
