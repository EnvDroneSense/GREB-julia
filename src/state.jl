# ── notebook cell 80ac789e-4fe7-4184-9946-b8d7c24b04ea  (orig lines 408-544) ──
begin
    """
            CirculationWorkspace

        Pre-allocated buffers for diffusion, advection, and circulation calculations.
        Reused across all time steps to eliminate allocations.
        """
    mutable struct CirculationWorkspace
        # Polar sub-stepping buffers
        T1h::Vector{Float64}      # polar sub-stepping
        dTxh::Vector{Float64}     # polar tendencies

        # Circulation work arrays
        X_work::Matrix{Float64}   # circulation work array
        dX_diff::Matrix{Float64}  # diffusion output
        dX_adv::Matrix{Float64}   # advection output
        dX_conv::Matrix{Float64}  # convection output
        dX_crcl::Matrix{Float64}  # circulation output

        # Tendency buffers
        temp_buf::Matrix{Float64}   # general workspace
        Q_sens_buf::Matrix{Float64} # Sensible heat flux buffer
        eva::Matrix{Float64}        # dq_eva
        rain::Matrix{Float64}       # dq_rain
        crcl::Matrix{Float64}       # dq_crcl

        # State buffers
        Ts0_buf::Matrix{Float64}    # Surface temperature output
        Ta0_buf::Matrix{Float64}    # Air temperature output
        To0_buf::Matrix{Float64}    # Ocean temperature output
        q0_buf::Matrix{Float64}     # Humidity output

        # LW radiation buffers
        e_co2_buf::Matrix{Float64}    # spatial CO₂ buffer
        e_vapor_buf::Matrix{Float64}  # spatial water vapor buffer
        em_buf::Matrix{Float64}       # spatial emissivity buffer
        LW_surf_buf::Matrix{Float64}  # Surface longwave
        LW_down_buf::Matrix{Float64}  # Downwelling longwave
        LW_up_buf::Matrix{Float64}    # Upwelling longwave

        # Hydrology buffers
        qs::Matrix{Float64}        # Saturation humidity buffer
        Tskin::Matrix{Float64}     # Skin temperature buffer
        rq::Matrix{Float64}        # Relative humidity buffer
        ws_base::Matrix{Float64}   # Base wind speed buffer
        # Hydrology
        Q_lat_buf::Matrix{Float64}
        Q_lat_air_buf::Matrix{Float64}
        dq_eva_buf::Matrix{Float64}
        dq_rain_buf::Matrix{Float64}
        cE_buf::Matrix{Float64}       # Surface exchange coefficient buffer

        # Deep_ocean
        dT_ocean_buf::Matrix{Float64}
        dTo_buf::Matrix{Float64}

        # Dedicated circulation output
        dTa_crcl::Matrix{Float64}   # temperature tendency
        dq_crcl::Matrix{Float64}    # humidity tendency

        # SWradiation
        ice_cover_buf::Matrix{Float64} # ice fraction
        a_surf_buf::Matrix{Float64}    # surface albedo
        albedo_buf::Matrix{Float64}    # combined albedo (surface + atmosphere)
        a_atmos_buf::Matrix{Float64}   # atmospheric albedo
        sw_buf::Matrix{Float64}        # net shortwave flux

        # time_loop
        precip_out::Matrix{Float64}   # precipitation output
        evap_out::Matrix{Float64}     # evaporation output
        qcrcl_out::Matrix{Float64}    # circulation moisture output
        term_north::Vector{Float64}   # northern boundary term
        term_south::Vector{Float64}   # southern boundary term
    end

    function CirculationWorkspace()
        CirculationWorkspace(
            zeros(Float64, xdim),# T1h
            zeros(Float64, xdim),# dTxh
            zeros(Float64, xdim, ydim),# X_work
            zeros(Float64, xdim, ydim),# dX_diff
            zeros(Float64, xdim, ydim),# dX_adv
            zeros(Float64, xdim, ydim),# dX_conv
            zeros(Float64, xdim, ydim),# dX_crcl
            zeros(Float64, xdim, ydim),# temp_buf
            zeros(Float64, xdim, ydim),# Q_sens_buf
            zeros(Float64, xdim, ydim),# rain
            zeros(Float64, xdim, ydim),# eva
            zeros(Float64, xdim, ydim),# crcl
            # State buffers
            zeros(Float64, xdim, ydim),# Ts0_buf
            zeros(Float64, xdim, ydim),# Ta0_buf
            zeros(Float64, xdim, ydim),# To0_buf
            zeros(Float64, xdim, ydim),# q0_buf
            # LWradiation buffers
            zeros(Float64, xdim, ydim),# e_co2_buf
            zeros(Float64, xdim, ydim),# e_vapor_buf
            zeros(Float64, xdim, ydim),# em_buf
            zeros(Float64, xdim, ydim),# LW_surf_buf
            zeros(Float64, xdim, ydim),# LW_down_buf
            zeros(Float64, xdim, ydim),# LW_up_buf
            # Hydrology buffers
            zeros(Float64, xdim, ydim),# qs
            zeros(Float64, xdim, ydim),# Tskin
            zeros(Float64, xdim, ydim),# rq
            zeros(Float64, xdim, ydim),# ws_base
            # Hydrology
            zeros(Float64, xdim, ydim),  # Q_lat_buf
            zeros(Float64, xdim, ydim),  # Q_lat_air_buf
            zeros(Float64, xdim, ydim),  # dq_eva_buf
            zeros(Float64, xdim, ydim),  # dq_rain_buf
            zeros(Float64, xdim, ydim),  # cE_buf
            # Deep_ocean
            zeros(Float64, xdim, ydim),  # dT_ocean_buf
            zeros(Float64, xdim, ydim),  # dTo_buf
            # Dedicated circulation output
            zeros(Float64, xdim, ydim),  # dTa_crcl
            zeros(Float64, xdim, ydim),  # dq_crcl
            # SWradiation buffers
            zeros(Float64, xdim, ydim),  # ice_cover_buf
            zeros(Float64, xdim, ydim),  # a_surf_buf
            zeros(Float64, xdim, ydim),  # albedo_buf
            zeros(Float64, xdim, ydim),  # a_atmos_buf
            zeros(Float64, xdim, ydim),  # sw_buf
            # time_loop
            zeros(Float64, xdim, ydim),  # precip_out
            zeros(Float64, xdim, ydim),  # evap_out
            zeros(Float64, xdim, ydim),  # qcrcl_out
            zeros(Float64, xdim),  # term_north
            zeros(Float64, xdim),  # term_south
        )
    end
end;

# ── notebook cell b9f7bde9-0aa4-4075-8fb9-14f84db0b3fa  (orig lines 545-623) ──
begin
    """
            MonthlyAccumulator

        Accumulates fields over a month for monthly-mean output.
        Reset after each month via `reset!`.
        """
    mutable struct MonthlyAccumulator
        Tmm::Matrix{Float64}          # Surface temperature accumulator
        Tamm::Matrix{Float64}         # Air temperature accumulator
        Tomm::Matrix{Float64}         # Ocean temperature accumulator
        qmm::Matrix{Float64}          # Humidity accumulator
        apmm::Matrix{Float64}         # Albedo accumulator
        icemm::Matrix{Float64}        # Ice fraction accumulator
        precipmm::Matrix{Float64}     # Precipitation accumulator
        evapmm::Matrix{Float64}       # Evaporation accumulator
        qcrclmm::Matrix{Float64}      # Circulation moisture accumulator
        swmm::Matrix{Float64}         # Shortwave radiation accumulator
        lwmm::Matrix{Float64}         # Longwave radiation accumulator
        qlatmm::Matrix{Float64}       # Latent heat accumulator
        qsensmm::Matrix{Float64}      # Sensible heat accumulator
        count::Int                    # Number of accumulations
    end

    function MonthlyAccumulator()
        MonthlyAccumulator(
            zeros(Float64, xdim, ydim),  # Tmm
            zeros(Float64, xdim, ydim),  # Tamm
            zeros(Float64, xdim, ydim),  # Tomm
            zeros(Float64, xdim, ydim),  # qmm
            zeros(Float64, xdim, ydim),  # apmm
            zeros(Float64, xdim, ydim),  # icemm
            zeros(Float64, xdim, ydim),  # precipmm
            zeros(Float64, xdim, ydim),  # evapmm
            zeros(Float64, xdim, ydim),  # qcrclmm
            zeros(Float64, xdim, ydim),  # swmm
            zeros(Float64, xdim, ydim),  # lwmm
            zeros(Float64, xdim, ydim),  # qlatmm
            zeros(Float64, xdim, ydim),  # qsensmm
            0
        )
    end

    function reset!(acc::MonthlyAccumulator)
        fill!(acc.Tmm, 0.0)
        fill!(acc.Tamm, 0.0)
        fill!(acc.Tomm, 0.0)
        fill!(acc.qmm, 0.0)
        fill!(acc.apmm, 0.0)
        fill!(acc.icemm, 0.0)
        fill!(acc.precipmm, 0.0)
        fill!(acc.evapmm, 0.0)
        fill!(acc.qcrclmm, 0.0)
        fill!(acc.swmm, 0.0)
        fill!(acc.lwmm, 0.0)
        fill!(acc.qlatmm, 0.0)
        fill!(acc.qsensmm, 0.0)
        acc.count = 0
    end

    function accumulate!(acc::MonthlyAccumulator, Ts, Ta, To, q, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens)
        @. acc.Tmm += Ts
        @. acc.Tamm += Ta
        @. acc.Tomm += To
        @. acc.qmm += q
        @. acc.apmm += albedo
        @. acc.icemm += ice
        @. acc.precipmm += precip
        @. acc.evapmm += evap
        @. acc.qcrclmm += qcrcl
        @. acc.swmm += sw
        @. acc.lwmm += lw
        @. acc.qlatmm += qlat
        @. acc.qsensmm += qsens
        acc.count += 1
    end
end;

# ── notebook cell f39520b7-a246-4980-b8b9-0215367d0b46  (orig lines 630-637) ──
begin
    # 🌍 Spatial CO₂ masking arrays ────────────────────────────────
    # Spatial fraction for regional CO₂ experiments
    co2_part = ones(Float64, xdim, ydim)    # Regional CO₂ mask (1.0 = full CO₂, 0.5 = half CO₂)
end;

# ── notebook cell 0dbfb663-46e7-4873-ac77-1e8e392fe69d  (orig lines 737-749, split: see constants.jl for the ΔT_AIR_FACTOR const from this cell) ──
begin
    # ☀️ Solar forcing storage
    global sw_solar_forcing_state = Ref(1.0)  # Runtime solar multiplier used by SWradiation!

    # 🔧 Flux correction arrays (initialised with zeros, overwritten if files exist)
    global TF_correct = zeros(Float64, xdim, ydim, nstep_yr)
    global qF_correct = zeros(Float64, xdim, ydim, nstep_yr)
    global ToF_correct = zeros(Float64, xdim, ydim, nstep_yr)
end;

# ── notebook cell 75f3b78f-924a-4a65-b2e3-a79c6f2082f9  (orig lines 814-824) ──
begin
    # 🗺️ 2D fields (xdim, ydim) ──────────────────────────────────────
    z_topo = zeros(Float64, xdim, ydim)  # topography [m] (<0: ocean)
    glacier = zeros(Float64, xdim, ydim)  # glacier mask (>0.5: glacier)
    z_ocean = zeros(Float64, xdim, ydim)  # derived ocean depth [m]
    cap_surf = zeros(Float64, xdim, ydim)  # surface heat capacity [J/K/m²]
    wz_air = zeros(Float64, xdim, ydim)  # exp(-z_topo / z_air)
    wz_vapor = zeros(Float64, xdim, ydim)  # exp(-z_topo / z_vapor)
end;

# ── notebook cell 0be9bc61-1a59-4dfe-84f9-bb1a27ca30fc  (orig lines 825-869) ──
begin
    # 🌡️ 3D climate fields (xdim, ydim, nstep_yr) ───────────────────
    Tclim = zeros(Float64, xdim, ydim, nstep_yr)   # surface temperature [K]
    uclim = zeros(Float64, xdim, ydim, nstep_yr)   # zonal wind [m/s]
    vclim = zeros(Float64, xdim, ydim, nstep_yr)   # meridional wind [m/s]
    qclim = zeros(Float64, xdim, ydim, nstep_yr)   # atmospheric humidity [kg/kg]
    mldclim = zeros(Float64, xdim, ydim, nstep_yr)   # mixed-layer depth [m]

    # additional climatology fields
    omegaclim = zeros(Float64, xdim, ydim, nstep_yr) # vertical velocity [Pa/s]
    omegastdclim = zeros(Float64, xdim, ydim, nstep_yr) # omega std deviation [Pa/s]
    wsclim = zeros(Float64, xdim, ydim, nstep_yr) # wind speed [m/s]

    # 📊 Anomaly Fields for ENSO/Climate Change Experiments ────────────────
    # ENSO anomaly fields
    Tclim_anom_enso = zeros(Float64, xdim, ydim, nstep_yr) # surface temperature [K]
    uclim_anom_enso = zeros(Float64, xdim, ydim, nstep_yr) # zonal wind [m/s]
    vclim_anom_enso = zeros(Float64, xdim, ydim, nstep_yr) # meridional wind [m/s]
    omegaclim_anom_enso = zeros(Float64, xdim, ydim, nstep_yr) # vertical velocity [Pa/s]
    wsclim_anom_enso = zeros(Float64, xdim, ydim, nstep_yr) # wind speed [m/s]

    # Climate change anomaly fields
    Tclim_anom_cc = zeros(Float64, xdim, ydim, nstep_yr) # surface temperature [K]
    uclim_anom_cc = zeros(Float64, xdim, ydim, nstep_yr) # zonal wind [m/s]
    vclim_anom_cc = zeros(Float64, xdim, ydim, nstep_yr) # meridional wind [m/s]
    omegaclim_anom_cc = zeros(Float64, xdim, ydim, nstep_yr) # vertical velocity [Pa/s]
    wsclim_anom_cc = zeros(Float64, xdim, ydim, nstep_yr) # wind speed [m/s]

    # 🌬️ Precomputed wind sign splits ──────────────────
    uclim_m = zeros(Float64, xdim, ydim, nstep_yr)   # negative u components
    uclim_p = zeros(Float64, xdim, ydim, nstep_yr)   # positive u components
    vclim_m = zeros(Float64, xdim, ydim, nstep_yr)   # negative v components
    vclim_p = zeros(Float64, xdim, ydim, nstep_yr)   # positive v components

    # Initialize wind component separation (CRITICAL: affects advection)
    @. uclim_m = ifelse(uclim >= 0.0, uclim, 0.0)  # positive winds only
    @. uclim_p = ifelse(uclim < 0.0, uclim, 0.0)   # negative winds only
    @. vclim_m = ifelse(vclim >= 0.0, vclim, 0.0)  # positive winds only
    @. vclim_p = ifelse(vclim < 0.0, vclim, 0.0)   # negative winds only
    Toclim = zeros(Float64, xdim, ydim, nstep_yr)   # deep ocean temperature [K]
    cldclim = zeros(Float64, xdim, ydim, nstep_yr)   # cloud cover fraction
    swetclim = zeros(Float64, xdim, ydim, nstep_yr)   # soil wetness [0-1]
end;

# ── notebook cell 0d200ce7-eadb-41de-8e31-45783a1faab9  (orig lines 870-876) ──
begin
    # ☀️ 2D solar field (ydim, nstep_yr) ─────────────────────────────
    sw_solar = zeros(Float64, ydim, nstep_yr) # 24hr mean solar radiation [W/m²]
    global dTrad = zeros(Float64, xdim, ydim, Int(nstep_yr))  # Tatmos-radiation offset
end;

# ── notebook cell 047b312f-8d6c-4732-aa0b-bea3de3e99e2  (orig lines 984-1002, split: see constants.jl for the calendar-lookup consts from this cell) ──
begin
    # 🕐 Time State Struct
    mutable struct TimeState
        jday::Int  # Current calendar day in year [1..365]
        ityr::Int  # Current timestep in year [1..730]
    end
end;

# ── notebook cell c06deaa7-2e6a-4143-a195-b473cbd84329  (orig lines 1083-1105) ──
begin
    # Annual-mean accumulators (xdim, ydim)
    Tsmn = zeros(Float64, xdim, ydim)   # surface temperature
    Tamn = zeros(Float64, xdim, ydim)   # air temperature
    Tomn = zeros(Float64, xdim, ydim)   # deep ocean temperature
    qmn = zeros(Float64, xdim, ydim)   # humidity
    amn = zeros(Float64, xdim, ydim)   # albedo

    swmn = zeros(Float64, xdim, ydim)   # shortwave radiation
    lwmn = zeros(Float64, xdim, ydim)   # longwave radiation
    qlatmn = zeros(Float64, xdim, ydim)   # latent heat flux
    qsensmn = zeros(Float64, xdim, ydim)   # sensible heat flux
    ftmn = zeros(Float64, xdim, ydim)   # temperature flux correction
    fqmn = zeros(Float64, xdim, ydim)   # humidity flux correction
end;

# ── notebook cell 7a06bf0d-a61c-4d28-b144-2725fe90ae62  (orig lines 1179-1182) ──
# Type alias for one monthly output record
const MonthlyRecord = NamedTuple{(:Ts, :Ta, :To, :q, :albedo, :ice, :precip, :evap, :qcrcl, :sw, :lw, :qlat, :qsens),NTuple{13,Matrix{Float64}}};
