begin
    """
            CirculationWorkspace

        Pre-allocated buffers for diffusion, advection, and circulation calculations.
        Reused across all time steps to eliminate allocations.
        """
    mutable struct CirculationWorkspace
        # Polar sub-stepping buffers.
        T1h::Vector{Float64}      # polar sub-stepping
        dTxh::Vector{Float64}     # polar increment (Jacobi scratch)

        # Circulation work arrays
        X_work::Matrix{Float64}   # circulation work array
        dX_diff::Matrix{Float64}  # diffusion output
        dX_adv::Matrix{Float64}   # advection output
        dX_conv::Matrix{Float64}  # convection output

        # Tendency buffers
        temp_buf::Matrix{Float64}   # general workspace (humidity-update scratch)
        Q_sens_buf::Matrix{Float64} # Sensible heat flux buffer
        crcl::Matrix{Float64}       # dq_crcl (zero stand-in when log_crcl_dmc is off)

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
            zeros(Float64, xdim, ydim),# temp_buf
            zeros(Float64, xdim, ydim),# Q_sens_buf
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

"""
    SurfaceState

A run's current surface/atmosphere state — `Ts`, `Ta`, `To`, `q` — passed as
one argument to [`diagnostics!`](@ref), [`output!`](@ref), [`time_loop!`](@ref),
and [`qflux_correction!`](@ref). A thin reference wrapper around
already-allocated arrays; construct once per run/call (like `ws`/`acc`),
never inside the per-timestep loop.
"""
struct SurfaceState
    Ts::Matrix{Float64}
    Ta::Matrix{Float64}
    To::Matrix{Float64}
    q::Matrix{Float64}
end

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
        Tmm = acc.Tmm; Tamm = acc.Tamm; Tomm = acc.Tomm; qmm = acc.qmm
        apmm = acc.apmm; icemm = acc.icemm
        precipmm = acc.precipmm; evapmm = acc.evapmm; qcrclmm = acc.qcrclmm
        swmm = acc.swmm; lwmm = acc.lwmm; qlatmm = acc.qlatmm; qsensmm = acc.qsensmm

        @turbo for j in 1:ydim
            for i in 1:xdim
                Tmm[i, j] += Ts[i, j]
                Tamm[i, j] += Ta[i, j]
                Tomm[i, j] += To[i, j]
                qmm[i, j] += q[i, j]
                apmm[i, j] += albedo[i, j]
                icemm[i, j] += ice[i, j]
                precipmm[i, j] += precip[i, j]
                evapmm[i, j] += evap[i, j]
                qcrclmm[i, j] += qcrcl[i, j]
                swmm[i, j] += sw[i, j]
                lwmm[i, j] += lw[i, j]
                qlatmm[i, j] += qlat[i, j]
                qsensmm[i, j] += qsens[i, j]
            end
        end
        acc.count += 1
    end
end;

"""
    ClimateFields

Loaded climatology, derived grid fields, flux corrections, and the
regional-CO2 mask/solar table — everything `load_greb_jld2!` fills in and
every physics function reads. One instance per `greb_model!` run; never
shared as global state.
"""
mutable struct ClimateFields
    # 2D fields (xdim, ydim)
    z_topo::Matrix{Float64}      # topography [m] (<0: ocean)
    glacier::Matrix{Float64}     # glacier mask (>0.5: glacier)
    z_ocean::Matrix{Float64}     # derived ocean depth [m]
    cap_surf::Matrix{Float64}    # surface heat capacity [J/K/m²]
    wz_air::Matrix{Float64}      # exp(-z_topo / z_air)
    wz_vapor::Matrix{Float64}    # exp(-z_topo / z_vapor)

    # 3D climate fields (xdim, ydim, nstep_yr)
    Tclim::Array{Float64,3}      # surface temperature [K]
    uclim::Array{Float64,3}      # zonal wind [m/s]
    vclim::Array{Float64,3}      # meridional wind [m/s]
    qclim::Array{Float64,3}      # atmospheric humidity [kg/kg]
    mldclim::Array{Float64,3}    # mixed-layer depth [m]
    omegaclim::Array{Float64,3}    # vertical velocity [Pa/s]
    omegastdclim::Array{Float64,3} # omega std deviation [Pa/s]
    wsclim::Array{Float64,3}       # wind speed [m/s]

    # Anomaly fields for ENSO/climate-change experiments
    Tclim_anom_enso::Array{Float64,3}
    uclim_anom_enso::Array{Float64,3}
    vclim_anom_enso::Array{Float64,3}
    omegaclim_anom_enso::Array{Float64,3}
    wsclim_anom_enso::Array{Float64,3}
    Tclim_anom_cc::Array{Float64,3}
    uclim_anom_cc::Array{Float64,3}
    vclim_anom_cc::Array{Float64,3}
    omegaclim_anom_cc::Array{Float64,3}
    wsclim_anom_cc::Array{Float64,3}

    # Precomputed wind sign splits
    uclim_m::Array{Float64,3}    # negative u components
    uclim_p::Array{Float64,3}    # positive u components
    vclim_m::Array{Float64,3}    # negative v components
    vclim_p::Array{Float64,3}    # positive v components

    Toclim::Array{Float64,3}     # deep ocean temperature [K]
    cldclim::Array{Float64,3}    # cloud cover fraction
    swetclim::Array{Float64,3}   # soil wetness [0-1]

    # Solar / radiation
    sw_solar::Matrix{Float64}    # 24hr mean solar radiation [W/m²] (ydim, nstep_yr)
    dTrad::Array{Float64,3}      # Tatmos-radiation offset

    # Flux correction arrays (zeros unless loaded from file)
    TF_correct::Array{Float64,3}
    qF_correct::Array{Float64,3}
    ToF_correct::Array{Float64,3}

    # Regional CO₂ mask (1.0 = full CO₂, 0.5 = half CO₂)
    co2_part::Matrix{Float64}
end

function ClimateFields()
    z2(args...) = zeros(Float64, args...)
    ClimateFields(
        z2(xdim, ydim), z2(xdim, ydim), z2(xdim, ydim), z2(xdim, ydim), z2(xdim, ydim), z2(xdim, ydim),
        z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr),
        z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr),
        z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr),
        z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr),
        z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr),
        z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr),
        z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr),
        z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr),
        z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr),
        z2(ydim, nstep_yr), z2(xdim, ydim, Int(nstep_yr)),
        z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr), z2(xdim, ydim, nstep_yr),
        ones(Float64, xdim, ydim),
    )
end

begin
    """
        TimeState

    Tracks the model's position within the current year: `jday` (calendar
    day, 1..365) and `ityr` (timestep-of-year, 1..`nstep_yr`). Mutated in
    place each timestep by [`time_loop!`](@ref)/[`qflux_correction!`](@ref).
    """
    mutable struct TimeState
        jday::Int  # Current calendar day in year [1..365]
        ityr::Int  # Current timestep in year [1..730]
    end
end;

"""
    ModelState

Per-run mutable state that isn't climatology: the runtime solar-forcing
multiplier (`SWradiation!` reads it) and the annual-mean diagnostic
accumulators (`diagnostics!` reads/writes them). One instance per
`greb_model!` run.
"""
mutable struct ModelState
    sw_solar_forcing::Float64   # runtime solar multiplier used by SWradiation!

    # Annual-mean accumulators (xdim, ydim)
    Tsmn::Matrix{Float64}    # surface temperature
    Tamn::Matrix{Float64}    # air temperature
    Tomn::Matrix{Float64}    # deep ocean temperature
    qmn::Matrix{Float64}     # humidity
    amn::Matrix{Float64}     # albedo
    swmn::Matrix{Float64}    # shortwave radiation
    lwmn::Matrix{Float64}    # longwave radiation
    qlatmn::Matrix{Float64}  # latent heat flux
    qsensmn::Matrix{Float64} # sensible heat flux
    ftmn::Matrix{Float64}    # temperature flux correction
    fqmn::Matrix{Float64}    # humidity flux correction
end

function ModelState()
    ModelState(1.0, (zeros(Float64, xdim, ydim) for _ in 1:11)...)
end

"""
    MonthlyRecord

One monthly-mean output record: a `NamedTuple` with fields `Ts`, `Ta`, `To`,
`q`, `albedo`, `ice`, `precip`, `evap`, `qcrcl`, `sw`, `lw`, `qlat`, `qsens`,
each an `(xdim, ydim)` `Matrix{Float64}`. Produced by [`output!`](@ref);
`greb_model!`'s `ctrl`/`scnr` results are `Vector{MonthlyRecord}`.
"""
const MonthlyRecord = NamedTuple{(:Ts, :Ta, :To, :q, :albedo, :ice, :precip, :evap, :qcrcl, :sw, :lw, :qlat, :qsens),NTuple{13,Matrix{Float64}}};
