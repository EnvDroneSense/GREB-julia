begin
    """
            CirculationWorkspace

        Pre-allocated buffers for diffusion, advection, and circulation calculations.
        Reused across all time steps to eliminate allocations.
        """
    mutable struct CirculationWorkspace
        # Polar sub-stepping buffers.
        T1h::Vector{Float32}      # polar sub-stepping
        dTxh::Vector{Float32}     # polar increment (Jacobi scratch)

        # Circulation work arrays
        X_work::Matrix{Float32}   # circulation work array
        dX_diff::Matrix{Float32}  # diffusion output
        dX_adv::Matrix{Float32}   # advection output
        dX_conv::Matrix{Float32}  # convection output

        # Tendency buffers
        temp_buf::Matrix{Float32}   # general workspace (humidity-update scratch)
        Q_sens_buf::Matrix{Float32} # Sensible heat flux buffer
        crcl::Matrix{Float32}       # dq_crcl (zero stand-in when log_crcl_dmc is off)

        # State buffers
        Ts0_buf::Matrix{Float32}    # Surface temperature output
        Ta0_buf::Matrix{Float32}    # Air temperature output
        To0_buf::Matrix{Float32}    # Ocean temperature output
        q0_buf::Matrix{Float32}     # Humidity output

        # LW radiation buffers
        e_co2_buf::Matrix{Float32}    # spatial CO₂ buffer
        e_vapor_buf::Matrix{Float32}  # spatial water vapor buffer
        em_buf::Matrix{Float32}       # spatial emissivity buffer
        LW_surf_buf::Matrix{Float32}  # Surface longwave
        LW_down_buf::Matrix{Float32}  # Downwelling longwave
        LW_up_buf::Matrix{Float32}    # Upwelling longwave

        # Hydrology buffers
        qs::Matrix{Float32}        # Saturation humidity buffer
        Tskin::Matrix{Float32}     # Skin temperature buffer
        rq::Matrix{Float32}        # Relative humidity buffer
        ws_base::Matrix{Float32}   # Base wind speed buffer
        # Hydrology
        Q_lat_buf::Matrix{Float32}
        Q_lat_air_buf::Matrix{Float32}
        dq_eva_buf::Matrix{Float32}
        dq_rain_buf::Matrix{Float32}
        cE_buf::Matrix{Float32}       # Surface exchange coefficient buffer

        # Deep_ocean
        dT_ocean_buf::Matrix{Float32}
        dTo_buf::Matrix{Float32}

        # Dedicated circulation output
        dTa_crcl::Matrix{Float32}   # temperature tendency
        dq_crcl::Matrix{Float32}    # humidity tendency

        # SWradiation
        ice_cover_buf::Matrix{Float32} # ice fraction
        a_surf_buf::Matrix{Float32}    # surface albedo
        albedo_buf::Matrix{Float32}    # combined albedo (surface + atmosphere)
        a_atmos_buf::Matrix{Float32}   # atmospheric albedo
        sw_buf::Matrix{Float32}        # net shortwave flux

        # time_loop
        precip_out::Matrix{Float32}   # precipitation output
        evap_out::Matrix{Float32}     # evaporation output
        qcrcl_out::Matrix{Float32}    # circulation moisture output
        term_north::Vector{Float32}   # northern boundary term
        term_south::Vector{Float32}   # southern boundary term
    end

    function CirculationWorkspace()
        CirculationWorkspace(
            zeros(Float32, xdim),# T1h
            zeros(Float32, xdim),# dTxh
            zeros(Float32, xdim, ydim),# X_work
            zeros(Float32, xdim, ydim),# dX_diff
            zeros(Float32, xdim, ydim),# dX_adv
            zeros(Float32, xdim, ydim),# dX_conv
            zeros(Float32, xdim, ydim),# temp_buf
            zeros(Float32, xdim, ydim),# Q_sens_buf
            zeros(Float32, xdim, ydim),# crcl
            # State buffers
            zeros(Float32, xdim, ydim),# Ts0_buf
            zeros(Float32, xdim, ydim),# Ta0_buf
            zeros(Float32, xdim, ydim),# To0_buf
            zeros(Float32, xdim, ydim),# q0_buf
            # LWradiation buffers
            zeros(Float32, xdim, ydim),# e_co2_buf
            zeros(Float32, xdim, ydim),# e_vapor_buf
            zeros(Float32, xdim, ydim),# em_buf
            zeros(Float32, xdim, ydim),# LW_surf_buf
            zeros(Float32, xdim, ydim),# LW_down_buf
            zeros(Float32, xdim, ydim),# LW_up_buf
            # Hydrology buffers
            zeros(Float32, xdim, ydim),# qs
            zeros(Float32, xdim, ydim),# Tskin
            zeros(Float32, xdim, ydim),# rq
            zeros(Float32, xdim, ydim),# ws_base
            # Hydrology
            zeros(Float32, xdim, ydim),  # Q_lat_buf
            zeros(Float32, xdim, ydim),  # Q_lat_air_buf
            zeros(Float32, xdim, ydim),  # dq_eva_buf
            zeros(Float32, xdim, ydim),  # dq_rain_buf
            zeros(Float32, xdim, ydim),  # cE_buf
            # Deep_ocean
            zeros(Float32, xdim, ydim),  # dT_ocean_buf
            zeros(Float32, xdim, ydim),  # dTo_buf
            # Dedicated circulation output
            zeros(Float32, xdim, ydim),  # dTa_crcl
            zeros(Float32, xdim, ydim),  # dq_crcl
            # SWradiation buffers
            zeros(Float32, xdim, ydim),  # ice_cover_buf
            zeros(Float32, xdim, ydim),  # a_surf_buf
            zeros(Float32, xdim, ydim),  # albedo_buf
            zeros(Float32, xdim, ydim),  # a_atmos_buf
            zeros(Float32, xdim, ydim),  # sw_buf
            # time_loop
            zeros(Float32, xdim, ydim),  # precip_out
            zeros(Float32, xdim, ydim),  # evap_out
            zeros(Float32, xdim, ydim),  # qcrcl_out
            zeros(Float32, xdim),  # term_north
            zeros(Float32, xdim),  # term_south
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
    Ts::Matrix{Float32}
    Ta::Matrix{Float32}
    To::Matrix{Float32}
    q::Matrix{Float32}
end

begin
    """
            MonthlyAccumulator

        Accumulates fields over a month for monthly-mean output.
        Reset after each month via `reset!`.
        """
    mutable struct MonthlyAccumulator
        Tmm::Matrix{Float32}          # Surface temperature accumulator
        Tamm::Matrix{Float32}         # Air temperature accumulator
        Tomm::Matrix{Float32}         # Ocean temperature accumulator
        qmm::Matrix{Float32}          # Humidity accumulator
        apmm::Matrix{Float32}         # Albedo accumulator
        icemm::Matrix{Float32}        # Ice fraction accumulator
        precipmm::Matrix{Float32}     # Precipitation accumulator
        evapmm::Matrix{Float32}       # Evaporation accumulator
        qcrclmm::Matrix{Float32}      # Circulation moisture accumulator
        swmm::Matrix{Float32}         # Shortwave radiation accumulator
        lwmm::Matrix{Float32}         # Longwave radiation accumulator
        qlatmm::Matrix{Float32}       # Latent heat accumulator
        qsensmm::Matrix{Float32}      # Sensible heat accumulator
        count::Int                    # Number of accumulations
    end

    function MonthlyAccumulator()
        MonthlyAccumulator(
            zeros(Float32, xdim, ydim),  # Tmm
            zeros(Float32, xdim, ydim),  # Tamm
            zeros(Float32, xdim, ydim),  # Tomm
            zeros(Float32, xdim, ydim),  # qmm
            zeros(Float32, xdim, ydim),  # apmm
            zeros(Float32, xdim, ydim),  # icemm
            zeros(Float32, xdim, ydim),  # precipmm
            zeros(Float32, xdim, ydim),  # evapmm
            zeros(Float32, xdim, ydim),  # qcrclmm
            zeros(Float32, xdim, ydim),  # swmm
            zeros(Float32, xdim, ydim),  # lwmm
            zeros(Float32, xdim, ydim),  # qlatmm
            zeros(Float32, xdim, ydim),  # qsensmm
            0
        )
    end

    function reset!(acc::MonthlyAccumulator)
        fill!(acc.Tmm, 0.0f0)
        fill!(acc.Tamm, 0.0f0)
        fill!(acc.Tomm, 0.0f0)
        fill!(acc.qmm, 0.0f0)
        fill!(acc.apmm, 0.0f0)
        fill!(acc.icemm, 0.0f0)
        fill!(acc.precipmm, 0.0f0)
        fill!(acc.evapmm, 0.0f0)
        fill!(acc.qcrclmm, 0.0f0)
        fill!(acc.swmm, 0.0f0)
        fill!(acc.lwmm, 0.0f0)
        fill!(acc.qlatmm, 0.0f0)
        fill!(acc.qsensmm, 0.0f0)
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

`ClimateFields()` builds an all-zero instance and leaves `loaded = false`;
`load_greb_jld2!` sets `loaded = true` once real climatology is in place.
[`greb_model!`](@ref) refuses to run unloaded fields unless explicitly
told to via `allow_uninitialized=true` - an all-zero climatology produces
a physically meaningless (~-40 °C) world rather than an error, so the
flag exists to keep that path opt-in.
"""
mutable struct ClimateFields
    # 2D fields (xdim, ydim)
    z_topo::Matrix{Float32}      # topography [m] (<0: ocean)
    glacier::Matrix{Float32}     # glacier mask (>0.5: glacier)
    z_ocean::Matrix{Float32}     # derived ocean depth [m]
    cap_surf::Matrix{Float32}    # surface heat capacity [J/K/m²]
    wz_air::Matrix{Float32}      # exp(-z_topo / z_air)
    wz_vapor::Matrix{Float32}    # exp(-z_topo / z_vapor)
    rain_limit::Matrix{Float32}  # -0.0015/(wz_vapor*r_qviwv*86400)
    # 3D climate fields (xdim, ydim, nstep_yr)
    Tclim::Array{Float32,3}      # surface temperature [K]
    uclim::Array{Float32,3}      # zonal wind [m/s]
    vclim::Array{Float32,3}      # meridional wind [m/s]
    qclim::Array{Float32,3}      # atmospheric humidity [kg/kg]
    mldclim::Array{Float32,3}    # mixed-layer depth [m]
    omegaclim::Array{Float32,3}    # vertical velocity [Pa/s]
    omegastdclim::Array{Float32,3} # omega std deviation [Pa/s]
    wsclim::Array{Float32,3}       # wind speed [m/s]

    # Anomaly fields for ENSO/climate-change experiments
    Tclim_anom_enso::Array{Float32,3}
    uclim_anom_enso::Array{Float32,3}
    vclim_anom_enso::Array{Float32,3}
    omegaclim_anom_enso::Array{Float32,3}
    wsclim_anom_enso::Array{Float32,3}
    Tclim_anom_cc::Array{Float32,3}
    uclim_anom_cc::Array{Float32,3}
    vclim_anom_cc::Array{Float32,3}
    omegaclim_anom_cc::Array{Float32,3}
    wsclim_anom_cc::Array{Float32,3}

    # Precomputed wind sign splits
    uclim_m::Array{Float32,3}    # negative u components
    uclim_p::Array{Float32,3}    # positive u components
    vclim_m::Array{Float32,3}    # negative v components
    vclim_p::Array{Float32,3}    # positive v components

    Toclim::Array{Float32,3}     # deep ocean temperature [K]
    cldclim::Array{Float32,3}    # cloud cover fraction
    swetclim::Array{Float32,3}   # soil wetness [0-1]

    # Solar / radiation
    sw_solar::Matrix{Float32}    # 24hr mean solar radiation [W/m²] (ydim, nstep_yr)
    dTrad::Array{Float32,3}      # Tatmos-radiation offset

    # Flux correction arrays (zeros unless loaded from file)
    TF_correct::Array{Float32,3}
    qF_correct::Array{Float32,3}
    ToF_correct::Array{Float32,3}

    # Regional CO₂ mask (1.0 = full CO₂, 0.5 = half CO₂)
    co2_part::Matrix{Float32}

    # false for a bare `ClimateFields()`; set by `load_greb_jld2!`. See the
    # docstring above and `greb_model!`'s `allow_uninitialized` keyword.
    loaded::Bool
end

function ClimateFields()
    z2(args...) = zeros(Float32, args...)
    ClimateFields(
        z2(xdim, ydim), z2(xdim, ydim), z2(xdim, ydim), z2(xdim, ydim), z2(xdim, ydim), z2(xdim, ydim), z2(xdim, ydim),
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
        ones(Float32, xdim, ydim),
        false,
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
    sw_solar_forcing::Float32   # runtime solar multiplier used by SWradiation!

    # Annual-mean accumulators (xdim, ydim)
    Tsmn::Matrix{Float32}    # surface temperature
    Tamn::Matrix{Float32}    # air temperature
    Tomn::Matrix{Float32}    # deep ocean temperature
    qmn::Matrix{Float32}     # humidity
    amn::Matrix{Float32}     # albedo
    swmn::Matrix{Float32}    # shortwave radiation
    lwmn::Matrix{Float32}    # longwave radiation
    qlatmn::Matrix{Float32}  # latent heat flux
    qsensmn::Matrix{Float32} # sensible heat flux
    ftmn::Matrix{Float32}    # temperature flux correction
    fqmn::Matrix{Float32}    # humidity flux correction
end

function ModelState()
    ModelState(1.0f0, (zeros(Float32, xdim, ydim) for _ in 1:11)...)
end

"""
    MonthlyRecord

One monthly-mean output record: a `NamedTuple` with fields `Ts`, `Ta`, `To`,
`q`, `albedo`, `ice`, `precip`, `evap`, `qcrcl`, `sw`, `lw`, `qlat`, `qsens`,
each an `(xdim, ydim)` `Matrix{Float32}`. Produced by [`output!`](@ref);
`greb_model!`'s `ctrl`/`scnr` results are `Vector{MonthlyRecord}`.
"""
const MonthlyRecord = NamedTuple{(:Ts, :Ta, :To, :q, :albedo, :ice, :precip, :evap, :qcrcl, :sw, :lw, :qlat, :qsens),NTuple{13,Matrix{Float32}}};
