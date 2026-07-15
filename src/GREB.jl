module GREB

# =============================================================================
# GREB — Globally Resolved Energy Balance model
#
# This module is a MECHANICAL extraction of the model definitions from the
# Pluto notebook `notebooks/GREB_julia.jl`. Function/struct/const bodies are
# copied VERBATIM (no logic changes). Only notebook scaffolding was removed:
#   * Pluto cell markers (`# ╔═╡ …`) and the cell-order footer
#   * markdown (`md"…"`) and `@bind` interactive-UI cells
#   * two notebook-only helpers: `current_physics_config` (reads @bind widget
#     globals) and `setup_benchmark` (BenchmarkTools/Profile) — reproduce their
#     role from the driver script `examples/run_greb.jl` instead.
# The original notebook is preserved unchanged under `notebooks/` for reference.
# =============================================================================

using Statistics
using NCDatasets
using StaticArrays        # static longitude indices
using LoopVectorization   # @turbo SIMD

export PhysicsConfig, CirculationWorkspace, MonthlyAccumulator, TimeState, MonthlyRecord
export read_jdal2, load_solar_forcing_jdal2, load_flux_corrections_jdal2!, load_greb_jdal2!
export create_experiment_config, set_hydrology_parameters!, init_model!
export SWradiation!, LWradiation!, hydro!, convergence!, seaice!, deep_ocean!
export diffusion!, advection!, circulation!, tendencies!, forcing
export diagnostics!, output!, time_loop!
export build_monthly_climatology, apply_scenario_anomalies, compute_annual_ice_climatology
export qflux_correction!, greb_model!
export xdim, ydim, nstep_yr

# ── notebook cell b303e4e9-49fa-45ad-967e-20f165fdf38c  (orig lines 73-112) ──
"""
    read_jdal2(filepath::String)

Read a JDAL2 Version 2 file into an Array{Float32}.

# Returns
- `(data, dim_names)` tuple where:
  - `data`: Array{Float32} with shape determined from file header
  - `dim_names`: Vector{String} of dimension names (e.g., ["lon", "lat", "time"])
"""
function read_jdal2(filepath::String)
    data, dims, dim_names = open(filepath, "r") do io
        # Magic bytes
        magic = UInt8[read(io, UInt8) for _ in 1:6]
        if magic[1:5] != UInt8[0x4A, 0x44, 0x41, 0x4C, 0x32]
            error("Not a valid JDAL2 file: $filepath")
        end
        version = magic[6]
        version == 0x02 || error("Only Version 2 supported, got $version")

        # Dimension names
        n_dim_names = read(io, Int32)
        dim_names = [String(read(io, read(io, Int32))) for _ in 1:n_dim_names]

        # Shape
        ndims = read(io, Int32)
        dims = [read(io, Int32) for _ in 1:ndims]

        # Data type
        read(io, UInt8) == 0x01 || error("Only Float32 supported")

        # Data
        data = Vector{Float32}(undef, prod(dims))
        read!(io, data)
        (data, dims, dim_names)
    end
    return (data=reshape(data, Tuple(dims)), dim_names=dim_names)
end

# ── notebook cell 8578d6aa-2782-4279-8f6b-78194b8ecc10  (orig lines 113-140) ──
function load_solar_forcing_jdal2(jdal2_dir::String, forcing_type::Symbol, index::Int=0)

    if forcing_type == :paleo
        filepath = joinpath(jdal2_dir, "solar_scenarios", "solar_paleo.jd2")
        result = read_jdal2(filepath)
        return result.data

    elseif forcing_type == :eccentricity
        filepath = joinpath(jdal2_dir, "solar_scenarios", "solar_eccentricity.jd2")
        result = read_jdal2(filepath)
        # result.data is (61, 48, 730) → index 0-based maps to position index+1
        @assert 0 <= index <= 60 "Eccentricity index must be 0-60"
        return result.data[index+1, :, :]

    elseif forcing_type == :obliquity
        filepath = joinpath(jdal2_dir, "solar_scenarios", "solar_obliquity.jd2")
        result = read_jdal2(filepath)
        indices = 0:25:230
        pos = findfirst(==(index), indices)
        @assert pos !== nothing "Obliquity index $index not found in $(indices)"
        return result.data[pos, :, :]

    else
        error("Unknown forcing type: $forcing_type. Use :paleo, :eccentricity, or :obliquity")
    end
end

# ── notebook cell 531589ab-c6b5-4048-9ba5-f9ad62ab00a6  (orig lines 202-227) ──
# Note: this should be in the script
begin
    # 📐 Grid dimensions ──────────────────────────────────
    const xdim = 96                                 # number of longitude grid points
    const ydim = 48                                 # number of latitude grid points
    const dlon = 360.0 / xdim                       # longitude spacing [degrees]
    const dlat = 180.0 / ydim                       # latitude spacing  [degrees]

    # ⏱️ Time stepping ────────────────────────────────────
    const ndays_yr = 365                            # days per year (no leap years)
    const Δt = 12.0 * 3600.0                        # main time step [s] (12 hours)
    const Δt_crcl = 1800.0                          # circulation sub-time step [s] (30 min)
    const ndt_days = Int(round(24 * 3600 / Δt))     # time steps per day
    const nstep_yr = Int(ndays_yr * ndt_days)       # time steps per year (= 730)
    const ntime = max(1, Int(round(Δt / Δt_crcl)))  # Number of sub-steps within one main time step

    # 📅 Calendar constants ───────────────────────────────
    const cjday_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    const jday_mon_cumsum = cumsum(cjday_mon)

    # 🔢 Physical Constants & Numerical Limits ────────────
    const min_T_K = 233.15              # 273.15 - 40°C, minimum allowed surface temperature [K]
    const max_humidity_change = 0.020   # Maximum humidity increment [kg/kg]
    const min_humidity_change = 0.9     # Fraction of humidity that can be removed
end;

# ── notebook cell 19f106e4-2b82-47e1-9284-799a105f30cb  (orig lines 235-291) ──
Base.@kwdef mutable struct PhysicsConfig
    # Mean Climate Switches
    log_clouds_dmc::Bool = true
    log_vapor_dmc::Bool = true
    log_ice_dmc::Bool = true
    log_crcl_dmc::Bool = true
    log_hydro_dmc::Bool = true
    log_atmos_dmc::Bool = true
    log_co2_dmc::Bool = true
    log_ocean_dmc::Bool = true
    log_qflux_dmc::Bool = true

    # CO₂ Response Switches
    log_clouds_drsp::Bool = true
    log_vapor_drsp::Bool = true
    log_ice_drsp::Bool = true
    log_crcl_drsp::Bool = true
    log_hydro_drsp::Bool = true
    log_topo_drsp::Bool = true
    log_humid_drsp::Bool = true
    log_ocean_drsp::Bool = true

    # Circulation Components
    log_ice::Bool = true
    log_hdif::Bool = true
    log_hadv::Bool = true
    log_vdif::Bool = true
    log_vadv::Bool = true
    log_conv::Bool = true

    # Hydrology Parameters
    log_rain::Int = 0
    log_eva::Int = -1
    log_clim::Int = 0

    # External Forcing
    log_tsurf_ext::Bool = false
    log_hwind_ext::Bool = false
    log_omega_ext::Bool = false

    # Experiment Type
    experiment::Symbol = :full_model  # :full_model, :constant_topo, :co2_double, etc.

    # CO₂ concentration for experiments (ppm)
    co2_concentration::Float64 = 340.0

    # Solar forcing multiplier
    solar_multiplier::Float64 = 1.0

    # Hydrology parameters (calculated by set_hydrology_parameters!)
    c_q::Float64 = 1.0
    c_rq::Float64 = 0.0
    c_omega::Float64 = 0.0
    c_omegastd::Float64 = 0.0
end

# ── notebook cell bb82d1dc-8f08-4351-9371-a8efed1dd9bc  (orig lines 292-356) ──
begin
    """
        create_experiment_config(experiment::Symbol) -> PhysicsConfig

    Create a pre-configured `PhysicsConfig` for common experiment types.

    # Experiments
    - `:full_model` - All processes active (default)
    - `:constant_topo` - Constant topography (log_topo_drsp = false)
    - `:co2_double` - 2×CO₂ (680 ppm)
    - `:co2_quadruple` - 4×CO₂ (1360 ppm)
    - `:solar_plus27` - +27 W/m² solar constant
    - `:elnino` - El Niño conditions
    - `:lanina` - La Niña conditions  
    - `:paleo_231kyr` - Paleoclimate (200 ppm CO₂)
    - `:rcp85` - RCP8.5 climate change scenario
    """
    function create_experiment_config(experiment::Symbol)::PhysicsConfig
        if experiment == :full_model
            return PhysicsConfig(experiment=:full_model)

        elseif experiment == :constant_topo
            cfg = PhysicsConfig(experiment=:constant_topo)
            cfg.log_topo_drsp = false  # Constant topography
            return cfg

        elseif experiment == :co2_double
            cfg = PhysicsConfig(experiment=:co2_double)
            cfg.co2_concentration = 680.0  # 2×CO₂
            return cfg

        elseif experiment == :co2_quadruple
            cfg = PhysicsConfig(experiment=:co2_quadruple)
            cfg.co2_concentration = 1360.0  # 4×CO₂
            return cfg

        elseif experiment == :solar_plus27
            cfg = PhysicsConfig(experiment=:solar_plus27)
            cfg.solar_multiplier = (1365.0 + 27.0) / 1365.0
            return cfg

        elseif experiment == :elnino
            cfg = PhysicsConfig(experiment=:elnino)
            return cfg

        elseif experiment == :lanina
            cfg = PhysicsConfig(experiment=:lanina)
            return cfg

        elseif experiment == :paleo_231kyr
            cfg = PhysicsConfig(experiment=:paleo_231kyr)
            cfg.co2_concentration = 200.0
            return cfg

        elseif experiment == :rcp85
            cfg = PhysicsConfig(experiment=:rcp85)
            return cfg

        else
            error("Unknown experiment: $experiment")
        end
    end
end;

# ── notebook cell 32b531ab-ee71-4685-af65-a2bae0b868f6  (orig lines 363-377) ──
begin
    # 💧 Optimized Hydrology Parameter Lookup Table ────────────────────────
    const HYDRO_PARAMS = (
        -1 => (1.0, 0.0, 0.0, 0.0),                     # Original GREB
        1 => (-1.391649, 3.018774, 0.0, 0.0),           # +Relative humidity
        2 => (0.862162, 0.0, -29.02096, 0.0),           # +Omega convergence
        3 => (-0.2685845, 1.4591853, -26.9858807, 0.0), # +RH & Omega
        0 => (-1.88, 2.25, -17.69, 59.07)               # Best GREB (ERA-Interim)
    )

    # 🎯 Cached Weight Arrays (avoid recomputation) ───
    global WZ_CACHE = Dict{Float64,Matrix{Float64}}()
end;

# ── notebook cell d0d74213-96ad-4a45-9a31-37f640a21a45  (orig lines 378-400) ──
# ⚙️ Optimized Parameter Initialization ──────────────────────────────────
"""
    set_hydrology_parameters!(cfg::PhysicsConfig)

Initialize precipitation parameters `c_q, c_rq, c_omega, c_omegastd` based on
`cfg.log_rain` and `cfg.log_clim` settings.
"""
function set_hydrology_parameters!(cfg::PhysicsConfig)
    global c_q, c_rq, c_omega, c_omegastd

    # Fast lookup instead of if-else chain
    params = get(HYDRO_PARAMS, cfg.log_rain, (1.0, 0.0, 0.0, 0.0))
    c_q, c_rq, c_omega, c_omegastd = params

    # NCEP parameter adjustment
    if cfg.log_rain == 0 && cfg.log_clim == 1
        c_q, c_rq, c_omega, c_omegastd = -1.27, 1.99, -16.54, 21.15
    end

    @info "⚙️ MSCM hydrology: log_rain=$(cfg.log_rain), log_clim=$(cfg.log_clim) → (c_q=$c_q, c_rq=$c_rq, c_omega=$c_omega, c_omegastd=$c_omegastd)"
end

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
    # Global instance (reused across all time steps)
    global circ_workspace = CirculationWorkspace()
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
    co2_part_scn = ones(Float64, xdim, ydim) # Scenario-specific spatial mask
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

# ── notebook cell 0dbfb663-46e7-4873-ac77-1e8e392fe69d  (orig lines 737-749) ──
begin
    # ☀️ Solar forcing storage
    global sw_solar_forcing_data = nothing  # Will hold (48, 730) array when loaded
    global sw_solar_forcing_state = Ref(1.0)  # Runtime solar multiplier used by SWradiation!

    # 🔧 Flux correction arrays (initialised with zeros, overwritten if files exist)
    global TF_correct = zeros(Float64, xdim, ydim, nstep_yr)
    global qF_correct = zeros(Float64, xdim, ydim, nstep_yr)
    global ToF_correct = zeros(Float64, xdim, ydim, nstep_yr)
    const ΔT_AIR_FACTOR = Δt / cap_air
end;

# ── notebook cell f578f25e-047e-4a7e-8483-d544c7b4bec3  (orig lines 750-772) ──
function load_flux_corrections_jdal2!(jdal2_dir::String)
    """Load flux corrections from JDAL2 files (zeros if missing)"""
    # Correction arrays are already defined as zeros in the global scope.
    correction_files = Dict(
        "Tsurf_flux_correction.jd2" => TF_correct,
        "vapour_flux_correction.jd2" => qF_correct,
        "Tocean_flux_correction.jd2" => ToF_correct
    )

    for (filename, array) in correction_files
        filepath = joinpath(jdal2_dir, "climatology", filename)
        if isfile(filepath)
            result = read_jdal2(filepath)
            array .= result.data
            println("✅ Loaded $filename")
        else
            fill!(array, 0.0)
            @warn "$filename not found, using zeros"
        end
    end
end

# ── notebook cell 898dd5aa-5273-4833-9a4a-0f2b94cc8d38  (orig lines 779-782) ──
const p_emi = [9.0721, 106.7252, 61.5562, 0.0179, 0.0028,
               0.0570, 0.3462, 2.3406, 0.7032, 1.0662]

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

# ── notebook cell 2bf0fe8e-5718-4c1e-863b-85db7b3ae7f3  (orig lines 877-977) ──
function load_greb_jdal2!(jdal2_dir::String; dataset::Symbol=:ncep)
    """Load all GREB input data from JDAL2 formatted files"""

    if !isdir(jdal2_dir)
        error("JDAL2 directory not found: $jdal2_dir")
    end

    # Static 2D files
    println("📂 Loading static fields...")
    topo_result = read_jdal2(joinpath(jdal2_dir, "static", "global.topography.jd2"))
    global z_topo .= topo_result.data

    glacier_result = read_jdal2(joinpath(jdal2_dir, "static", "greb.glaciers.jd2"))
    global glacier .= glacier_result.data

    # Dataset-specific file mapping
    file_map = Dict(
        :ncep => Dict(
            "Tclim" => "ncep.tsurf.1948-2007.clim.jd2",
            "uclim" => "ncep.zonal_wind.850hpa.clim.jd2",
            "vclim" => "ncep.meridional_wind.850hpa.clim.jd2",
            "qclim" => "ncep.atmospheric_humidity.clim.jd2",
            "swetclim" => "ncep.soil_moisture.clim.jd2"
        ),
        :era => Dict(
            "Tclim" => "erainterim.tsurf.1979-2015.clim.jd2",
            "uclim" => "erainterim.zonal_wind.850hpa.clim.jd2",
            "vclim" => "erainterim.meridional_wind.850hpa.clim.jd2",
            "qclim" => "erainterim.atmospheric_humidity.clim.jd2",
            "swetclim" => "ncep.soil_moisture.clim.jd2"
        )
    )

    # Use mixed dataset as fallback
    files = get(file_map, dataset, file_map[:ncep])

    println("📂 Loading 3D climatology ($dataset dataset)...")
    climatology_dir = joinpath(jdal2_dir, "climatology")

    # Load each variable individually with unique result names
    tsurf_result = read_jdal2(joinpath(climatology_dir, files["Tclim"]))
    global Tclim .= tsurf_result.data

    uwind_result = read_jdal2(joinpath(climatology_dir, files["uclim"]))
    global uclim .= uwind_result.data

    vwind_result = read_jdal2(joinpath(climatology_dir, files["vclim"]))
    global vclim .= vwind_result.data

    humid_result = read_jdal2(joinpath(climatology_dir, files["qclim"]))
    global qclim .= humid_result.data

    swet_result = read_jdal2(joinpath(climatology_dir, files["swetclim"]))
    global swetclim .= swet_result.data

    # Common climatology files
    println("📂 Loading common climatology fields...")

    cld_result = read_jdal2(joinpath(climatology_dir, "isccp.cloud_cover.clim.jd2"))
    global cldclim .= cld_result.data

    mld_result = read_jdal2(joinpath(climatology_dir, "woce.ocean_mixed_layer_depth.clim.jd2"))
    global mldclim .= mld_result.data

    tocean_result = read_jdal2(joinpath(climatology_dir, "Tocean.clim.jd2"))
    global Toclim .= tocean_result.data

    omega_result = read_jdal2(joinpath(climatology_dir, "erainterim.omega.vertmean.clim.jd2"))
    global omegaclim .= omega_result.data

    omegastd_result = read_jdal2(joinpath(climatology_dir, "erainterim.omega_std.vertmean.clim.jd2"))
    global omegastdclim .= omegastd_result.data

    ws_result = read_jdal2(joinpath(climatology_dir, "erainterim.windspeed.850hpa.clim.jd2"))
    global wsclim .= ws_result.data

    # Solar radiation (special: lat × time)
    println("📂 Loading solar radiation...")
    solar_path = joinpath(jdal2_dir, "solar", "solar_radiation.clim.jd2")
    if isfile(solar_path)
        solar_result = read_jdal2(solar_path)
        @assert size(solar_result.data) == (ydim, nstep_yr) "Wrong solar dimensions"
        global sw_solar .= solar_result.data
    else
        error("Solar radiation file not found: $solar_path")
    end

    # Optional: Load flux corrections
    println("📂 Loading flux corrections...")
    load_flux_corrections_jdal2!(jdal2_dir)

    # Update wind sign splits
    @. uclim_m = ifelse(uclim >= 0.0, uclim, 0.0)
    @. uclim_p = ifelse(uclim < 0.0, uclim, 0.0)
    @. vclim_m = ifelse(vclim >= 0.0, vclim, 0.0)
    @. vclim_p = ifelse(vclim < 0.0, vclim, 0.0)

    println("✅ All GREB data loaded successfully from JDAL2")
end

# ── notebook cell 047b312f-8d6c-4732-aa0b-bea3de3e99e2  (orig lines 984-1002) ──
begin
    # 🕐 Time State Struct
    mutable struct TimeState
        jday::Int  # Current calendar day in year [1..365]
        ityr::Int  # Current timestep in year [1..730]
    end

    # 📅 Precomputed Calendar Lookup
    const max_timesteps = 200 * nstep_yr  # Support up to 200-year runs
    const calendar_lookup = [(
        day=mod((it - 1) ÷ ndt_days, ndays_yr) + 1,
        step=mod(it - 1, nstep_yr) + 1
    ) for it in 1:max_timesteps]

    const polar_treshold = 2.5e5  # 250 km in meters
    const IS_POLAR = [dxlat_grid[k] <= polar_treshold for k in 1:ydim]
end;

# ── notebook cell c06deaa7-2e6a-4143-a195-b473cbd84329  (orig lines 1083-1105) ──
begin
    # Annual-mean accumulators (xdim, ydim) 
    Tsmn = zeros(Float64, xdim, ydim)   # surface temperature
    Tamn = zeros(Float64, xdim, ydim)   # air temperature
    Tomn = zeros(Float64, xdim, ydim)   # deep ocean temperature
    qmn = zeros(Float64, xdim, ydim)   # humidity
    amn = zeros(Float64, xdim, ydim)   # albedo

    # additional accumulators
    icmn = zeros(Float64, xdim, ydim)   # ice cover fraction
    prmn = zeros(Float64, xdim, ydim)   # precipitation tendency
    evmn = zeros(Float64, xdim, ydim)   # evaporation tendency
    qcrclmn = zeros(Float64, xdim, ydim)   # circulation tendency

    swmn = zeros(Float64, xdim, ydim)   # shortwave radiation
    lwmn = zeros(Float64, xdim, ydim)   # longwave radiation
    qlatmn = zeros(Float64, xdim, ydim)   # latent heat flux
    qsensmn = zeros(Float64, xdim, ydim)   # sensible heat flux
    ftmn = zeros(Float64, xdim, ydim)   # temperature flux correction
    fqmn = zeros(Float64, xdim, ydim)   # humidity flux correction
end;

# ── notebook cell a8b5fa01-0526-46f4-9aa0-31b52398a2bd  (orig lines 1124-1143) ──
begin
    # Monthly-mean buffers (xdim, ydim) 
    Tmm = zeros(Float64, xdim, ydim)   # surface temperature
    Tamm = zeros(Float64, xdim, ydim)   # air temperature
    Tomm = zeros(Float64, xdim, ydim)   # deep ocean temperature
    qmm = zeros(Float64, xdim, ydim)   # humidity
    apmm = zeros(Float64, xdim, ydim)   # albedo

    # additional monthly buffers
    icmm = zeros(Float64, xdim, ydim)   # ice cover fraction
    prmm = zeros(Float64, xdim, ydim)   # precipitation tendency
    evmm = zeros(Float64, xdim, ydim)   # evaporation tendency
    qcrclmm = zeros(Float64, xdim, ydim)   # circulation tendency
    swmm = zeros(Float64, xdim, ydim)   # shortwave radiation
    lwmm = zeros(Float64, xdim, ydim)   # longwave radiation
    qlatmm = zeros(Float64, xdim, ydim)   # latent heat flux
    qsensmm = zeros(Float64, xdim, ydim)   # sensible heat flux
end;

# ── notebook cell 1bdeea50-a73b-47a5-b941-d9d3374f54cf  (orig lines 1144-1156) ──
begin
    # Control run monthly means (for anomaly calculation)
    Tmn_ctrl = zeros(Float64, xdim, ydim, 12)  # surface temperature
    Tamn_ctrl = zeros(Float64, xdim, ydim, 12)  # air temperature
    Tomn_ctrl = zeros(Float64, xdim, ydim, 12)  # deep ocean temperature
    qmn_ctrl = zeros(Float64, xdim, ydim, 12)  # humidity
    icmn_ctrl = zeros(Float64, xdim, ydim, 12)  # ice cover
    prmn_ctrl = zeros(Float64, xdim, ydim, 12)  # precipitation
    evamn_ctrl = zeros(Float64, xdim, ydim, 12)  # evaporation
    qcrclmn_ctrl = zeros(Float64, xdim, ydim, 12) # circulation
end;

# ── notebook cell 7a06bf0d-a61c-4d28-b144-2725fe90ae62  (orig lines 1179-1182) ──
# Type alias for one monthly output record
const MonthlyRecord = NamedTuple{(:Ts, :Ta, :To, :q, :albedo, :ice, :precip, :evap, :qcrcl, :sw, :lw, :qlat, :qsens),NTuple{13,Matrix{Float64}}};

# ── notebook cell d404043f-8080-4262-9ab6-d9bb13eee504  (orig lines 1200-1304) ──
function init_model!(cfg::PhysicsConfig)

    # ── Hydrology Parameter Initialization ────────────
    set_hydrology_parameters!(cfg)

    # ── dTrad: offset between T_atm and radiation temperature ────
    @. dTrad = -0.16 * Tclim - 5.0

    # ── z_ocean: 3× maximum mixed-layer depth over the year ──────
    z_ocean .= 3.0 .* dropdims(maximum(mldclim; dims=3); dims=3)

    # ── Sensitivity experiment overrides ─────────────────────────
    if !cfg.log_topo_drsp
        @. z_topo = min(z_topo, 1.0)       # constant topography
    end

    # Apply cloud deconstruction switch
    if !cfg.log_clouds_dmc
        cldclim .= 0.0  # zero cloud climatology
    end

    # Apply flux correction conditional zeroing (MSCM feature)
    if !cfg.log_topo_drsp && !cfg.log_qflux_dmc
        TF_correct .= 0.0
        qF_correct .= 0.0
        ToF_correct .= 0.0
    end

    # Climatology modifications
    if !cfg.log_hydro_dmc
        qclim .= 0.0  # zero out humidity climatology
    end

    if !cfg.log_vapor_dmc
        qclim .= 0.0052          # constant water vapor
    end

    if !cfg.log_ocean_dmc
        mldclim .= d_ocean       # no deep ocean
    end

    # ── Experiment Handler ───────────────────────────────────
    # Apply advanced experiment forcing
    if cfg.experiment == :rcp85
        @info "Applying CMIP5 RCP8.5 climate change forcing"
        Tclim .+= Tclim_anom_cc
        uclim .+= uclim_anom_cc
        vclim .+= vclim_anom_cc
        omegaclim .+= omegaclim_anom_cc
        wsclim .+= wsclim_anom_cc
    elseif cfg.experiment == :elnino
        @info "Applying ERA-Interim El Niño forcing"
        Tclim .+= Tclim_anom_enso
        uclim .+= uclim_anom_enso
        vclim .+= vclim_anom_enso
        omegaclim .+= omegaclim_anom_enso
        wsclim .+= wsclim_anom_enso
    elseif cfg.experiment == :lanina
        @info "Applying ERA-Interim La Niña forcing"
        Tclim .-= Tclim_anom_enso
        uclim .-= uclim_anom_enso
        vclim .-= vclim_anom_enso
        omegaclim .-= omegaclim_anom_enso
        wsclim .-= wsclim_anom_enso
    end

    # ── Topography pressure weights ─────
    @. wz_air = exp(-z_topo / z_air)
    @. wz_vapor = exp(-z_topo / z_vapor)

    # ── Surface heat capacity ────────────────────────────────────
    @inbounds for j in 1:ydim
        for i in 1:xdim
            if z_topo[i, j] > 0.0
                cap_surf[i, j] = cap_land
            else
                cap_surf[i, j] = cfg.log_ocean_dmc ? cap_ocean * mldclim[i, j, 1] : cap_land
            end
        end
    end

    # ── Initial conditions from last time step of climatology ────
    Ts_ini = Tclim[:, :, nstep_yr] |> copy   # surface temperature
    Ta_ini = copy(Ts_ini)                      # air temperature = Tsurf
    To_ini = Toclim[:, :, nstep_yr] |> copy   # deep ocean temperature
    q_ini = qclim[:, :, nstep_yr] |> copy   # atmospheric water vapor

    # ── Control CO₂ level ───────────────────────────────────────
    CO2_ctrl = cfg.co2_concentration

    if cfg.experiment == :a1b_scenario
        CO2_ctrl = 298.0  # A1B scenario baseline
    elseif cfg.experiment in (:a1b_enhanced, :rcp26, :rcp45, :rcp60, :rcp85, :custom_co2)
        CO2_ctrl = 280.0  # IPCC scenarios baseline
    end

    if !cfg.log_co2_dmc
        CO2_ctrl = 0.0  # Zero CO2 for deconstruction experiments
    end

    return (Ts_ini=Ts_ini, Ta_ini=Ta_ini, To_ini=To_ini,
        q_ini=q_ini, CO2_ctrl=CO2_ctrl)
end

# ── notebook cell 1df2b91b-be14-427a-87b3-95cdef26ce00  (orig lines 1359-1419) ──
function SWradiation!(Ts, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    # Reuse workspace buffers
    ice_cover = ws.ice_cover_buf # output: ice fraction
    a_surf = ws.a_surf_buf    # surface albedo
    albedo = ws.albedo_buf    # output: combined albedo (surface + atmosphere) [hernoemen]
    a_atmos = ws.a_atmos_buf   # atmospheric albedo
    sw = ws.sw_buf        # output: net shortwave flux

    # 1. Ice cover fraction – branch‑free with ifelse, vectorized
    @turbo for i in 1:xdim, j in 1:ydim
        T = Ts[i, j]
        land_expr = ifelse(T <= Tl_ice1, 1.0,
            ifelse(T < Tl_ice2,
                1.0 - (T - Tl_ice1) * inv_Tl_ice_range,
                0.0))
        ocean_expr = ifelse(T <= To_ice1, 1.0,
            ifelse(T < To_ice2,
                1.0 - (T - To_ice1) * inv_To_ice_range,
                0.0))
        ice_cover[i, j] = ifelse(z_topo[i, j] >= 0.0, land_expr, ocean_expr)
    end

    # 2. Atmospheric albedo – simple multiplication, use broadcasting
    cld = @view cldclim[:, :, timestate.ityr]
    @. a_atmos = cld * a_cloud

    # 3. Surface albedo – conditional logic, @turbo beneficial
    if cfg.log_ice
        @turbo for i in 1:xdim, j in 1:ydim
            T = Ts[i, j]
            # Land albedo expression
            land_alb = ifelse(T <= Tl_ice1, a_no_ice + da_ice,
                ifelse(T >= Tl_ice2, a_no_ice,
                    a_no_ice + da_ice * (1.0 - (T - Tl_ice1) * inv_Tl_ice_range)))
            # Ocean albedo expression
            ocean_alb = ifelse(T <= To_ice1, a_no_ice + da_ice,
                ifelse(T >= To_ice2, a_no_ice,
                    a_no_ice + da_ice * (1.0 - (T - To_ice1) * inv_To_ice_range)))
            # Choose based on topography
            a_surf[i, j] = ifelse(z_topo[i, j] >= 0.0, land_alb, ocean_alb)
            # Glacier override: if glacier mask > 0.5, set to ice albedo
            a_surf[i, j] = ifelse(glacier[i, j] > 0.5, a_no_ice + da_ice, a_surf[i, j])
        end
    else
        @. a_surf = a_no_ice
    end

    # 4. Combined albedo
    @. albedo = a_surf + a_atmos - a_surf * a_atmos

    # 5. Shortwave flux
    multiplier = sw_solar_forcing_state[] * 0.01 * S0_var
    for j in 1:ydim
        sf = sw_solar[j, timestate.ityr] * multiplier
        @. sw[:, j] = sf * (1.0 - albedo[:, j])
    end

    return (SW=sw, albedo=albedo, ice_cover=ice_cover)
end

# ── notebook cell c0c40037-4169-4d38-bebe-2086cebc24f2  (orig lines 1437-1476) ──
function LWradiation!(Ts, Ta, q, CO2, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    # Extract workspace buffers
    e_co2 = ws.e_co2_buf    # CO₂ [ppm scaled by pressure]
    e_vapor = ws.e_vapor_buf  # water vapour [kg/m²]
    em = ws.em_buf       # emissivity ε_atmos
    LW_surf = ws.LW_surf_buf  # surface long-wave flux [W/m²]
    LW_down = ws.LW_down_buf  # downward long-wave flux [W/m²]
    LW_up = ws.LW_up_buf    # upward long-wave flux [W/m²]

    # Current cloud cover (climatology, 3D array)
    e_cloud = @view cldclim[:, :, timestate.ityr]

    ## 1. Effective columns (topography scaling via wz_air, precomputed)
    @. e_vapor = wz_air * r_qviwv * q
    @. e_co2 = wz_air * CO2 * co2_part

    # ── Emissivity (log-regression with 10 parameters) ─────────
    @. em = p_emi[4] * log(p_emi[1] * e_co2 + p_emi[2] * e_vapor + p_emi[3]) +
            p_emi[7] +
            p_emi[5] * log(p_emi[1] * e_co2 + p_emi[3]) +
            p_emi[6] * log(p_emi[2] * e_vapor + p_emi[3])

    # Cloud adjustment
    @. em = (p_emi[8] - e_cloud) / p_emi[9] * (em - p_emi[10]) + p_emi[10]

    # 4. Radiation temperature (precomputed offset dTrad = -0.16*Tclim - 5 K)
    dTr = @view dTrad[:, :, timestate.ityr]
    @. LW_surf = -σ * Ts^4
    @. LW_down = -em * σ * (Ta + dTr)^4
    LW_up .= LW_down

    if !cfg.log_atmos_dmc
        LW_down .= 0.0
        LW_up .= 0.0
    end

    return (LW_surf=LW_surf, LW_up=LW_up, LW_down=LW_down, em=em)
end

# ── notebook cell 606032a2-b2ca-4fd8-9930-afd83aecec7a  (orig lines 1495-1608) ──
function hydro!(Ts, q, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    c_q = cfg.c_q
    c_rq = cfg.c_rq
    c_omega = cfg.c_omega
    c_omegastd = cfg.c_omegastd

    fill!(ws.Q_lat_buf, 0.0)
    fill!(ws.Q_lat_air_buf, 0.0)
    fill!(ws.dq_eva_buf, 0.0)
    fill!(ws.dq_rain_buf, 0.0)

    if !cfg.log_atmos_dmc || !cfg.log_hydro_dmc || !cfg.log_hydro_drsp
        return (Q_lat=ws.Q_lat_buf, Q_lat_air=ws.Q_lat_air_buf,
            dq_eva=ws.dq_eva_buf, dq_rain=ws.dq_rain_buf)
    end

    u = @view uclim[:, :, timestate.ityr]
    v = @view vclim[:, :, timestate.ityr]
    swet = @view swetclim[:, :, timestate.ityr]
    omega = @view omegaclim[:, :, timestate.ityr]
    omegastd = @view omegastdclim[:, :, timestate.ityr]

    const_factor1 = 3.75e-3
    const_factor2 = 17.08085
    const_factor3 = 234.175
    gust_land = 4.0
    gust_ocean = 9.0
    cE_land = 0.25 * ce
    cE_ocean = 0.58 * ce
    const_latent = cq_latent * ρ_air * ce

    # Saturation humidity
    @turbo for j in 1:ydim
        for i in 1:xdim
            T = Ts[i, j] - 273.15
            ws.qs[i, j] = const_factor1 * exp(const_factor2 * T / (T + const_factor3)) * wz_air[i, j]
            ws.qs[i, j] = max(ws.qs[i, j], 1e-8)
        end
    end

    # Evaporation
    if cfg.log_eva == -1
        @turbo for j in 1:ydim
            for i in 1:xdim
                u_val = u[i, j];
                v_val = v[i, j]
                wind = sqrt(u_val*u_val + v_val*v_val)
                wind = sqrt(wind*wind + ifelse(z_topo[i, j] > 0.0, gust_land, gust_ocean))
                ws.Q_lat_buf[i, j] = (q[i, j] - ws.qs[i, j]) * wind * const_latent * swet[i, j]
            end
        end
    elseif cfg.log_eva == 0
        ws_view = @view wsclim[:, :, timestate.ityr]
        @turbo for j in 1:ydim
            for i in 1:xdim
                ws.Tskin[i, j] = ifelse(z_topo[i, j] > 0.0, Ts[i, j] + 5.0, Ts[i, j] + 1.0)
                ws.Tskin[i, j] = ifelse(ws.Tskin[i, j] < 200.0, 200.0, ws.Tskin[i, j])
                T = ws.Tskin[i, j] - 273.15
                qs_val = const_factor1 * exp(const_factor2 * T / (T + const_factor3)) * wz_air[i, j]
                ws.qs[i, j] = qs_val

                ws.ws_base[i, j] = ws_view[i, j]
                gust = ifelse(z_topo[i, j] > 0.0, 132.25, 29.16)
                wind = sqrt(ws.ws_base[i, j]*ws.ws_base[i, j] + gust)

                ws.cE[i, j] = ifelse(z_topo[i, j] > 0.0, cE_land, cE_ocean)
                ws.Q_lat_buf[i, j] = ws.cE[i, j] * wind * ρ_air * cq_latent * (q[i, j] - qs_val) * swet[i, j]
            end
        end
    else
        @turbo for j in 1:ydim
            for i in 1:xdim
                u_val = u[i, j];
                v_val = v[i, j]
                wind = sqrt(u_val*u_val + v_val*v_val)
                wind = sqrt(wind*wind + ifelse(z_topo[i, j] > 0.0, gust_land, gust_ocean))
                ws.Q_lat_buf[i, j] = (q[i, j] - ws.qs[i, j]) * wind * const_latent * swet[i, j]
            end
        end
    end

    # Precipitation - use ws.rq buffer
    @turbo for j in 1:ydim
        for i in 1:xdim
            ws.rq[i, j] = q[i, j] / max(ws.qs[i, j], 1e-8)
            ws.dq_rain_buf[i, j] = (c_q + c_rq * ws.rq[i, j] + c_omega * omega[i, j] + c_omegastd * omegastd[i, j]) * cq_rain * q[i, j]
        end
    end

    # Apply rain limit
    if cfg.log_rain == 1
        limit_val = -0.0015 / (wz_vapor[1, 1] * r_qviwv * 86400.0)
        @turbo for j in 1:ydim
            for i in 1:xdim
                ws.dq_rain_buf[i, j] = ifelse(ws.dq_rain_buf[i, j] >= limit_val, limit_val, ws.dq_rain_buf[i, j])
            end
        end
    end

    # Water vapor tendencies
    @turbo for j in 1:ydim
        for i in 1:xdim
            ws.dq_eva_buf[i, j] = -ws.Q_lat_buf[i, j] / cq_latent / r_qviwv
            min_dq = -0.9 * q[i, j] / Δt
            ws.dq_rain_buf[i, j] = ifelse(ws.dq_rain_buf[i, j] < min_dq, min_dq, ws.dq_rain_buf[i, j])
            ws.Q_lat_air_buf[i, j] = -ws.dq_rain_buf[i, j] * cq_latent * r_qviwv
        end
    end

    return (Q_lat=ws.Q_lat_buf,
        Q_lat_air=ws.Q_lat_air_buf,
        dq_eva=ws.dq_eva_buf,
        dq_rain=ws.dq_rain_buf)
end

# ── notebook cell c6a0d656-8289-4f74-b8d8-f94c236e541d  (orig lines 1620-1637) ──
function convergence!(T1, omegaclim, timestate, ws::CirculationWorkspace)
    """Calculate moisture flux convergence using omega vertical velocity.
    	
    	Implements Eq. 18 from Stassen et al 2019.
    	
    	Args:
    		T1: Input field (typically specific humidity) [kg/kg]
    		omegaclim: Vertical velocity climatology [Pa/s] 
    		timestate: Time state object
    		ws: Workspace for temporary arrays
    	"""
    omega = @view omegaclim[:, :, timestate.ityr]

    @. ws.dX_conv = -T1 * omega * const_factor
    return nothing
end

# ── notebook cell 9bff59c5-4631-4091-8230-989a835788e5  (orig lines 1650-1688) ──
function seaice!(Ts0, timestate, cfg::PhysicsConfig)
    mld = @view mldclim[:, :, timestate.ityr]

    if !cfg.log_ocean_dmc
        return     # No ice feedback: skip sea ice calculation
    end

    # Compute ice‑dependent heat capacity for ocean points
    @turbo for i in 1:xdim, j in 1:ydim
        is_ocean = z_topo[i, j] < 0.0
        T = Ts0[i, j]
        mld_val = mld[i, j]
        cap_open = cap_ocean * mld_val

        # Ice fraction (0 = no ice, 1 = full ice)
        ice_frac = ifelse(T <= To_ice1, 1.0,
            ifelse(T >= To_ice2, 0.0,
                1.0 - (T - To_ice1) * inv_To_ice_range))

        # Blend between land (ice) and open ocean capacities
        cap_with_ice = cap_land * ice_frac + cap_open * (1.0 - ice_frac)

        # Apply only to ocean points; keep land points unchanged
        cap_surf[i, j] = ifelse(is_ocean, cap_with_ice, cap_surf[i, j])
    end

    # Override for experiments without ice‑albedo feedback
    if !cfg.log_ice
        @turbo for i in 1:xdim, j in 1:ydim
            cap_surf[i, j] = ifelse(z_topo[i, j] > 0.0, cap_land, cap_ocean * mld[i, j])
        end
        return
    end

    # Glacier override: ice sheets have land heat capacity
    @. cap_surf = ifelse(glacier > 0.5, cap_land, cap_surf)
end

# ── notebook cell 625089e2-ef77-4821-a6d6-d0a0f88207f2  (orig lines 1703-1750) ──
function deep_ocean!(Ts, To, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    # Use pre-allocated zero buffers
    dT_ocean = ws.dT_ocean_buf
    dTo = ws.dTo_buf

    # no deep-ocean coupling
    if !cfg.log_ocean_dmc || !cfg.log_ocean_drsp
        fill!(dT_ocean, 0.0)
        fill!(dTo, 0.0)
        return (dT_ocean=dT_ocean, dTo=dTo)
    end

    # ── Change in mixed-layer depth ─────────────────────────
    mld_now = @view mldclim[:, :, timestate.ityr]
    mld_prev = timestate.ityr > 1 ? @view(mldclim[:, :, timestate.ityr-1]) : @view(mldclim[:, :, nstep_yr])

    # Zero buffers first (one pass, cheap)
    fill!(dT_ocean, 0.0)
    fill!(dTo, 0.0)

    # ── Entrainment & detrainment & turbulent mixing ──────
    @turbo for i in 1:xdim, j in 1:ydim
        active = (z_topo[i, j] < 0.0) & (Ts[i, j] >= To_ice2)
        h_now = mld_now[i, j]
        h_prev = mld_prev[i, j]
        dh = h_now - h_prev
        z_deep = z_ocean[i, j]
        z_rem = z_deep - h_now

        # Entrainment/detrainment contributions (only when active)
        dTo_entr = ifelse(active & (dh < 0.0), c_effmix * (-dh / z_rem) *
                                               (Ts[i, j] - To[i, j]), 0.0)
        dT_ocean_entr = ifelse(active & (dh > 0.0), c_effmix * (dh / h_now) *
                                                    (To[i, j] - Ts[i, j]), 0.0)

        # Turbulent mixing (only when active)
        Tx = ifelse(Ts[i, j] > To_ice2, Ts[i, j], To_ice2)
        dTo_turb = ifelse(active, turb_coeff * (Tx - To[i, j]) / z_rem, 0.0)
        dT_ocean_turb = ifelse(active, turb_coeff * (To[i, j] - Tx) / h_now, 0.0)

        # Combine (buffer was zeroed before loop)
        dTo[i, j] = dTo_entr + dTo_turb
        dT_ocean[i, j] = dT_ocean_entr + dT_ocean_turb
    end
    return (dT_ocean=dT_ocean, dTo=dTo)
end

# ── notebook cell ba96178d-77d4-4f26-a94f-5ad43c5242db  (orig lines 1777-1888) ──
function diffusion!(T1, h_scl, ws::CirculationWorkspace, timestate)
    # Zero output buffer (we will accumulate into it)
    fill!(ws.dX_diff, 0.0)

    # Topographic scaling (choose based on scale height)
    wz = if h_scl == z_air
        wz_air
    elseif h_scl == z_vapor
        wz_vapor
    else
        error("Invalid h_scl = $h_scl (must be z_air or z_vapor)")
    end

    # Precomputed geometry/coefficients
    dxlat = dxlat_grid
    ccy = ccy_diff
    ccx = ccx_diff

    # Pre‑cached neighbour indices
    jm1, jp1 = lon_jm1, lon_jp1
    jm2, jp2 = lon_jm2, lon_jp2
    jm3, jp3 = lon_jm3, lon_jp3

    # ----- Precompute k‑independent terms for the poles -----
    # For k == 1 (North Pole)
    @view @. ws.term_north = ccy * wz[:, 2] * (T1[:, 2] - T1[:, 1])
    # For k == ydim (South Pole)
    @view @. ws.term_south = ccy * wz[:, ydim-1] * (T1[:, ydim-1] - T1[:, ydim])

    for k in 1:ydim
        # ----- Meridional diffusion -----
        if k == 1
            @turbo for i in 1:xdim
                ws.dX_diff[i, k] += wz[i, k] * ws.term_north[i]
            end
        elseif k == ydim
            @turbo for i in 1:xdim
                ws.dX_diff[i, k] += wz[i, k] * ws.term_south[i]
            end
        else
            # Mid‑latitudes: no precomputation possible (depends on k‑1, k+1)
            @view @. ws.dX_diff[:, k] += wz[:, k] * ccy * (
                wz[:, k-1] * (T1[:, k-1] - T1[:, k]) +
                wz[:, k+1] * (T1[:, k+1] - T1[:, k])
            )
        end

        # ----- Zonal diffusion -----
        if dxlat[k] > 2.5e5   # mid‑latitudes, normal time step
            # Complex stencil – keep @turbo
            @turbo for j in 1:xdim
                jm1v = jm1[j];
                jp1v = jp1[j]
                jm2v = jm2[j];
                jp2v = jp2[j]
                jm3v = jm3[j];
                jp3v = jp3[j]

                dTx = ccx[k] * 0.05 * (
                    10.0 * (wz[jm1v, k] * (T1[jm1v, k] - T1[j, k]) +
                            wz[jp1v, k] * (T1[jp1v, k] - T1[j, k])) +
                    4.0 * (wz[jm2v, k] * (T1[jm2v, k] - T1[jm1v, k]) +
                           wz[jm1v, k] * (T1[j, k] - T1[jm1v, k])) +
                    4.0 * (wz[jp1v, k] * (T1[j, k] - T1[jp1v, k]) +
                           wz[jp2v, k] * (T1[jp2v, k] - T1[jp1v, k])) +
                    1.0 * (wz[jm3v, k] * (T1[jm3v, k] - T1[jm2v, k]) +
                           wz[jm2v, k] * (T1[jm1v, k] - T1[jm2v, k])) +
                    1.0 * (wz[jp2v, k] * (T1[jp1v, k] - T1[jp2v, k]) +
                           wz[jp3v, k] * (T1[jp3v, k] - T1[jp2v, k]))
                )
                ws.dX_diff[j, k] += wz[j, k] * dTx
            end
        else   # polar regions – sub‑timestepping (complex, keep @turbo)
            # Number of sub‑steps for stability
            dd = max(1, round(Int, Δt_crcl / (dxlat[k]^2 / κ)))
            dtdff2 = Δt_crcl / dd
            time2 = max(1, round(Int, Δt_crcl / dtdff2))
            ccx2 = κ * dtdff2 / dxlat[k]^2

            # Copy current row into temporary buffer
            ws.T1h .= @view T1[:, k]

            for _ in 1:time2
                @turbo for j in 1:xdim
                    jm1v = jm1[j];
                    jp1v = jp1[j]
                    jm2v = jm2[j];
                    jp2v = jp2[j]
                    jm3v = jm3[j];
                    jp3v = jp3[j]

                    dq = ccx2 * 0.05 * (
                        10.0 * (wz[jm1v, k] * (ws.T1h[jm1v] - ws.T1h[j]) +
                                wz[jp1v, k] * (ws.T1h[jp1v] - ws.T1h[j])) +
                        4.0 * (wz[jm2v, k] * (ws.T1h[jm2v] - ws.T1h[jm1v]) +
                               wz[jm1v, k] * (ws.T1h[j] - ws.T1h[jm1v])) +
                        4.0 * (wz[jp1v, k] * (ws.T1h[j] - ws.T1h[jp1v]) +
                               wz[jp2v, k] * (ws.T1h[jp2v] - ws.T1h[jp1v])) +
                        1.0 * (wz[jm3v, k] * (ws.T1h[jm3v] - ws.T1h[jm2v]) +
                               wz[jm2v, k] * (ws.T1h[jm1v] - ws.T1h[jm2v])) +
                        1.0 * (wz[jp2v, k] * (ws.T1h[jp1v] - ws.T1h[jp2v]) +
                               wz[jp3v, k] * (ws.T1h[jp3v] - ws.T1h[jp2v]))
                    )
                    # Stability clamp
                    dq = ifelse(dq <= -ws.T1h[j], -0.9 * ws.T1h[j], dq)
                    ws.T1h[j] += dq
                end
            end

            # Add total change (scaled by outer wz) to output buffer
            @view @. ws.dX_diff[:, k] += wz[:, k] * (ws.T1h - T1[:, k])
        end
    end

    return nothing
end

# ── notebook cell 2bab06b9-ca98-4142-99cb-d2ad4f1cde93  (orig lines 1903-2046) ──
function advection!(T1, h_scl, ws::CirculationWorkspace, timestate, cfg::PhysicsConfig)
    # Disable advection for water vapour or heat according to switches
    if (h_scl == z_vapor && !cfg.log_vadv) || (h_scl == z_air && !cfg.log_hadv)
        fill!(ws.dX_adv, 0.0)
        return nothing
    end

    # Pre-zero the output buffer (we will accumulate into it)
    fill!(ws.dX_adv, 0.0)

    # Extract 2D views for current time step
    vclim_p_t = @view vclim_p[:, :, timestate.ityr]
    vclim_m_t = @view vclim_m[:, :, timestate.ityr]
    uclim_p_t = @view uclim_p[:, :, timestate.ityr]
    uclim_m_t = @view uclim_m[:, :, timestate.ityr]

    # Topographic scaling (choose based on scale height)
    wz = if h_scl == z_air
        wz_air
    elseif h_scl == z_vapor
        wz_vapor
    else
        error("Invalid h_scl = $h_scl (must be z_air or z_vapor)")
    end

    # Precomputed constants
    dxlat = dxlat_grid
    ccy = ccy_adv
    ccx = ccx_adv
    is_polar = IS_POLAR

    @inbounds for k in 1:ydim
        # ----- Meridional (v) advection -----
        if k == 1          # North Pole
            @turbo for j in 1:xdim
                v_p = vclim_p_t[j, k]
                ws.dX_adv[j, k] += ccy * v_p * (
                    wz[j, 2] * (T1[j, 1] - T1[j, 2]) +
                    wz[j, 3] * (T1[j, 1] - T1[j, 3])
                ) / 3.0
            end
        elseif k == 2
            @turbo for j in 1:xdim
                v_m = vclim_m_t[j, k]
                v_p = vclim_p_t[j, k]
                ws.dX_adv[j, k] += ccy * (
                    -v_m * wz[j, 1] * (T1[j, 2] - T1[j, 1]) +
                    v_p * (wz[j, 3] * (T1[j, 2] - T1[j, 3]) +
                           wz[j, 4] * (T1[j, 2] - T1[j, 4])) / 3.0
                )
            end
        elseif k >= 3 && k <= ydim-2
            @turbo for j in 1:xdim
                km1, km2 = k-1, k-2
                kp1, kp2 = k+1, k+2
                v_m = vclim_m_t[j, k]
                v_p = vclim_p_t[j, k]
                ws.dX_adv[j, k] += ccy * (
                    -v_m * (wz[j, km1] * (T1[j, k] - T1[j, km1]) +
                            wz[j, km2] * (T1[j, k] - T1[j, km2])) +
                    v_p * (wz[j, kp1] * (T1[j, k] - T1[j, kp1]) +
                           wz[j, kp2] * (T1[j, k] - T1[j, kp2]))
                ) / 3.0
            end
        elseif k == ydim-1
            @turbo for j in 1:xdim
                km1, km2 = k-1, k-2
                kp1 = k+1
                v_m = vclim_m_t[j, k]
                v_p = vclim_p_t[j, k]
                ws.dX_adv[j, k] += ccy * (
                    -v_m * (wz[j, km1] * (T1[j, k] - T1[j, km1]) +
                            wz[j, km2] * (T1[j, k] - T1[j, km2])) / 3.0 +
                    v_p * wz[j, kp1] * (T1[j, k] - T1[j, kp1])
                )
            end
        else               # k == ydim (South Pole)
            @turbo for j in 1:xdim
                km1, km2 = k-1, k-2
                v_m = vclim_m_t[j, k]
                ws.dX_adv[j, k] += ccy * (
                    -v_m * (wz[j, km1] * (T1[j, k] - T1[j, km1]) +
                            wz[j, km2] * (T1[j, k] - T1[j, km2]))
                ) / 3.0
            end
        end

        # ----- Zonal (u) advection -----
        if !is_polar[k]   # mid‑latitudes, normal timestep
            @turbo for j in 1:xdim
                jm1, jp1 = lon_jm1[j], lon_jp1[j]
                jm2, jp2 = lon_jm2[j], lon_jp2[j]
                u_m = uclim_m_t[j, k]
                u_p = uclim_p_t[j, k]
                ws.dX_adv[j, k] += ccx[k] * (
                    -u_m * (wz[jm1, k] * (T1[j, k] - T1[jm1, k]) +
                            wz[jm2, k] * (T1[j, k] - T1[jm2, k])) +
                    u_p * (wz[jp1, k] * (T1[j, k] - T1[jp1, k]) +
                           wz[jp2, k] * (T1[j, k] - T1[jp2, k]))
                ) / 3.0
            end
        else               # polar regions – sub‑timestepping
            # Number of sub‑steps (CFL stability)
            dd = max(1, round(Int, Δt_crcl / (dxlat[k] / 10.0)))
            dtdff2 = Δt_crcl / dd
            time2 = max(1, round(Int, Δt_crcl / dtdff2))
            ccx2 = dtdff2 / dxlat[k] / 2.0

            # Copy current row into temporary buffer
            ws.T1h .= @view T1[:, k]

            for _ in 1:time2
                # One fused loop: compute increment and update in place
                @turbo for j in 1:xdim
                    jm1, jp1 = lon_jm1[j], lon_jp1[j]
                    jm2, jp2 = lon_jm2[j], lon_jp2[j]
                    jm3, jp3 = lon_jm3[j], lon_jp3[j]   # precomputed!
                    u_m = uclim_m_t[j, k]
                    u_p = uclim_p_t[j, k]

                    dq = ccx2 * (
                        -u_m * (10.0 * wz[jm1, k] * (ws.T1h[j] - ws.T1h[jm1]) +
                                4.0 * wz[jm2, k] * (ws.T1h[jm1] - ws.T1h[jm2]) +
                                1.0 * wz[jm3, k] * (ws.T1h[jm2] - ws.T1h[jm3])) +
                        u_p * (10.0 * wz[jp1, k] * (ws.T1h[j] - ws.T1h[jp1]) +
                               4.0 * wz[jp2, k] * (ws.T1h[jp1] - ws.T1h[jp2]) +
                               1.0 * wz[jp3, k] * (ws.T1h[jp2] - ws.T1h[jp3]))
                    ) / 20.0

                    # Stability clamp (avoid negative water vapour)
                    dq = ifelse(dq <= -ws.T1h[j], -0.9 * ws.T1h[j], dq)
                    ws.T1h[j] += dq
                end
            end

            # Add total change to the output buffer
            @view  @. ws.dX_adv[:, k] += ws.T1h - T1[:, k]
        end
    end

    return nothing
end

# ── notebook cell db75ea52-9387-4c15-bde5-61777ac9b570  (orig lines 2047-2079) ──
function circulation!(X_in, h_scl, dX_out, ws::CirculationWorkspace, timestate, cfg::PhysicsConfig)
    # Early exit if atmospheric processes disabled
    if (!cfg.log_atmos_dmc || !cfg.log_crcl_dmc || !cfg.log_crcl_drsp)
        fill!(dX_out, 0.0)
        return nothing
    end

    # Precompute flags (hoist conditionals)
    do_diff_v = cfg.log_vdif == 1 && h_scl == z_vapor
    do_diff_h = cfg.log_hdif == 1 && h_scl == z_air
    do_adv_v = cfg.log_vadv == 1 && h_scl == z_vapor
    do_adv_h = cfg.log_hadv == 1 && h_scl == z_air
    do_conv = cfg.log_conv == 0 && h_scl == z_vapor

    copyto!(ws.X_work, X_in)

    for _tt in 1:ntime
        do_diff_v && diffusion!(ws.X_work, h_scl, ws, timestate)
        do_diff_h && diffusion!(ws.X_work, h_scl, ws, timestate)
        do_adv_v && advection!(ws.X_work, h_scl, ws, timestate, cfg)
        do_adv_h && advection!(ws.X_work, h_scl, ws, timestate, cfg)
        do_conv && convergence!(ws.X_work, omegaclim, timestate, ws)

        @. ws.X_work += ws.dX_diff + ws.dX_adv + ws.dX_conv
    end

    # Final difference
    @. dX_out = ws.X_work - X_in

    return nothing
end

# ── notebook cell e493fae7-239a-494c-9a59-728446d70f7a  (orig lines 2080-2127) ──
function tendencies!(CO2, Ts, Ta, To, q, ws::CirculationWorkspace,
    timestate, cfg::PhysicsConfig)

    # Short-wave radiation → albedo, SW flux
    sw_out = SWradiation!(Ts, timestate, cfg, ws)

    # Long-wave radiation → LW_surf, LW_up, LW_down, emissivity
    lw_out = LWradiation!(Ts, Ta, q, CO2, timestate, cfg, ws)

    # Sensible heat flux
    Q_sens = ws.Q_sens_buf
    if cfg.log_atmos_dmc
        @. Q_sens = ct_sens * (Ta - Ts)
    else
        fill!(Q_sens, 0.0)
    end

    # Hydrological cycle → latent heat + evaporation/rain tendencies
    hy_out = hydro!(Ts, q, timestate, cfg, ws)

    # Atmospheric circulation — temperature diffusion/advection 
    circulation!(Ta, z_air, ws.dTa_crcl, ws, timestate, cfg)

    # Atmospheric circulation — water-vapour diffusion/advection 
    circulation!(q, z_vapor, ws.dq_crcl, ws, timestate, cfg)

    # Deep ocean coupling
    do_out = deep_ocean!(Ts, To, timestate, cfg, ws)

    return (albedo=sw_out.albedo,
        SW=sw_out.SW,
        ice_cover=sw_out.ice_cover,
        LW_surf=lw_out.LW_surf,
        Q_lat=hy_out.Q_lat,
        Q_sens=Q_sens,
        Q_lat_air=hy_out.Q_lat_air,
        dq_eva=hy_out.dq_eva,
        dq_rain=hy_out.dq_rain,
        dq_crcl=ws.dq_crcl,
        dTa_crcl=ws.dTa_crcl,
        dT_ocean=do_out.dT_ocean,
        dTo=do_out.dTo,
        LW_down=lw_out.LW_down,
        LW_up=lw_out.LW_up,
        em=lw_out.em)
end

# ── notebook cell 1894ad94-cdf8-4e79-a0e5-b72088db31be  (orig lines 2162-2343) ──
function forcing(it, year, cfg::PhysicsConfig, icmn_ctrl; nstep_yr=nstep_yr)
    # Default CO₂ concentration
    CO2 = cfg.co2_concentration
    sw_solar_forcing = 1.0

    # 📜 Legacy experiments ───────────
    if cfg.experiment == :constant_topo
        CO2 = 550.0  # 550 ppm CO₂ steady state

    elseif cfg.experiment == :a1b_scenario
        CO2_1950 = 310.0;
        CO2_2000 = 370.0;
        CO2_2050 = 520.0
        if year <= 2000
            CO2 = CO2_1950 + 60.0 / 50.0 * (year - 1950)
        elseif year <= 2050
            CO2 = CO2_2000 + 150.0 / 50.0 * (year - 2000)
        elseif year <= 2100
            CO2 = CO2_2050 + 180.0 / 50.0 * (year - 2050)
        end

        # 💨 CO₂ scaling experiments ──────────────────────────────────────────────
    elseif cfg.experiment == :co2_double
        CO2 = 680.0  # 2×CO₂ (already set, but explicit)

    elseif cfg.experiment == :co2_quadruple
        CO2 = 1360.0  # 4×CO₂

    elseif cfg.experiment == :co2_10x
        CO2 = 3400.0  # 10×CO₂

    elseif cfg.experiment == :co2_half
        CO2 = 170.0  # 0.5×CO₂

    elseif cfg.experiment == :co2_zero
        CO2 = 0.0  # 0×CO₂ (no greenhouse effect)

    # ☀️ Solar forcing experiments ───────────────────────────────────────────
    elseif cfg.experiment == :solar_plus27
        CO2 = 340.0
        sw_solar_forcing = (1365.0 + 27.0) / 1365.0

    elseif cfg.experiment == :solar_cycle_11yr
        CO2 = 340.0
        sw_solar_forcing = (1365.0 + 1.0 * sin(2π * year / 11.0)) / 1365.0

        # 📈 Enhanced A1B scenario ──────────────────────────────────────────────
    elseif cfg.experiment == :a1b_enhanced
        CO2_1950 = 310.0;
        CO2_2000 = 370.0;
        CO2_2050 = 520.0
        if year <= 2000
            CO2 = CO2_1950 + 60.0 / 50.0 * (year - 1950)
        elseif year <= 2050
            CO2 = CO2_2000 + 150.0 / 50.0 * (year - 2000)
        elseif year <= 2100
            CO2 = CO2_2050 + 180.0 / 50.0 * (year - 2050)
        end

        # ── Time-varying CO₂ experiments ────────────
    elseif cfg.experiment == :co2_sine_wave
        CO2 = 340.0 + 170.0 + 170.0 * cos(2π * (year - 13.0) / 30.0)

    elseif cfg.experiment == :co2_step
        CO2 = year >= 1980 ? 340.0 : 680.0

        # ── Paleoclimate experiments ────────────────────   
    elseif cfg.experiment == :paleo_231kyr
        CO2 = 200.0

    elseif cfg.experiment == :paleo_solar_modern_co2
        CO2 = 340.0

    elseif cfg.experiment == :modern_solar_paleo_co2
        CO2 = 200.0

        # ── Orbital forcing experiments ─────────────────
    elseif cfg.experiment == :obliquity
        CO2 = 340.0     # Solar forcing loaded externally

    elseif cfg.experiment == :eccentricity
        CO2 = 340.0     # Solar forcing loaded externally

    elseif cfg.experiment == :earth_sun_distance
        CO2 = 340.0     # Solar constant varies with Earth-Sun distance

    # 📂 File I/O dependent experiments (placeholders) ───────────────────────
    elseif cfg.experiment == :rcp26
        error("RCP2.6 scenario requires external CO₂ data file. Not yet implemented.")

    elseif cfg.experiment == :rcp45
        error("RCP4.5 scenario requires external CO₂ data file. Not yet implemented.")

    elseif cfg.experiment == :rcp60
        error("RCP6.0 scenario requires external CO₂ data file. Not yet implemented.")

    elseif cfg.experiment == :rcp85
        CO2 = 340.0  # Handled by boundary conditions

    elseif cfg.experiment == :custom_co2
        error("Custom CO₂ scenario requires external trajectory file. Not yet implemented.")

        # 🌍 Regional/partial CO₂ experiments ────────────────────────────────────
    elseif startswith(string(cfg.experiment), "regional_co2_")
        # Reset co2_part to full CO₂ first
        co2_part .= 1.0

        if cfg.experiment == :regional_co2_nh
            # 2×CO₂ Northern Hemisphere only
            CO2 = 680.0
            co2_part[:, 1:24] .= 0.5

        elseif cfg.experiment == :regional_co2_sh
            # 2×CO₂ Southern Hemisphere only
            CO2 = 680.0
            co2_part[:, 25:48] .= 0.5

        elseif cfg.experiment == :regional_co2_tropics
            # 2×CO₂ Tropics only
            CO2 = 680.0
            co2_part[:, 1:15] .= 0.5
            co2_part[:, 33:48] .= 0.5
            for i in 4:4:96
                co2_part[i, 33] = 1.0
                co2_part[i, 15] = 1.0
            end

        elseif cfg.experiment == :regional_co2_extratropics
            # 2×CO₂ Extratropics only
            CO2 = 680.0
            co2_part[:, 16:32] .= 0.5
            for i in 4:4:96
                co2_part[i, 32] = 1.0
                co2_part[i, 16] = 1.0
            end

        elseif cfg.experiment == :regional_co2_ocean
            # 2×CO₂ Ocean only
            CO2 = 680.0
            for j in 1:ydim, i in 1:xdim
                if z_topo[i, j] > 0.0
                    co2_part[i, j] = 0.5
                end
            end
            icmn_ctrl1 = @view icmn_ctrl[:, :, 1]
            for j in 1:ydim, i in 1:xdim
                if icmn_ctrl1[i, j] >= 0.5
                    co2_part[i, j] = 0.5
                end
            end

        elseif cfg.experiment == :regional_co2_land_ice
            # 2×CO₂ Land/Ice only
            CO2 = 680.0
            for j in 1:ydim, i in 1:xdim
                if z_topo[i, j] <= 0.0
                    co2_part[i, j] = 0.5
                end
            end
            icmn_ctrl1 = @view icmn_ctrl[:, :, 1]
            for j in 1:ydim, i in 1:xdim
                if icmn_ctrl1[i, j] >= 0.5
                    co2_part[i, j] = 1.0
                end
            end

        elseif cfg.experiment == :regional_co2_winter
            # 2×CO₂ Boreal Winter only
            ityr_step = mod(it - 1, nstep_yr) + 1
            CO2 = (ityr_step <= 181 || ityr_step >= 547) ? 680.0 : 340.0

        elseif cfg.experiment == :regional_co2_summer
            # 2×CO₂ Boreal Summer only
            ityr_step = mod(it - 1, nstep_yr) + 1
            CO2 = (ityr_step <= 181 || ityr_step >= 547) ? 340.0 : 680.0
        end

        # 📊 Forced boundary condition experiments (handled in scenario loop) ───── 
    elseif cfg.experiment == :elnino || cfg.experiment == :lanina || cfg.experiment == :rcp85
        CO2 = 340.0
    end

    return (CO2=CO2, sw_solar_forcing=sw_solar_forcing)
end

# ── notebook cell cc0e682c-8767-498e-8178-2de4e796b3a8  (orig lines 2355-2392) ──
function diagnostics!(it, year, CO2, Ts0, Ta0, To0, q0, albedo, sw, lw_surf, q_lat, q_sens, timestate)
    # Accumulate
    Tsmn .+= Ts0;
    Tamn .+= Ta0;
    Tomn .+= To0
    qmn .+= q0;
    amn .+= albedo
    swmn .+= sw;
    lwmn .+= lw_surf
    qlatmn .+= q_lat;
    qsensmn .+= q_sens
    ftmn .+= @view TF_correct[:, :, timestate.ityr]
    fqmn .+= @view qF_correct[:, :, timestate.ityr]

    if timestate.ityr == nstep_yr
        # Compute annual means
        n = nstep_yr
        Tsmn ./= n;
        Tamn ./= n;
        Tomn ./= n
        qmn ./= n;
        amn ./= n
        swmn ./= n;
        lwmn ./= n
        qlatmn ./= n;
        qsensmn ./= n
        ftmn ./= n;
        fqmn ./= n

        # Global mean and sample points (°C)
        global_mean = sum(Tsmn) / (xdim * ydim) - 273.15
        point1 = Tsmn[48, 27] - 273.15   # Tropical Pacific
        point2 = Tsmn[16, 38] - 273.15   # Hamburg/North Europe

        println(year, "  ", round(global_mean, digits=2),
            "  ", round(point1, digits=2),
            "  ", round(point2, digits=2))

        # Reset accumulators
        fill!(Tsmn, 0.0);
        fill!(Tamn, 0.0);
        fill!(Tomn, 0.0)
        fill!(qmn, 0.0);
        fill!(amn, 0.0)
        fill!(swmn, 0.0);
        fill!(lwmn, 0.0)
        fill!(qlatmn, 0.0);
        fill!(qsensmn, 0.0)
        fill!(ftmn, 0.0);
        fill!(fqmn, 0.0)
    end
    return nothing
end

# ── notebook cell dfdde9f1-b226-4af0-9ac2-36f1b01622fa  (orig lines 2405-2437) ──
function output!(it, irec, mon, Ts0, Ta0, To0, q0, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens,
    output_buf::Vector{MonthlyRecord}, acc::MonthlyAccumulator, timestate)
    # ----- SAFETY: clamp month to 1..12 -----
    mon = clamp(mon, 1, 12)

    accumulate!(acc, Ts0, Ta0, To0, q0, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens)

    # ----- Check end of month -----
    if timestate.jday == jday_mon_cumsum[mon] && (it % ndt_days == 0)
        ndm = jday_mon[mon] * ndt_days
        irec += 1
        push!(output_buf, (
            Ts=copy(acc.Tmm ./ ndm),
            Ta=copy(acc.Tamm ./ ndm),
            To=copy(acc.Tomm ./ ndm),
            q=copy(acc.qmm ./ ndm),
            albedo=copy(acc.apmm ./ ndm),
            ice=copy(acc.icemm ./ ndm),
            precip=copy(acc.precipmm ./ ndm),
            evap=copy(acc.evapmm ./ ndm),
            qcrcl=copy(acc.qcrclmm ./ ndm),
            sw=copy(acc.swmm ./ ndm),
            lw=copy(acc.lwmm ./ ndm),
            qlat=copy(acc.qlatmm ./ ndm),
            qsens=copy(acc.qsensmm ./ ndm)
        ))
        reset!(acc)
        mon = mon == 12 ? 1 : mon + 1
    end
    return (irec=irec, mon=mon)
end

# ── notebook cell 33fa7b1f-938b-481c-bf25-eca8d7fb33a7  (orig lines 2438-2500) ──
function time_loop!(it, year, CO2, mon, irec, Ts, Ta, q, To, output_buf,
    ws::CirculationWorkspace, acc::MonthlyAccumulator,
    timestate, cfg::PhysicsConfig)
    # Calendar lookup
    cal = it <= max_timesteps ? calendar_lookup[it] : (
        day=mod((it - 1) ÷ ndt_days, ndays_yr) + 1,
        step=mod(it - 1, nstep_yr) + 1
    )
    timestate.jday = cal.day
    timestate.ityr = cal.step
    ityr = timestate.ityr

    # Compute tendencies
    tend = tendencies!(CO2, Ts, Ta, To, q, ws, timestate, cfg)

    # Correction views
    TF_corr = @view TF_correct[:, :, ityr]
    qF_corr = @view qF_correct[:, :, ityr]
    ToF_corr = @view ToF_correct[:, :, ityr]

    # Surface temperature
    @. Ts = Ts + tend.dT_ocean + Δt * (tend.SW + tend.LW_surf - tend.LW_down +
                                       tend.Q_lat + tend.Q_sens + TF_corr) / cap_surf

    # Air temperature
    @. Ta = Ta + tend.dTa_crcl + Δt * (tend.LW_up + tend.LW_down - tend.em *
                                                                   tend.LW_surf + tend.Q_lat_air - tend.Q_sens) / cap_air

    # Clamps
    @. Ts = max(Ts, min_T_K)
    @. Ta = max(Ta, min_T_K)

    # Deep ocean
    @. To = To + tend.dTo + ToF_corr

    # Humidity (with clamp)
    dq_eva_use = cfg.log_hydro_dmc ? tend.dq_eva : ws.eva
    dq_rain_use = cfg.log_hydro_dmc ? tend.dq_rain : ws.rain
    dq_crcl_use = cfg.log_crcl_dmc ? tend.dq_crcl : ws.crcl
    @. q = q + clamp(Δt * (dq_eva_use + dq_rain_use) + dq_crcl_use + qF_corr,
        -min_humidity_change * q, max_humidity_change)

    # Sea ice heat capacity
    seaice!(Ts, timestate, cfg)

    # Conversion to mm/day (analysis units)
    @. ws.precip_out = (-dq_rain_use) * wz_vapor * conv_factor
    @. ws.evap_out = dq_eva_use * wz_vapor * conv_factor
    @. ws.qcrcl_out = dq_crcl_use


    # Output and diagnostics
    (mon, irec) = output!(it, irec, mon, Ts, Ta, To, q, tend.albedo,
        tend.ice_cover, ws.precip_out, ws.evap_out, ws.qcrcl_out,
        tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens,
        output_buf, acc, timestate)
    diagnostics!(it, year, CO2, Ts, Ta, To, q, tend.albedo,
        tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens, timestate)

    return (mon=mon, irec=irec)
end

# ── notebook cell 4f97badf-a501-4c70-a943-d3c86b48f8a1  (orig lines 2501-2537) ──
function build_monthly_climatology(records::Vector{MonthlyRecord})::Vector{MonthlyRecord}
    isempty(records) && return MonthlyRecord[]

    fields = propertynames(records[1])
    nmonths = 12
    counts = zeros(Int, nmonths)
    clim_acc = [Dict{Symbol,Matrix{Float64}}() for _ in 1:nmonths]

    # Accumulate by mont
    for (idx, rec) in enumerate(records)
        mon = mod(idx - 1, nmonths) + 1
        counts[mon] += 1
        for fld in fields
            arr = getfield(rec, fld)
            if haskey(clim_acc[mon], fld)
                clim_acc[mon][fld] .+= arr
            else
                clim_acc[mon][fld] = copy(arr)
            end
        end
    end

    # Build climatology
    clim = MonthlyRecord[]
    for mon in 1:nmonths
        if counts[mon] == 0
            push!(clim, records[1])
        else
            push!(clim, NamedTuple{fields}(
                Tuple(clim_acc[mon][fld] ./ counts[mon] for fld in fields)
            ))
        end
    end
    return clim
end

# ── notebook cell bf5ccbb7-fff9-42bc-8f2f-2e628b448b3a  (orig lines 2538-2554) ──
function apply_scenario_anomalies(scnr_records::Vector{MonthlyRecord}, ctrl_clim::Vector{MonthlyRecord})::Vector{MonthlyRecord}
    isempty(scnr_records) && return scnr_records
    isempty(ctrl_clim) && return scnr_records

    fields = propertynames(scnr_records[1])
    anom = MonthlyRecord[]

    for (idx, rec) in enumerate(scnr_records)
        mon = mod(idx - 1, 12) + 1
        ref = ctrl_clim[mon]
        push!(anom, NamedTuple{fields}(
            Tuple(getfield(rec, fld) .- getfield(ref, fld) for fld in fields)))
    end
    return anom
end

# ── notebook cell 867b193e-3390-4b3b-b5fe-5c399a250660  (orig lines 2587-2606) ──
function compute_annual_ice_climatology(ctrl_output::Vector{MonthlyRecord})
    ice_months = zeros(Float64, xdim, ydim, 12)
    count = zeros(Int, 12)

    for (idx, rec) in enumerate(ctrl_output)
        mon = mod1(idx, 12)          # month 1..12 based on record index
        @. ice_months[:, :, mon] += rec.ice   # note: field name is `ice`, not `ice_cover`
        count[mon] += 1
    end

    for mon in 1:12
        cnt = count[mon]
        cnt > 0 || continue
        @. ice_months[:, :, mon] *= 1.0 / cnt
    end

    return ice_months
end

# ── notebook cell 584e767e-4dc5-4821-af63-d6a825326d9e  (orig lines 2869-2928) ──
function qflux_correction!(CO2_ctrl, Ts, Ta, q, To, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace, time_flux)
    for it in 1:(time_flux*ndt_days*ndays_yr)
        timestate.jday = mod((it - 1) ÷ ndt_days, ndays_yr) + 1
        timestate.ityr = mod(it - 1, nstep_yr) + 1
        ityr = timestate.ityr

        tend = tendencies!(CO2_ctrl, Ts, Ta, To, q, ws, timestate, cfg)

        # Views into climatology & correction fields
        Tc = @view Tclim[:, :, ityr]
        Toc = @view Toclim[:, :, ityr]
        qc = @view qclim[:, :, ityr]
        TFc = @view TF_correct[:, :, ityr]
        ToFc = @view ToF_correct[:, :, ityr]
        qFc = @view qF_correct[:, :, ityr]

        # ── Surface temperature ──────────────────────────────
        # Uncorrected state (store in workspace buffer)
        @. ws.Ts0_buf = Ts + tend.dT_ocean + Δt * (
            tend.SW + tend.LW_surf - tend.LW_down +
            tend.Q_lat + tend.Q_sens
        ) / cap_surf

        # Correction and corrected state
        @. TFc = (Tc - ws.Ts0_buf) * cap_surf / Δt
        @. ws.Ts0_buf = ws.Ts0_buf + TFc * Δt / cap_surf

        # ── Air temperature ──────────────────────────────────
        @. ws.Ta0_buf = Ta + tend.dTa_crcl + ΔT_AIR_FACTOR * (
            tend.LW_up + tend.LW_down - tend.em * tend.LW_surf +
            tend.Q_lat_air - tend.Q_sens
        )

        # ── Deep ocean ───────────────────────────────────────
        @. ws.To0_buf = To + tend.dTo
        @. ToFc = Toc - ws.To0_buf
        @. ws.To0_buf = ws.To0_buf + ToFc

        # ── Humidity ─────────────────────────────────────────
        @. ws.q0_buf = q + tend.dq_crcl + Δt * (tend.dq_eva + tend.dq_rain)
        @. qFc = qc - ws.q0_buf
        @. ws.q0_buf = ws.q0_buf + qFc

        # Sea ice (updates cap_surf in place)
        seaice!(ws.Ts0_buf, timestate, cfg)

        # Diagnostics
        diagnostics!(it, 0.0, CO2_ctrl, ws.Ts0_buf, ws.Ta0_buf, ws.To0_buf, ws.q0_buf,
            tend.albedo, tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens, timestate)

        # Advance state (broadcast – same as `.=`)
        @. Ts = ws.Ts0_buf
        @. Ta = ws.Ta0_buf
        @. q = ws.q0_buf
        @. To = ws.To0_buf
    end
    return nothing
end

# ── notebook cell 92c3bd68-bd07-4381-9c04-e6611650cd1e  (orig lines 2929-3043) ──
function greb_model!(time_flux, time_ctrl, time_scnr, cfg::PhysicsConfig; jdal2_dir::AbstractString="")

    # ── 1. Initialisation ───────────────────────────────────────
    ini = init_model!(cfg)
    Ts_ini = ini.Ts_ini;
    Ta_ini = ini.Ta_ini
    To_ini = ini.To_ini;
    q_ini = ini.q_ini
    CO2_ctrl = ini.CO2_ctrl

    # Determine experiment type from cfg
    is_orbital_exp = cfg.experiment in (:obliquity, :eccentricity, :earth_sun_distance)
    is_forced_boundary = cfg.experiment in (:rcp85, :elnino, :lanina)
    is_sst_plus1 = cfg.experiment == :sst_plus1

    # Workspace and accumulator
    ws = CirculationWorkspace()
    acc = MonthlyAccumulator()

    # Initialize time state
    timestate = TimeState(1, 1)

    # ── 2. Flux-correction spin-up ──────────────────────────────
    if cfg.log_topo_drsp || cfg.log_qflux_dmc
        if !cfg.log_topo_drsp && cfg.log_qflux_dmc
            println("% loading flux correction fields...")
            load_flux_corrections_jdal2!(jdal2_dir)
        end
        println("% flux correction  CO2 = ", CO2_ctrl)
        qflux_correction!(CO2_ctrl, Ts_ini, Ta_ini, q_ini, To_ini, timestate, cfg, ws, time_flux)
    else
        println("Flux correction skipped")
    end

    # Reset accumulators after spin-up
    reset!(acc)

    # ── 3. Control run ──────────────────────────────────────────
    println("CONTROL RUN: CO2 = ", CO2_ctrl, " time = ", time_ctrl, " yr")

    # Initialize state arrays
    Ts = copy(Ts_ini);
    Ta = copy(Ta_ini)
    To = copy(To_ini);
    q = copy(q_ini)
    sw_solar_forcing_state[] = 1.0
    mon = 1;
    year = 1970;
    irec = 0

    ctrl_output = MonthlyRecord[]
    sizehint!(ctrl_output, time_ctrl * 12)  # Pre-allocate for all months
    timestate = TimeState(1, 1)  # Initialize time state

    for it in 1:(time_ctrl*nstep_yr)
        (mon, irec) = time_loop!(it, year, CO2_ctrl, mon, irec,
            Ts, Ta, q, To, ctrl_output, ws, acc, timestate, cfg)
        if mod(it, nstep_yr) == 0
            year += 1
        end
    end

    # ── Build ice climatology from control output ───────────────
    # Compute annual mean ice cover from the stored control monthly means
    # (Assuming control output contains monthly ice cover; adjust if needed)
    ice_forcing = compute_annual_ice_climatology(ctrl_output)

    # ── 4. Scenario run ─────────────────────────────────────────
    println("SCENARIO: ", cfg.experiment, "  time = ", time_scnr, " yr")

    # Reset state to initial conditions
    Ts .= Ts_ini;
    Ta .= Ta_ini
    q .= q_ini;
    To .= To_ini
    year = is_orbital_exp ? 1 : 1950
    CO2 = 340.0;
    mon = 1;
    irec = 0

    sw_solar_forcing_state[] = 1.0
    reset!(acc)  # Use accumulator reset

    scnr_output = MonthlyRecord[]
    if time_scnr > 0
        sizehint!(scnr_output, time_scnr * 12)
    end

    for it in 1:(time_scnr*nstep_yr)
        # Obtain forcing (CO2 and solar multiplier)
        forcing_result = forcing(it, year, cfg, ice_forcing; nstep_yr=nstep_yr)
        CO2 = forcing_result.CO2
        sw_solar_forcing_state[] = forcing_result.sw_solar_forcing

        # Forced‑boundary experiments: overwrite Ts with climatology
        if is_forced_boundary
            ityr_now = mod(it - 1, nstep_yr) + 1
            Ts .= @view Tclim[:, :, ityr_now]
        end

        # SST+1 K experiment
        if is_sst_plus1
            CO2 = CO2_ctrl
            ityr_now = mod(it - 1, nstep_yr) + 1
            @. Ts = ifelse(z_topo < 0.0, Tclim[:, :, ityr_now] + 1.0, Ts)
        end

        (mon, irec) = time_loop!(it, year, CO2, mon, irec,
            Ts, Ta, q, To, scnr_output, ws, acc, timestate, cfg)

        if mod(it, nstep_yr) == 0
            year += 1
        end
    end

    # Post‑processing: anomalies for non‑orbital experiments
    if !is_orbital_exp && !isempty(ctrl_output) && !isempty(scnr_output)
        ctrl_clim = build_monthly_climatology(ctrl_output)
        scnr_output = apply_scenario_anomalies(scnr_output, ctrl_clim)
    end

    return (ctrl=ctrl_output, scnr=scnr_output)
end

end # module GREB
