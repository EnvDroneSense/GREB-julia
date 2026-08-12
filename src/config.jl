"""
    PhysicsConfig

All model switches and parameters: mean-climate/CO₂-response toggles,
circulation components, hydrology parameters, external-forcing flags, and
the current experiment type. Passed explicitly to every physics function;
build one with [`create_experiment_config`](@ref) rather than the bare
keyword constructor for anything beyond `:full_model`.
"""
Base.@kwdef mutable struct PhysicsConfig
    # Mean Climate Switches
    log_clouds_dmc::Bool = true
    log_vapor_dmc::Bool = true
    log_crcl_dmc::Bool = true
    log_hydro_dmc::Bool = true
    log_atmos_dmc::Bool = true
    log_co2_dmc::Bool = true
    log_ocean_dmc::Bool = true
    log_qflux_dmc::Bool = true

    # CO₂ Response Switches
    log_clouds_drsp::Bool = true
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

    # External anomaly forcing gates
    log_tsurf_ext::Bool = false
    log_hwind_ext::Bool = false
    log_omega_ext::Bool = false

    # Experiment Type
    experiment::Symbol = :full_model  # :full_model, :constant_topo, :co2_double, etc.

    # CO₂ concentration for experiments (ppm)
    co2_concentration::Float64 = 340.0

    # Orbital-forcing experiments: which solar_scenarios table row to load
    # (:eccentricity / :obliquity, see load_solar_forcing_jld2)
    orbital_index::Int = 0

    # Earth-Sun distance experiment: percent change in orbital radius
    earth_sun_distance_pct::Float64 = 0.0

    # Hydrology parameters (calculated by set_hydrology_parameters!)
    c_q::Float64 = 1.0
    c_rq::Float64 = 0.0
    c_omega::Float64 = 0.0
    c_omegastd::Float64 = 0.0

    # IPCC scenario CO₂ lookup table (year => ppm)
    co2_scenario::Dict{Int,Float64} = Dict{Int,Float64}()

    # :custom_co2 experiment: path to a user-supplied "year CO2" text file
    # (same format as the IPCC scenario files), loaded into `co2_scenario`
    # at scenario start.
    custom_co2_path::String = ""
end

"""
    RunSpec

Run durations (in years) for [`greb_model!`](@ref): `flux` (flux-correction
spin-up), `ctrl` (control run), `scnr` (scenario run). A keyword struct
instead of three bare positional ints, whose order was easy to swap by
mistake.
"""
Base.@kwdef struct RunSpec
    flux::Int = 0
    ctrl::Int = 1
    scnr::Int = 1
end

begin
    """
        create_experiment_config(experiment::Symbol; co2_path="",
            log_clouds_dmc=true, log_ocean_dmc=true, log_atmos_dmc=true,
            log_co2_dmc=true, log_hydro_dmc=true, log_qflux_dmc=true,
            log_topo_drsp=true, log_clouds_drsp=true, log_humid_drsp=true,
            log_ocean_drsp=true, log_hydro_drsp=true, log_ice=true, log_hdif=true,
            log_hadv=true, log_vdif=true, log_vadv=true) -> PhysicsConfig

    Create a pre-configured `PhysicsConfig` for common experiment types.
    The `log_*` keywords are only used by `:decon_mean_climate`/`:decon_2xco2`
    (see below); `co2_path` is only used by `:custom_co2`.

    # Experiments
    - `:full_model` - All processes active (default)
    - `:constant_topo` - Constant topography (log_topo_drsp = false)
    - `:co2_double` - 2×CO₂ (680 ppm)
    - `:co2_quadruple` - 4×CO₂ (1360 ppm)
    - `:solar_plus27` - +27 W/m² solar constant
    - `:elnino` - El Niño conditions
    - `:lanina` - La Niña conditions
    - `:paleo_231kyr` - Paleoclimate (200 ppm CO₂)
    - `:rcp26`/`:rcp45`/`:rcp60`/`:rcp85` - IPCC RCP climate change scenarios
    - `:custom_co2` - user-supplied CO₂ trajectory, `co2_path` keyword gives
      the "year CO2" text file path (see [`load_custom_co2_scenario`](@ref))
    - `:ssp119`/`:ssp126`/`:ssp245`/`:ssp460`/`:ssp585`
    - `:historical_co2` - Observed CO₂ 1850–2017 (year starts at 1850, not 1950)
    - `:decon_mean_climate` - deconstruct-mean-state experiment 
    - `:decon_2xco2` - deconstruct-2×CO₂-response experiment 
    """
    function create_experiment_config(experiment::Symbol; co2_path::AbstractString="",
        log_clouds_dmc::Bool=true, log_ocean_dmc::Bool=true, log_atmos_dmc::Bool=true,
        log_co2_dmc::Bool=true, log_hydro_dmc::Bool=true, log_qflux_dmc::Bool=true,
        log_topo_drsp::Bool=true, log_clouds_drsp::Bool=true, log_humid_drsp::Bool=true,
        log_ocean_drsp::Bool=true, log_hydro_drsp::Bool=true, log_ice::Bool=true, log_hdif::Bool=true,
        log_hadv::Bool=true, log_vdif::Bool=true, log_vadv::Bool=true)::PhysicsConfig
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
            return PhysicsConfig(experiment=:solar_plus27)

        elseif experiment == :elnino
            cfg = PhysicsConfig(experiment=:elnino)
            cfg.log_tsurf_ext = cfg.log_hwind_ext = cfg.log_omega_ext = true
            return cfg

        elseif experiment == :lanina
            cfg = PhysicsConfig(experiment=:lanina)
            cfg.log_tsurf_ext = cfg.log_hwind_ext = cfg.log_omega_ext = true
            return cfg

        elseif experiment == :paleo_231kyr
            cfg = PhysicsConfig(experiment=:paleo_231kyr)
            cfg.co2_concentration = 200.0
            return cfg

        elseif experiment == :rcp85
            cfg = PhysicsConfig(experiment=:rcp85)
            cfg.log_tsurf_ext = cfg.log_hwind_ext = cfg.log_omega_ext = true
            return cfg

        elseif experiment in (:rcp26, :rcp45, :rcp60, :ssp119, :ssp126, :ssp245, :ssp460, :ssp585)
            cfg = PhysicsConfig(experiment=experiment)
            return cfg

        elseif experiment == :historical_co2
            cfg = PhysicsConfig(experiment=:historical_co2)
            return cfg

        elseif experiment == :custom_co2
            cfg = PhysicsConfig(experiment=:custom_co2)
            cfg.custom_co2_path = co2_path
            return cfg

        elseif experiment == :decon_mean_climate
            cfg = PhysicsConfig(experiment=:decon_mean_climate)
            cfg.log_clouds_dmc = log_clouds_dmc
            cfg.log_ocean_dmc = log_ocean_dmc
            cfg.log_atmos_dmc = log_atmos_dmc
            cfg.log_co2_dmc = log_co2_dmc
            cfg.log_hydro_dmc = log_hydro_dmc
            cfg.log_qflux_dmc = log_qflux_dmc
            cfg.log_ice = log_ice
            cfg.log_hdif = log_hdif
            cfg.log_hadv = log_hadv
            cfg.log_vdif = log_vdif
            cfg.log_vadv = log_vadv
            return cfg

        elseif experiment == :decon_2xco2
            cfg = PhysicsConfig(experiment=:decon_2xco2, co2_concentration=680.0)
            cfg.log_topo_drsp = log_topo_drsp
            cfg.log_clouds_drsp = log_clouds_drsp
            cfg.log_humid_drsp = log_humid_drsp
            cfg.log_ocean_drsp = log_ocean_drsp
            cfg.log_hydro_drsp = log_hydro_drsp
            cfg.log_ice = log_ice
            cfg.log_hdif = log_hdif
            cfg.log_hadv = log_hadv
            cfg.log_vdif = log_vdif
            cfg.log_vadv = log_vadv
            return cfg

        else
            error("Unknown experiment: $experiment")
        end
    end
end;

# - Optimized Parameter Initialization ──────────────────────────────────
"""
    set_hydrology_parameters!(cfg::PhysicsConfig)

Initialize precipitation parameters `c_q, c_rq, c_omega, c_omegastd` based on
`cfg.log_rain` and `cfg.log_clim` settings.
"""
function set_hydrology_parameters!(cfg::PhysicsConfig)
    # Fast lookup instead of if-else chain
    haskey(HYDRO_PARAMS, cfg.log_rain) ||
        error("Unknown log_rain value: $(cfg.log_rain). Valid values: $(sort(collect(keys(HYDRO_PARAMS))))")
    params = HYDRO_PARAMS[cfg.log_rain]
    cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd = params

    # NCEP parameter adjustment
    if cfg.log_rain == 0 && cfg.log_clim == 1
        cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd = -1.27, 1.99, -16.54, 21.15
    end

    @info "⚙️ MSCM hydrology: log_rain=$(cfg.log_rain), log_clim=$(cfg.log_clim) → (c_q=$(cfg.c_q), c_rq=$(cfg.c_rq), c_omega=$(cfg.c_omega), c_omegastd=$(cfg.c_omegastd))"
end
