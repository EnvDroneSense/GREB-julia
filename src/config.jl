# ── notebook cell 19f106e4-2b82-47e1-9284-799a105f30cb  (orig lines 235-291) ──
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

    # Orbital-forcing experiments: which solar_scenarios table row to load
    # (:eccentricity / :obliquity, see load_solar_forcing_jld2)
    orbital_index::Int = 0

    # Earth-Sun distance experiment: percent change in orbital radius
    # (Fortran's `dradius`; rS0 = (1/(1+0.01*earth_sun_distance_pct))^2)
    earth_sun_distance_pct::Float64 = 0.0

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

# ── notebook cell d0d74213-96ad-4a45-9a31-37f640a21a45  (orig lines 378-400) ──
# ⚙️ Optimized Parameter Initialization ──────────────────────────────────
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
