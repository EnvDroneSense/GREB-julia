module GREBClimate

# =============================================================================
# GREB - Globally Resolved Energy Balance model
#
# A global climate model that steps surface/air/ocean temperature and
# humidity forward under shortwave/longwave radiation, hydrology, sea ice,
# deep-ocean coupling, and atmospheric circulation (diffusion, advection,
# moisture convergence). An interactive Pluto version of the same model
# lives in `notebooks/GREB_julia.jl`; see the package docs or README for
# usage.
#
# Files below are included in dependency order: constants → config → data
# → state → io → physics/{radiation,hydrology,ocean} → circulation →
# tendencies → output → postprocess → model.
# =============================================================================

using LoopVectorization   # @turbo SIMD
using JLD2
using DataDeps: DataDeps, DataDep, register, unpack, @datadep_str

export PhysicsConfig, RunSpec, CirculationWorkspace, MonthlyAccumulator, TimeState, MonthlyRecord
export ClimateFields, ModelState, SurfaceState
export greb_data_dir
export read_jld2, load_solar_forcing_jld2, load_flux_corrections_jld2!, load_greb_jld2!
export load_co2_scenario_jld2, load_custom_co2_scenario, load_cc_anomaly_jld2!, load_enso_anomaly_jld2!
export create_experiment_config, set_hydrology_parameters!, init_model!
export SWradiation!, LWradiation!, hydro!, convergence!, seaice!, deep_ocean!
export diffusion!, advection!, circulation!, tendencies!, forcing
export diagnostics!, output!, time_loop!
export build_monthly_climatology, apply_scenario_anomalies, compute_annual_ice_climatology
export qflux_correction!, greb_model!
export xdim, ydim, nstep_yr

include("constants.jl")
include("config.jl")
include("data.jl")
include("state.jl")
include("io.jl")
include("physics/radiation.jl")
include("physics/hydrology.jl")
include("physics/ocean.jl")
include("circulation.jl")
include("tendencies.jl")
include("output.jl")
include("postprocess.jl")
include("model.jl")

function __init__()
    # Registration only: nothing is downloaded until `greb_data_dir()` has to
    # fall through to the DataDep, which cannot happen during precompilation.
    register_greb_datadep()
    return nothing
end

using PrecompileTools: @compile_workload
using Logging: with_logger, NullLogger

@compile_workload begin
    with_logger(NullLogger()) do
        redirect_stdout(devnull) do
            cfg = create_experiment_config(:full_model)
            greb_model!(RunSpec(scnr=0), cfg; jld2_dir="", allow_uninitialized=true)

            cfg_eva0 = create_experiment_config(:full_model)
            cfg_eva0.log_eva = 0
            greb_model!(RunSpec(scnr=0), cfg_eva0; jld2_dir="", allow_uninitialized=true)
        end
    end
end

end # module GREBClimate
