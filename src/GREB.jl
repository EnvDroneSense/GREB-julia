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
#
# The module used to be a single 2245-line file; it's since been split into
# topical files below (see claude/IMPROVEMENTS.md §1.2). The split is purely
# mechanical — code was moved, not rewritten — and `include()` order follows
# each file's actual dependencies (constants before everything; state before
# the io/physics files that read its globals; physics/circulation before
# tendencies; tendencies+postprocess before model).
# =============================================================================

using LoopVectorization   # @turbo SIMD
using JLD2

export PhysicsConfig, CirculationWorkspace, MonthlyAccumulator, TimeState, MonthlyRecord
export read_jld2, load_solar_forcing_jld2, load_flux_corrections_jld2!, load_greb_jld2!
export create_experiment_config, set_hydrology_parameters!, init_model!
export SWradiation!, LWradiation!, hydro!, convergence!, seaice!, deep_ocean!
export diffusion!, advection!, circulation!, tendencies!, forcing
export diagnostics!, output!, time_loop!
export build_monthly_climatology, apply_scenario_anomalies, compute_annual_ice_climatology
export qflux_correction!, greb_model!
export xdim, ydim, nstep_yr

include("constants.jl")
include("config.jl")
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

using PrecompileTools: @compile_workload

@compile_workload begin
    redirect_stdout(devnull) do
        cfg = create_experiment_config(:full_model)
        greb_model!(0, 1, 0, cfg; jld2_dir="")

        cfg_eva0 = create_experiment_config(:full_model)
        cfg_eva0.log_eva = 0
        greb_model!(0, 1, 0, cfg_eva0; jld2_dir="")
    end
end

end # module GREB
