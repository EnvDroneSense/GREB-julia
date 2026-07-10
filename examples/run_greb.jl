# =============================================================================
# run_greb.jl — plain-Julia driver reproducing the GREB_julia.jl notebook flow
# on top of the extracted `GREB` package (no Pluto / no @bind widgets).
#
# Mirrors the notebook's non-interactive path:
#   1. load JDAL2 input data      (notebook data-loading cell)
#   2. build a PhysicsConfig       (replaces `current_physics_config()` widgets)
#   3. run the model               (replaces the `run_toggle` execute cell)
#   4. print a summary             (the notebook's diagnostics cell)
#   5. plot global-mean Ts          (the notebook's plotting cell)
#
# Two ways to run:
#   * as a script:  julia --project=. examples/run_greb.jl [path/to/greb_dataset_jdal2]
#   * from a REPL:   include("examples/run_greb.jl"); run_greb("path/to/greb_dataset_jdal2")
#
# NOTE: the JDAL2 input data is large and not committed (see DATA_README.md);
# supply its directory as the argument or via the GREB_DATA env var. This script
# never calls `exit()`, so `include`-ing it from a REPL will not kill the session.
# =============================================================================

using GREB
using Statistics

"""
    run_greb(jdal2_dir; time_flux=0, time_ctrl=1, time_scnr=1)

Load the JDAL2 dataset from `jdal2_dir`, run a GREB control+scenario simulation,
print a summary and (if Plots.jl is available) save a global-mean Ts plot.
Returns the result NamedTuple, or `nothing` if the data directory is missing.
"""
function run_greb(jdal2_dir::AbstractString;
                  time_flux::Int=0, time_ctrl::Int=1, time_scnr::Int=1)

    # ── 1. locate + load input data ─────────────────────────────────────────
    if !isdir(jdal2_dir)
        @warn """
        JDAL2 data directory not found: $jdal2_dir
        Pass it as an argument:  julia --project=. examples/run_greb.jl <dir>
        or set the GREB_DATA environment variable. See DATA_README.md.
        """
        return nothing
    end

    println("Loading GREB dataset from: ", jdal2_dir)
    load_greb_jdal2!(jdal2_dir; dataset=:ncep)

    # ── 2. configure the experiment (replaces the interactive widgets) ──────
    cfg = create_experiment_config(:full_model)

    # ── 3. run the model ────────────────────────────────────────────────────
    println("Running GREB (flux=$time_flux, ctrl=$time_ctrl, scnr=$time_scnr years)...")
    result = greb_model!(time_flux, time_ctrl, time_scnr, cfg; jdal2_dir=jdal2_dir)
    println("Run complete. control months: ", length(result.ctrl),
            ", scenario months: ", length(result.scnr))

    # ── 4. summary of the first control month ───────────────────────────────
    if !isempty(result.ctrl)
        rec = result.ctrl[1]
        println("\n" * "="^50)
        println("GREB MODEL OUTPUT (control month 1)")
        println("="^50)
        println("🌡️  Ts (K):   mean=$(round(mean(rec.Ts), digits=1))  " *
                "min=$(round(minimum(rec.Ts), digits=1))  max=$(round(maximum(rec.Ts), digits=1))")
        println("💧  precip:   mean=$(round(mean(rec.precip), digits=2))")
        println("☀️  SW (W/m²): mean=$(round(mean(rec.sw), digits=1))")
        println("❄️  ice:      mean=$(round(mean(rec.ice), digits=2))")
        println("✅ all finite: $(all(isfinite, rec.Ts))")
    end

    # ── 5. plot global-mean surface temperature (optional; needs Plots) ─────
    if !isempty(result.ctrl)
        Ts_global_mean = [mean(rec.Ts) for rec in result.ctrl]
        try
            @eval using Plots
            plt = Base.invokelatest(plot, Ts_global_mean;
                                    xlabel="Month", ylabel="Global Mean Ts [K]",
                                    legend=false, title="GREB control run")
            Base.invokelatest(savefig, plt, joinpath(@__DIR__, "greb_global_mean_Ts.png"))
            println("\nSaved plot to examples/greb_global_mean_Ts.png")
        catch err
            println("\n(Plots.jl not available — skipping plot; series has ",
                    length(Ts_global_mean), " points)")
        end
    end

    return result
end

# Default data directory: first CLI arg, else $GREB_DATA, else repo-root dataset.
const DEFAULT_JDAL2_DIR = get(ENV, "GREB_DATA",
    !isempty(ARGS) ? ARGS[1] : joinpath(@__DIR__, "..", "greb_dataset_jdal2"))

# Run automatically when executed as a script (`julia run_greb.jl`), but NOT when
# `include`-d into an interactive session — so a REPL is never terminated.
if abspath(PROGRAM_FILE) == @__FILE__
    run_greb(DEFAULT_JDAL2_DIR)
end
