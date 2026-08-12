# =============================================================================
# run_benchmarks.jl - simple wall-clock timing harness for GREB.jl
#
# Times a real 1-simulated-year `:full_model` control run (post-JIT-warmup,
# against the real dataset) for comparing performance across code, config,
# or thread-count changes. No dependencies beyond GREB itself.
#
# Usage:
#   julia --project=. -t 3 benchmark/run_benchmarks.jl [path/to/greb_dataset_jld2]
#   (or from a REPL: include("benchmark/run_benchmarks.jl"); time_1yr(dir))
#
# `deepcopy(fields)` (needed because `greb_model!` mutates `fields` in place)
# is deliberately kept OUTSIDE every timed region - it's a large,
# thread-count-independent constant that would otherwise swamp the real
# per-run signal. Each run is repeated a few times because wall-clock timing
# on a shared/dev machine has real run-to-run noise; report the mean and
# range, not a single reading.
# =============================================================================

using GREB

"""
    time_1yr(jld2_dir; cfg=create_experiment_config(:full_model), reps=3)

Runs a real 1-simulated-year `:full_model` control run `reps` times (after
one untimed warm-up call to exclude JIT compilation) and prints each timing
plus the mean/min/max. Returns the vector of timings in seconds.
"""
function time_1yr(jld2_dir::AbstractString; cfg=create_experiment_config(:full_model), reps::Int=3)
    if !isdir(jld2_dir)
        @warn "JLD2 data directory not found: $jld2_dir - see DATA_README.md"
        return Float64[]
    end

    println("Threads.nthreads() = ", Threads.nthreads())
    fields = load_greb_jld2!(jld2_dir; dataset=:ncep)

    # Warm-up: excludes JIT compilation from every timed run below.
    redirect_stdout(devnull) do
        greb_model!(RunSpec(scnr=0), cfg; jld2_dir=jld2_dir, fields=deepcopy(fields))
    end

    times = Float64[]
    for r in 1:reps
        fields_r = deepcopy(fields)  # outside the timed region, see header note
        t = @elapsed redirect_stdout(devnull) do
            greb_model!(RunSpec(scnr=0), cfg; jld2_dir=jld2_dir, fields=fields_r)
        end
        push!(times, t)
        println("  run $r: ", round(t, digits=3), " s")
    end

    println("mean: ", round(sum(times) / length(times), digits=3), " s  ",
        "(min ", round(minimum(times), digits=3), "s, max ", round(maximum(times), digits=3), "s)")
    return times
end

const DEFAULT_JLD2_DIR = get(ENV, "GREB_DATA",
    !isempty(ARGS) ? ARGS[1] : joinpath(@__DIR__, "..", "greb_dataset_jld2"))

# Run automatically when executed as a script, but NOT when `include`-d into
# an interactive session - so a REPL is never terminated.
if abspath(PROGRAM_FILE) == @__FILE__
    time_1yr(DEFAULT_JLD2_DIR)
end
