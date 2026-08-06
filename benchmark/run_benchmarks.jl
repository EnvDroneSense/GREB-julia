# =============================================================================
# benchmark/run_benchmarks.jl — per-kernel BenchmarkTools micro-benchmarks
#
# Run with:
#   julia --project=benchmark benchmark/run_benchmarks.jl ["tag for this run"]
#
# Benchmarks GREB's hot per-timestep kernels against synthetic (or, if
# greb_dataset_jld2/ is present, real) climatology data, then APPENDS a dated,
# tagged section to claude/BENCHMARKS.md — the log is never overwritten, so
# re-running this after each optimization pass grows a before/after history.
# =============================================================================

using GREB
using BenchmarkTools
using Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))
const LOG_PATH = joinpath(ROOT, "claude", "BENCHMARKS.md")
const DATA_DIR = joinpath(ROOT, "greb_dataset_jld2")

tag = isempty(ARGS) ? "unlabeled run" : ARGS[1]

# ── Load real data if available, else benchmark against GREB's zero fields ──
using_real_data = isdir(DATA_DIR)
if using_real_data
    println("Loading real GREB dataset from $DATA_DIR ...")
    fields = load_greb_jld2!(DATA_DIR; dataset=:ncep)
else
    println("No greb_dataset_jld2/ found — benchmarking against zero-initialized fields.")
    fields = ClimateFields()
end

cfg = create_experiment_config(:full_model)
ws = CirculationWorkspace()
state = ModelState()
timestate = TimeState(1, 1)

# Physically-plausible synthetic state (not all-zero) so the kernels' ifelse
# branches (ice/no-ice, land/ocean, ...) actually get exercised either way.
T_test  = 288.0 .+ 5.0 .* randn(GREB.xdim, GREB.ydim)
Ta_test = copy(T_test)
To_test = fill(285.0, GREB.xdim, GREB.ydim)
q_test  = fill(0.006, GREB.xdim, GREB.ydim)

bench_data = (; T_test, Ta_test, To_test, q_test, fields, state, ws, timestate, cfg)

println("Running kernel benchmarks (this takes a minute or two)...")

results = Pair{String,BenchmarkTools.Trial}[]

push!(results, "convergence!" => @benchmark convergence!($bench_data.T_test, $bench_data.fields, $bench_data.timestate, $bench_data.ws) seconds=3)
push!(results, "diffusion!" => @benchmark diffusion!($bench_data.T_test, $(GREB.z_air), $bench_data.fields, $bench_data.ws, $bench_data.timestate) seconds=3)
push!(results, "advection!" => @benchmark advection!($bench_data.T_test, $(GREB.z_air), $bench_data.fields, $bench_data.ws, $bench_data.timestate, $bench_data.cfg) seconds=3)
push!(results, "SWradiation!" => @benchmark SWradiation!($bench_data.T_test, $bench_data.fields, $bench_data.state, $bench_data.timestate, $bench_data.cfg, $bench_data.ws) seconds=3)
push!(results, "LWradiation!" => @benchmark LWradiation!($bench_data.T_test, $bench_data.Ta_test, $bench_data.q_test, 340.0, $bench_data.fields, $bench_data.timestate, $bench_data.cfg, $bench_data.ws) seconds=3)
push!(results, "hydro!" => @benchmark hydro!($bench_data.T_test, $bench_data.q_test, $bench_data.fields, $bench_data.timestate, $bench_data.cfg, $bench_data.ws) seconds=3)
push!(results, "seaice!" => @benchmark seaice!($bench_data.T_test, $bench_data.fields, $bench_data.timestate, $bench_data.cfg) seconds=3)
push!(results, "deep_ocean!" => @benchmark deep_ocean!($bench_data.T_test, $bench_data.To_test, $bench_data.fields, $bench_data.timestate, $bench_data.cfg, $bench_data.ws) seconds=3)
push!(results, "circulation!" => @benchmark circulation!($bench_data.Ta_test, $(GREB.z_air), $bench_data.ws.dTa_crcl, $bench_data.fields, $bench_data.ws, $bench_data.timestate, $bench_data.cfg) seconds=3)
push!(results, "tendencies!" => @benchmark tendencies!(340.0, $bench_data.T_test, $bench_data.Ta_test, $bench_data.To_test, $bench_data.q_test, $bench_data.fields, $bench_data.state, $bench_data.ws, $bench_data.timestate, $bench_data.cfg) seconds=5)

for (name, trial) in results
    println(name, ": ", BenchmarkTools.prettytime(median(trial).time),
             " (min ", BenchmarkTools.prettytime(minimum(trial).time), "), ",
             BenchmarkTools.prettymemory(trial.memory), ", ", trial.allocs, " allocs")
end

# ── Append a dated, tagged section to claude/BENCHMARKS.md ─────────────────
function git_field(args...)
    try
        readchomp(Cmd([["git", "-C", ROOT]; collect(args)]))
    catch
        "unknown"
    end
end

branch = git_field("rev-parse", "--abbrev-ref", "HEAD")
sha = git_field("rev-parse", "--short", "HEAD")
dirty = !isempty(git_field("status", "--porcelain"))

mkpath(dirname(LOG_PATH))
is_new_file = !isfile(LOG_PATH)

open(LOG_PATH, "a") do io
    if is_new_file
        println(io, "# GREB.jl — Benchmark Log")
        println(io)
        println(io, "Running log of `BenchmarkTools.@benchmark` results for GREB's hot ",
                     "per-timestep kernels, appended by `benchmark/run_benchmarks.jl` ",
                     "each time it's run. Compare entries across runs to see the effect ",
                     "of each optimization pass.")
        println(io)
        println(io, "Reproduce: `julia --project=benchmark benchmark/run_benchmarks.jl \"<tag>\"`")
        println(io)
    end
    println(io, "---")
    println(io)
    println(io, "## ", Dates.format(now(), "yyyy-mm-dd HH:MM"), " — ", tag)
    println(io)
    println(io, "- branch: `", branch, "`  commit: `", sha, "`", dirty ? " (dirty working tree)" : "")
    println(io, "- data: ", using_real_data ? "real `greb_dataset_jld2/` (NCEP)" : "zero-initialized fields (no dataset present)")
    println(io, "- julia: `", VERSION, "`")
    println(io)
    println(io, "| kernel | min | median | memory | allocs |")
    println(io, "|---|---|---|---|---|")
    for (name, trial) in results
        t_min = BenchmarkTools.prettytime(minimum(trial).time)
        t_med = BenchmarkTools.prettytime(median(trial).time)
        mem = BenchmarkTools.prettymemory(trial.memory)
        println(io, "| `", name, "` | ", t_min, " | ", t_med, " | ", mem, " | ", trial.allocs, " |")
    end
    println(io)
end

println("\nAppended results to ", LOG_PATH)
