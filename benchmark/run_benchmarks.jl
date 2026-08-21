# =============================================================================
# run_benchmarks.jl - dependency-free timing harness for GREBClimate.jl
#
# Four benchmarking modes, all against the real dataset, post-JIT-warmup:
#
#   year     - wall-clock time for a real 1-simulated-year `:full_model`
#              control run (the original/default benchmark).
#   stages   - per-timestep breakdown: times each of tendencies!'s component
#              stages (circulation Ta/q, SW/LW radiation, hydro, deep ocean)
#              individually and reports each one's share of the total.
#   threads  - runs `year` in separate `-t N` subprocesses (thread count is
#              fixed at Julia startup, so this can't be done in-process) and
#              reports relative speedup vs. `-t 1`.
#   alloc    - reports bytes allocated by one `tendencies!` call, a cheap
#              regression guard for the zero-allocation hot path.
#
# Usage:
#   julia --project=. -t 2 benchmark/run_benchmarks.jl [mode] [jld2_dir] [reps]
#   (mode defaults to "year"; jld2_dir defaults to GREB_DATA or
#   ../greb_input_data; reps only applies to "year", default 3)
#
#   -t 2, not -t 3, is the recommended thread count for `year` - see the
#   `sweep_threads` docstring below for why.
#
#   julia --project=. -t 2 benchmark/run_benchmarks.jl                     # year, default dir
#   julia --project=. -t 2 benchmark/run_benchmarks.jl greb_input_data   # year (legacy call form)
#   julia --project=. -t 2 benchmark/run_benchmarks.jl stages
#   julia --project=. -t 1 benchmark/run_benchmarks.jl threads
#   julia --project=. -t 1 benchmark/run_benchmarks.jl alloc
#
#   (or from a REPL: include("benchmark/run_benchmarks.jl"); time_1yr(dir))
#
# `deepcopy(fields)` (needed because `greb_model!` mutates `fields` in place)
# is deliberately kept OUTSIDE every timed region in `time_1yr` - it's a
# large, thread-count-independent constant that would otherwise swamp the
# real per-run signal. Runs are repeated because wall-clock timing on a
# shared/dev machine has real run-to-run noise; report the mean and range,
# not a single reading - see `.claude/skills/benchmark/SKILL.md` for the
# known noise sources on this machine.
# =============================================================================

using GREBClimate

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

"""
    time_stages(jld2_dir; cfg=create_experiment_config(:full_model), reps=2000)

Breaks one timestep's `tendencies!` pipeline into its component stages
(circulation for `Ta` and `q`, SW/LW radiation, `hydro!`, `deep_ocean!`) and
times each individually - `reps` back-to-back calls after a warm-up,
averaged - reporting per-call time and each stage's share of the pipeline
total. This is the same per-stage breakdown methodology
`claude/IMPROVEMENTS.md` §2.10 used to establish circulation's ~65-75%
share, now scripted instead of ad hoc, so it can be re-run after any change
that might shift the balance (e.g. a kernel-specific optimization).

Threads are not spawned here (single-workspace calls only) even if
`Threads.nthreads() > 1` - this measures each stage's serial cost, which is
what `tendencies!`'s parallel path is actually racing against.
"""
function time_stages(jld2_dir::AbstractString; cfg=create_experiment_config(:full_model), reps::Int=2000)
    if !isdir(jld2_dir)
        @warn "JLD2 data directory not found: $jld2_dir - see DATA_README.md"
        return nothing
    end

    fields = load_greb_jld2!(jld2_dir; dataset=:ncep)
    init_model!(cfg, fields)
    state = ModelState()
    ws = CirculationWorkspace()
    timestate = TimeState(1, 1)
    CO2 = cfg.co2_concentration

    Ts = copy(fields.Tclim[:, :, end])
    Ta = copy(Ts)
    To = copy(fields.Toclim[:, :, end])
    q = copy(fields.qclim[:, :, end])

    stages = [
        ("circulation!(Ta)", () -> circulation!(Ta, GREBClimate.z_air, ws.dTa_crcl, fields, ws, timestate, cfg)),
        ("circulation!(q)", () -> circulation!(q, GREBClimate.z_vapor, ws.dq_crcl, fields, ws, timestate, cfg)),
        ("SWradiation!", () -> SWradiation!(Ts, fields, state, timestate, cfg, ws)),
        ("LWradiation!", () -> LWradiation!(Ts, Ta, q, CO2, fields, timestate, cfg, ws)),
        ("hydro!", () -> hydro!(Ts, q, fields, timestate, cfg, ws)),
        ("deep_ocean!", () -> deep_ocean!(Ts, To, fields, timestate, cfg, ws)),
    ]

    # Warm-up: excludes JIT compilation from every timed loop below.
    for (_, f) in stages
        f()
    end

    println("Per-stage timing (", reps, " calls each, single workspace, ",
        Threads.nthreads(), " thread(s) available but unused here):")
    results = Tuple{String,Float64}[]
    for (name, f) in stages
        t = @elapsed for _ in 1:reps
            f()
        end
        push!(results, (name, t / reps))
    end

    total = sum(last(r) for r in results)
    for (name, per_call) in results
        share = 100 * per_call / total
        println("  ", rpad(name, 18), round(per_call * 1e6, digits=2), " µs/call   (",
            round(share, digits=1), "% of pipeline total)")
    end
    println("  pipeline total (sum of stages): ", round(total * 1e6, digits=1), " µs")
    println("  (excludes convergence!/seaice!/output!/diagnostics! - see time_1yr for whole-model cost)")
    return results
end

"""
    sweep_threads(jld2_dir; thread_counts=(1,2,3,4), reps=3)

Runs `time_1yr` in separate `julia -t N` subprocesses for each `N` in
`thread_counts` - thread count is fixed at Julia startup, so this can't be
swept in-process - and reports relative speedup vs. `-t 1`.

`tendencies!` spawns a 3-lane split (`circulation!(Ta)` ‖ `circulation!(q)`
‖ everything else), but since the `Float32` conversion
(`claude/IMPROVEMENTS.md` §2.3) circulation is ~98% of per-timestep cost and
"everything else" is ~2% - too little work to justify a dedicated third
thread. `-t 2` already gets nearly all the real parallelism (the main task
finishes its ~30µs of synchronous work and blocks on `wait()`, at which
point Julia's work-stealing scheduler runs whichever spawned circulation
task hasn't started yet, so the two big tasks still overlap). Measured
across 5 independent sweeps on this machine: `-t 2` gave a consistent
~1.5× (range 1.24-1.62×, low variance); `-t 3`/`-t 4` ranged ~1.0-1.65×
(high variance, no longer a dependable edge over `-t 2`) - more active
threads means more exposure to background system load when there's this
little real work to parallelize. Run this a few times before trusting a
single `-t 3`/`-t 4` reading either way; `-t 2` vs `-t 1` is the comparison
that reliably holds up in one run.
"""
function sweep_threads(jld2_dir::AbstractString; thread_counts=(1, 2, 3, 4), reps::Int=3)
    project_dir = joinpath(@__DIR__, "..")
    julia_bin = joinpath(Sys.BINDIR, Base.julia_exename())
    script = @__FILE__

    results = Dict{Int,Float64}()
    for n in thread_counts
        println("--- -t $n ---")
        cmd = `$julia_bin --project=$project_dir -t $n $script year $jld2_dir $reps`
        output = read(cmd, String)
        print(output)
        m = match(r"mean:\s*([\d.]+)\s*s", output)
        if m === nothing
            @warn "Could not parse mean time from -t $n run"
        else
            results[n] = parse(Float64, m.captures[1])
        end
    end

    base_n = first(thread_counts)
    if haskey(results, base_n)
        base = results[base_n]
        println("\nSpeedup vs -t $base_n:")
        for n in thread_counts
            haskey(results, n) || continue
            println("  -t $n: ", round(results[n], digits=3), "s  (", round(base / results[n], digits=2), "x)")
        end
    end
    return results
end

"""
    check_allocations(jld2_dir)

Reports bytes allocated by one `tendencies!` call (after a warm-up call to
exclude JIT/compilation allocations) - a cheap regression guard for the
zero-allocation-per-timestep invariant established in
`claude/IMPROVEMENTS.md` §2.2/§2.5/§2.7. Uses the default `ws_a=ws_q=ws`
(synchronous, non-`Threads.@spawn`) path, since spawning tasks allocates
regardless of the physics code's own allocation behavior.
"""
function check_allocations(jld2_dir::AbstractString)
    if !isdir(jld2_dir)
        @warn "JLD2 data directory not found: $jld2_dir - see DATA_README.md"
        return nothing
    end

    cfg = create_experiment_config(:full_model)
    fields = load_greb_jld2!(jld2_dir; dataset=:ncep)
    init_model!(cfg, fields)
    state = ModelState()
    ws = CirculationWorkspace()
    timestate = TimeState(1, 1)
    CO2 = cfg.co2_concentration

    Ts = copy(fields.Tclim[:, :, end])
    Ta = copy(Ts)
    To = copy(fields.Toclim[:, :, end])
    q = copy(fields.qclim[:, :, end])

    tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, timestate, cfg)  # warm-up
    bytes = @allocated tendencies!(CO2, Ts, Ta, To, q, fields, state, ws, timestate, cfg)

    println("tendencies! allocations (single-workspace path): ", bytes, " bytes")
    return bytes
end

const DEFAULT_DATA_DIR = get(ENV, "GREB_DATA", joinpath(@__DIR__, "..", "greb_input_data"))

const _MODES = ("year", "stages", "threads", "alloc")

# Run automatically when executed as a script, but NOT when `include`-d into
# an interactive session - so a REPL is never terminated. Backward-compatible
# with the original single-argument call form (`... run_benchmarks.jl <dir>`)
# - ARGS[1] is only treated as a mode if it's one of `_MODES`.
if abspath(PROGRAM_FILE) == @__FILE__
    mode, rest = !isempty(ARGS) && ARGS[1] in _MODES ? (ARGS[1], ARGS[2:end]) : ("year", ARGS)
    jld2_dir = !isempty(rest) ? rest[1] : DEFAULT_DATA_DIR

    if mode == "year"
        reps = length(rest) >= 2 ? parse(Int, rest[2]) : 3
        time_1yr(jld2_dir; reps=reps)
    elseif mode == "stages"
        time_stages(jld2_dir)
    elseif mode == "threads"
        sweep_threads(jld2_dir)
    elseif mode == "alloc"
        check_allocations(jld2_dir)
    end
end
