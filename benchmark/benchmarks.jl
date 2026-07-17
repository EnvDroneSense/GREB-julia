# =============================================================================
# benchmarks.jl — kernel-level performance suite for GREB
#
# Mirrors the notebook's `setup_benchmark()` cell: it times each physics kernel
# individually plus one full `greb_model!` control year, so performance work can
# be tracked systematically (change code → re-run → compare medians).
#
# Two ways to run:
#   # one-time setup of the benchmark environment
#   julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
#   # standalone — prints a table of median times / allocations
#   julia --project=benchmark benchmark/benchmarks.jl
#
#   # or via PkgBenchmark (defines the global `SUITE`):
#   using PkgBenchmark; benchmarkpkg("GREB")
#
# NOTE: kernels are exercised on representative *synthetic* fields (as in the
# notebook). Module-level climatology globals are zero-initialised (no JDAL2 data
# needed) — fine for *relative* comparisons, which is what perf tuning needs.
# Data-dependent branches (e.g. seaice!) may take an unrepresentative path, but
# consistently so across runs. See claude/IMPROVEMENTS.md for context.
# =============================================================================

using GREB
using BenchmarkTools

const SUITE = BenchmarkGroup()

# ── shared, read-only fixtures (built once) ──────────────────────────────────
const WS        = CirculationWorkspace()
const TIMESTATE = TimeState(1, 1)
const CFG       = PhysicsConfig()
const OMEGA     = rand(xdim, ydim, nstep_yr) .* 0.1   # cf. dummy omegaclim
const Z_AIR     = GREB.z_air
const Z_VAPOR   = GREB.z_vapor
const CO2       = 340.0

if size(WS.dX_crcl) != (xdim, ydim)
    WS.dX_crcl = zeros(Float64, xdim, ydim)
end

# fresh synthetic 2-D fields per sample (cheap: 96×48). Ranges mirror the notebook.
Ts_field()    = rand(xdim, ydim) .* 80  .+ 233.15   # 233–313 K
To_field()    = rand(xdim, ydim) .* 30  .+ 270.0    # 270–300 K
q_field()     = rand(xdim, ydim) .* 0.02            # 0–0.02 kg/kg
unit_field()  = rand(xdim, ydim)                    # generic tracer 0–1

# ── kernel benchmarks (one per physics routine) ──────────────────────────────
k = SUITE["kernels"] = BenchmarkGroup()

k["SWradiation!"] = @benchmarkable SWradiation!(Ts, $TIMESTATE, $CFG, $WS) setup = (Ts = Ts_field())

k["LWradiation!"] = @benchmarkable LWradiation!(Ts, Ta, q, $CO2, $TIMESTATE, $CFG, $WS) setup =
    (Ts = Ts_field(); Ta = Ts_field(); q = q_field())

k["hydro!"] = @benchmarkable hydro!(Ts, q, $TIMESTATE, $CFG, $WS) setup =
    (Ts = Ts_field(); q = q_field())

k["seaice!"] = @benchmarkable seaice!(Ts, $TIMESTATE, $CFG) setup = (Ts = Ts_field())

k["deep_ocean!"] = @benchmarkable deep_ocean!(Ts, To, $TIMESTATE, $CFG, $WS) setup =
    (Ts = Ts_field(); To = To_field())

k["diffusion!"] = @benchmarkable diffusion!(T, $Z_AIR, $WS, $TIMESTATE) setup = (T = unit_field())

k["advection!"] = @benchmarkable advection!(T, $Z_AIR, $WS, $TIMESTATE, $CFG) setup = (T = unit_field())

k["convergence!"] = @benchmarkable convergence!(T, $OMEGA, $TIMESTATE, $WS) setup = (T = unit_field())

k["circulation! (heat)"] = @benchmarkable circulation!(T, $Z_AIR, $WS.dX_crcl, $WS, $TIMESTATE, $CFG) setup =
    (T = Ts_field())

k["circulation! (vapor)"] = @benchmarkable circulation!(T, $Z_VAPOR, $WS.dX_crcl, $WS, $TIMESTATE, $CFG) setup =
    (T = q_field())

k["tendencies!"] = @benchmarkable tendencies!($CO2, Ts, Ta, To, q, $WS, $TIMESTATE, $CFG) setup =
    (Ts = Ts_field(); Ta = Ts_field(); To = To_field(); q = q_field())

# ── full-model benchmark (one control year = 730 steps) ──────────────────────
m = SUITE["model"] = BenchmarkGroup()
m["greb_model! (ctrl=1yr)"] = @benchmarkable(
    (redirect_stdout(devnull) do; greb_model!(0, 1, 0, cfg); end);
    setup = (cfg = create_experiment_config(:full_model)),
    evals = 1, samples = 5, seconds = 60,
)

# ── standalone runner: pretty-print the results ──────────────────────────────
function main()
    println("Running GREB benchmark suite (median of samples)...\n")
    results = run(SUITE; verbose = true)
    println("\n", "="^64)
    println(rpad("benchmark", 30), rpad("median", 14), rpad("allocs", 10), "memory")
    println("="^64)
    for gname in ("kernels", "model")
        haskey(results, gname) || continue
        for (name, t) in sort(collect(results[gname]); by = first)
            tr = median(t)
            println(rpad(name, 30),
                    rpad(BenchmarkTools.prettytime(tr.time), 14),
                    rpad(string(tr.allocs), 10),
                    BenchmarkTools.prettymemory(tr.memory))
        end
    end
    println("="^64)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
