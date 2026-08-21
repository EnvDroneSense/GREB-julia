# Tutorial

This walks through the same flow as [`examples/run_greb.jl`](https://github.com/EnvDroneSense/GREBClimate.jl/blob/main/examples/run_greb.jl):
load input data, configure an experiment, run the model, and inspect the
result.

## 1. Load input data

```julia
using GREBClimate

jld2_dir = "greb_input_data"   # see DATA_README.md for the expected layout
fields = load_greb_jld2!(jld2_dir; dataset = :ncep)   # or :era
```

`fields` is a [`ClimateFields`](@ref) — climatology, grid geometry, flux
corrections, and the regional-CO₂ mask/solar table. Every physics function
takes it as an explicit argument; nothing is shared as module-global state,
so you can hold several independent `fields` instances (e.g. for parameter
sweeps) in the same session.

## 2. Configure the experiment

[`create_experiment_config`](@ref) returns a [`PhysicsConfig`](@ref) preset
for a named experiment:

```julia
cfg = create_experiment_config(:full_model)   # or :co2_double, :elnino, :rcp85, ...
```

`cfg` is a mutable struct — override individual switches after construction,
e.g. `cfg.log_rain = 1` to pick a different hydrology parameterization.
See the [Physics Switches](@ref) page for the full list of switches and
what each one controls.

## 3. Run the model

[`greb_model!`](@ref) takes a [`RunSpec`](@ref) (how many years of
flux-correction spin-up, control, and scenario to run) and the config:

```julia
run = RunSpec(flux = 0, ctrl = 1, scnr = 1)
result = greb_model!(run, cfg; jld2_dir = jld2_dir, fields = fields)
```

This runs, in order: an optional flux-correction spin-up (nudges toward
climatology), a control run at fixed CO₂, and a scenario run under
time-varying forcing (e.g. a CO₂ ramp).

## 4. Inspect results

```julia
result.ctrl    # Vector{MonthlyRecord}, one per control-run month
result.scnr    # Vector{MonthlyRecord}, one per scenario-run month
```

Each [`MonthlyRecord`](@ref) is a `NamedTuple` with fields
`Ts, Ta, To, q, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens` — each
a `(96, 48)` matrix of that month's mean.

```julia
using Statistics
Ts_global_mean = [mean(rec.Ts) for rec in result.ctrl]
```

## Next steps

- The [API Reference](@ref) lists every exported function and type.
- `benchmark/run_benchmarks.jl` micro-benchmarks the per-timestep physics
  kernels; `claude/BENCHMARKS.md` in the repository tracks results across
  optimization passes.
- `test/runtests.jl` doubles as executable documentation for individual
  kernels' behavior under different config switches.
