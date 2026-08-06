# GREB.jl

A Julia translation of the **Globally Resolved Energy Balance (GREB)**
climate model, originally developed by Dietmar Dommenget and colleagues at
Monash University.

The model simulates the global energy balance on a **3.75° × 3.75°** grid
(96 longitudes × 48 latitudes), stepping every 12 hours with 30-minute
sub-steps for atmospheric circulation (730 time steps per year).

## Features

- Global grid resolution: 96×48 (longitude × latitude)
- 12-hour main time steps with 30-minute sub-steps for circulation
- Support for multiple climate datasets (NCEP, ERA-Interim)
- Future climate scenarios (RCP 2.6, 4.5, 6.0, 8.5)
- Orbital forcing and paleoclimate experiments
- SIMD-vectorized physics kernels (`LoopVectorization.@turbo`)

## Installation

Requires Julia 1.9 or later.

```julia
using Pkg
Pkg.develop(url = "https://github.com/EnvDroneSense/GREB-julia")
```

Or, working from a clone:

```julia
julia --project=.
using Pkg; Pkg.instantiate()
```

## Input data

The model reads **JLD2**-formatted climatology, flux-correction, and solar
forcing files. These are too large to ship with the package; see
[`DATA_README.md`](https://github.com/EnvDroneSense/GREB-julia/blob/main/DATA_README.md)
in the repository for the expected layout and the conversion script that
builds a `greb_dataset_jld2/` directory from the original GREB `.bin` files.

See the [Tutorial](@ref) for a runnable end-to-end example, or the
[API Reference](@ref) for every exported function and type.

## An interactive alternative

The package also ships an interactive [Pluto.jl](https://github.com/fonsp/Pluto.jl)
notebook (`notebooks/GREB_julia.jl`) with widget-driven experiment
configuration — useful for exploration, though this documentation covers the
plain-Julia API.
