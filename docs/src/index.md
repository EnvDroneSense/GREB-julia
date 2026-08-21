# GREBClimate.jl

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

Requires Julia 1.10 (the current LTS) or later, matching the `julia = "1.10"`
compat bound in `Project.toml`.

```julia
using Pkg
Pkg.develop(url = "https://github.com/EnvDroneSense/GREBClimate.jl")
```

Or, working from a clone:

```julia
julia --project=.
using Pkg; Pkg.instantiate()
```

## Input data

The model reads a **JLD2** dataset - climatology, flux corrections, solar
forcing and scenario tables - laid out under a single directory. At ~580 MB it
is not shipped with the package: request the prepared bundle from the
maintainers, unpack it anywhere, and pass its path in.

```julia
fields = load_greb_jld2!("/path/to/greb_input_data"; dataset = :ncep)
```

`load_greb_jld2!` *returns* the loaded state; it does not populate globals. The
returned value must be handed to [`greb_model!`](@ref) as `fields = ...`, which
refuses to run on an unloaded [`ClimateFields`](@ref) rather than silently
producing a meaningless ≈233 K world. The [Tutorial](@ref) walks through this.

!!! note "Planned: automatic download"
    Distribution via [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl) is
    planned, after which the dataset will be fetched and cached on first use.
    Until then the manual path above is the only one.

Regenerating the dataset from the original GREB `.bin` files is a maintainer
task, not a prerequisite for using the package; see
[`DATA_README.md`](https://github.com/EnvDroneSense/GREBClimate.jl/blob/main/DATA_README.md)
and `scripts/convert_greb_to_jld2.jl` in the repository.

See the [Tutorial](@ref) for a runnable end-to-end example, or the
[API Reference](@ref) for every exported function and type.

## An interactive alternative

The repository also ships an interactive [Pluto.jl](https://github.com/fonsp/Pluto.jl)
notebook (`notebooks/GREB_julia.jl`) with widget-driven experiment
configuration - useful for exploration. It is a front-end onto the same
package; this documentation covers the plain-Julia API it calls.
