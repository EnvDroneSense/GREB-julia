# GREB Climate Model - Julia Implementation

[![Julia](https://img.shields.io/badge/Julia-1.10+-9558B2?logo=julia)](https://julialang.org/)
[![Pluto](https://img.shields.io/badge/Pluto-Interactive-purple)](https://github.com/fonsp/Pluto.jl)
[![docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://EnvDroneSense.github.io/GREB-julia/)

A high-performance Julia translation of the **Globally Resolved Energy Balance (GREB)** climate model, originally developed by Dietmar Dommenget and colleagues at Monash University. This implementation runs in an interactive [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebook with process isolation capabilities for decomposition experiments.

---
> **Repository layout (v0.1):** GREB is now organized as a standard Julia package.
> The model code lives under `src/` (a `module GREB`, originally extracted
> **verbatim** from the notebook and since split into topical files - see
> `src/GREB.jl` for the include order), tests in `test/`, a [Documenter.jl](https://EnvDroneSense.github.io/GREB-julia/)
> site in `docs/`, a plain-Julia driver in `examples/run_greb.jl`, and the
> original interactive Pluto notebook - unchanged - in `notebooks/GREB_julia.jl`.
> See [Project Structure](#project-structure) for the full layout.
>
> ```julia
> julia --project=.                 # activate the package env
> using GREB
> cfg = create_experiment_config(:full_model)
> load_greb_jld2!("greb_dataset_jld2"; dataset=:ncep)
> result = greb_model!(RunSpec(), cfg)   # flux=0, ctrl=1, scnr=1 years
> ```
> Run the tests with `julia --project=. -e 'using Pkg; Pkg.test()'`, or the full
> driver with `julia --project=. examples/run_greb.jl <path/to/greb_dataset_jld2>`.

## 📖 Table of Contents

- [About the Model](#about-the-model)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Input Data](#input-data)
- [Quick Start](#quick-start)
- [Running the Model](#running-the-model)
- [Project Structure](#project-structure)
- [Key Model Components](#key-model-components)
- [References](#references)
- [Contributing](#reporting-issues)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## About the Model

The GREB model is a conceptual climate model that simulates the global energy balance on a **3.75° × 3.75°** grid (96 longitudes × 48 latitudes). It uses a **12-hour main time step** with **30-minute sub-steps** for atmospheric circulation (730 time steps per year).

This implementation has been translated from Fortran90 to Julia with a focus on:

- **Performance optimizations** using `@turbo` (SIMD vectorization)
- **Interactive visualization** through Pluto.jl
- **Multiple climate scenarios** (e.g., IPCC RCP scenarios)
- **Flexible experiment configurations**

## Features

- 🌍 Global grid resolution: 96×48 (longitude × latitude)
- ⏱️ 12-hour main time steps with 30-minute sub-steps for circulation
- 📊 Real-time visualization of climate variables
- 🔬 Support for multiple climate datasets (NCEP, ERA-Interim)
- 🌡️ Future climate scenarios (RCP 2.6, 4.5, 6.0, 8.5)
- ☀️ Orbital forcing and paleoclimate experiments

## 🚀 Quick Start

### Prerequisites

Requires **Julia 1.10** (the current LTS) or later. Download from [julialang.org](https://julialang.org/downloads/).

```bash
git clone https://github.com/EnvDroneSense/GREB-julia
cd GREB-julia
```

### Installation

Open Julia and run:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This installs all dependencies from `Project.toml`:

| Package | Purpose |
|:--------|:--------|
| `JLD2` | Reading/writing the model's `.jld2` input data |
| `LoopVectorization` | SIMD performance |
| `PrecompileTools` | Precompiles hot kernels at build time (faster first run) |

The Pluto notebook environment (`notebooks/`) separately depends on `PlutoUI`
for its interactive controls. The [Documenter.jl](https://EnvDroneSense.github.io/GREB-julia/)
site under `docs/` has its own environment too.

### Launch Pluto (optional)

The interactive notebook is one way to run the model - see
[Running the Model](#running-the-model) below for the plain-Julia path.

```julia
using Pluto
Pluto.run()
```

Open `GREB_julia.jl` from the Pluto interface.

## 📂 Input Data

The model reads **JLD2** formatted files ([JuliaIO/JLD2.jl](https://github.com/JuliaIO/JLD2.jl)) - a standard Julia data container. Each field file stores plain Julia values under the keys `"data"` (an `Array{Float32}`), `"dim_names"`, and optionally `"coords"` (physical coordinate values, e.g. an orbital scenario's index) and `"ctl"` (the original GrADS `.ctl` metadata text).

In the original model these were all separate `.bin` files. `scripts/convert_greb_to_jld2.jl` converts the raw GREB `.bin` input files (see [DATA_README.md](DATA_README.md) for their expected layout, normally under `Data/input/`) into this `.jld2` layout:

```bash
julia --project=. scripts/convert_greb_to_jld2.jl [input_dir] [output_dir]
# defaults: input_dir=Data/input, output_dir=greb_dataset_jld2
```

These data files are too large to upload to GitHub but can be made available on request, or regenerated from the raw `.bin` files with the converter script.

### Directory Structure

```
greb_dataset_jld2/
├── static/
│   ├── global.topography.jld2      # 2D (96×48)
│   └── greb.glaciers.jld2          # 2D (96×48)
├── climatology/
│   ├── ncep.tsurf.1948-2007.clim.jld2       # 3D (96×48×730)
│   ├── ncep.zonal_wind.850hpa.clim.jld2
│   ├── ncep.meridional_wind.850hpa.clim.jld2
│   ├── ncep.atmospheric_humidity.clim.jld2
│   ├── ncep.soil_moisture.clim.jld2
│   ├── isccp.cloud_cover.clim.jld2
│   ├── woce.ocean_mixed_layer_depth.clim.jld2
│   ├── Tocean.clim.jld2
│   ├── erainterim.omega.vertmean.clim.jld2
│   ├── erainterim.omega_std.vertmean.clim.jld2
│   ├── erainterim.windspeed.850hpa.clim.jld2
│   └── [flux_correction files]
├── solar/
│   └── solar_radiation.clim.jld2   # 2D (48×730)
└── solar_scenarios/                # Optional
    ├── solar_paleo.jld2
    ├── solar_eccentricity.jld2    # (ecc_index, lat, time), coords[1] = actual eccentricity values
    └── solar_obliquity.jld2       # (obl_index, lat, time), coords[1] = actual obliquity angles
```

### Loading Data

In the notebook, set the `jld2_dir` variable and run:

```julia
load_greb_jld2!(jld2_dir; dataset=:ncep)   # or :era
```

---
## 🎮 Running the Model

### 1. Load Data

```julia
jld2_dir = joinpath(@__DIR__, "greb_dataset_jld2")
fields = load_greb_jld2!(jld2_dir; dataset=:ncep)
```

### 2. Configure the Experiment

```julia
cfg = create_experiment_config(:full_model)   # or :co2_double, :elnino, :rcp85, ...
cfg.log_rain = 1                              # override any switch after construction
```

### 3. Run the Model

```julia
run = RunSpec(flux=0, ctrl=1, scnr=1)   # flux-correction, control, scenario years
result = greb_model!(run, cfg; jld2_dir=jld2_dir, fields=fields)
```

This runs, in order: an optional flux-correction spin-up (nudges toward
climatology), a control run at fixed CO₂, and a scenario run under
time-varying forcing (e.g. a CO₂ ramp).

### 4. Access Results

```julia
result.ctrl    # Vector{MonthlyRecord} (control)
result.scnr    # Vector{MonthlyRecord} (scenario)
```

Each `MonthlyRecord` is a `NamedTuple` with fields:
`Ts, Ta, To, q, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens`

See the [Tutorial](https://EnvDroneSense.github.io/GREB-julia/tutorial/) or
[`examples/run_greb.jl`](examples/run_greb.jl) for the full runnable version
of the above, including a global-mean summary and an optional plot.

### Or, interactively

The Pluto notebook (`notebooks/GREB_julia.jl`) exposes the same options as
widgets instead of code:

| Control | Description |
|:--------|:------------|
| **Experiment** | Preset experiments (2×CO₂, El Niño, RCP8.5, etc.) |
| **Configuration Preset** | Full physics, no feedbacks, MSCM, custom |
| **Mean Climate Switches** | Toggle clouds, vapor, ice, circulation, etc. |
| **CO₂ Response Switches** | Process-specific response toggles |
| **Circulation Components** | Diffusion, advection, convergence |
| **Hydrology Parameters** | Rain/EVA modes, climatology dataset |
| **Run Duration** | Flux correction, control, and scenario years |

Toggle the **Execute Model** checkbox to run; results land in `last_run`
(same `.ctrl`/`.scnr` shape as above).

## 📁 Project Structure

```
GREB-julia/
├── src/                        # the GREB package (module GREB)
│   ├── GREB.jl                 # module shell + include order
│   ├── constants.jl            # grid/physical constants
│   ├── config.jl               # PhysicsConfig, RunSpec, experiment presets
│   ├── state.jl                # ClimateFields, ModelState, workspaces
│   ├── io.jl                   # JLD2 loaders
│   ├── physics/                # radiation.jl, hydrology.jl, ocean.jl
│   ├── circulation.jl          # diffusion/advection/convergence
│   ├── tendencies.jl           # per-timestep physics pipeline, forcing()
│   ├── output.jl               # diagnostics!/output!/time_loop!
│   ├── postprocess.jl          # monthly climatology/anomalies
│   └── model.jl                # init_model!/qflux_correction!/greb_model!
├── test/runtests.jl            # unit, integration, and golden-regression tests
├── docs/                       # Documenter.jl site (API reference + tutorial)
├── examples/run_greb.jl        # plain-Julia driver (no Pluto)
├── notebooks/GREB_julia.jl     # original interactive Pluto notebook (unchanged)
├── scripts/convert_greb_to_jld2.jl  # raw .bin -> JLD2 converter
└── claude/                     # dev notes: IMPROVEMENTS.md
```

## 🔬 Key Model Components

### Energy Balance
- Shortwave radiation absorption
- Longwave radiation emission
- Surface energy fluxes

### Hydrological Cycle (MSCM)
- Precipitation calculation
- Evaporation processes
- Soil moisture dynamics

### Ocean Model
- Mixed-layer temperature evolution
- Deep ocean heat exchange
- Sea ice formation and melting

### Atmosphere
- Atmospheric heat transport
- Moisture transport
- Simplified circulation patterns

---
## ⚠️ Known Issues


### Reporting Issues

If you encounter these or other problems:
1. Check that all input data files are correctly formatted and located
2. Verify Julia and package versions match requirements
3. Try restarting the Pluto notebook
4. Open an issue on GitHub with:
   - Julia version (`versioninfo()`)
   - Error messages or unexpected behavior description
   - Steps to reproduce

Contributions to fix these issues are welcome! Open a pull request or an issue on GitHub.

---
## 🔭 Future Plans

- **NetCDF output** - optional direct‑write of monthly means 
- **Parallelisation** - multi‑threading for longer runs  
- **Visualisation dashboard** - embedded interactive maps and time series (similar to the [interactive database](https://mscm.dkrz.de/GREB_model.html?locale=EN) )
- **Physics guide** - the [Documenter.jl site](https://EnvDroneSense.github.io/GREB-julia/) now covers the API and a runnable tutorial; a deeper physics-derivation guide is still open
- **Package registration** - formally register GREB.jl with the Julia General Registry for easy installation

---
## 📚 References

### Primary Publications

1. **Dommenget, D., and Flöter, J. (2011)**. Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model. *Journal of Climate Dynamics*, 37: 2143. [doi:10.1007/s00382-011-1026-0](https://doi.org/10.1007/s00382-011-1026-0)

2. **Stassen, C., Dommenget, D., and Loveday, N. (2019)**. A hydrological cycle model for the Globally Resolved Energy Balance (GREB) model v1.0. *Geoscientific Model Development*, 12, 425-440. [doi:10.5194/gmd-12-425-2019](https://doi.org/10.5194/gmd-12-425-2019)

3. **Dommenget, D., Nice, K., Bayr, T., Kasang, D., Stassen, C., and Rezny, M.** The Monash Simple Climate Model Experiments: An interactive database of the mean climate, climate change and scenarios simulations. *Geoscientific Model Development*, 12, 2155-2179. [doi:10.5194/gmd-12-2155-2019](https://doi.org/10.5194/gmd-12-2155-2019)

### Original GREB Model
- [Monash University GREB Homepage](http://www.monash.edu/science/research/climate)

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Original GREB Model**: Dietmar Dommenget, Janine Flöter, Tobias Bayr, Christian Stassen (Monash University)
- **Julia Translation and Optimization**: Thomas Struys (UGent)
- **Pluto.jl**: For the interactive notebook environment
- **Julia Community**: For excellent scientific computing tools

---

**Contact**: For questions or support, please open an issue on GitHub.
