# GREB Climate Model - Julia Implementation

[![Julia](https://img.shields.io/badge/Julia-1.9+-9558B2?logo=julia)](https://julialang.org/)
[![Pluto](https://img.shields.io/badge/Pluto-Interactive-purple)](https://github.com/fonsp/Pluto.jl)

A high-performance Julia translation of the **Globally Resolved Energy Balance (GREB)** climate model, originally developed by Dietmar Dommenget and colleagues at Monash University. This implementation runs in an interactive [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebook with process isolation capabilities for decomposition experiments.

---
> **Repository layout (v0.1):** GREB is now organized as a standard Julia package.
> The model code lives in `src/GREB.jl` (a `module GREB`, extracted **verbatim**
> from the notebook), tests in `test/`, a plain-Julia driver in
> `examples/run_greb.jl`, and the original interactive Pluto notebook — unchanged —
> in `notebooks/GREB_julia.jl`.
>
> ```julia
> julia --project=.                 # activate the package env
> using GREB
> cfg = create_experiment_config(:full_model)
> load_greb_jdal2!("greb_dataset_jdal2"; dataset=:ncep)
> result = greb_model!(0, 1, 1, cfg)   # (flux, ctrl, scenario) years
> ```
> Run the tests with `julia --project=. -e 'using Pkg; Pkg.test()'`, or the full
> driver with `julia --project=. examples/run_greb.jl <path/to/greb_dataset_jdal2>`.

## 📖 Table of Contents

- [About the Model](#About-the-Model)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Input Data](#input-data)
- [Quick Start](#quick-start)
- [Running the Model](#running-the-model)
- [Project Structure](#project-structure)
- [Key Model Components](#key-model-components)
- [References](#references)
- [Contributing](#contributing)
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

### 1. Clone the Repository

```bash
git clone https://github.com/EnvDroneSense/GREB-julia
cd GREB_julia
```

### 2. Install Julia

Requires **Julia 1.9** or later. Download from [julialang.org](https://julialang.org/downloads/).

### 3. Activate the Environment

Open Julia and run:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This installs all dependencies from `Project.toml`:

| Package | Purpose |
|:--------|:--------|
| `PlutoUI` | Interactive controls |
| `NCDatasets` | NetCDF I/O (optional) |
| `LoopVectorization` | SIMD performance |
| `StaticArrays` | Optimized array operations |
| `BenchmarkTools`, `Profile` | Performance analysis |
| `Statistics` | Statistical functions |

### 4. Launch Pluto

```julia
using Pluto
Pluto.run()
```

Open `GREB_julia.jl` from the Pluto interface.

## 📂 Input Data

The model reads **JDAL2** formatted files. JDAL2 is a self-describing binary format with embedded dimensions.

In the original model these were al seperate BIN files, but to improve loading efficiency and folder clarity these have been converted to JDAL2. these data files were to large to upload to github but can be made available on request.

### Directory Structure

```
greb_dataset_jdal2/
├── static/
│   ├── global.topography.jd2      # 2D (96×48)
│   └── greb.glaciers.jd2          # 2D (96×48)
├── climatology/
│   ├── ncep.tsurf.1948-2007.clim.jd2       # 3D (96×48×730)
│   ├── ncep.zonal_wind.850hpa.clim.jd2
│   ├── ncep.meridional_wind.850hpa.clim.jd2
│   ├── ncep.atmospheric_humidity.clim.jd2
│   ├── ncep.soil_moisture.clim.jd2
│   ├── isccp.cloud_cover.clim.jd2
│   ├── woce.ocean_mixed_layer_depth.clim.jd2
│   ├── Tocean.clim.jd2
│   ├── erainterim.omega.vertmean.clim.jd2
│   ├── erainterim.omega_std.vertmean.clim.jd2
│   ├── erainterim.windspeed.850hpa.clim.jd2
│   └── [flux_correction files]
├── solar/
│   └── solar_radiation.clim.jd2   # 2D (48×730)
└── solar_scenarios/                # Optional
    ├── solar_paleo.jd2
    ├── solar_eccentricity.jd2
    └── solar_obliquity.jd2
```

### Loading Data

In the notebook, set the `jdal2_dir` variable and run:

```julia
load_greb_jdal2!(jdal2_dir; dataset=:ncep)   # or :era
```

---
## 🎮 Quick Start

### 1. Load Data

```julia
jdal2_dir = joinpath(@__DIR__, "greb_dataset_jdal2")
load_greb_jdal2!(jdal2_dir; dataset=:ncep)
```

### 2. Configure the Experiment

Use the interactive widgets in the notebook:

| Control | Description |
|:--------|:------------|
| **Experiment** | Preset experiments (2×CO₂, El Niño, RCP8.5, etc.) |
| **Configuration Preset** | Full physics, no feedbacks, MSCM, custom |
| **Mean Climate Switches** | Toggle clouds, vapor, ice, circulation, etc. |
| **CO₂ Response Switches** | Process-specific response toggles |
| **Circulation Components** | Diffusion, advection, convergence |
| **Hydrology Parameters** | Rain/EVA modes, climatology dataset |
| **Run Duration** | Flux correction, control, and scenario years |

### 3. Run the Model

Toggle the **Execute Model** checkbox. The model runs three phases:

1. **Flux Correction** (optional) - computes correction fields to nudge toward climatology
2. **Control Run** - steady-state at fixed CO₂
3. **Scenario Run** - time-varying forcing (e.g., CO₂ ramp, solar changes)

### 4. Access Results

Results are stored in `last_run`:

```julia
ctrl = last_run.ctrl    # Vector of MonthlyRecord (control)
scnr = last_run.scnr    # Vector of MonthlyRecord (scenario)
```

Each `MonthlyRecord` is a `NamedTuple` with fields:  
`Ts, Ta, To, q, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens`

## 🎛️ Interactive Controls

| Section              | Controls                                                                   |
| :------------------- | :------------------------------------------------------------------------- |
| **Experiment**       | Dropdown: full_model, co2_double, elnino, rcp85, etc.                      |
| **Physics Preset**   | Full / No Feedbacks / MSCM / Sensitivity / Custom                          |
| **Mean Climate**     | Clouds, Vapor, Ice, Circulation, Hydrology, Atmosphere, CO₂, Ocean, Q-Flux |
| **CO₂ Response**     | Clouds, Vapor, Circulation, Hydrology, Topography, Humidity                |
| **Circulation**      | Ice albedo, Horizontal/Vertical diffusion & advection, Convergence         |
| **Hydrology**        | Rain mode (-1..3), Evaporation mode (-1..2), Climatology (ERA/NCEP)        |
| **External Forcing** | Surface temperature, Horizontal wind, Vertical velocity                    |
| **Run Duration**     | Flux correction, Control, Scenario years (0-100 each)                      |
| **Execute**          | Run checkbox                                                               |
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

The following issues are currently being worked on:

### Qflux correction
I'm in the process of debugging the qflux_correction module as I am not sure that it is fully working.

### Reporting Issues

If you encounter these or other problems:
1. Check that all input data files are correctly formatted and located
2. Verify Julia and package versions match requirements
3. Try restarting the Pluto notebook
4. Open an issue on GitHub with:
   - Julia version (`versioninfo()`)
   - Error messages or unexpected behavior description
   - Steps to reproduce

Contributions to fix these issues are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md).

---
## 🔭 Future Plans

- **NetCDF output** - optional direct‑write of monthly means 
- **Parallelisation** - multi‑threading for longer runs  
- **Visualisation dashboard** - embedded interactive maps and time series (similar to the [interactive database](https://mscm.dkrz.de/GREB_model.html?locale=EN) )
- **Improved documentation** - detailed physics guide and tutorial notebooks  

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
