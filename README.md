# GREB Climate Model - Julia Implementation

[![Julia](https://img.shields.io/badge/Julia-1.10+-9558B2?logo=julia)](https://julialang.org/)
[![Pluto](https://img.shields.io/badge/Pluto-Interactive-purple)](https://github.com/fonsp/Pluto.jl)
[![docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://EnvDroneSense.github.io/GREBClimate.jl/)

A high-performance Julia translation of the **Globally Resolved Energy Balance (GREB)** climate model, originally developed by Dietmar Dommenget and colleagues at Monash University.

GREBClimate is a **Julia package**: you call it from a script or the REPL, and nothing is held as module-global state. An interactive [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebook ships alongside it (`notebooks/GREB_julia.jl`) for widget-driven exploration and process-isolation decomposition experiments, but it is one front-end onto the package rather than the model itself.

---
> **Repository layout:** GREB is now organized as a standard Julia package,
> `GREBClimate.jl`. The model code lives under `src/` (a `module GREBClimate`,
> originally extracted **verbatim** from the notebook and since split into
> topical files - see `src/GREBClimate.jl` for the include order), tests in
> `test/`, a [Documenter.jl](https://EnvDroneSense.github.io/GREBClimate.jl/)
> site in `docs/`, a plain-Julia driver in `examples/run_greb.jl`, and the
> original interactive Pluto notebook - unchanged - in `notebooks/GREB_julia.jl`.
> See [Project Structure](#project-structure) for the full layout.
>
> ```julia
> julia --project=.                 # activate the package env
> using GREBClimate
> cfg    = create_experiment_config(:full_model)
> fields = load_greb_jld2!("greb_input_data"; dataset=:ncep)
> result = greb_model!(RunSpec(), cfg;    # flux=0, ctrl=1, scnr=1 years
>                      jld2_dir="greb_input_data", fields=fields)
> ```
> `load_greb_jld2!` **returns** the loaded climatology - it does not populate
> global state - so the result must be handed to `greb_model!` via `fields=`.
> Omitting it is refused rather than silently run; see
> [Uninitialized fields](#uninitialized-fields) below.
>
> Run the tests with `julia --project=. -e 'using Pkg; Pkg.test()'`, or the full
> driver with `julia --project=. examples/run_greb.jl <path/to/greb_input_data>`.

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

- **Performance optimizations** using `@turbo` (SIMD vectorization), multi-threaded circulation, and a `Float32` compute path throughout
- **Interactive visualization** through Pluto.jl
- **Multiple climate scenarios** (e.g., IPCC RCP and SSP scenarios)
- **Flexible experiment configurations**

## Features

- 🌍 Global grid resolution: 96×48 (longitude × latitude)
- ⏱️ 12-hour main time steps with 30-minute sub-steps for circulation
- 📊 Real-time visualization of climate variables
- 🔬 Support for multiple climate datasets (NCEP, ERA-Interim)
- 🌡️ IPCC climate scenarios - all four RCPs and five SSPs are implemented (`:rcp26`/`:rcp45`/`:rcp60`/`:rcp85`, `:ssp119`/`:ssp126`/`:ssp245`/`:ssp460`/`:ssp585`), plus a historical CO2 (1850–2017) hindcast (`:historical_co2`) and a user-supplied CO2 trajectory (`:custom_co2`)
- 🧩 "Deconstruct" experiments (`:decon_mean_climate`, `:decon_2xco2`) - toggle individual mean-climate/2×CO₂-response feedback processes on or off via `log_*_dmc`/`log_*_drsp` keywords passed to `create_experiment_config`
- ⚡ Multi-threaded physics - the temperature and humidity `circulation!` calls run concurrently via `Threads.@spawn` when Julia is started with 2+ threads (e.g. `julia --project=. -t 2`, the recommended count - see `.claude/skills/benchmark/SKILL.md` for why more threads no longer reliably helps further)
- 🎯 `Float32` throughout - climatology, workspace buffers, and model state all compute in single precision (matching the `Float32` JLD2 input data natively, no upconversion), for a further ~1.6× wall-clock speedup on top of threading, with output validated against the previous `Float64` path to well under 0.01 K
- ☀️ Orbital forcing and paleoclimate experiments

## 🚀 Quick Start

### Prerequisites

Requires **Julia 1.10** (the current LTS) or later. Download from [julialang.org](https://julialang.org/downloads/).

```bash
git clone https://github.com/EnvDroneSense/GREBClimate.jl
cd GREBClimate.jl
```

### Installation

GREBClimate is not yet registered in the Julia General Registry — once it
is, `Pkg.add("GREBClimate")` will work directly. Until then, install from
the git clone above:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This installs all dependencies from `Project.toml`:

| Package | Purpose |
|:--------|:--------|
| `JLD2` | Reading/writing the model's `.jld2` input data |
| `DataDeps` | Fetching and caching the input dataset on first use |
| `LoopVectorization` | SIMD performance |
| `PrecompileTools` | Precompiles hot kernels at build time (faster first run) |

The Pluto notebook environment (`notebooks/`) separately depends on `PlutoUI`
for its interactive controls. The [Documenter.jl](https://EnvDroneSense.github.io/GREBClimate.jl/)
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

The model reads **JLD2** formatted files ([JuliaIO/JLD2.jl](https://github.com/JuliaIO/JLD2.jl)) - a standard Julia data container. Each field file stores plain Julia values under the keys `"data"` (an `Array{Float32}`), `"dim_names"`, and optionally `"coords"` (physical coordinate values, e.g. an orbital scenario's index) and `"ctl"` (the original GrADS `.ctl` metadata text). `load_greb_jld2!` loads this data directly into `Float32` `ClimateFields` - the same precision the model computes in throughout, so no promotion happens at load time.

### Getting the data

The dataset is ~439 MB unpacked, so it is **not** committed to the repository
(`greb_input_data/` is gitignored). It is fetched on demand via
[DataDeps.jl](https://github.com/oxinabox/DataDeps.jl):

```julia
using GREBClimate
dir    = greb_data_dir()        # prompts, downloads (~353 MB) and caches on first use
fields = load_greb_jld2!(dir; dataset=:ncep)
```

`greb_data_dir()` resolves in this order, and only the last step touches the
network:

| | Source |
|--:|:-------|
| 1 | an explicit path — `greb_data_dir("/path/to/greb_input_data")` |
| 2 | `$GREB_DATA` |
| 3 | `greb_input_data/` beside the repository |
| 4 | the `GREB-input-data` DataDep — downloaded and cached |

So if you already have the dataset, nothing is downloaded. Point at it directly,
or set the environment variable once:

```bash
export GREB_DATA=/path/to/greb_input_data
julia --project=. examples/run_greb.jl
```

The download is verified against a recorded SHA256 and cached under
`~/.julia/scratchspaces/.../datadeps/GREB-input-data`, so it happens once per
machine rather than once per project.

Two environment variables are worth knowing:

| Variable | Effect |
|:---------|:-------|
| `DATADEPS_ALWAYS_ACCEPT=true` | Skip the download prompt. **Required** in CI or any non-interactive session, which otherwise blocks waiting on stdin. |
| `DATADEPS_DISABLE_DOWNLOAD=true` | Make a would-be download throw instead. Useful on metered connections. |

Running the tests or the benchmarks never downloads anything — both resolve with
`allow_download=false` and skip data-dependent work when no local dataset is
found.

### Regenerating the data from raw `.bin` files (maintainers)

You do **not** need this to run the model - it is how the `.jld2` bundle above
is produced in the first place. The raw GREB `.bin` inputs are collated from
several upstream sources and are not redistributed here;
[DATA_README.md](DATA_README.md) documents the files and their layout.

```bash
julia --project=. tools/convert_greb_to_jld2.jl <input_dir> [output_dir]
# output_dir defaults to greb_input_data/
```

To publish a regenerated dataset, build the distributable archive and its
checksum, then attach it to the `data-v1` release:

```bash
julia --project=. tools/package_dataset.jl greb_input_data greb_input_data-v1.tar.gz
```

That script validates the tree against the converter's allowlist first (so a
stray field cannot silently enlarge every user's download), builds a
reproducible `tar.gz`, and prints the SHA256 to paste into `DATA_SHA256` in
[`src/data.jl`](src/data.jl).

### Directory Structure

```
greb_input_data/                    # 39 files, ~439 MB
├── static/
│   ├── global.topography.jld2      # 2D (96×48)
│   └── greb.glaciers.jld2          # 2D (96×48)
├── climatology/                    # 31 files; all 3D (96×48×730) unless noted
│   │   # dataset=:ncep
│   ├── ncep.tsurf.1948-2007.clim.jld2
│   ├── ncep.zonal_wind.850hpa.clim.jld2
│   ├── ncep.meridional_wind.850hpa.clim.jld2
│   ├── ncep.atmospheric_humidity.clim.jld2
│   ├── ncep.soil_moisture.clim.jld2        # also used by dataset=:era
│   │   # dataset=:era (alternative to the ncep.* four above)
│   ├── erainterim.tsurf.1979-2015.clim.jld2
│   ├── erainterim.zonal_wind.850hpa.clim.jld2
│   ├── erainterim.meridional_wind.850hpa.clim.jld2
│   ├── erainterim.atmospheric_humidity.clim.jld2
│   │   # common to both datasets
│   ├── isccp.cloud_cover.clim.jld2
│   ├── woce.ocean_mixed_layer_depth.clim.jld2
│   ├── Tocean.clim.jld2
│   ├── erainterim.omega.vertmean.clim.jld2
│   ├── erainterim.omega_std.vertmean.clim.jld2
│   ├── erainterim.windspeed.850hpa.clim.jld2
│   ├── flux_corrections.jld2       # Tsurf/vapour/Tocean corrections, combined (always loaded together)
│   │   # CMIP5 RCP8.5 anomalies - climate-change experiments only
│   ├── cmip5.{tsurf,zonal.wind,meridional.wind,windspeed,omega}.rcp85.ensmean.forcing.jld2
│   │   # ENSO anomalies - :elnino / :lanina only (10 files)
│   └── erainterim.{tsurf,zonal.wind,meridional.wind,windspeed,omega}.{elnino,lanina}.forcing.jld2
├── solar/
│   └── solar_radiation.clim.jld2   # 2D (48×730)
├── solar_scenarios/                # Optional
│   ├── solar_paleo.jld2
│   ├── solar_eccentricity.jld2    # (ecc_index, lat, time), coords[1] = actual eccentricity values
│   └── solar_obliquity.jld2       # (obl_index, lat, time), coords[1] = actual obliquity angles
└── scenario/                        # Optional
    ├── ipcc_scenarios.jld2        # Dict{String,Dict{Int,Float64}}, keyed "rcp85"/"ssp585"/"hist"/...
    │                              # ("hist" backs :historical_co2, which starts at year 1850 not 1950)
    └── historical_emissions_population.jld2   # year => (co2_emissions_gt_co2_yr, population_billions)
                                    # written by the converter; not read by the model
                                    # or the tests. Kept for analysis use (11 KB).
```

### Loading Data

```julia
fields = load_greb_jld2!(jld2_dir; dataset=:ncep)   # or :era
```

`load_greb_jld2!` returns a [`ClimateFields`](https://EnvDroneSense.github.io/GREBClimate.jl/)
holding the climatology, derived grid geometry, flux corrections and solar
table. It is a value, not global state: hold several independent instances in
one session (e.g. for parameter sweeps), and pass the one you want into
`greb_model!` with `fields=`.

#### Uninitialized fields

A bare `ClimateFields()` is all zeros. Stepping the model on a zero
climatology runs to completion but yields a physically meaningless world
(global-mean Ts ≈ 233 K / −40 °C), so `greb_model!` **refuses** it rather than
returning plausible-looking nonsense:

```julia
julia> greb_model!(RunSpec(), cfg)          # forgot fields=
ERROR: greb_model! was given an uninitialized ClimateFields (all-zero climatology).
```

Data-free runs are legitimate for config- and scenario-plumbing tests (and for
package precompilation, which must not require the dataset). Those opt in
explicitly:

```julia
greb_model!(RunSpec(scnr=0), cfg; jld2_dir="", allow_uninitialized=true)
```

---
## 🎮 Running the Model

### 1. Load Data

```julia
jld2_dir = joinpath(@__DIR__, "greb_input_data")
fields = load_greb_jld2!(jld2_dir; dataset=:ncep)
```

### 2. Configure the Experiment

```julia
cfg = create_experiment_config(:full_model)   # or :co2_double, :elnino, :rcp85, :ssp585, :historical_co2, ...
cfg.log_rain = 1                              # override any switch after construction
```

Advanced experiment types take extra keywords:

```julia
# User-supplied CO2 trajectory (a "year CO2" text file, same format as the
# IPCC scenario files)
cfg = create_experiment_config(:custom_co2; co2_path="my_co2_trajectory.txt")

# Deconstruct experiments: toggle individual feedback processes off
cfg = create_experiment_config(:decon_mean_climate; log_ocean_dmc=false)
cfg = create_experiment_config(:decon_2xco2; log_clouds_drsp=false)
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

Each `MonthlyRecord` is a `NamedTuple` of `(xdim, ydim)` `Matrix{Float32}` fields:
`Ts, Ta, To, q, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens`

See the [Tutorial](https://EnvDroneSense.github.io/GREBClimate.jl/tutorial/) or
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
GREBClimate.jl/
├── src/                        # the GREBClimate package (module GREBClimate)
│   ├── GREBClimate.jl          # module shell + include order
│   ├── constants.jl            # grid/physical constants
│   ├── config.jl               # PhysicsConfig, RunSpec, experiment presets
│   ├── state.jl                # ClimateFields, ModelState, workspaces
│   ├── data.jl                 # greb_data_dir(): dataset location + DataDep
│   ├── io.jl                   # JLD2 loaders
│   ├── physics/                # radiation.jl, hydrology.jl, ocean.jl
│   ├── circulation.jl          # diffusion/advection/convergence
│   ├── tendencies.jl           # per-timestep physics pipeline, forcing()
│   ├── output.jl               # diagnostics!/output!/time_loop!
│   ├── postprocess.jl          # monthly climatology/anomalies
│   └── model.jl                # init_model!/qflux_correction!/greb_model!
├── test/runtests.jl            # unit, integration, and golden-regression tests
├── benchmark/run_benchmarks.jl # timing/allocation suite for the physics kernels
├── docs/                       # Documenter.jl site (index, tutorial, switches, API)
├── examples/run_greb.jl        # plain-Julia driver (no Pluto) — start here
├── notebooks/
│   ├── GREB_julia.jl           # interactive Pluto notebook (own Project.toml)
│   ├── PultoUI.jl              # work-in-progress UI experiments — not wired up
│   └── launch_pluto.jl         # convenience launcher
├── tools/
│   ├── convert_greb_to_jld2.jl  # raw .bin -> JLD2 converter (maintainers only)
│   └── package_dataset.jl       # build the published dataset archive + SHA256
├── DATA_README.md              # raw .bin input inventory (maintainers only)
├── CHANGELOG.md                # user-facing changelog
└── .claude/skills/             # agent task playbooks (benchmarking, docs checks)
```

Not committed, but expected at runtime or used by maintainers (all gitignored):

```
greb_input_data/                # the .jld2 dataset the model reads (~439 MB, auto-downloaded)
Data/                           # raw GREB .bin inputs, only to regenerate the above (~581 MB)
.claude/notes/                  # maintainers' working notes
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

### 🐛 Reporting Issues

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
- **Visualisation dashboard** - embedded interactive maps and time series (similar to the [interactive database](https://mscm.dkrz.de/GREB_model.html?locale=EN) )
- **Physics guide** - the [Documenter.jl site](https://EnvDroneSense.github.io/GREBClimate.jl/) now covers the API and a runnable tutorial; a deeper physics-derivation guide is still open
- **Package registration** - formally register GREBClimate.jl with the Julia General Registry for easy installation

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
