# GREB Climate Model - Julia Implementation

[![Julia](https://img.shields.io/badge/Julia-1.7+-9558B2?logo=julia)](https://julialang.org/)
[![Pluto](https://img.shields.io/badge/Pluto-Interactive-purple)](https://github.com/fonsp/Pluto.jl)

A Julia translation and optimization of the **Globally Resolved Energy Balance (GREB)** climate model, originally developed by Dietmar Dommenget and colleagues at Monash University. This implementation provides an interactive Pluto.jl notebook interface for running climate simulations and exploring climate dynamics.

## 📖 Table of Contents

- [About the Model](#about-the-model)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Data](#input-data)
- [Running the Model](#running-the-model)
- [Project Structure](#project-structure)
- [Key Model Components](#key-model-components)
- [References](#references)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## 📚 About the Model

The GREB model is a conceptual climate model that simulates the global energy balance and climate dynamics on a simplified grid. It includes:

- **Energy balance** computations for surface temperature
- **Hydrological cycle** with precipitation and evaporation
- **Ocean mixed-layer dynamics**
- **Sea ice** representation
- **Atmospheric circulation** patterns
- **Cloud and albedo** feedbacks

This Julia implementation features:
- **Performance optimizations** using `@turbo` (SIMD vectorization)
- **Interactive visualization** through Pluto.jl
- **Multiple climate scenarios** (e.g., IPCC RCP scenarios)
- **Flexible experiment configurations**

## ✨ Features

- 🌍 Global grid resolution: 96×48 (longitude × latitude)
- ⏱️ 12-hour main time steps with 30-minute sub-steps for circulation
- 📊 Real-time visualization of climate variables
- 🔬 Support for multiple climate datasets (NCEP, ERA-Interim)
- 🌡️ Future climate scenarios (RCP 2.6, 4.5, 6.0, 8.5)
- ☀️ Orbital forcing and paleoclimate experiments

## 🔧 Prerequisites

### Required Software

- **Julia 1.7 or later** - [Download Julia](https://julialang.org/downloads/)
- **Git** (for cloning the repository)

### System Requirements

- **RAM**: Minimum 4 GB (8 GB recommended)
- **Storage**: ~2-5 GB for input data files
- **OS**: Windows, macOS, or Linux

## 📥 Installation

### 1. Clone the Repository

```bash
git clone https://github.com/EnvDroneSense/GREB-julia
cd greb-julia
```

### 2. Set Up Julia Environment

Open Julia and activate the project environment:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This will install all required packages listed in `Project.toml`:

- **Plots** - Visualization
- **PlutoUI** - Interactive interface components
- **NCDatasets** - NetCDF file I/O
- **LoopVectorization** - Performance optimization
- **StaticArrays** - Efficient array operations
- **Statistics** - Statistical functions
- And other dependencies...

### 3. Download Input Data

The GREB model requires climate input data files in GrADS binary format. These files should be placed in a directory structure as follows:

```
ClimaModel/
├── Data/
│   └── input/
│       ├── global.topography.bin
│       ├── greb.glaciers.bin
│       ├── ncep.tsurf.1948-2007.clim.bin
│       ├── ncep.zonal_wind.850hpa.clim.bin
│       ├── ncep.atmospheric_humidity.clim.bin
│       ├── isccp.cloud_cover.clim.bin
│       ├── ncep.soil_moisture.clim.bin
│       ├── woce.ocean_mixed_layer_depth.clim.bin
│       ├── Tocean.clim.bin
│       ├── Tsurf_flux_correction.bin
│       ├── vapour_flux_correction.bin
│       ├── Tocean_flux_correction.bin
│       ├── solar_radiation.clim.bin
│       ├── cmip5.omega.rcp85.ensmean.forcing.new.bin
│       ├── cmip5.omegastd.rcp85.ensmean.forcing.new.bin
│       ├── ipcc.scenario.rcp26.forcing.txt
│       ├── ipcc.scenario.rcp45.forcing.txt
│       ├── ipcc.scenario.rcp6.forcing.txt
│       ├── ipcc.scenario.rcp85.forcing.txt
│       └── solar_forcing_scenarios/
│           └── [various solar forcing files]
```

**Note**: Input data files are available from the original GREB model repository or can be obtained by contacting the model developers. See the [Input Data](#input-data) section for more details.

## 🚀 Quick Start

### Launch the Interactive Notebook

1. **Start Pluto**:

```julia
using Pluto
Pluto.run()
```

2. **Open the notebook**:
   - In the Pluto interface that opens in your browser
   - Navigate to and open `GREB_julia.jl`

3. **Configure the model**:
   - Set the input data directory path in the notebook
   - Choose a climate scenario (e.g., RCP 8.5)
   - Adjust simulation parameters if needed

4. **Run the simulation**:
   - Click "Run Model" button
   - Watch the simulation progress in real-time
   - Explore interactive visualizations

### Basic Command Line Usage

Alternatively, you can run the model from a Julia script:

```julia
include("GREB_julia.jl")

# Load input data
input_dir = "Data/input"
load_greb_input_data!(input_dir; dataset=:ncep)

# Configure and run simulation
# [See notebook for detailed configuration]
```

## 📂 Input Data

### Required Input Files

The model requires the following binary data files:

#### Static 2D Fields
- `global.topography.bin` - Global topography (96×48)
- `greb.glaciers.bin` - Glacier mask (96×48)

#### 3D Climatology Fields (96×48×730)
- Surface temperature climatology
- Atmospheric winds (zonal and meridional at 850 hPa)
- Atmospheric humidity
- Cloud cover
- Soil moisture
- Ocean mixed-layer depth
- Deep ocean temperature

#### Flux Corrections
- `Tsurf_flux_correction.bin`
- `vapour_flux_correction.bin`
- `Tocean_flux_correction.bin`

#### Solar and Forcing
- `solar_radiation.clim.bin`
- IPCC scenario forcing files (`.txt`)
- Optional: Solar forcing scenarios for paleoclimate experiments

### Climate Datasets

The model supports three climate dataset configurations:

1. **`:ncep`** (default) - NCEP/NCAR Reanalysis (1948-2007)
2. **`:era`** - ERA-Interim (1979-2015)
3. **`:era_ncep`** - Mixed dataset

Specify the dataset when loading data:

```julia
load_greb_input_data!(input_dir; dataset=:ncep)
```

## 🎮 Running the Model

### Simulation Parameters

Key parameters you can adjust:

- **Simulation length**: Number of years to simulate
- **Climate scenario**: Control CO₂ levels or use IPCC scenarios
- **Output frequency**: How often to save results
- **Physics options**: Enable/disable specific processes

### Output

The model generates:

- **NetCDF files** with climate variables (optional)
- **CSV files** for specific variables
- **Interactive plots** within the Pluto notebook

Output variables include:
- Surface temperature (`Ts`)
- Ocean temperature (`To`)
- Atmospheric temperature (`Ta`)
- Specific humidity (`q`)
- Precipitation
- Evaporation
- Sea ice fraction
- Surface albedo

## 📁 Project Structure

```
ClimaModel/
├── GREB_julia.jl              # Main Pluto notebook
├── GREB_julia.jl.bak          # Backup copy
├── Hydrologisch_model.jl      # Separate hydrological model
├── Project.toml               # Julia project dependencies
├── Manifest.toml              # Locked dependency versions
├── README.md                  # This file
├── LICENSE                    # License information
├── Data/                      # Input data directory
│   └── input/                 # Climate input files (.bin)
├── greb_output_*.csv          # Model output files
└── greb_output.nc             # NetCDF output (optional)
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

## 📚 References

### Primary Publications

1. **Dommenget, D., and Flöter, J. (2011)**. Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model. *Journal of Climate Dynamics*, 37: 2143. [doi:10.1007/s00382-011-1026-0](https://doi.org/10.1007/s00382-011-1026-0)

2. **Stassen, C., Dommenget, D., and Loveday, N. (2019)**. A hydrological cycle model for the Globally Resolved Energy Balance (GREB) model v1.0. *Geoscientific Model Development*, 12, 425-440. [doi:10.5194/gmd-12-425-2019](https://doi.org/10.5194/gmd-12-425-2019)

3. **Dommenget, D., Nice, K., Bayr, T., Kasang, D., Stassen, C., and Rezny, M.** The Monash Simple Climate Model Experiments: An interactive database of the mean climate, climate change and scenarios simulations. *Geoscientific Model Development* (submitted)

### Original GREB Model
- [Monash University GREB Homepage](http://www.monash.edu/science/research/climate)

## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

### Development Workflow

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Guidelines

- Follow Julia style conventions
- Add documentation for new functions
- Include tests for new features
- Update README as needed

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Original GREB Model**: Dietmar Dommenget, Janine Flöter, Tobias Bayr, Christian Stassen (Monash University)
- **Julia Translation and Optimization**: Thomas Struys (UGent)
- **Pluto.jl**: For the interactive notebook environment
- **Julia Community**: For excellent scientific computing tools

---

**Contact**: For questions or support, please open an issue on GitHub.

**Last Updated**: March 2026

