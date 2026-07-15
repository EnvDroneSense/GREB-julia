### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 1165eeb4-10d8-4fb5-859d-7fe410189608
begin
	using Statistics
	using PlutoUI
	using NCDatasets
	using StaticArrays  # For optimization: static longitude indices
	using LoopVectorization  # For SIMD vectorization with @turbo
end

# ╔═╡ 3f63ce64-effd-49c7-9906-92cda40d59f0
md"""
# The Globally Resolved Energy Balance (GREB) Model 

Authors: Dietmar Dommenget, Janine Flöter, Tobias Bayr and Christian Stassen

Translated + optimized for Julia: Thomas Struys

**References:**

- Dommenget, D., and Flöter, J. (2011). Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model. *Journal of Climate Dynamics*, 37: 2143. doi:[10.1007/s00382-011-1026-0](https://doi.org/10.1007/s00382-011-1026-0)

- Stassen, C., Dommenget, D., and Loveday, N. (2019). A hydrological cycle model for the Globally Resolved Energy Balance (GREB) model v1.0. *Geoscientific Model Development*, 12, 425-440. doi:[10.5194/gmd-12-425-2019](https://doi.org/10.5194/gmd-12-425-2019)

- Dommenget, D., Nice, K., Bayr, T., Kasang, D., Stassen, C., and Rezny, M. The Monash Simple Climate Model Experiments: An interactive database of the mean climate, climate change and scenarios simulations. *Geoscientific Model Development* (submitted)


input fields: The GREB model needs the following fields to be specified before

| Variable | Dimensions | Description |
|:---------|:-----------|:------------|
| `z_topo(xdim,ydim)` | 2D | topography (<0 are ocean points) [m] |
| `glacier(xdim,ydim)` | 2D | glacier mask (>0.5 are glacier points) |
| `Tclim(xdim,ydim,nstep_yr)` | 3D | mean Tsurf [K] |
| `uclim(xdim,ydim,nstep_yr)` | 3D | mean zonal wind speed [m/s] |
| `vclim(xdim,ydim,nstep_yr)` | 3D | mean meridional wind speed [m/s] |
| `qclim(xdim,ydim,nstep_yr)` | 3D | mean atmospheric humidity [kg/kg] |
| `mldclim(xdim,ydim,nstep_yr)` | 3D | mean ocean mixed layer depth [m] |
| `Toclim(xdim,ydim,nstep_yr)` | 3D | mean deep ocean temperature [K] |
| `swetclim(xdim,ydim,nstep_yr)` | 3D | soil wetness, fraction of total [0-1] |
| `sw_solar(ydim,nstep_yr)` | 2D | 24hrs mean solar radiation [W/m²] |
"""

# ╔═╡ 0780e492-cdcd-43fb-af85-f9a4bab37450
PlutoUI.TableOfContents(title="📋 Table of Contents", depth=2)

# ╔═╡ d0ebd34f-a83a-4dc1-9c66-24a3878b681d
md"""
## 💾 Input/Output

### Data Loading Module

Functions to read JDAL2 formatted input files (self-describing binary format with embedded dimensions).
"""

# ╔═╡ b303e4e9-49fa-45ad-967e-20f165fdf38c
"""
    read_jdal2(filepath::String)

Read a JDAL2 Version 2 file into an Array{Float32}.

# Returns
- `(data, dim_names)` tuple where:
  - `data`: Array{Float32} with shape determined from file header
  - `dim_names`: Vector{String} of dimension names (e.g., ["lon", "lat", "time"])
"""
function read_jdal2(filepath::String)
    data, dims, dim_names = open(filepath, "r") do io
        # Magic bytes
        magic = UInt8[read(io, UInt8) for _ in 1:6]
        if magic[1:5] != UInt8[0x4A, 0x44, 0x41, 0x4C, 0x32]
            error("Not a valid JDAL2 file: $filepath")
        end
        version = magic[6]
        version == 0x02 || error("Only Version 2 supported, got $version")
        
        # Dimension names
        n_dim_names = read(io, Int32)
        dim_names = [String(read(io, read(io, Int32))) for _ in 1:n_dim_names]
        
        # Shape
        ndims = read(io, Int32)
        dims = [read(io, Int32) for _ in 1:ndims]
        
        # Data type
        read(io, UInt8) == 0x01 || error("Only Float32 supported")
        
        # Data
        data = Vector{Float32}(undef, prod(dims))
        read!(io, data)
        (data, dims, dim_names)
    end
    return (data = reshape(data, Tuple(dims)), dim_names = dim_names)
end

# ╔═╡ 8578d6aa-2782-4279-8f6b-78194b8ecc10
function load_solar_forcing_jdal2(jdal2_dir::String, forcing_type::Symbol, index::Int=0)
    
    if forcing_type == :paleo
        filepath = joinpath(jdal2_dir, "solar_scenarios", "solar_paleo.jd2")
        result = read_jdal2(filepath)
        return result.data
        
    elseif forcing_type == :eccentricity
        filepath = joinpath(jdal2_dir, "solar_scenarios", "solar_eccentricity.jd2")
        result = read_jdal2(filepath)
        # result.data is (61, 48, 730) → index 0-based maps to position index+1
        @assert 0 <= index <= 60 "Eccentricity index must be 0-60"
        return result.data[index+1, :, :]
        
    elseif forcing_type == :obliquity
        filepath = joinpath(jdal2_dir, "solar_scenarios", "solar_obliquity.jd2")
        result = read_jdal2(filepath)
        indices = 0:25:230
        pos = findfirst(==(index), indices)
        @assert pos !== nothing "Obliquity index $index not found in $(indices)"
        return result.data[pos, :, :]
        
    else
        error("Unknown forcing type: $forcing_type. Use :paleo, :eccentricity, or :obliquity")
    end
end

# ╔═╡ ebca5877-1305-40f1-ac52-66356dc17661
md"""
### Solar Forcing Storage

Global variable for orbital/paleoclimate solar forcing. Loaded on-demand using grouped JDAL2 files.

**Usage:**
```julia
# Load eccentricity forcing (index 0-60)
global sw_solar_forcing_data = load_solar_forcing_jdal2(jdal2_dir, :eccentricity, 36)

# Load obliquity forcing (index from 0:25:230)
global sw_solar_forcing_data = load_solar_forcing_jdal2(jdal2_dir, :obliquity, 230)

# Load paleoclimate forcing
global sw_solar_forcing_data = load_solar_forcing_jdal2(jdal2_dir, :paleo)
```
"""

# ╔═╡ 8165ed3f-ede0-45c9-83ac-e4e49262457c
md"""
## 🕹️ Advanced Physics Control (MSCM)

**Process Isolation Capabilities**

The MSCM physics switches enable decomposition experiments by selectively enabling/disabling physical processes:

### Mean Climate Switches (`*_dmc`)
- `log_clouds_dmc`: Cloud feedback in mean climate
- `log_vapor_dmc`: Water vapor feedback in mean climate  
- `log_ice_dmc`: Ice-albedo feedback in mean climate
- `log_crcl_dmc`: Atmospheric circulation in mean climate
- `log_hydro_dmc`: Hydrological cycle in mean climate

### 2×CO₂ Response Switches (`*_drsp`) 
- `log_clouds_drsp`: Cloud feedback in CO₂ response
- `log_vapor_drsp`: Water vapor feedback in CO₂ response
- `log_crcl_drsp`: Atmospheric circulation in CO₂ response
- `log_hydro_drsp`: Hydrological cycle in CO₂ response  

**Implementation Patterns**:
1. **Early returns**: Functions exit early when switches are disabled
2. **Climatology control**: Mean climate switches modify initialization
3. **Tendency zeroing**: Computed changes are zeroed based on switches

This enables precise process attribution studies and decomposition experiments.
"""

# ╔═╡ 07a88c94-93bb-4564-8737-980144c4af43
md"""
## ⚙️ Configuration

### 📐 Grid Dimensions & Time-Stepping

The model uses a **3.75° × 3.75°** latitude-longitude grid with **96** longitude points × **48** latitude points.

Time-stepping uses a **12-hour** main time step (`Δt`) and a **30-minute** sub-time step for circulation (`Δt_crcl`), giving **730** time steps per year.

Calendar arrays track days per month for output functions.
"""

# ╔═╡ 531589ab-c6b5-4048-9ba5-f9ad62ab00a6
begin
	# 📐 Grid dimensions ──────────────────────────────────
	xdim = 96                                 # number of longitude grid points
	ydim = 48                                 # number of latitude grid points
	dlon = 360.0 / xdim                       # longitude spacing [degrees]
	dlat = 180.0 / ydim                       # latitude spacing  [degrees]

	# ⏱️ Time stepping ────────────────────────────────────
	ndays_yr = 365                            # days per year (no leap years)
	Δt = 12 * 3600                            # main time step [s] (12 hours)
	Δt_crcl = round(Int, 0.5 * 3600)          # circulation sub-time step [s] (30 min)
	ndt_days = 24 * 3600 / Δt                 # time steps per day
	nstep_yr = Int(ndays_yr * ndt_days)            # time steps per year (= 730)
	ntime = max(1, round(Int, Δt / Δt_crcl))  # Number of sub-steps within one main time step

	# 📅 Calendar constants ───────────────────────────────
	jday_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
	jday_mon_cumsum = cumsum(jday_mon)
	
	# 🔢 Physical Constants & Numerical Limits ────────────
	min_T_K = 233.15       # # 273.15 - 40°C, minimum allowed surface temperature [K]
	max_humidity_change = 0.020      # Maximum humidity increment [kg/kg]
	min_humidity_change = 0.9      # Fraction of humidity that can be removed
end;

# ╔═╡ 813d59f7-eeea-402a-96e8-e5cbe2fd4582
md"""
## 𝚿 `mo_physics` - Physical Parameters & Climate Fields

This module defines all physical constants, model parameters, emissivity coefficients, and declares the climate field arrays.
"""

# ╔═╡ 19f106e4-2b82-47e1-9284-799a105f30cb
Base.@kwdef mutable struct PhysicsConfig
	# Mean Climate Switches
    log_clouds_dmc::Bool = true
    log_vapor_dmc::Bool = true
    log_ice_dmc::Bool = true
    log_crcl_dmc::Bool = true
    log_hydro_dmc::Bool = true
    log_atmos_dmc::Bool = true
    log_co2_dmc::Bool = true
    log_ocean_dmc::Bool = true
    log_qflux_dmc::Bool = true
    
    # CO₂ Response Switches
    log_clouds_drsp::Bool = true
    log_vapor_drsp::Bool = true
    log_ice_drsp::Bool = true
    log_crcl_drsp::Bool = true
    log_hydro_drsp::Bool = true
    log_topo_drsp::Bool = true
    log_humid_drsp::Bool = true
	log_ocean_drsp::Bool = true
    
    # Circulation Components
    log_ice::Bool = true
    log_hdif::Bool = true
    log_hadv::Bool = true
    log_vdif::Bool = true
    log_vadv::Bool = true
    log_conv::Bool = true
    
    # Hydrology Parameters
    log_rain::Int = 0
    log_eva::Int = -1
    log_clim::Int = 0
    
    # External Forcing
    log_tsurf_ext::Bool = false
    log_hwind_ext::Bool = false
    log_omega_ext::Bool = false
    
    # Experiment Type
    experiment::Symbol = :full_model  # :full_model, :constant_topo, :co2_double, etc.
    
    # CO₂ concentration for experiments (ppm)
    co2_concentration::Float64 = 340.0
    
    # Solar forcing multiplier
    solar_multiplier::Float64 = 1.0

	# Hydrology parameters (calculated by set_hydrology_parameters!)
    c_q::Float64 = 1.0     
    c_rq::Float64 = 0.0
    c_omega::Float64 = 0.0
    c_omegastd::Float64 = 0.0
end

# ╔═╡ bb82d1dc-8f08-4351-9371-a8efed1dd9bc
begin
"""
    create_experiment_config(experiment::Symbol) -> PhysicsConfig

Create a pre-configured `PhysicsConfig` for common experiment types.

# Experiments
- `:full_model` - All processes active (default)
- `:constant_topo` - Constant topography (log_topo_drsp = false)
- `:co2_double` - 2×CO₂ (680 ppm)
- `:co2_quadruple` - 4×CO₂ (1360 ppm)
- `:solar_plus27` - +27 W/m² solar constant
- `:elnino` - El Niño conditions
- `:lanina` - La Niña conditions  
- `:paleo_231kyr` - Paleoclimate (200 ppm CO₂)
- `:rcp85` - RCP8.5 climate change scenario
"""
function create_experiment_config(experiment::Symbol)::PhysicsConfig
    if experiment == :full_model
        return PhysicsConfig(experiment=:full_model)
        
    elseif experiment == :constant_topo
        cfg = PhysicsConfig(experiment=:constant_topo)
        cfg.log_topo_drsp = false  # Constant topography
        return cfg
        
    elseif experiment == :co2_double
        cfg = PhysicsConfig(experiment=:co2_double)
        cfg.co2_concentration = 680.0  # 2×CO₂
        return cfg
        
    elseif experiment == :co2_quadruple
        cfg = PhysicsConfig(experiment=:co2_quadruple)
        cfg.co2_concentration = 1360.0  # 4×CO₂
        return cfg
        
    elseif experiment == :solar_plus27
        cfg = PhysicsConfig(experiment=:solar_plus27)
        cfg.solar_multiplier = (1365.0 + 27.0) / 1365.0
        return cfg
        
    elseif experiment == :elnino
        cfg = PhysicsConfig(experiment=:elnino)
        return cfg
        
    elseif experiment == :lanina
        cfg = PhysicsConfig(experiment=:lanina)
        return cfg
        
    elseif experiment == :paleo_231kyr
        cfg = PhysicsConfig(experiment=:paleo_231kyr)
        cfg.co2_concentration = 200.0
        return cfg
        
    elseif experiment == :rcp85
        cfg = PhysicsConfig(experiment=:rcp85)
        return cfg
        
    else
        error("Unknown experiment: $experiment")
    end
end
end;

# ╔═╡ 54e63724-44ca-42b9-9246-7afb271b8154
md"""
### Hydrology Parameters
Optimized parameter system with lookup table for different parameterization modes and climatology datasets.
"""

# ╔═╡ 32b531ab-ee71-4685-af65-a2bae0b868f6
begin
	# 💧 Optimized Hydrology Parameter Lookup Table ────────────────────────
	const HYDRO_PARAMS = (
		-1 => (1.0, 0.0, 0.0, 0.0),                      # Original GREB
		 1 => (-1.391649, 3.018774, 0.0, 0.0),        	 # +Relative humidity
		 2 => (0.862162, 0.0, -29.02096, 0.0),        	 # +Omega convergence
		 3 => (-0.2685845, 1.4591853, -26.9858807, 0.0), # +RH & Omega
		 0 => (-1.88, 2.25, -17.69, 59.07)            	 # Best GREB (ERA-Interim)
	)     
	
	# 🎯 Cached Weight Arrays (avoid recomputation) ───
	global WZ_CACHE = Dict{Float64, Matrix{Float64}}()
end;

# ╔═╡ d0d74213-96ad-4a45-9a31-37f640a21a45
# ⚙️ Optimized Parameter Initialization ──────────────────────────────────
"""
    set_hydrology_parameters!(cfg::PhysicsConfig)

Initialize precipitation parameters `c_q, c_rq, c_omega, c_omegastd` based on
`cfg.log_rain` and `cfg.log_clim` settings.
"""
function set_hydrology_parameters!(cfg::PhysicsConfig)
	global c_q, c_rq, c_omega, c_omegastd
	
	# Fast lookup instead of if-else chain
	params = get(HYDRO_PARAMS, cfg.log_rain, (1.0, 0.0, 0.0, 0.0))
	c_q, c_rq, c_omega, c_omegastd = params
	
	# NCEP parameter adjustment
	if cfg.log_rain == 0 && cfg.log_clim == 1
		c_q, c_rq, c_omega, c_omegastd = -1.27, 1.99, -16.54, 21.15
	end
	
	@info "⚙️ MSCM hydrology: log_rain=$(cfg.log_rain), log_clim=$(cfg.log_clim) → (c_q=$c_q, c_rq=$c_rq, c_omega=$c_omega, c_omegastd=$c_omegastd)"
end

# ╔═╡ caec8065-9874-4d8c-8e7e-598a97506f06
md"""
## 🔧 Workspace Structs

Pre-allocated structures to eliminate runtime allocations.
"""

# ╔═╡ 80ac789e-4fe7-4184-9946-b8d7c24b04ea
begin
	"""
        CirculationWorkspace
    
    Pre-allocated buffers for diffusion, advection, and circulation calculations.
    Reused across all time steps to eliminate allocations.
    """
	mutable struct CirculationWorkspace
		# Polar sub-stepping buffers
		T1h::Vector{Float64}      # polar sub-stepping
		dTxh::Vector{Float64}     # polar tendencies
		
		# Circulation work arrays
		X_work::Matrix{Float64}   # circulation work array
		dX_diff::Matrix{Float64}  # diffusion output
		dX_adv::Matrix{Float64}   # advection output
		dX_conv::Matrix{Float64}  # convection output
		dX_crcl::Matrix{Float64}  # circulation output
		
		# Tendency buffers 
		temp_buf::Matrix{Float64}   # general workspace
        Q_sens_buf::Matrix{Float64} # Sensible heat flux buffer
		eva::Matrix{Float64}        # dq_eva
		rain::Matrix{Float64}       # dq_rain
		crcl::Matrix{Float64}     	# dq_crcl
		
		# State buffers
		Ts0_buf::Matrix{Float64}    # Surface temperature output
		Ta0_buf::Matrix{Float64}    # Air temperature output
		To0_buf::Matrix{Float64}    # Ocean temperature output
		q0_buf::Matrix{Float64}     # Humidity output
				
		# LW radiation buffers
		e_co2_buf::Matrix{Float64}    # spatial CO₂ buffer
		e_vapor_buf::Matrix{Float64}  # spatial water vapor buffer
		em_buf::Matrix{Float64}       # spatial emissivity buffer
		LW_surf_buf::Matrix{Float64}  # Surface longwave
    	LW_down_buf::Matrix{Float64}  # Downwelling longwave
    	LW_up_buf::Matrix{Float64}    # Upwelling longwave
		
		# Hydrology buffers
		qs::Matrix{Float64}        # Saturation humidity buffer
		Tskin::Matrix{Float64}     # Skin temperature buffer
		rq::Matrix{Float64}        # Relative humidity buffer
		ws_base::Matrix{Float64}   # Base wind speed buffer
		# Hydrology
    	Q_lat_buf::Matrix{Float64}
    	Q_lat_air_buf::Matrix{Float64}
    	dq_eva_buf::Matrix{Float64}
    	dq_rain_buf::Matrix{Float64}
        cE_buf::Matrix{Float64}       # Surface exchange coefficient buffer

    	# Deep_ocean
    	dT_ocean_buf::Matrix{Float64}
    	dTo_buf::Matrix{Float64}

		# Dedicated circulation output
		dTa_crcl::Matrix{Float64}   # temperature tendency
		dq_crcl::Matrix{Float64}    # humidity tendency

        # SWradiation
        ice_cover_buf::Matrix{Float64} # ice fraction
        a_surf_buf::Matrix{Float64}    # surface albedo
        albedo_buf::Matrix{Float64}    # combined albedo (surface + atmosphere)
        a_atmos_buf::Matrix{Float64}   # atmospheric albedo
        sw_buf::Matrix{Float64}        # net shortwave flux

		# time_loop
		precip_out::Matrix{Float64}   # precipitation output
		evap_out::Matrix{Float64}     # evaporation output
        qcrcl_out::Matrix{Float64}    # circulation moisture output
		term_north::Vector{Float64}      # northern boundary term
        term_south::Vector{Float64}      # southern boundary term
	end
	
	function CirculationWorkspace()
		CirculationWorkspace(
			zeros(Float64, xdim),		# T1h
			zeros(Float64, xdim),		# dTxh
			zeros(Float64, xdim, ydim),	# X_work
			zeros(Float64, xdim, ydim),	# dX_diff
			zeros(Float64, xdim, ydim),	# dX_adv
			zeros(Float64, xdim, ydim),	# dX_conv
			zeros(Float64, xdim, ydim),	# dX_crcl
			zeros(Float64, xdim, ydim),	# temp_buf
            zeros(Float64, xdim, ydim), # Q_sens_buf
			zeros(Float64, xdim, ydim), # rain	
			zeros(Float64, xdim, ydim),	# eva
			zeros(Float64, xdim, ydim),	# crcl
			# State buffers
			zeros(Float64, xdim, ydim),	# Ts0_buf
			zeros(Float64, xdim, ydim),	# Ta0_buf
			zeros(Float64, xdim, ydim),	# To0_buf
			zeros(Float64, xdim, ydim),	# q0_buf
			# LWradiation buffers
			zeros(Float64, xdim, ydim),	# e_co2_buf
			zeros(Float64, xdim, ydim),	# e_vapor_buf
			zeros(Float64, xdim, ydim),	# em_buf
			zeros(Float64, xdim, ydim),	# LW_surf_buf
			zeros(Float64, xdim, ydim),	# LW_down_buf
			zeros(Float64, xdim, ydim),	# LW_up_buf
			# Hydrology buffers
			zeros(Float64, xdim, ydim),	# qs
			zeros(Float64, xdim, ydim),	# Tskin
			zeros(Float64, xdim, ydim),	# rq
			zeros(Float64, xdim, ydim),	# ws_base
			# Hydrology
        	zeros(Float64, xdim, ydim),  # Q_lat_buf
        	zeros(Float64, xdim, ydim),  # Q_lat_air_buf
        	zeros(Float64, xdim, ydim),  # dq_eva_buf
        	zeros(Float64, xdim, ydim),  # dq_rain_buf
            zeros(Float64, xdim, ydim),  # cE_buf
        	# Deep_ocean
        	zeros(Float64, xdim, ydim),  # dT_ocean_buf
        	zeros(Float64, xdim, ydim),  # dTo_buf
			# Dedicated circulation output
			zeros(Float64, xdim, ydim),  # dTa_crcl
        	zeros(Float64, xdim, ydim),  # dq_crcl
            # SWradiation buffers
            zeros(Float64, xdim, ydim),  # ice_cover_buf
            zeros(Float64, xdim, ydim),  # a_surf_buf
            zeros(Float64, xdim, ydim),  # albedo_buf
            zeros(Float64, xdim, ydim),  # a_atmos_buf
            zeros(Float64, xdim, ydim),   # sw_buf
			# time_loop
			zeros(Float64, xdim, ydim),  # precip_out
			zeros(Float64, xdim, ydim),  # evap_out
        	zeros(Float64, xdim, ydim),  # qcrcl_out
			zeros(Float64, xdim),  # term_north
        	zeros(Float64, xdim),  # term_south
		)
	end
	# Global instance (reused across all time steps)
	global circ_workspace = CirculationWorkspace()
end;

# ╔═╡ b9f7bde9-0aa4-4075-8fb9-14f84db0b3fa
begin
	"""
        MonthlyAccumulator
    
    Accumulates fields over a month for monthly-mean output.
    Reset after each month via `reset!`.
    """
	mutable struct MonthlyAccumulator
		Tmm::Matrix{Float64}          # Surface temperature accumulator
		Tamm::Matrix{Float64}         # Air temperature accumulator
		Tomm::Matrix{Float64}         # Ocean temperature accumulator
		qmm::Matrix{Float64}          # Humidity accumulator
		apmm::Matrix{Float64}         # Albedo accumulator
		icemm::Matrix{Float64}        # Ice fraction accumulator
		precipmm::Matrix{Float64}     # Precipitation accumulator
		evapmm::Matrix{Float64}       # Evaporation accumulator
		qcrclmm::Matrix{Float64}      # Circulation moisture accumulator
		swmm::Matrix{Float64}         # Shortwave radiation accumulator
		lwmm::Matrix{Float64}         # Longwave radiation accumulator
		qlatmm::Matrix{Float64}       # Latent heat accumulator
		qsensmm::Matrix{Float64}      # Sensible heat accumulator
		count::Int                    # Number of accumulations
	end
	
	function MonthlyAccumulator()
		MonthlyAccumulator(
			zeros(Float64, xdim, ydim),  # Tmm
			zeros(Float64, xdim, ydim),  # Tamm
			zeros(Float64, xdim, ydim),  # Tomm
			zeros(Float64, xdim, ydim),  # qmm
			zeros(Float64, xdim, ydim),  # apmm
			zeros(Float64, xdim, ydim),  # icemm
			zeros(Float64, xdim, ydim),  # precipmm
			zeros(Float64, xdim, ydim),  # evapmm
			zeros(Float64, xdim, ydim),  # qcrclmm
			zeros(Float64, xdim, ydim),  # swmm
			zeros(Float64, xdim, ydim),  # lwmm
			zeros(Float64, xdim, ydim),  # qlatmm
			zeros(Float64, xdim, ydim),  # qsensmm
			0
		)
	end
	
	function reset!(acc::MonthlyAccumulator)
		fill!(acc.Tmm, 0.0)
		fill!(acc.Tamm, 0.0)
		fill!(acc.Tomm, 0.0)
		fill!(acc.qmm, 0.0)
		fill!(acc.apmm, 0.0)
		fill!(acc.icemm, 0.0)
		fill!(acc.precipmm, 0.0)
		fill!(acc.evapmm, 0.0)
		fill!(acc.qcrclmm, 0.0)
		fill!(acc.swmm, 0.0)
		fill!(acc.lwmm, 0.0)
		fill!(acc.qlatmm, 0.0)
		fill!(acc.qsensmm, 0.0)
		acc.count = 0
	end
	
	function accumulate!(acc::MonthlyAccumulator, Ts, Ta, To, q, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens)
		@. acc.Tmm  += Ts
		@. acc.Tamm += Ta
		@. acc.Tomm += To
		@. acc.qmm  += q
		@. acc.apmm += albedo
		@. acc.icemm += ice
		@. acc.precipmm += precip
		@. acc.evapmm += evap
		@. acc.qcrclmm += qcrcl
		@. acc.swmm += sw
		@. acc.lwmm += lw
		@. acc.qlatmm += qlat
		@. acc.qsensmm += qsens
		acc.count += 1
	end
end;

# ╔═╡ 12a89062-bc27-4e36-8b5e-1670b4bd0b28
md"""
### Spatial CO₂ Fields
Spatial masking arrays for regional/partial CO₂ experiments.
"""

# ╔═╡ f39520b7-a246-4980-b8b9-0215367d0b46
begin
	# 🌍 Spatial CO₂ masking arrays ────────────────────────────────
	# Spatial fraction for regional CO₂ experiments
	co2_part = ones(Float64, xdim, ydim)    # Regional CO₂ mask (1.0 = full CO₂, 0.5 = half CO₂)
	co2_part_scn = ones(Float64, xdim, ydim) # Scenario-specific spatial mask
end;

# ╔═╡ 1879fd40-0117-478b-acdd-18112738b81f
md"""
### Natural Constants
Fundamental physical constants used throughout the model.
"""

# ╔═╡ f8a2c2de-5045-4ab6-a6fa-7bca502afc9b
begin
	# ── Natural constants ────────────────────────────────────────────
	const_pi   = pi               # π (model precision)
	σ          = 5.6704e-8        # Stefan-Boltzmann constant [W/m²/K⁴]
	ρ_ocean    = 999.1            # density of water at T=15°C [kg/m³]
	ρ_land     = 2600.0           # density of solid rock [kg/m³]
	ρ_air      = 1.2              # density of air at 20°C at sea 
	grav       = 9.80665          # gravitational acceleration [m/s²]
	cp_ocean   = 4186.0           # specific heat of water at T=15°C [J/kg/K]
	cp_land    = cp_ocean / 4.5   # specific heat of dry land [J/kg/K]
	cp_air     = 1005.0           # specific heat of air [J/kg/K]
	ε          = 1.0              # emissivity for IR
end;

# ╔═╡ b23bc922-e9bc-4012-9782-1258e3cc8e7b
begin
	# ── Column depths [m] ───────────────────────────────────────────
	d_ocean   = 50.0                        # ocean column
	d_land    = 2.0                         # land column 
	d_air     = 5000.0                      # air column 

	# ── Heat capacities [J/K/m²] ────────────────────────────────────
	cap_ocean = cp_ocean * ρ_ocean          # 1m ocean
	cap_land  = cp_land * ρ_land * d_land   # land column
	cap_air   = cp_air * ρ_air * d_air      # air column 

	# ── Sensible heat [W/K/m²] ──────────────────────────────────────
	ct_sens   = 22.5                        # sensible heat coupling 

	# ── Albedo parameters ───────────────────────────────────────────
	da_ice    = 0.25                       # albedo increase for ice-cover
	a_no_ice  = 0.1                        # albedo ice-free
	a_cloud   = 0.35                       # cloud albedo

	# ── Ice/snow temperature thresholds [K] ─────────────────────────
	Tl_ice1 = 273.15 - 10.0                 # land: full ice albedo
    Tl_ice2 = 273.15                        # land: no ice albedo
    To_ice1 = 273.15 - 7.0                  # ocean: full ice
    To_ice2 = 273.15 - 1.7                  # ocean: no ic
	
	# Precomputed inverse ranges (avoids division in hot loops)
	inv_To_ice_range = 1.0 / (To_ice2 - To_ice1)
	inv_Tl_ice_range = 1.0 / (Tl_ice2 - Tl_ice1)

	# ── Deep ocean ──────────────────────────────────────────────────
	co_turb   = 5.0                        # turbulent mixing coefficient [W/K/m²]
	c_effmix  = 0.5						   # mixing efficiency
	turb_coeff = Δt * co_turb / cap_ocean  # precomputed mixing coefficient

	# ── Atmospheric transport ───────────────────────────────────────
	κ = 8e5                        # diffusion coefficient [m²/s]

	# ── Latent heat / hydrology ─────────────────────────────────────
	ce        = 2e-3                   # latent heat transfer coefficient
	cq_latent = 2.257e6                # latent heat of evaporation [J/kg]
	cq_rain   = -0.1 / 24.0 / 3600.0   # rain-related vapor decrease [1/s]

	# ── Scaling heights [m] ─────────────────────────────────────────
	z_air     = 8400.0                 # heat & CO₂ scaling height
	z_vapor   = 5000.0                 # water vapor scaling height
	const_factor = Δt_crcl / z_vapor * 2.5 / (ρ_air * grav) # precompute

	# ── Regression factor [kg/m³] ───────────────────────────────────
	r_qviwv   = 2.6736e3               # VIWV ↔ q_air regression factor
	conv_factor = r_qviwv * 86400.0	   # kg/kg → mm/day conversion

	# ── solar factor [%] ────────────────────────────────────────────
	S0_var = 100.0  		     	   # default 100%

	# ── transport geometry & coefficients ───────────────────────────
	deg_grid  = 2.0 * const_pi * 6.371e6 / 360.0
	dyy_grid  = dlat * deg_grid
	lat_grid  = [dlat * k - dlat / 2.0 - 90.0 for k in 1:ydim]
	dxlat_grid = [dlon * deg_grid * cos(2.0 * const_pi / 360.0 * lat_grid[k]) for k in 1:ydim]

	# ── Diffusion coefficients ──────────────────────────────────────
	ccy_diff  = κ * Δt_crcl / dyy_grid^2
	ccx_diff  = [κ * Δt_crcl / dxlat_grid[k]^2 for k in 1:ydim]
	
	# ── Advection coefficients ──────────────────────────────────────
	ccy_adv   = Δt_crcl / dyy_grid / 2.0
	ccx_adv   = [Δt_crcl / dxlat_grid[k] / 2.0 for k in 1:ydim]

	# ── periodic longitude neighbour indices ───────────────────────
	const lon_jm1 = [mod1(i-1, xdim) for i in 1:xdim]
	const lon_jp1 = [mod1(i+1, xdim) for i in 1:xdim]
	const lon_jm2 = [mod1(i-2, xdim) for i in 1:xdim]
	const lon_jp2 = [mod1(i+2, xdim) for i in 1:xdim]
	const lon_jm3 = [mod1(i-3, xdim) for i in 1:xdim]
	const lon_jp3 = [mod1(i+3, xdim) for i in 1:xdim]
end

# ╔═╡ 0dbfb663-46e7-4873-ac77-1e8e392fe69d
begin
    # ☀️ Solar forcing storage
    global sw_solar_forcing_data = nothing  # Will hold (48, 730) array when loaded
	global sw_solar_forcing_state = Ref(1.0)  # Runtime solar multiplier used by SWradiation!

	# 🔧 Flux correction arrays (initialised with zeros, overwritten if files exist)
	global TF_correct  = zeros(Float64, xdim, ydim, nstep_yr)
	global qF_correct  = zeros(Float64, xdim, ydim, nstep_yr)
	global ToF_correct = zeros(Float64, xdim, ydim, nstep_yr)
	const ΔT_AIR_FACTOR = Δt / cap_air
end;

# ╔═╡ f578f25e-047e-4a7e-8483-d544c7b4bec3
function load_flux_corrections_jdal2!(jdal2_dir::String)
    """Load flux corrections from JDAL2 files (zeros if missing)"""
    # Correction arrays are already defined as zeros in the global scope.
    correction_files = Dict(
        "Tsurf_flux_correction.jd2" => TF_correct,
        "vapour_flux_correction.jd2" => qF_correct,
        "Tocean_flux_correction.jd2" => ToF_correct
    )

    for (filename, array) in correction_files
        filepath = joinpath(jdal2_dir, "climatology", filename)
        if isfile(filepath)
            result = read_jdal2(filepath)
            array .= result.data
            println("✅ Loaded $filename")
        else
            fill!(array, 0.0)
            @warn "$filename not found, using zeros"
        end
    end
end

# ╔═╡ 0308940f-16bc-4b03-8321-9ea8c91a21c7
md"""
### Emissivity Regression Parameters
10-element vector used in the LW radiation subroutine to compute atmospheric emissivity from CO₂, water vapor, and cloud cover.
"""

# ╔═╡ 898dd5aa-5273-4833-9a4a-0f2b94cc8d38
p_emi = [9.0721, 106.7252, 61.5562, 0.0179, 0.0028,
         0.0570,   0.3462,  2.3406, 0.7032, 1.0662]

# ╔═╡ 78605b5f-b03d-4851-b339-aea63bad3688
md"""
### Climate Field Arrays (Placeholders)

These arrays will be populated with real data once available. For now we initialize them with zeros so the notebook is runnable.

**2D fields** `(xdim, ydim)` = (96, 48):
- `z_topo` — topography [m] (< 0 = ocean)
- `glacier` — glacier mask (> 0.5 = glaciated)
- `z_ocean` — derived ocean depth (3× max MLD)
- `cap_surf` — surface heat capacity [J/K/m²]

**3D fields** `(xdim, ydim, nstep_yr)` = (96, 48, 730):
- `Tclim` — climatological surface temperature [K]
- `uclim`, `vclim` — zonal/meridional wind [m/s]
- `qclim` — atmospheric humidity [kg/kg]
- `mldclim` — mixed-layer depth [m]
- `Toclim` — deep ocean temperature [K]
- `cldclim` — cloud cover fraction
- `swetclim` — soil wetness [0–1]
- `omegaclim` — vertical velocity [Pa/s] (MSCM)
- `omegastdclim` — omega standard deviation [Pa/s] (MSCM)
- `wsclim` — wind speed climatology [m/s] (MSCM)
- `uclim_m/p`, `vclim_m/p` — negative/positive wind components (MSCM)
- `TF_correct`, `qF_correct`, `ToF_correct` — flux corrections
- `dTrad` — Tatmos–radiation offset

**2D solar** `(ydim, nstep_yr)` = (48, 730):
- `sw_solar` — 24-hr mean solar radiation [W/m²]
"""

# ╔═╡ 75f3b78f-924a-4a65-b2e3-a79c6f2082f9
begin
	# 🗺️ 2D fields (xdim, ydim) ──────────────────────────────────────
	z_topo   = zeros(Float64, xdim, ydim)  # topography [m] (<0: ocean)
	glacier  = zeros(Float64, xdim, ydim)  # glacier mask (>0.5: glacier)
	z_ocean  = zeros(Float64, xdim, ydim)  # derived ocean depth [m]
	cap_surf = zeros(Float64, xdim, ydim)  # surface heat capacity [J/K/m²]
	wz_air   = zeros(Float64, xdim, ydim)  # exp(-z_topo / z_air)
	wz_vapor = zeros(Float64, xdim, ydim)  # exp(-z_topo / z_vapor)
end;

# ╔═╡ 0be9bc61-1a59-4dfe-84f9-bb1a27ca30fc
begin
	# 🌡️ 3D climate fields (xdim, ydim, nstep_yr) ───────────────────
	Tclim   = zeros(Float64, xdim, ydim, nstep_yr)   # surface temperature [K]
	uclim   = zeros(Float64, xdim, ydim, nstep_yr)   # zonal wind [m/s]
	vclim   = zeros(Float64, xdim, ydim, nstep_yr)   # meridional wind [m/s]
	qclim   = zeros(Float64, xdim, ydim, nstep_yr)   # atmospheric humidity [kg/kg]
	mldclim = zeros(Float64, xdim, ydim, nstep_yr)   # mixed-layer depth [m]
	
	# additional climatology fields
	omegaclim    = zeros(Float64, xdim, ydim, nstep_yr) # vertical velocity [Pa/s]
	omegastdclim = zeros(Float64, xdim, ydim, nstep_yr) # omega std deviation [Pa/s]
	wsclim       = zeros(Float64, xdim, ydim, nstep_yr) # wind speed [m/s]

	# 📊 Anomaly Fields for ENSO/Climate Change Experiments ────────────────
	# ENSO anomaly fields
	Tclim_anom_enso  = zeros(Float64, xdim, ydim, nstep_yr) # surface temperature [K]
	uclim_anom_enso  = zeros(Float64, xdim, ydim, nstep_yr) # zonal wind [m/s]
	vclim_anom_enso  = zeros(Float64, xdim, ydim, nstep_yr) # meridional wind [m/s]
	omegaclim_anom_enso = zeros(Float64, xdim, ydim, nstep_yr) # vertical velocity [Pa/s]
	wsclim_anom_enso = zeros(Float64, xdim, ydim, nstep_yr) # wind speed [m/s]
	
	# Climate change anomaly fields
	Tclim_anom_cc   = zeros(Float64, xdim, ydim, nstep_yr) # surface temperature [K]
	uclim_anom_cc   = zeros(Float64, xdim, ydim, nstep_yr) # zonal wind [m/s] 
	vclim_anom_cc   = zeros(Float64, xdim, ydim, nstep_yr) # meridional wind [m/s]
	omegaclim_anom_cc = zeros(Float64, xdim, ydim, nstep_yr) # vertical velocity [Pa/s]
	wsclim_anom_cc  = zeros(Float64, xdim, ydim, nstep_yr) # wind speed [m/s]
		
	# 🌬️ Precomputed wind sign splits ──────────────────
	uclim_m  = zeros(Float64, xdim, ydim, nstep_yr)   # negative u components
	uclim_p  = zeros(Float64, xdim, ydim, nstep_yr)   # positive u components  
	vclim_m  = zeros(Float64, xdim, ydim, nstep_yr)   # negative v components
	vclim_p  = zeros(Float64, xdim, ydim, nstep_yr)   # positive v components
	
	# Initialize wind component separation (CRITICAL: affects advection)
	@. uclim_m = ifelse(uclim >= 0.0, uclim, 0.0)  # positive winds only
	@. uclim_p = ifelse(uclim < 0.0, uclim, 0.0)   # negative winds only
	@. vclim_m = ifelse(vclim >= 0.0, vclim, 0.0)  # positive winds only
	@. vclim_p = ifelse(vclim < 0.0, vclim, 0.0)   # negative winds only
	Toclim   = zeros(Float64, xdim, ydim, nstep_yr)   # deep ocean temperature [K]
	cldclim  = zeros(Float64, xdim, ydim, nstep_yr)   # cloud cover fraction
	swetclim = zeros(Float64, xdim, ydim, nstep_yr)   # soil wetness [0-1]
end;

# ╔═╡ 0d200ce7-eadb-41de-8e31-45783a1faab9
begin
	# ☀️ 2D solar field (ydim, nstep_yr) ─────────────────────────────
	sw_solar = zeros(Float64, ydim, nstep_yr) # 24hr mean solar radiation [W/m²]
	global dTrad = zeros(Float64, xdim, ydim, Int(nstep_yr))  # Tatmos-radiation offset
end;

# ╔═╡ 2bf0fe8e-5718-4c1e-863b-85db7b3ae7f3
function load_greb_jdal2!(jdal2_dir::String; dataset::Symbol=:ncep)
"""Load all GREB input data from JDAL2 formatted files"""

	if !isdir(jdal2_dir)
		error("JDAL2 directory not found: $jdal2_dir")
	end

	# Static 2D files
    println("📂 Loading static fields...")
    topo_result = read_jdal2(joinpath(jdal2_dir, "static", "global.topography.jd2"))
    global z_topo .= topo_result.data
    
    glacier_result = read_jdal2(joinpath(jdal2_dir, "static", "greb.glaciers.jd2"))
    global glacier .= glacier_result.data

	# Dataset-specific file mapping
    file_map = Dict(
        :ncep => Dict(
            "Tclim" => "ncep.tsurf.1948-2007.clim.jd2",
            "uclim" => "ncep.zonal_wind.850hpa.clim.jd2",
            "vclim" => "ncep.meridional_wind.850hpa.clim.jd2",
            "qclim" => "ncep.atmospheric_humidity.clim.jd2",
            "swetclim" => "ncep.soil_moisture.clim.jd2"
        ),
        :era => Dict(
            "Tclim" => "erainterim.tsurf.1979-2015.clim.jd2",
            "uclim" => "erainterim.zonal_wind.850hpa.clim.jd2",
            "vclim" => "erainterim.meridional_wind.850hpa.clim.jd2",
            "qclim" => "erainterim.atmospheric_humidity.clim.jd2",
            "swetclim" => "ncep.soil_moisture.clim.jd2"
        )
    )

	# Use mixed dataset as fallback
    files = get(file_map, dataset, file_map[:ncep])

	println("📂 Loading 3D climatology ($dataset dataset)...")
    climatology_dir = joinpath(jdal2_dir, "climatology")
    
    # Load each variable individually with unique result names
    tsurf_result = read_jdal2(joinpath(climatology_dir, files["Tclim"]))
    global Tclim .= tsurf_result.data
    
    uwind_result = read_jdal2(joinpath(climatology_dir, files["uclim"]))
    global uclim .= uwind_result.data
    
    vwind_result = read_jdal2(joinpath(climatology_dir, files["vclim"]))
    global vclim .= vwind_result.data
    
    humid_result = read_jdal2(joinpath(climatology_dir, files["qclim"]))
    global qclim .= humid_result.data
    
    swet_result = read_jdal2(joinpath(climatology_dir, files["swetclim"]))
    global swetclim .= swet_result.data

	# Common climatology files
    println("📂 Loading common climatology fields...")
    
    cld_result = read_jdal2(joinpath(climatology_dir, "isccp.cloud_cover.clim.jd2"))
    global cldclim .= cld_result.data
    
    mld_result = read_jdal2(joinpath(climatology_dir, "woce.ocean_mixed_layer_depth.clim.jd2"))
    global mldclim .= mld_result.data
    
    tocean_result = read_jdal2(joinpath(climatology_dir, "Tocean.clim.jd2"))
    global Toclim .= tocean_result.data
    
    omega_result = read_jdal2(joinpath(climatology_dir, "erainterim.omega.vertmean.clim.jd2"))
    global omegaclim .= omega_result.data
    
    omegastd_result = read_jdal2(joinpath(climatology_dir, "erainterim.omega_std.vertmean.clim.jd2"))
    global omegastdclim .= omegastd_result.data
    
    ws_result = read_jdal2(joinpath(climatology_dir, "erainterim.windspeed.850hpa.clim.jd2"))
    global wsclim .= ws_result.data
	
    # Solar radiation (special: lat × time)
    println("📂 Loading solar radiation...")
    solar_path = joinpath(jdal2_dir, "solar", "solar_radiation.clim.jd2")
    if isfile(solar_path)
        solar_result = read_jdal2(solar_path)
        @assert size(solar_result.data) == (ydim, nstep_yr) "Wrong solar dimensions"
        global sw_solar .= solar_result.data
    else
        error("Solar radiation file not found: $solar_path")
    end
	
    # Optional: Load flux corrections
    println("📂 Loading flux corrections...")
    load_flux_corrections_jdal2!(jdal2_dir)
    
    # Update wind sign splits
    @. uclim_m = ifelse(uclim >= 0.0, uclim, 0.0)
    @. uclim_p = ifelse(uclim < 0.0, uclim, 0.0)
    @. vclim_m = ifelse(vclim >= 0.0, vclim, 0.0)
    @. vclim_p = ifelse(vclim < 0.0, vclim, 0.0)
    
    println("✅ All GREB data loaded successfully from JDAL2")
end

# ╔═╡ 22d5e751-f095-42d7-9c24-78e546b3ce37
md"""
### Time-Step Indices
`jday` and `ityr` track the current calendar day and time-step within the year during integration. In Fortran these were module-level integers; here they start at 1.
"""

# ╔═╡ 047b312f-8d6c-4732-aa0b-bea3de3e99e2
begin
	# 🕐 Time State Struct
	mutable struct TimeState
		jday::Int  # Current calendar day in year [1..365]
		ityr::Int  # Current timestep in year [1..730]
	end

	# 📅 Precomputed Calendar Lookup
	const max_timesteps = 200 * nstep_yr  # Support up to 200-year runs
	const calendar_lookup = [(
		day = mod((it - 1) ÷ ndt_days, ndays_yr) + 1,
		step = mod(it - 1, nstep_yr) + 1
	) for it in 1:max_timesteps]

	const polar_treshold = 2.5e5  # 250 km in meters
	const IS_POLAR = [dxlat_grid[k] <= polar_treshold for k in 1:ydim]
end;

# ╔═╡ 698a990b-cd78-43a1-a747-e79102767d62
begin
	using BenchmarkTools, Profile, ProfileSVG
# Create test data for benchmarking
function setup_benchmark()
    xdim_local = 96
    ydim_local = 48
    
    # Create workspace
    ws = CirculationWorkspace()
    
    # Ensure dX_crcl is properly sized
    if size(ws.dX_crcl) != (xdim_local, ydim_local)
        ws.dX_crcl = zeros(Float64, xdim_local, ydim_local)
    end
    
    # Time state
    timestate = TimeState(1, 1)
    
    # Configuration
    cfg = PhysicsConfig()
    
    # Climate fields
    Ts = rand(Float64, xdim_local, ydim_local) .* 80 .+ 233.15
    Ta = copy(Ts)
    To = rand(Float64, xdim_local, ydim_local) .* 30 .+ 270.0
    q = rand(Float64, xdim_local, ydim_local) .* 0.02
    CO2 = 340.0
    
    # For advection and convergence tests
    T_test = rand(Float64, xdim_local, ydim_local)
    
    # Create dummy omegaclim if not available
    if !@isdefined(omegaclim)
        omegaclim_local = rand(Float64, xdim_local, ydim_local, 730) .* 0.1
    else
        omegaclim_local = omegaclim
    end
    
    return (ws=ws, timestate=timestate, cfg=cfg, 
            Ts=Ts, Ta=Ta, To=To, q=q, CO2=CO2,
            z_air=z_air, z_vapor=z_vapor,
            T_test=T_test, omegaclim=omegaclim_local)
end

# Setup once
bench_data = setup_benchmark()
end;

# ╔═╡ d403b01e-7d9f-4778-813d-bbbb0bfdbfb9
md"""
## ⚙️ `mo_diagnostics` - Diagnostic & Output Fields

This module declares the accumulator arrays for annual-mean diagnostics and monthly-mean output.
"""

# ╔═╡ c94558a5-0ae2-4563-93d7-0609ea3cce41
md"""
### Annual-Mean Accumulators
These fields accumulate values over a year and are divided by `nstep_yr` at year-end in the `diagnostics` subroutine, then reset to zero.

| Variable | Description |
|:---------|:------------|
| `Tsmn` | Surface temperature |
| `Tamn` | Air temperature |
| `Tomn` | Deep ocean temperature |
| `qmn` | Atmospheric humidity |
| `amn` | Albedo |
| `icmn` | Ice-cover fraction (MSCM) |
| `prmn` | Precipitation tendency (MSCM) |
| `evmn` | Evaporation tendency (MSCM) |
| `qcrclmn` | Humidity circulation tendency (MSCM) |
| `swmn` | Shortwave radiation |
| `lwmn` | Longwave radiation |
| `qlatmn` | Latent heat flux |
| `qsensmn` | Sensible heat flux |
| `ftmn` | Temperature flux correction |
| `fqmn` | Humidity flux correction |
"""

# ╔═╡ c06deaa7-2e6a-4143-a195-b473cbd84329
begin
	# Annual-mean accumulators (xdim, ydim) 
	Tsmn    = zeros(Float64, xdim, ydim)   # surface temperature
	Tamn    = zeros(Float64, xdim, ydim)   # air temperature
	Tomn    = zeros(Float64, xdim, ydim)   # deep ocean temperature
	qmn     = zeros(Float64, xdim, ydim)   # humidity
	amn     = zeros(Float64, xdim, ydim)   # albedo
	
	# additional accumulators
	icmn    = zeros(Float64, xdim, ydim)   # ice cover fraction
	prmn    = zeros(Float64, xdim, ydim)   # precipitation tendency
	evmn    = zeros(Float64, xdim, ydim)   # evaporation tendency
	qcrclmn = zeros(Float64, xdim, ydim)   # circulation tendency
	
	swmn    = zeros(Float64, xdim, ydim)   # shortwave radiation
	lwmn    = zeros(Float64, xdim, ydim)   # longwave radiation
	qlatmn  = zeros(Float64, xdim, ydim)   # latent heat flux
	qsensmn = zeros(Float64, xdim, ydim)   # sensible heat flux
	ftmn    = zeros(Float64, xdim, ydim)   # temperature flux correction
	fqmn    = zeros(Float64, xdim, ydim)   # humidity flux correction
end;

# ╔═╡ 2d02d36f-efee-411e-befb-52d61b3fb16e
md"""
### Monthly-Mean Output Buffers
Accumulated within each calendar month, divided by month length at month-end.

| Variable | Description |
|:---------|:------------|
| `Tmm` | Surface temperature |
| `Tamm` | Air temperature |
| `Tomm` | Deep ocean temperature |
| `qmm` | Atmospheric humidity |
| `apmm` | Albedo |
| `icmm` | Ice-cover fraction (MSCM) |
| `prmm` | Precipitation tendency (MSCM) |
| `evmm` | Evaporation tendency (MSCM) |
| `qcrclmm` | Humidity circulation tendency (MSCM) |
"""

# ╔═╡ a8b5fa01-0526-46f4-9aa0-31b52398a2bd
begin
	# Monthly-mean buffers (xdim, ydim) 
	Tmm     = zeros(Float64, xdim, ydim)   # surface temperature
	Tamm    = zeros(Float64, xdim, ydim)   # air temperature
	Tomm    = zeros(Float64, xdim, ydim)   # deep ocean temperature
	qmm     = zeros(Float64, xdim, ydim)   # humidity
	apmm    = zeros(Float64, xdim, ydim)   # albedo
		
	# additional monthly buffers
	icmm     = zeros(Float64, xdim, ydim)   # ice cover fraction
	prmm     = zeros(Float64, xdim, ydim)   # precipitation tendency
	evmm     = zeros(Float64, xdim, ydim)   # evaporation tendency
	qcrclmm  = zeros(Float64, xdim, ydim)   # circulation tendency
	swmm     = zeros(Float64, xdim, ydim)   # shortwave radiation
	lwmm     = zeros(Float64, xdim, ydim)   # longwave radiation
	qlatmm   = zeros(Float64, xdim, ydim)   # latent heat flux
	qsensmm  = zeros(Float64, xdim, ydim)   # sensible heat flux
end;

# ╔═╡ 1bdeea50-a73b-47a5-b941-d9d3374f54cf
begin
	# Control run monthly means (for anomaly calculation)
	Tmn_ctrl    = zeros(Float64, xdim, ydim, 12)  # surface temperature
	Tamn_ctrl   = zeros(Float64, xdim, ydim, 12)  # air temperature
	Tomn_ctrl   = zeros(Float64, xdim, ydim, 12)  # deep ocean temperature
	qmn_ctrl    = zeros(Float64, xdim, ydim, 12)  # humidity
	icmn_ctrl   = zeros(Float64, xdim, ydim, 12)  # ice cover
	prmn_ctrl   = zeros(Float64, xdim, ydim, 12)  # precipitation
	evamn_ctrl  = zeros(Float64, xdim, ydim, 12)  # evaporation
	qcrclmn_ctrl = zeros(Float64, xdim, ydim, 12) # circulation
end;

# ╔═╡ 1b38b24d-bb82-4d52-a26d-850d142e25d5
md"""
## 🎮 `greb_model` - Main Driver

The driver orchestrates four phases:

1. **Initialization** - Derived fields and sensitivity overrides
2. **Flux correction** - Spin-up to compute correction fields
3. **Control run** - Integration at fixed CO₂
4. **Scenario run** - Integration with time-varying CO₂

Output is stored as in-memory `Vector` of `MonthlyRecord` named tuples.
"""

# ╔═╡ 0ff4be7c-580c-4db1-a874-175badc8c11a
md"""
### Monthly Output Record Type
Each monthly snapshot stores nine 2-D fields as a `NamedTuple`: 
`Ts`, `Ta`, `To`, `q`, `albedo`, `ice`, `precip`, `evap`, `qcrcl`.
This replaces the Fortran direct-access binary writes to units 21/22.
"""

# ╔═╡ 7a06bf0d-a61c-4d28-b144-2725fe90ae62
# Type alias for one monthly output record
const MonthlyRecord = NamedTuple{(:Ts, :Ta, :To, :q, :albedo, :ice, :precip, :evap, :qcrcl, :sw, :lw, :qlat, :qsens), NTuple{13, Matrix{Float64}}};

# ╔═╡ 3beff4df-89da-4f25-a0e4-be03048c5f2c
md"""
## `init_model!`

Initializes the model state from climatology and experiment configuration:

- Sets hydrology parameters (`c_q`, `c_rq`, `c_omega`, `c_omegastd`) via `set_hydrology_parameters!`
- Computes derived fields: `dTrad` (radiation offset), `z_ocean` (3× max mixed-layer depth), `wz_air`/`wz_vapor` (pressure weights)
- Applies experiment overrides based on `cfg` switches: zeroes clouds, fixes topography, adjusts humidity/ocean climatology
- Handles special experiments (El Niño, La Niña, RCP8.5) by adding/subtracting anomaly fields
- Sets `cap_surf` per grid cell (land vs. ocean heat capacity, with ice-aware ocean capacity)
- Extracts initial conditions from the last time step of climatology: `Ts`, `Ta` (= `Ts`), `To`, `q`
- Determines control CO₂ level based on experiment type (default 340 ppm, IPCC scenarios 280 ppm, or 0 for deconstruction)

Returns `(Ts_ini, Ta_ini, To_ini, q_ini, CO2_ctrl)`.
"""

# ╔═╡ d404043f-8080-4262-9ab6-d9bb13eee504
function init_model!(cfg::PhysicsConfig)

	# ── Hydrology Parameter Initialization ────────────
	set_hydrology_parameters!(cfg)
	
	# ── dTrad: offset between T_atm and radiation temperature ────
	@. dTrad = -0.16 * Tclim - 5.0

	# ── z_ocean: 3× maximum mixed-layer depth over the year ──────
	z_ocean .= 3.0 .* dropdims(maximum(mldclim; dims=3); dims=3) 

	# ── Sensitivity experiment overrides ─────────────────────────
	if !cfg.log_topo_drsp
		@. z_topo = min(z_topo, 1.0)       # constant topography
	end

	# Apply cloud deconstruction switch
	if !cfg.log_clouds_dmc
		cldclim .= 0.0  # zero cloud climatology
	end
	
	# Apply flux correction conditional zeroing (MSCM feature)
	if !cfg.log_topo_drsp && !cfg.log_qflux_dmc
		TF_correct .= 0.0
		qF_correct .= 0.0
		ToF_correct .= 0.0
	end

	# Climatology modifications
	if !cfg.log_hydro_dmc
		qclim .= 0.0  # zero out humidity climatology
	end
	
	if !cfg.log_vapor_dmc
		qclim .= 0.0052          # constant water vapor
	end
	
	if !cfg.log_ocean_dmc
		mldclim .= d_ocean       # no deep ocean
	end

	# ── Experiment Handler ───────────────────────────────────
	# Apply advanced experiment forcing
	if cfg.experiment == :rcp85
		@info "Applying CMIP5 RCP8.5 climate change forcing"
		Tclim      .+= Tclim_anom_cc
		uclim      .+= uclim_anom_cc
		vclim      .+= vclim_anom_cc
		omegaclim  .+= omegaclim_anom_cc
		wsclim     .+= wsclim_anom_cc
	elseif cfg.experiment == :elnino
        @info "Applying ERA-Interim El Niño forcing"
        Tclim .+= Tclim_anom_enso
        uclim .+= uclim_anom_enso
        vclim .+= vclim_anom_enso
        omegaclim .+= omegaclim_anom_enso
        wsclim .+= wsclim_anom_enso
    elseif cfg.experiment == :lanina
        @info "Applying ERA-Interim La Niña forcing"
        Tclim .-= Tclim_anom_enso
        uclim .-= uclim_anom_enso
        vclim .-= vclim_anom_enso
        omegaclim .-= omegaclim_anom_enso
        wsclim .-= wsclim_anom_enso
	end
	
	# ── Topography pressure weights ─────
	@. wz_air = exp(-z_topo / z_air)
	@. wz_vapor = exp(-z_topo / z_vapor)

	# ── Surface heat capacity ────────────────────────────────────
	@inbounds for j in 1:ydim
    	for i in 1:xdim
        	if z_topo[i,j] > 0.0
            	cap_surf[i,j] = cap_land
        	else
            	cap_surf[i,j] = cfg.log_ocean_dmc ? cap_ocean * mldclim[i,j,1] : cap_land
        	end
    	end
	end

	# ── Initial conditions from last time step of climatology ────
	Ts_ini  = Tclim[:, :, nstep_yr]   |> copy   # surface temperature
	Ta_ini  = copy(Ts_ini)                      # air temperature = Tsurf
	To_ini  = Toclim[:, :, nstep_yr]  |> copy   # deep ocean temperature
	q_ini   = qclim[:, :, nstep_yr]   |> copy   # atmospheric water vapor

	# ── Control CO₂ level ───────────────────────────────────────
	CO2_ctrl = cfg.co2_concentration

	if cfg.experiment == :a1b_scenario  
        CO2_ctrl = 298.0  # A1B scenario baseline
    elseif cfg.experiment in (:a1b_enhanced, :rcp26, :rcp45, :rcp60, :rcp85, :custom_co2)  
        CO2_ctrl = 280.0  # IPCC scenarios baseline
    end
    
    if !cfg.log_co2_dmc
        CO2_ctrl = 0.0  # Zero CO2 for deconstruction experiments
    end

	return (Ts_ini = Ts_ini, Ta_ini = Ta_ini, To_ini = To_ini,
	        q_ini  = q_ini,  CO2_ctrl = CO2_ctrl)
end

# ╔═╡ fd193b13-336f-42f8-8f93-3ca3bea2e620
md"""
## 🔄 `time_loop` - Single Time-Step Integration

Advances the model by one main time step (12 hours):

1. Looks up calendar day and time-step index from iteration counter
2. Computes all physics tendencies via `tendencies!`
3. Updates `Ts` (surface) and `Ta` (air) from energy budgets, with clamps at 233.15 K
4. Updates `To` (deep ocean) with mixing tendencies and flux correction
5. Updates `q` (humidity) with evaporation, rain, circulation, and flux correction (clamped)
6. Adjusts ocean heat capacity via `seaice!` based on current ice cover
7. Converts moisture tendencies to mm/day for output
8. Calls `output!` (monthly accumulation) and `diagnostics!` (annual accumulation)

Returns updated `(mon, irec)`.
"""

# ╔═╡ 241496f6-580b-4e81-acbb-3295696fe40e
md"""
## ⚔️ `tendencies` - Physics Aggregation

Calls all physical process modules in sequence for one time step:

| Step | Process | Key Output |
|:-----|:--------|:-----------|
| 1 | Shortwave radiation | `albedo`, `SW`, `ice_cover` |
| 2 | Longwave radiation | `LW_surf`, `LW_down`, `LW_up`, `em` |
| 3 | Sensible heat flux | `Q_sens` (zero if atmosphere disabled) |
| 4 | Hydrological cycle | `Q_lat`, `dq_eva`, `dq_rain` |
| 5 | Circulation (heat) | `dTa_crcl` |
| 6 | Circulation (moisture) | `dq_crcl` |
| 7 | Deep ocean coupling | `dT_ocean`, `dTo` |

Processes respect MSCM isolation switches - disabled processes return zero tendencies.

Returns a NamedTuple with all tendencies and diagnostics.
"""

# ╔═╡ a4ab662c-5876-4861-aecf-2f695d15e15d
md"""
## ☀️ `SWradiation` - Shortwave Radiation

Computes ice cover, surface albedo, and net shortwave flux.

- **Ice fraction**: linear ramp between ice-free and fully frozen thresholds, separate for land (`Tl_ice1` → `Tl_ice2`) and ocean (`To_ice1` → `To_ice2`); glaciers always fully iced
- **Atmospheric albedo**: `a_cloud × cloud_fraction`
- **Surface albedo**: ramps from `a_no_ice` to `a_no_ice + da_ice` based on ice fraction; fixed to `a_no_ice` if ice feedback disabled (`log_ice = false`)
- **Combined albedo**: `α = α_surf + α_atmos − α_surf × α_atmos`
- **Net SW**: `SW = solar_flux × (1 − α)` with optional solar multiplier

Returns `(SW, albedo, ice_cover)`.
"""

# ╔═╡ 1df2b91b-be14-427a-87b3-95cdef26ce00
function SWradiation!(Ts, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    # Reuse workspace buffers
    ice_cover = ws.ice_cover_buf # output: ice fraction
    a_surf    = ws.a_surf_buf    # surface albedo
    albedo    = ws.albedo_buf    # output: combined albedo (surface + atmosphere) [hernoemen]
    a_atmos   = ws.a_atmos_buf   # atmospheric albedo
    sw        = ws.sw_buf        # output: net shortwave flux
    
    # 1. Ice cover fraction – branch‑free with ifelse, vectorized
    @turbo for i in 1:xdim, j in 1:ydim
        T = Ts[i, j]
        land_expr = ifelse(T <= Tl_ice1, 1.0,
                    ifelse(T < Tl_ice2,
                           1.0 - (T - Tl_ice1) * inv_Tl_ice_range,
                           0.0))
        ocean_expr = ifelse(T <= To_ice1, 1.0,
                    ifelse(T < To_ice2,
                            1.0 - (T - To_ice1) * inv_To_ice_range,
                            0.0))
        ice_cover[i, j] = ifelse(z_topo[i, j] >= 0.0, land_expr, ocean_expr)
    end

    # 2. Atmospheric albedo – simple multiplication, use broadcasting
    cld = @view cldclim[:, :, timestate.ityr]
    @. a_atmos = cld * a_cloud

    # 3. Surface albedo – conditional logic, @turbo beneficial
    if cfg.log_ice
        @turbo for i in 1:xdim, j in 1:ydim
            T = Ts[i, j]
            # Land albedo expression
            land_alb = ifelse(T <= Tl_ice1, a_no_ice + da_ice,
                       ifelse(T >= Tl_ice2, a_no_ice,
                                a_no_ice + da_ice * (1.0 - (T - Tl_ice1) * inv_Tl_ice_range)))
            # Ocean albedo expression
            ocean_alb = ifelse(T <= To_ice1, a_no_ice + da_ice,
                            ifelse(T >= To_ice2, a_no_ice,
                                a_no_ice + da_ice * (1.0 - (T - To_ice1) * inv_To_ice_range)))
            # Choose based on topography
            a_surf[i, j] = ifelse(z_topo[i, j] >= 0.0, land_alb, ocean_alb)
            # Glacier override: if glacier mask > 0.5, set to ice albedo
            a_surf[i, j] = ifelse(glacier[i, j] > 0.5, a_no_ice + da_ice, a_surf[i, j])
        end
    else
        @. a_surf = a_no_ice
    end

    # 4. Combined albedo
    @. albedo = a_surf + a_atmos - a_surf * a_atmos
	
    # 5. Shortwave flux
    multiplier = sw_solar_forcing_state[] * 0.01 * S0_var
    for j in 1:ydim
        sf = sw_solar[j, timestate.ityr] * multiplier
        @. sw[:, j] = sf * (1.0 - albedo[:, j])
    end
    
    return (SW = sw, albedo = albedo, ice_cover = ice_cover)
end

# ╔═╡ 2a859d43-5768-4b98-88fa-978b12203513
md"""
## 🌙 `LWradiation` - Long-Wave Radiation

Computes atmospheric emissivity and longwave fluxes.

- **Effective columns**: CO₂ and water vapor scaled by pressure weight `wz_air`; CO₂ also scaled by spatial mask `co2_part`
- **Emissivity**: log-regression with 10 coefficients combining CO₂, vapor, and cloud effects
- **Fluxes**:  
  `LW_surf = −σ·Ts⁴`  
  `LW_down = −ε·σ·(Ta + dTrad)⁴`  
  `LW_up = LW_down`
- Zeroed if atmosphere disabled (`log_atmos_dmc = false`)

Returns `(LW_surf, LW_up, LW_down, em)`.
"""

# ╔═╡ c0c40037-4169-4d38-bebe-2086cebc24f2
function LWradiation!(Ts, Ta, q, CO2, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
	# Extract workspace buffers
	e_co2   = ws.e_co2_buf    # CO₂ [ppm scaled by pressure]
	e_vapor = ws.e_vapor_buf  # water vapour [kg/m²]
    em      = ws.em_buf       # emissivity ε_atmos
	LW_surf = ws.LW_surf_buf  # surface long-wave flux [W/m²]
	LW_down = ws.LW_down_buf  # downward long-wave flux [W/m²]
	LW_up   = ws.LW_up_buf    # upward long-wave flux [W/m²]

    # Current cloud cover (climatology, 3D array)
	e_cloud = @view cldclim[:, :, timestate.ityr]  
	
	## 1. Effective columns (topography scaling via wz_air, precomputed)
	@. e_vapor = wz_air * r_qviwv * q         
	@. e_co2   = wz_air * CO2 * co2_part

	# ── Emissivity (log-regression with 10 parameters) ─────────
	@. em = p_emi[4] * log(p_emi[1] * e_co2 + p_emi[2] * e_vapor + p_emi[3]) +
		       p_emi[7] +
		       p_emi[5] * log(p_emi[1] * e_co2 + p_emi[3]) +
		       p_emi[6] * log(p_emi[2] * e_vapor + p_emi[3])

	# Cloud adjustment
	@. em = (p_emi[8] - e_cloud) / p_emi[9] * (em - p_emi[10]) + p_emi[10]

	# 4. Radiation temperature (precomputed offset dTrad = -0.16*Tclim - 5 K)
    dTr = @view dTrad[:, :, timestate.ityr]
	@. LW_surf = -σ * Ts^4
	@. LW_down = -em * σ * (Ta + dTr)^4
	LW_up .= LW_down

	if !cfg.log_atmos_dmc
		LW_down .= 0.0
        LW_up   .= 0.0
	end

	return (LW_surf = LW_surf, LW_up = LW_up, LW_down = LW_down, em = em)
end

# ╔═╡ 0febe534-237e-4921-b39c-3828dbae9d19
md"""
## 💧 `hydro` - Hydrological Cycle

Computes latent heat flux and water vapor tendencies (evaporation, precipitation).

- **Saturation humidity**: Clausius-Clapeyron approximation with pressure scaling `wz_air`
- **Evaporation**: bulk aerodynamic formula with gustiness; three modes selectable via `log_eva`:
  - `-1` - original GREB (wind from `uclim`, `vclim`)
  - `0`  - skin temperature method with wind speed climatology
  - `1`  - enhanced
- **Precipitation**: `dq_rain = (c_q + c_rq·RH + c_omega·ω + c_omegastd·σ_ω) · cq_rain · q`
- **Vapor tendencies**: `dq_eva = −Q_lat / cq_latent / r_qviwv`, rain limited to `−0.9·q/Δt`
- Returns zero if hydrology or atmosphere disabled

Returns `(Q_lat, Q_lat_air, dq_eva, dq_rain)`.
"""

# ╔═╡ 606032a2-b2ca-4fd8-9930-afd83aecec7a
function hydro!(Ts, q, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    c_q = cfg.c_q
    c_rq = cfg.c_rq
    c_omega = cfg.c_omega
    c_omegastd = cfg.c_omegastd

    fill!(ws.Q_lat_buf, 0.0)
    fill!(ws.Q_lat_air_buf, 0.0)
    fill!(ws.dq_eva_buf, 0.0)
    fill!(ws.dq_rain_buf, 0.0)

    if !cfg.log_atmos_dmc || !cfg.log_hydro_dmc || !cfg.log_hydro_drsp
        return (Q_lat = ws.Q_lat_buf, Q_lat_air = ws.Q_lat_air_buf,
                dq_eva = ws.dq_eva_buf, dq_rain = ws.dq_rain_buf)
    end

    u = @view uclim[:, :, timestate.ityr]
    v = @view vclim[:, :, timestate.ityr]
    swet = @view swetclim[:, :, timestate.ityr]
    omega = @view omegaclim[:, :, timestate.ityr]
    omegastd = @view omegastdclim[:, :, timestate.ityr]

    const_factor1 = 3.75e-3
    const_factor2 = 17.08085
    const_factor3 = 234.175
    gust_land = 4.0
    gust_ocean = 9.0
    cE_land = 0.25 * ce
    cE_ocean = 0.58 * ce
    const_latent = cq_latent * ρ_air * ce

    # Saturation humidity
    @turbo for j in 1:ydim
        for i in 1:xdim
            T = Ts[i,j] - 273.15
            ws.qs[i,j] = const_factor1 * exp(const_factor2 * T / (T + const_factor3)) * wz_air[i,j]
            ws.qs[i,j] = max(ws.qs[i,j], 1e-8)
        end
    end

    # Evaporation
    if cfg.log_eva == -1
        @turbo for j in 1:ydim
            for i in 1:xdim
                u_val = u[i,j]; v_val = v[i,j]
                wind = sqrt(u_val*u_val + v_val*v_val)
                wind = sqrt(wind*wind + ifelse(z_topo[i,j] > 0.0, gust_land, gust_ocean))
                ws.Q_lat_buf[i,j] = (q[i,j] - ws.qs[i,j]) * wind * const_latent * swet[i,j]
            end
        end
    elseif cfg.log_eva == 0
        ws_view = @view wsclim[:, :, timestate.ityr]
        @turbo for j in 1:ydim
            for i in 1:xdim
                ws.Tskin[i,j] = ifelse(z_topo[i,j] > 0.0, Ts[i,j] + 5.0, Ts[i,j] + 1.0)
                ws.Tskin[i,j] = ifelse(ws.Tskin[i,j] < 200.0, 200.0, ws.Tskin[i,j])
                T = ws.Tskin[i,j] - 273.15
                qs_val = const_factor1 * exp(const_factor2 * T / (T + const_factor3)) * wz_air[i,j]
                ws.qs[i,j] = qs_val
                
                ws.ws_base[i,j] = ws_view[i,j]
                gust = ifelse(z_topo[i,j] > 0.0, 132.25, 29.16)
                wind = sqrt(ws.ws_base[i,j]*ws.ws_base[i,j] + gust)
                
                ws.cE[i,j] = ifelse(z_topo[i,j] > 0.0, cE_land, cE_ocean)
                ws.Q_lat_buf[i,j] = ws.cE[i,j] * wind * ρ_air * cq_latent * (q[i,j] - qs_val) * swet[i,j]
            end
        end
    else
        @turbo for j in 1:ydim
            for i in 1:xdim
                u_val = u[i,j]; v_val = v[i,j]
                wind = sqrt(u_val*u_val + v_val*v_val)
                wind = sqrt(wind*wind + ifelse(z_topo[i,j] > 0.0, gust_land, gust_ocean))
                ws.Q_lat_buf[i,j] = (q[i,j] - ws.qs[i,j]) * wind * const_latent * swet[i,j]
            end
        end
    end

    # Precipitation - use ws.rq buffer
    @turbo for j in 1:ydim
        for i in 1:xdim
            ws.rq[i,j] = q[i,j] / max(ws.qs[i,j], 1e-8)
            ws.dq_rain_buf[i,j] = (c_q + c_rq * ws.rq[i,j] + c_omega * omega[i,j] + c_omegastd * omegastd[i,j]) * cq_rain * q[i,j]
        end
    end

    # Apply rain limit
    if cfg.log_rain == 1
        limit_val = -0.0015 / (wz_vapor[1,1] * r_qviwv * 86400.0)
        @turbo for j in 1:ydim
            for i in 1:xdim
                ws.dq_rain_buf[i,j] = ifelse(ws.dq_rain_buf[i,j] >= limit_val, limit_val, ws.dq_rain_buf[i,j])
            end
        end
    end

    # Water vapor tendencies
    @turbo for j in 1:ydim
        for i in 1:xdim
            ws.dq_eva_buf[i,j] = -ws.Q_lat_buf[i,j] / cq_latent / r_qviwv
            min_dq = -0.9 * q[i,j] / Δt
            ws.dq_rain_buf[i,j] = ifelse(ws.dq_rain_buf[i,j] < min_dq, min_dq, ws.dq_rain_buf[i,j])
            ws.Q_lat_air_buf[i,j] = -ws.dq_rain_buf[i,j] * cq_latent * r_qviwv
        end
    end

    return (Q_lat = ws.Q_lat_buf,
            Q_lat_air = ws.Q_lat_air_buf,
            dq_eva = ws.dq_eva_buf,
            dq_rain = ws.dq_rain_buf)
end

# ╔═╡ a1c04c52-f7c9-430a-8791-7fea15650b2c
md"""
## 🌫️ `convergence` - Moisture Flux Convergence

Computes moisture tendency from vertical velocity (omega) convergence:

`dX_conv = −T1 · ω · (Δt_crcl / z_vapor · 2.5 / (ρ_air·g))`

Implements Eq. 18 from Stassen et al. (2019). Called during circulation sub-stepping for water vapor only.
"""

# ╔═╡ c6a0d656-8289-4f74-b8d8-f94c236e541d
function convergence!(T1, omegaclim, timestate, ws::CirculationWorkspace)
	"""Calculate moisture flux convergence using omega vertical velocity.
	
	Implements Eq. 18 from Stassen et al 2019.
	
	Args:
		T1: Input field (typically specific humidity) [kg/kg]
		omegaclim: Vertical velocity climatology [Pa/s] 
		timestate: Time state object
		ws: Workspace for temporary arrays
	"""
	omega = @view omegaclim[:, :, timestate.ityr]

	@. ws.dX_conv = -T1 * omega * const_factor
	return nothing
end

# ╔═╡ 34fb479a-9834-4354-8f36-2680bedec798
md"""
## ❄️ `seaice` - Sea-Ice Heat Capacity

Updates `cap_surf` in-place for ocean points based on ice cover fraction.

- Ocean heat capacity blends linearly between `cap_land` (fully frozen, `Ts ≤ To_ice1`) and `cap_ocean × mld` (ice-free, `Ts ≥ To_ice2`)
- Land and glacier points keep `cap_land`
- Skipped entirely if ocean disabled (`log_ocean_dmc = false`)
- If ice feedback disabled (`log_ice = false`): land → `cap_land`, ocean → `cap_ocean × mld`
"""

# ╔═╡ 9bff59c5-4631-4091-8230-989a835788e5
function seaice!(Ts0, timestate, cfg::PhysicsConfig)
	mld = @view mldclim[:, :, timestate.ityr]

	if !cfg.log_ocean_dmc 
		return     # No ice feedback: skip sea ice calculation
	end
	
	# Compute ice‑dependent heat capacity for ocean points
    @turbo for i in 1:xdim, j in 1:ydim
        is_ocean = z_topo[i, j] < 0.0
        T = Ts0[i, j]
        mld_val = mld[i, j]
        cap_open = cap_ocean * mld_val

        # Ice fraction (0 = no ice, 1 = full ice)
        ice_frac = ifelse(T <= To_ice1, 1.0,
                   ifelse(T >= To_ice2, 0.0,
                          1.0 - (T - To_ice1) * inv_To_ice_range))

        # Blend between land (ice) and open ocean capacities
        cap_with_ice = cap_land * ice_frac + cap_open * (1.0 - ice_frac)

        # Apply only to ocean points; keep land points unchanged
        cap_surf[i, j] = ifelse(is_ocean, cap_with_ice, cap_surf[i, j])
    end

    # Override for experiments without ice‑albedo feedback
    if !cfg.log_ice
        @turbo for i in 1:xdim, j in 1:ydim
            cap_surf[i, j] = ifelse(z_topo[i, j] > 0.0, cap_land, cap_ocean * mld[i, j])
        end
        return
    end

    # Glacier override: ice sheets have land heat capacity
    @. cap_surf = ifelse(glacier > 0.5, cap_land, cap_surf)
end

# ╔═╡ ae38f814-fe9c-443f-ae3b-42fa7a7d199a
md"""
## 🌊 `deep_ocean` - Deep Ocean Coupling

Computes heat exchange between mixed layer and deep ocean.

- **Entrainment** (`mld` shallowing): deep ocean warms, surface cools — scaled by `c_effmix`
- **Detrainment** (`mld` deepening): surface cools, deep ocean warms
- **Turbulent mixing**: proportional to `co_turb × (max(Ts, To_ice2) − To)`
- Active only for ocean points with `Ts ≥ To_ice2`; returns zeros if ocean disabled

Returns `(dT_ocean, dTo)`.
"""

# ╔═╡ 625089e2-ef77-4821-a6d6-d0a0f88207f2
function deep_ocean!(Ts, To, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
	# Use pre-allocated zero buffers
	dT_ocean = ws.dT_ocean_buf
    dTo      = ws.dTo_buf

	# no deep-ocean coupling
	if !cfg.log_ocean_dmc || !cfg.log_ocean_drsp
        fill!(dT_ocean, 0.0)
        fill!(dTo, 0.0)
        return (dT_ocean = dT_ocean, dTo = dTo)
    end
	
	# ── Change in mixed-layer depth ─────────────────────────
	mld_now  = @view mldclim[:, :, timestate.ityr]
    mld_prev = timestate.ityr > 1 ? @view(mldclim[:, :, timestate.ityr-1]) : @view(mldclim[:, :, nstep_yr])	
	
    # Zero buffers first (one pass, cheap)
    fill!(dT_ocean, 0.0)
    fill!(dTo, 0.0)

	# ── Entrainment & detrainment & turbulent mixing ──────
	@turbo for i in 1:xdim, j in 1:ydim
        active = (z_topo[i,j] < 0.0) & (Ts[i,j] >= To_ice2)
        h_now  = mld_now[i,j]
        h_prev = mld_prev[i,j]
        dh     = h_now - h_prev
        z_deep = z_ocean[i,j]
        z_rem  = z_deep - h_now

        # Entrainment/detrainment contributions (only when active)
        dTo_entr      = ifelse(active & (dh < 0.0), c_effmix * (-dh / z_rem) *
                                        (Ts[i,j] - To[i,j]), 0.0)
        dT_ocean_entr = ifelse(active & (dh > 0.0), c_effmix * (dh / h_now) * 
                                        (To[i,j] - Ts[i,j]), 0.0)

        # Turbulent mixing (only when active)
        Tx = ifelse(Ts[i,j] > To_ice2, Ts[i,j], To_ice2)
        dTo_turb      = ifelse(active, turb_coeff * (Tx - To[i,j]) / z_rem, 0.0)
        dT_ocean_turb = ifelse(active, turb_coeff * (To[i,j] - Tx) / h_now, 0.0)

        # Combine (buffer was zeroed before loop)
        dTo[i,j]      = dTo_entr + dTo_turb
        dT_ocean[i,j] = dT_ocean_entr + dT_ocean_turb
    end
    return (dT_ocean = dT_ocean, dTo = dTo)
end

# ╔═╡ a4cd8d40-eedd-4a4b-ac17-51e68334328c
md"""
## 🌀 `circulation` - Atmospheric Circulation

Runs sub‑stepped diffusion + advection at the finer time step `Δt_crcl`.

- Loops `ntime` sub‑steps, each calling selected processes: horizontal/vertical diffusion, horizontal/vertical advection, moisture convergence
- Process selection controlled by switches: `log_hdif`, `log_vdif`, `log_hadv`, `log_vadv`, `log_conv`
- Early return with zeros if atmosphere or circulation disabled
- Returns `dX_crcl = X_final − X_in`
"""

# ╔═╡ cc7db228-0a06-4812-8b32-6541988b2115
md"""
## 🌀 `diffusion` - 3rd-Order Diffusion

Topography‑weighted horizontal diffusion on a periodic lat‑lon grid.

- **Meridional**: 2nd‑order centred differences; zero‑flux at poles
- **Zonal (mid‑latitudes)**: 3rd‑order stencil with ±1, ±2, ±3 neighbours (weights 10/4/1, normalised by 20)
- **Zonal (polar, `dxlat ≤ 250 km`)**: sub‑time‑stepping for stability, with clamp `dX ≥ −0.9·X`
- Scaled by pressure weight `exp(−z_topo / h_scl)`

Writes directly to `ws.dX_diff`.
"""

# ╔═╡ ba96178d-77d4-4f26-a94f-5ad43c5242db
function diffusion!(T1, h_scl, ws::CirculationWorkspace, timestate)
    # Zero output buffer (we will accumulate into it)
    fill!(ws.dX_diff, 0.0)

    # Topographic scaling (choose based on scale height)
    wz = if h_scl == z_air
        wz_air
    elseif h_scl == z_vapor
        wz_vapor
    else
        error("Invalid h_scl = $h_scl (must be z_air or z_vapor)")
    end

    # Precomputed geometry/coefficients
    dxlat = dxlat_grid
    ccy   = ccy_diff
    ccx   = ccx_diff

    # Pre‑cached neighbour indices
    jm1, jp1 = lon_jm1, lon_jp1
    jm2, jp2 = lon_jm2, lon_jp2
    jm3, jp3 = lon_jm3, lon_jp3

    # ----- Precompute k‑independent terms for the poles -----
    # For k == 1 (North Pole)
    @. ws.term_north = ccy * wz[:, 2] * (T1[:, 2] - T1[:, 1])
    # For k == ydim (South Pole)
    @. ws.term_south = ccy * wz[:, ydim-1] * (T1[:, ydim-1] - T1[:, ydim])

    for k in 1:ydim
        # ----- Meridional diffusion -----
        if k == 1
            @turbo for i in 1:xdim
            ws.dX_diff[i, k] += wz[i, k] * ws.term_north[i]
        end
        elseif k == ydim
            @turbo for i in 1:xdim
            ws.dX_diff[i, k] += wz[i, k] * ws.term_south[i]
        end
        else
            # Mid‑latitudes: no precomputation possible (depends on k‑1, k+1)
            @. ws.dX_diff[:, k] += wz[:, k] * ccy * (
                wz[:, k-1] * (T1[:, k-1] - T1[:, k]) +
                wz[:, k+1] * (T1[:, k+1] - T1[:, k])
            )
        end

        # ----- Zonal diffusion -----
        if dxlat[k] > 2.5e5   # mid‑latitudes, normal time step
            # Complex stencil – keep @turbo
            @turbo for j in 1:xdim
                jm1v = jm1[j]; jp1v = jp1[j]
                jm2v = jm2[j]; jp2v = jp2[j]
                jm3v = jm3[j]; jp3v = jp3[j]

                dTx = ccx[k] * 0.05 * (
                    10.0 * (wz[jm1v, k] * (T1[jm1v, k] - T1[j, k]) +
                          wz[jp1v, k] * (T1[jp1v, k] - T1[j, k])) +
                     4.0 * (wz[jm2v, k] * (T1[jm2v, k] - T1[jm1v, k]) +
                          wz[jm1v, k] * (T1[j, k] - T1[jm1v, k])) +
                     4.0 * (wz[jp1v, k] * (T1[j, k] - T1[jp1v, k]) +
                          wz[jp2v, k] * (T1[jp2v, k] - T1[jp1v, k])) +
                     1.0 * (wz[jm3v, k] * (T1[jm3v, k] - T1[jm2v, k]) +
                          wz[jm2v, k] * (T1[jm1v, k] - T1[jm2v, k])) +
                     1.0 * (wz[jp2v, k] * (T1[jp1v, k] - T1[jp2v, k]) +
                          wz[jp3v, k] * (T1[jp3v, k] - T1[jp2v, k]))
                )
                ws.dX_diff[j, k] += wz[j, k] * dTx
            end
        else   # polar regions – sub‑timestepping (complex, keep @turbo)
            # Number of sub‑steps for stability
            dd = max(1, round(Int, Δt_crcl / (dxlat[k]^2 / κ)))
            dtdff2 = Δt_crcl / dd
            time2 = max(1, round(Int, Δt_crcl / dtdff2))
            ccx2 = κ * dtdff2 / dxlat[k]^2

            # Copy current row into temporary buffer
            ws.T1h .= @view T1[:, k]

            for _ in 1:time2
                @turbo for j in 1:xdim
                    jm1v = jm1[j]; jp1v = jp1[j]
                    jm2v = jm2[j]; jp2v = jp2[j]
                    jm3v = jm3[j]; jp3v = jp3[j]

                    dq = ccx2 * 0.05 * (
                        10.0 * (wz[jm1v, k] * (ws.T1h[jm1v] - ws.T1h[j]) +
                              wz[jp1v, k] * (ws.T1h[jp1v] - ws.T1h[j])) +
                         4.0 * (wz[jm2v, k] * (ws.T1h[jm2v] - ws.T1h[jm1v]) +
                              wz[jm1v, k] * (ws.T1h[j] - ws.T1h[jm1v])) +
                         4.0 * (wz[jp1v, k] * (ws.T1h[j] - ws.T1h[jp1v]) +
                              wz[jp2v, k] * (ws.T1h[jp2v] - ws.T1h[jp1v])) +
                         1.0 * (wz[jm3v, k] * (ws.T1h[jm3v] - ws.T1h[jm2v]) +
                              wz[jm2v, k] * (ws.T1h[jm1v] - ws.T1h[jm2v])) +
                         1.0 * (wz[jp2v, k] * (ws.T1h[jp1v] - ws.T1h[jp2v]) +
                              wz[jp3v, k] * (ws.T1h[jp3v] - ws.T1h[jp2v]))
                    )
                    # Stability clamp
                    dq = ifelse(dq <= -ws.T1h[j], -0.9 * ws.T1h[j], dq)
                    ws.T1h[j] += dq
                end
            end

            # Add total change (scaled by outer wz) to output buffer – use broadcast
            @. ws.dX_diff[:, k] += wz[:, k] * (ws.T1h - T1[:, k])
        end
    end

    return nothing
end

# ╔═╡ 99f0e9ad-438f-4fa5-be9d-36d0fa78d89c
md"""
## ➡️ `advection` - Upwind Advection

Upwind transport driven by climatological wind fields.

- **Meridional**: 2nd‑order upwind with ±1/±2 neighbours; one‑sided at pole boundaries
- **Zonal (mid‑latitudes)**: 2nd‑order upwind (±1/±2)
- **Zonal (polar)**: sub‑time‑stepping with 3rd‑order stencil (weights 10/4/1, normalised by 20), clamp `dX ≥ −0.9·X`
- Returns zeros if advection disabled for the transported quantity (`log_hadv` for heat, `log_vadv` for vapor)

Writes directly to `ws.dX_adv`.
"""

# ╔═╡ 2bab06b9-ca98-4142-99cb-d2ad4f1cde93
function advection!(T1, h_scl, ws::CirculationWorkspace, timestate, cfg::PhysicsConfig)
    # Disable advection for water vapour or heat according to switches
    if (h_scl == z_vapor && !cfg.log_vadv) || (h_scl == z_air && !cfg.log_hadv)
        fill!(ws.dX_adv, 0.0)
        return nothing
    end
    
    # Pre-zero the output buffer (we will accumulate into it)
    fill!(ws.dX_adv, 0.0)

    # Extract 2D views for current time step
    vclim_p_t = @view vclim_p[:, :, timestate.ityr]
    vclim_m_t = @view vclim_m[:, :, timestate.ityr]
    uclim_p_t = @view uclim_p[:, :, timestate.ityr]
    uclim_m_t = @view uclim_m[:, :, timestate.ityr]

    # Topographic scaling (choose based on scale height)
    wz = if h_scl == z_air
        wz_air
    elseif h_scl == z_vapor
        wz_vapor
    else
        error("Invalid h_scl = $h_scl (must be z_air or z_vapor)")
    end

    # Precomputed constants
    dxlat = dxlat_grid
    ccy = ccy_adv
    ccx = ccx_adv
    is_polar = IS_POLAR

    @inbounds for k in 1:ydim
        # ----- Meridional (v) advection -----
        if k == 1          # North Pole
            @turbo for j in 1:xdim
                v_p = vclim_p_t[j, k]
                ws.dX_adv[j, k] += ccy * v_p * (
                    wz[j, 2] * (T1[j, 1] - T1[j, 2]) +
                    wz[j, 3] * (T1[j, 1] - T1[j, 3])
                ) / 3.0
            end
        elseif k == 2
            @turbo for j in 1:xdim
                v_m = vclim_m_t[j, k]
                v_p = vclim_p_t[j, k]
                ws.dX_adv[j, k] += ccy * (
                    -v_m * wz[j, 1] * (T1[j, 2] - T1[j, 1]) +
                    v_p * (wz[j, 3] * (T1[j, 2] - T1[j, 3]) +
                           wz[j, 4] * (T1[j, 2] - T1[j, 4])) / 3.0
                )
            end
        elseif k >= 3 && k <= ydim-2
            @turbo for j in 1:xdim
                km1, km2 = k-1, k-2
                kp1, kp2 = k+1, k+2
                v_m = vclim_m_t[j, k]
                v_p = vclim_p_t[j, k]
                ws.dX_adv[j, k] += ccy * (
                    -v_m * (wz[j, km1] * (T1[j, k] - T1[j, km1]) +
                            wz[j, km2] * (T1[j, k] - T1[j, km2])) +
                    v_p * (wz[j, kp1] * (T1[j, k] - T1[j, kp1]) +
                           wz[j, kp2] * (T1[j, k] - T1[j, kp2]))
                ) / 3.0
            end
        elseif k == ydim-1
            @turbo for j in 1:xdim
                km1, km2 = k-1, k-2
                kp1 = k+1
                v_m = vclim_m_t[j, k]
                v_p = vclim_p_t[j, k]
                ws.dX_adv[j, k] += ccy * (
                    -v_m * (wz[j, km1] * (T1[j, k] - T1[j, km1]) +
                            wz[j, km2] * (T1[j, k] - T1[j, km2])) / 3.0 +
                    v_p * wz[j, kp1] * (T1[j, k] - T1[j, kp1])
                )
            end
        else               # k == ydim (South Pole)
            @turbo for j in 1:xdim
                km1, km2 = k-1, k-2
                v_m = vclim_m_t[j, k]
                ws.dX_adv[j, k] += ccy * (
                    -v_m * (wz[j, km1] * (T1[j, k] - T1[j, km1]) +
                            wz[j, km2] * (T1[j, k] - T1[j, km2]))
                ) / 3.0
            end
        end

        # ----- Zonal (u) advection -----
        if !is_polar[k]   # mid‑latitudes, normal timestep
            @turbo for j in 1:xdim
                jm1, jp1 = lon_jm1[j], lon_jp1[j]
                jm2, jp2 = lon_jm2[j], lon_jp2[j]
                u_m = uclim_m_t[j, k]
                u_p = uclim_p_t[j, k]
                ws.dX_adv[j, k] += ccx[k] * (
                    -u_m * (wz[jm1, k] * (T1[j, k] - T1[jm1, k]) +
                            wz[jm2, k] * (T1[j, k] - T1[jm2, k])) +
                    u_p * (wz[jp1, k] * (T1[j, k] - T1[jp1, k]) +
                           wz[jp2, k] * (T1[j, k] - T1[jp2, k]))
                ) / 3.0
            end
        else               # polar regions – sub‑timestepping
            # Number of sub‑steps (CFL stability)
            dd = max(1, round(Int, Δt_crcl / (dxlat[k] / 10.0)))
            dtdff2 = Δt_crcl / dd
            time2 = max(1, round(Int, Δt_crcl / dtdff2))
            ccx2 = dtdff2 / dxlat[k] / 2.0

            # Copy current row into temporary buffer
            ws.T1h .= @view T1[:,k]

            for _ in 1:time2
                # One fused loop: compute increment and update in place
                @turbo for j in 1:xdim
                    jm1, jp1 = lon_jm1[j], lon_jp1[j]
                    jm2, jp2 = lon_jm2[j], lon_jp2[j]
                    jm3, jp3 = lon_jm3[j], lon_jp3[j]   # precomputed!
                    u_m = uclim_m_t[j, k]
                    u_p = uclim_p_t[j, k]

                    dq = ccx2 * (
                        -u_m * (10.0 * wz[jm1, k] * (ws.T1h[j] - ws.T1h[jm1]) +
                                 4.0 * wz[jm2, k] * (ws.T1h[jm1] - ws.T1h[jm2]) +
                                 1.0 * wz[jm3, k] * (ws.T1h[jm2] - ws.T1h[jm3])) +
                         u_p * (10.0 * wz[jp1, k] * (ws.T1h[j] - ws.T1h[jp1]) +
                                 4.0 * wz[jp2, k] * (ws.T1h[jp1] - ws.T1h[jp2]) +
                                 1.0 * wz[jp3, k] * (ws.T1h[jp2] - ws.T1h[jp3]))
                    ) / 20.0

                    # Stability clamp (avoid negative water vapour)
                    dq = ifelse(dq <= -ws.T1h[j], -0.9 * ws.T1h[j], dq)
                    ws.T1h[j] += dq
                end
            end

            # Add total change to the output buffer
            @. ws.dX_adv[:, k] += ws.T1h - T1[:, k]
        end
    end

    return nothing
end

# ╔═╡ db75ea52-9387-4c15-bde5-61777ac9b570
function circulation!(X_in, h_scl, dX_out, ws::CirculationWorkspace, timestate, cfg::PhysicsConfig)
    # Early exit if atmospheric processes disabled
    if (!cfg.log_atmos_dmc || !cfg.log_crcl_dmc || !cfg.log_crcl_drsp)
        fill!(dX_out, 0.0)
        return nothing
    end

    # Precompute flags (hoist conditionals)
    do_diff_v = cfg.log_vdif == 1 && h_scl == z_vapor
    do_diff_h = cfg.log_hdif == 1 && h_scl == z_air
    do_adv_v  = cfg.log_vadv == 1 && h_scl == z_vapor
    do_adv_h  = cfg.log_hadv == 1 && h_scl == z_air
    do_conv   = cfg.log_conv == 0 && h_scl == z_vapor

    copyto!(ws.X_work, X_in)

    for _tt in 1:ntime
        do_diff_v && diffusion!(ws.X_work, h_scl, ws, timestate)
        do_diff_h && diffusion!(ws.X_work, h_scl, ws, timestate)
        do_adv_v  && advection!(ws.X_work, h_scl, ws, timestate, cfg)
        do_adv_h  && advection!(ws.X_work, h_scl, ws, timestate, cfg)
        do_conv   && convergence!(ws.X_work, omegaclim, timestate, ws)

        @. ws.X_work += ws.dX_diff + ws.dX_adv + ws.dX_conv
    end

    # Final difference
    @. dX_out = ws.X_work - X_in

    return nothing
end

# ╔═╡ e493fae7-239a-494c-9a59-728446d70f7a
function tendencies!(CO2, Ts, Ta, To, q, ws::CirculationWorkspace,
					 timestate, cfg::PhysicsConfig)
	
	# Short-wave radiation → albedo, SW flux
	sw_out   = SWradiation!(Ts, timestate, cfg, ws)
	
	# Long-wave radiation → LW_surf, LW_up, LW_down, emissivity
	lw_out   = LWradiation!(Ts, Ta, q, CO2, timestate, cfg, ws)
	
	# Sensible heat flux
    Q_sens = ws.Q_sens_buf
    if cfg.log_atmos_dmc
        @. Q_sens = ct_sens * (Ta - Ts)
    else
        fill!(Q_sens, 0.0)
    end
	
	# Hydrological cycle → latent heat + evaporation/rain tendencies
	hy_out   = hydro!(Ts, q, timestate, cfg, ws)
	
	# Atmospheric circulation — temperature diffusion/advection 
	circulation!(Ta, z_air, ws.dTa_crcl, ws, timestate, cfg)
	
	# Atmospheric circulation — water-vapour diffusion/advection 
	circulation!(q,  z_vapor, ws.dq_crcl, ws, timestate, cfg)
	
	# Deep ocean coupling
	do_out   = deep_ocean!(Ts, To, timestate, cfg, ws)

	return (albedo     = sw_out.albedo,
	        SW         = sw_out.SW,
			ice_cover  = sw_out.ice_cover,
	        LW_surf    = lw_out.LW_surf,
	        Q_lat      = hy_out.Q_lat,
	        Q_sens     = Q_sens,
	        Q_lat_air  = hy_out.Q_lat_air,
	        dq_eva     = hy_out.dq_eva,
	        dq_rain    = hy_out.dq_rain,
	        dq_crcl    = ws.dq_crcl,
        	dTa_crcl   = ws.dTa_crcl,
	        dT_ocean   = do_out.dT_ocean,
	        dTo        = do_out.dTo,
	        LW_down = lw_out.LW_down,
	        LW_up   = lw_out.LW_up,
	        em         = lw_out.em)
end

# ╔═╡ 8fdcfefd-0490-434a-a2bb-d171557b6ae7
md"""
## 🔄 `qflux_correction` - Flux Correction Spin-Up

Runs a `time_flux`-year simulation to compute correction fields that keep the model close to observed climatology.

Each time step:
1. Computes all physics tendencies
2. Calculates uncorrected state for `Ts`, `Ta`, `To`, `q`
3. Derives corrections: `TF_correct = (Tclim − Ts_uncorrected) × cap_surf / Δt` (and similarly for ocean and humidity)
4. Applies corrections to get corrected state
5. Updates sea ice and accumulates diagnostics

Corrections are stored in `TF_correct`, `ToF_correct`, `qF_correct` (3D arrays: lon × lat × time‑step‑in‑year) and applied during the main simulation to nudge the model toward climatology.
"""

# ╔═╡ 7e514484-9d81-4c5c-83ab-191d68c13043
md"""
---
## 🧪 `forcing` — Experiment Forcing

Returns `(CO2, sw_solar_forcing)` for the current year/time step based on experiment type.

| Category | Experiments | Description |
|:---------|:------------|:------------|
| **CO₂ scaling** | `:co2_double`, `:co2_quadruple`, `:co2_10x`, `:co2_half`, `:co2_zero` | Fixed CO₂ multipliers |
| **Solar** | `:solar_plus27`, `:solar_cycle_11yr` | Solar constant +27 W/m² or 11-year cycle |
| **Time‑varying CO₂** | `:co2_sine_wave`, `:co2_step`, `:a1b_scenario`, `:a1b_enhanced` | CO₂ as function of year |
| **Paleoclimate** | `:paleo_231kyr`, `:paleo_solar_modern_co2`, `:modern_solar_paleo_co2` | Paleo solar + CO₂ combinations |
| **Orbital** | `:obliquity`, `:eccentricity`, `:earth_sun_distance` | Solar forcing loaded externally |
| **Regional CO₂** | `:regional_co2_nh/sh/tropics/extratropics/ocean/land_ice/winter/summer` | Spatial masks via `co2_part` |
| **Forced boundary** | `:elnino`, `:lanina`, `:rcp85` | CO₂ = 340 ppm; anomalies applied in `init_model!` |
"""

# ╔═╡ 1894ad94-cdf8-4e79-a0e5-b72088db31be
function forcing(it, year, cfg::PhysicsConfig, icmn_ctrl=zeros(xdim,ydim,12); nstep_yr=nstep_yr)
	# Default CO₂ concentration
	CO2 = cfg.co2_concentration
	sw_solar_forcing = 1.0
	
	# 📜 Legacy experiments ───────────
	if cfg.experiment == :constant_topo
        CO2 = 550.0  # 550 ppm CO₂ steady state
        
    elseif cfg.experiment == :a1b_scenario
        CO2_1950 = 310.0; CO2_2000 = 370.0; CO2_2050 = 520.0
        if year <= 2000
            CO2 = CO2_1950 + 60.0 / 50.0 * (year - 1950)
        elseif year <= 2050
            CO2 = CO2_2000 + 150.0 / 50.0 * (year - 2000)
        elseif year <= 2100
            CO2 = CO2_2050 + 180.0 / 50.0 * (year - 2050)
        end
	
	# 💨 CO₂ scaling experiments ──────────────────────────────────────────────
	elseif cfg.experiment == :co2_double
        CO2 = 680.0  # 2×CO₂ (already set, but explicit)
        
    elseif cfg.experiment == :co2_quadruple
        CO2 = 1360.0  # 4×CO₂
        
    elseif cfg.experiment == :co2_10x
        CO2 = 3400.0  # 10×CO₂
        
    elseif cfg.experiment == :co2_half
        CO2 = 170.0  # 0.5×CO₂
        
    elseif cfg.experiment == :co2_zero
        CO2 = 0.0  # 0×CO₂ (no greenhouse effect)
	
	# ☀️ Solar forcing experiments ───────────────────────────────────────────
	elseif cfg.experiment == :solar_plus27
        CO2 = 340.0
        sw_solar_forcing = (1365.0 + 27.0) / 1365.0
        
    elseif cfg.experiment == :solar_cycle_11yr
        CO2 = 340.0
        sw_solar_forcing = (1365.0 + 1.0 * sin(2π * year / 11.0)) / 1365.0
	
	# 📈 Enhanced A1B scenario ──────────────────────────────────────────────
	elseif cfg.experiment == :a1b_enhanced
        CO2_1950 = 310.0; CO2_2000 = 370.0; CO2_2050 = 520.0
        if year <= 2000
            CO2 = CO2_1950 + 60.0 / 50.0 * (year - 1950)
        elseif year <= 2050
            CO2 = CO2_2000 + 150.0 / 50.0 * (year - 2000)
        elseif year <= 2100
            CO2 = CO2_2050 + 180.0 / 50.0 * (year - 2050)
        end

	# ── Time-varying CO₂ experiments ────────────
	elseif cfg.experiment == :co2_sine_wave
        CO2 = 340.0 + 170.0 + 170.0 * cos(2π * (year - 13.0) / 30.0)
        
    elseif cfg.experiment == :co2_step
        CO2 = year >= 1980 ? 340.0 : 680.0
	
	# ── Paleoclimate experiments ────────────────────   
	elseif cfg.experiment == :paleo_231kyr
        CO2 = 200.0
        
    elseif cfg.experiment == :paleo_solar_modern_co2
        CO2 = 340.0
        
    elseif cfg.experiment == :modern_solar_paleo_co2
        CO2 = 200.0
	
	# ── Orbital forcing experiments ─────────────────
	elseif cfg.experiment == :obliquity
        CO2 = 340.0     # Solar forcing loaded externally
        
    elseif cfg.experiment == :eccentricity
        CO2 = 340.0     # Solar forcing loaded externally
        
    elseif cfg.experiment == :earth_sun_distance
        CO2 = 340.0     # Solar constant varies with Earth-Sun distance
	
	# 📂 File I/O dependent experiments (placeholders) ───────────────────────
	elseif cfg.experiment == :rcp26
        error("RCP2.6 scenario requires external CO₂ data file. Not yet implemented.")
        
    elseif cfg.experiment == :rcp45
        error("RCP4.5 scenario requires external CO₂ data file. Not yet implemented.")
        
    elseif cfg.experiment == :rcp60
        error("RCP6.0 scenario requires external CO₂ data file. Not yet implemented.")
        
    elseif cfg.experiment == :rcp85
        CO2 = 340.0  # Handled by boundary conditions
        
    elseif cfg.experiment == :custom_co2
        error("Custom CO₂ scenario requires external trajectory file. Not yet implemented.")

	# 🌍 Regional/partial CO₂ experiments ────────────────────────────────────
	elseif startswith(string(cfg.experiment), "regional_co2_")
        # Reset co2_part to full CO₂ first
        co2_part .= 1.0
        
        if cfg.experiment == :regional_co2_nh
            # 2×CO₂ Northern Hemisphere only
            CO2 = 680.0
            co2_part[:, 1:24] .= 0.5
            
        elseif cfg.experiment == :regional_co2_sh
            # 2×CO₂ Southern Hemisphere only
            CO2 = 680.0
            co2_part[:, 25:48] .= 0.5
            
        elseif cfg.experiment == :regional_co2_tropics
            # 2×CO₂ Tropics only
            CO2 = 680.0
            co2_part[:, 1:15] .= 0.5
            co2_part[:, 33:48] .= 0.5
            for i in 4:4:96
                co2_part[i, 33] = 1.0
                co2_part[i, 15] = 1.0
            end
            
        elseif cfg.experiment == :regional_co2_extratropics
            # 2×CO₂ Extratropics only
            CO2 = 680.0
            co2_part[:, 16:32] .= 0.5
            for i in 4:4:96
                co2_part[i, 32] = 1.0
                co2_part[i, 16] = 1.0
            end
            
        elseif cfg.experiment == :regional_co2_ocean
            # 2×CO₂ Ocean only
            CO2 = 680.0
            for j in 1:ydim, i in 1:xdim
                if z_topo[i, j] > 0.0
                    co2_part[i, j] = 0.5
                end
            end
            icmn_ctrl1 = @view icmn_ctrl[:, :, 1]
            for j in 1:ydim, i in 1:xdim
                if icmn_ctrl1[i, j] >= 0.5
                    co2_part[i, j] = 0.5
                end
            end
            
        elseif cfg.experiment == :regional_co2_land_ice
            # 2×CO₂ Land/Ice only
            CO2 = 680.0
            for j in 1:ydim, i in 1:xdim
                if z_topo[i, j] <= 0.0
                    co2_part[i, j] = 0.5
                end
            end
            icmn_ctrl1 = @view icmn_ctrl[:, :, 1]
            for j in 1:ydim, i in 1:xdim
                if icmn_ctrl1[i, j] >= 0.5
                    co2_part[i, j] = 1.0
                end
            end
            
        elseif cfg.experiment == :regional_co2_winter
            # 2×CO₂ Boreal Winter only
            ityr_step = mod(it - 1, nstep_yr) + 1
            CO2 = (ityr_step <= 181 || ityr_step >= 547) ? 680.0 : 340.0
            
        elseif cfg.experiment == :regional_co2_summer
            # 2×CO₂ Boreal Summer only
            ityr_step = mod(it - 1, nstep_yr) + 1
            CO2 = (ityr_step <= 181 || ityr_step >= 547) ? 340.0 : 680.0
        end
	
	# 📊 Forced boundary condition experiments (handled in scenario loop) ───── 
	elseif cfg.experiment == :elnino || cfg.experiment == :lanina || cfg.experiment == :rcp85
		CO2 = 340.0
	end
	
	return (CO2 = CO2, sw_solar_forcing = sw_solar_forcing)
end

# ╔═╡ 7b8e4130-7d8e-492a-a949-bd7fb6808dd6
md"""
---
## 📊 `diagnostics` - Annual Diagnostics

Accumulates annual‑mean fields and prints a summary at year end.

- Accumulates 12 fields every time step: `Ts`, `Ta`, `To`, `q`, `albedo`, `SW`, `LW`, `Q_lat`, `Q_sens`, and flux corrections `TF_correct`, `qF_correct`
- At `ityr == nstep_yr`: divides by number of steps, prints global mean temperature and two sample points (°C), then resets all accumulators to zero
"""

# ╔═╡ cc0e682c-8767-498e-8178-2de4e796b3a8
function diagnostics!(it, year, CO2, Ts0, Ta0, To0, q0, albedo, sw, lw_surf, q_lat, q_sens, timestate)
	# Accumulate
	Tsmn    .+= Ts0;   Tamn    .+= Ta0;   Tomn    .+= To0
	qmn     .+= q0;    amn     .+= albedo
	swmn    .+= sw;    lwmn    .+= lw_surf
	qlatmn  .+= q_lat; qsensmn .+= q_sens
	ftmn    .+= @view TF_correct[:, :, timestate.ityr]
	fqmn    .+= @view qF_correct[:, :, timestate.ityr]

	if timestate.ityr == nstep_yr
		# Compute annual means
		n = nstep_yr
		Tsmn ./= n; 	 Tamn ./= n; 	Tomn ./= n
        qmn ./= n; 		  amn ./= n
        swmn ./= n; 	 lwmn ./= n
        qlatmn ./= n; qsensmn ./= n
        ftmn ./= n; 	 fqmn ./= n

		# Global mean and sample points (°C)
		global_mean = sum(Tsmn) / (xdim * ydim) - 273.15
        point1 = Tsmn[48, 27] - 273.15   # Tropical Pacific
        point2 = Tsmn[16, 38] - 273.15   # Hamburg/North Europe
        
        println(year, "  ", round(global_mean, digits=2),
                "  ", round(point1, digits=2),
                "  ", round(point2, digits=2))

		# Reset accumulators
		fill!(Tsmn, 0.0); 	fill!(Tamn, 0.0); 	fill!(Tomn, 0.0)
        fill!(qmn, 0.0); 	fill!(amn, 0.0)
        fill!(swmn, 0.0); 	fill!(lwmn, 0.0)
        fill!(qlatmn, 0.0); fill!(qsensmn, 0.0)
        fill!(ftmn, 0.0); 	fill!(fqmn, 0.0)
	end
	return nothing
end

# ╔═╡ 7281dc60-e0e2-4a34-b4b6-c3f63a97f60f
md"""
---
## 💾 `output` - Monthly Output

Accumulates monthly‑mean fields and pushes a `MonthlyRecord` to the output buffer at each month boundary.

- Accumulates 13 fields every time step via `MonthlyAccumulator`
- At the last time step of each calendar month: averages by `days_in_month × steps_per_day`, pushes a copy as a NamedTuple, resets accumulator, advances month counter
- Returns updated `(irec, mon)` for bookkeeping
"""

# ╔═╡ dfdde9f1-b226-4af0-9ac2-36f1b01622fa
function output!(it, irec, mon, Ts0, Ta0, To0, q0, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens,
                 output_buf::Vector{MonthlyRecord}, acc::MonthlyAccumulator, timestate)
    # ----- SAFETY: clamp month to 1..12 -----
    mon = clamp(mon, 1, 12)
    
    accumulate!(acc, Ts0, Ta0, To0, q0, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens)

    # ----- Check end of month -----
    if timestate.jday == jday_mon_cumsum[mon] && (it % ndt_days == 0)
        ndm = jday_mon[mon] * ndt_days
        irec += 1
        push!(output_buf, (
            Ts     = copy(acc.Tmm    ./ ndm),
            Ta     = copy(acc.Tamm   ./ ndm),
            To     = copy(acc.Tomm   ./ ndm),
            q      = copy(acc.qmm    ./ ndm),
            albedo = copy(acc.apmm   ./ ndm),
            ice    = copy(acc.icemm  ./ ndm),
            precip = copy(acc.precipmm ./ ndm),
            evap   = copy(acc.evapmm ./ ndm),
            qcrcl  = copy(acc.qcrclmm ./ ndm),
            sw     = copy(acc.swmm   ./ ndm),
            lw     = copy(acc.lwmm   ./ ndm),
            qlat   = copy(acc.qlatmm ./ ndm),
            qsens  = copy(acc.qsensmm ./ ndm)
        ))
        reset!(acc)
        mon = mon == 12 ? 1 : mon + 1
    end
    return (irec = irec, mon = mon)
end

# ╔═╡ 33fa7b1f-938b-481c-bf25-eca8d7fb33a7
function time_loop!(it, year, CO2, mon, irec, Ts, Ta, q, To, output_buf,
                    ws::CirculationWorkspace, acc::MonthlyAccumulator,
                    timestate, cfg::PhysicsConfig)
    # Calendar lookup
    cal = it <= max_timesteps ? calendar_lookup[it] : (
        day = mod((it - 1) ÷ ndt_days, ndays_yr) + 1,
        step = mod(it - 1, nstep_yr) + 1
    )
    timestate.jday = cal.day
    timestate.ityr = cal.step
    ityr = timestate.ityr

    # Compute tendencies
    tend = tendencies!(CO2, Ts, Ta, To, q, ws, timestate, cfg)

    # Correction views
    TF_corr = @view TF_correct[:, :, ityr]
    qF_corr = @view qF_correct[:, :, ityr]
    ToF_corr = @view ToF_correct[:, :, ityr]

    # Surface temperature
    @. Ts = Ts + tend.dT_ocean + Δt * (tend.SW + tend.LW_surf - tend.LW_down +
                 tend.Q_lat + tend.Q_sens + TF_corr) / cap_surf

    # Air temperature
    @. Ta = Ta + tend.dTa_crcl + Δt * (tend.LW_up + tend.LW_down - tend.em * 
                 tend.LW_surf + tend.Q_lat_air - tend.Q_sens) / cap_air

    # Clamps
    @. Ts = max(Ts, min_T_K)
    @. Ta = max(Ta, min_T_K)

    # Deep ocean
    @. To = To + tend.dTo + ToF_corr

    # Humidity (with clamp)
    dq_eva_use = cfg.log_hydro_dmc ? tend.dq_eva : ws.eva
    dq_rain_use = cfg.log_hydro_dmc ? tend.dq_rain : ws.rain
    dq_crcl_use = cfg.log_crcl_dmc ? tend.dq_crcl : ws.crcl
    @. q = q + clamp(Δt * (dq_eva_use + dq_rain_use) + dq_crcl_use + qF_corr,
                     -min_humidity_change * q, max_humidity_change)

    # Sea ice heat capacity
    seaice!(Ts, timestate, cfg)

    # Conversion to mm/day (analysis units)
    @. ws.precip_out = (-dq_rain_use) * wz_vapor * conv_factor
    @. ws.evap_out  = dq_eva_use * wz_vapor * conv_factor
    @. ws.qcrcl_out = dq_crcl_use
     

    # Output and diagnostics
    (mon, irec) = output!(it, irec, mon, Ts, Ta, To, q, tend.albedo,
                      tend.ice_cover, ws.precip_out, ws.evap_out, ws.qcrcl_out,
                      tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens,
                      output_buf, acc, timestate)
    diagnostics!(it, year, CO2, Ts, Ta, To, q, tend.albedo,
                 tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens, timestate)

    return (mon = mon, irec = irec)
end

# ╔═╡ 4f97badf-a501-4c70-a943-d3c86b48f8a1
function build_monthly_climatology(records::Vector{MonthlyRecord})::Vector{MonthlyRecord}
	isempty(records) && return MonthlyRecord[]

	fields = propertynames(records[1])
	nmonths = 12
	counts = zeros(Int, nmonths)
	clim_acc = [Dict{Symbol, Matrix{Float64}}() for _ in 1:nmonths]

	# Accumulate by mont
	for (idx, rec) in enumerate(records)
		mon = mod(idx - 1, nmonths) + 1
		counts[mon] += 1
		for fld in fields
			arr = getfield(rec, fld)
			if haskey(clim_acc[mon], fld)
				clim_acc[mon][fld] .+= arr
			else
				clim_acc[mon][fld] = copy(arr)
			end
		end
	end

	# Build climatology
	clim = MonthlyRecord[]
	for mon in 1:nmonths
        if counts[mon] == 0
            push!(clim, records[1])
        else
            push!(clim, NamedTuple{fields}(
                Tuple(clim_acc[mon][fld] ./ counts[mon] for fld in fields)
            ))
        end
    end
    return clim
end

# ╔═╡ bf5ccbb7-fff9-42bc-8f2f-2e628b448b3a
function apply_scenario_anomalies(scnr_records::Vector{MonthlyRecord}, ctrl_clim::Vector{MonthlyRecord})::Vector{MonthlyRecord}
	isempty(scnr_records) && return scnr_records
	isempty(ctrl_clim) && return scnr_records

	fields = propertynames(scnr_records[1])
	anom = MonthlyRecord[]

	for (idx, rec) in enumerate(scnr_records)
        mon = mod(idx - 1, 12) + 1
        ref = ctrl_clim[mon]
        push!(anom, NamedTuple{fields}(
          Tuple(getfield(rec, fld) .- getfield(ref, fld) for fld in fields)))
    end
	return anom
end

# ╔═╡ db95cf44-4d37-463c-8fa7-8f084c2a20e8
md"""
---
## 🔄 Data Loading

Before running the model, load the climate input data from the `input/` directory:
"""

# ╔═╡ 4995d3d8-1f95-41d8-be6c-50663edbce10
begin
	# Path to input data directory (adjust if needed)
	jdal2_dir = joinpath(@__DIR__, "greb_dataset_jdal2")
	
	# Load the data
	load_greb_jdal2!(jdal2_dir, dataset=:ncep)
end

# ╔═╡ 83e54812-2291-44e6-9b6e-0c0be57865d3
md"""
---
## `greb_model` - Main Integration Driver

Orchestrates the three simulation phases:

1. **Initialisation** - calls `init_model!` to set up derived fields, experiment overrides, initial conditions, and control CO₂
2. **Flux correction spin‑up** (`time_flux` years) - computes `TF_correct`, `ToF_correct`, `qF_correct` by nudging the model toward climatology; skipped if topography and q‑flux corrections are both disabled
3. **Control run** (`time_ctrl` years) - integrates at constant CO₂ (or constant forcing for orbital experiments); stores monthly output; builds ice climatology from control results
4. **Scenario run** (`time_scnr` years) - integrates with time‑varying CO₂ from `forcing()`; for forced‑boundary experiments (`:elnino`, `:lanina`, `:rcp85`), `Ts` is overwritten with climatology each step; for orbital experiments, starts at year 1 instead of 1950

Returns `(ctrl, scnr)` - vectors of `MonthlyRecord`. Scenario output is converted to anomalies relative to control climatology (except for orbital experiments).
"""

# ╔═╡ 867b193e-3390-4b3b-b5fe-5c399a250660
function compute_annual_ice_climatology(ctrl_output::Vector{MonthlyRecord})
    ice_months = zeros(Float64, xdim, ydim, 12)
    count = zeros(Int, 12)
    
    for (idx, rec) in enumerate(ctrl_output)
        mon = mod1(idx, 12)          # month 1..12 based on record index
        @. ice_months[:, :, mon] += rec.ice   # note: field name is `ice`, not `ice_cover`
        count[mon] += 1
    end
    
    for mon in 1:12
        cnt = count[mon]
        cnt > 0 || continue
        @. ice_months[:, :, mon] *= 1.0 / cnt
    end
    
    return ice_months
end

# ╔═╡ d17c570f-ad56-4d43-98cb-ccd03a8fcb4a
md"""
---
## 🎛️ Interactive Control Panel

Configure and run GREB experiments interactively via Pluto widgets.

**Controls:**
- **Experiment** - preset experiment type (full model, 2×CO₂, El Niño, etc.)
- **Configuration preset** - process isolation presets (full physics, no feedbacks, MSCM, sensitivity test, custom)
- **Mean climate / CO₂ response / circulation panels** - individual process toggles for custom configurations
- **Hydrology** - rain parameterization mode, evaporation mode, climatology dataset
- **External forcing** - toggles for forced surface temperature, wind, and vertical velocity

**Run controls:**
- Flux correction years, control run years, scenario run years (sliders)
- Run toggle checkbox to execute the model
- Results stored in `last_run` (control and scenario `MonthlyRecord` vectors)
"""

# ╔═╡ 71676720-4685-4216-a55e-5918829d258f
md"""
### Experiment Selection
"""

# ╔═╡ b7cef546-f0f4-4f42-a0ed-00ffb92a422e
begin
@bind experiment_type Select([
    :full_model => "Full Model (default)",
    :constant_topo => "Constant Topography",
    :co2_double => "2×CO₂",
    :co2_quadruple => "4×CO₂",
    :solar_plus27 => "Solar +27W/m²",
    :elnino => "El Niño",
    :lanina => "La Niña",
    :paleo_231kyr => "Paleo 231kyr",
    :rcp85 => "RCP8.5 Climate Change"
], default=:full_model)
end

# ╔═╡ 7969e897-121e-4fe0-9b75-36d21931357f
md"""
#### 🎛️ Configuration Preset
"""

# ╔═╡ dec9de54-eede-45ae-94cd-5bcc477298bf
@bind config_preset Select([
	"full" => "Full Physics (All processes active)",
	"no_feedbacks" => "No Feedbacks (Fixed albedo, clouds, vapor)",
	"mscm" => "MSCM Original (Mean State Climate Model)",
	"sensitivity" => "Sensitivity Test (Minimal processes)",
	"custom" => "Custom Configuration"
], default="full")

# ╔═╡ 2d2aa153-1b51-4972-89a0-de0f5d803838
md"""
#### 🌡️ Mean Climate Processes
"""

# ╔═╡ 66527a77-1926-41f7-b86f-84e0d65d9f30
@bind mean_climate_panel PlutoUI.combine() do Child
	md"""
	$(Child("clouds", CheckBox(default=true))) Clouds  
	$(Child("vapor", CheckBox(default=true))) Water Vapor  
	$(Child("ice", CheckBox(default=true))) Ice-Albedo Feedback  
	$(Child("crcl", CheckBox(default=true))) Circulation  
	$(Child("hydro", CheckBox(default=true))) Hydrology  
	$(Child("atmos", CheckBox(default=true))) Atmosphere  
	$(Child("co2", CheckBox(default=true))) CO₂  
	$(Child("ocean", CheckBox(default=true))) Ocean  
	$(Child("qflux", CheckBox(default=true))) Q-Flux Corrections
	"""
end

# ╔═╡ 393a8bae-c5ee-4a46-85f3-2c6a7cead357
md"""
#### 🔥 CO₂ Response Processes
"""

# ╔═╡ 506e51a4-d054-4e5d-8741-c4d36a9f2714
@bind co2_response_panel PlutoUI.combine() do Child
	md"""
	$(Child("clouds", CheckBox(default=true))) Clouds  
	$(Child("vapor", CheckBox(default=true))) Water Vapor  
	$(Child("crcl", CheckBox(default=true))) Circulation  
	$(Child("hydro", CheckBox(default=true))) Hydrology  
	$(Child("topo", CheckBox(default=true))) Topography  
	$(Child("humid", CheckBox(default=true))) Humidity Climatology
	"""
end

# ╔═╡ 3ac049a9-50e2-4db6-acb9-aacc20586ca4
md"""
#### ⚡ Circulation Components
"""

# ╔═╡ 26b4d172-72ec-4dd6-83c0-0779904634ae
@bind circulation_panel PlutoUI.combine() do Child
	md"""
	$(Child("ice", CheckBox(default=true))) Ice-Albedo Feedback  
	$(Child("hdif", CheckBox(default=true))) Horizontal Diffusion  
	$(Child("hadv", CheckBox(default=true))) Horizontal Advection  
	$(Child("vdif", CheckBox(default=true))) Vertical Diffusion  
	$(Child("vadv", CheckBox(default=true))) Vertical Advection  
	$(Child("conv", CheckBox(default=true))) Moisture Convergence
	"""
end

# ╔═╡ 2d4f533e-a240-459d-9575-f9982c03184b
md"""
#### 💧 Hydrology Parameterizations
"""

# ╔═╡ 8ef38e96-42f4-4bb2-9063-88db9842fb47
@bind hydro_rain_mode Select([
	-1 => "Original GREB",
	0 => "Best Performance",
	1 => "+Relative Humidity",
	2 => "+Omega Convergence",
	3 => "+Both RH & Omega"
], default=0)

# ╔═╡ d5331962-6bb1-485d-b538-ff35221e14e6
@bind hydro_eva_mode Select([
	-1 => "Original GREB",
	0 => "Skin Temperature",
	1 => "Enhanced",
	2 => "Wind Speed Climatology"
], default=-1)

# ╔═╡ 54e19b14-4bc1-4086-9015-ff2cd8f4afbf
@bind hydro_clim_dataset Select([
	0 => "ERA-Interim",
	1 => "NCEP"
], default=0)

# ╔═╡ efb23f98-fb75-42e6-8681-c5147a82a09c
md"""
#### 📡 External Forcing
"""

# ╔═╡ d74882fa-3a38-4cc7-b476-e52cb5a7db31
@bind external_forcing_panel PlutoUI.combine() do Child
	md"""
	$(Child("tsurf", CheckBox(default=false))) Surface Temperature  
	$(Child("hwind", CheckBox(default=false))) Horizontal Wind  
	$(Child("omega", CheckBox(default=false))) Vertical Velocity
	"""
end

# ╔═╡ dc0f1a5e-1e7b-4a9d-816e-731d2747ac4e
begin
	# Apply preset configurations or use custom settings
	if config_preset == "full"
		# Full Physics: All processes active
		log_clouds_dmc = log_vapor_dmc = log_ice_dmc = true
        log_crcl_dmc = log_hydro_dmc = true
        log_atmos_dmc = log_co2_dmc = log_ocean_dmc = log_qflux_dmc = true
        
        log_clouds_drsp = log_vapor_drsp = true
        log_crcl_drsp = log_hydro_drsp = true
        log_topo_drsp = log_humid_drsp = true
        
        log_ice = log_hdif = log_hadv = true
        log_vdif = log_vadv = log_conv = true
        
        log_rain, log_eva, log_clim = 0, -1, 0
        log_tsurf_ext = log_hwind_ext = log_omega_ext = false
		
	elseif config_preset == "no_feedbacks"
		# No Feedbacks: Fixed albedo, clouds, vapor
		log_clouds_dmc = log_vapor_dmc = log_ice_dmc = false
        log_crcl_dmc = log_hydro_dmc = true
        log_atmos_dmc = log_co2_dmc = log_ocean_dmc = log_qflux_dmc = true
        
        log_clouds_drsp = log_vapor_drsp = false
        log_crcl_drsp = log_hydro_drsp = true
        log_topo_drsp = log_humid_drsp = true
        
        log_ice = false
        log_hdif = log_hadv = true
        log_vdif = log_vadv = log_conv = true
        
        log_rain, log_eva, log_clim = 0, -1, 0
        log_tsurf_ext = log_hwind_ext = log_omega_ext = false
		
	elseif config_preset == "mscm"
		# MSCM Original: Mean State Climate Model configuration
		log_clouds_dmc = log_vapor_dmc = log_ice_dmc = true
        log_crcl_dmc = log_hydro_dmc = true
        log_atmos_dmc = log_co2_dmc = log_ocean_dmc = log_qflux_dmc = true
        
        log_clouds_drsp = log_vapor_drsp = true
        log_crcl_drsp = log_hydro_drsp = true
        log_topo_drsp = log_humid_drsp = true
        
        log_ice = log_hdif = log_hadv = true
        log_vdif = log_vadv = log_conv = true
        
        log_rain, log_eva, log_clim = 0, -1, 0
        log_tsurf_ext = log_hwind_ext = log_omega_ext = false
		
	elseif config_preset == "sensitivity"
		# Sensitivity Test: Minimal processes for testing
		log_clouds_dmc = log_vapor_dmc = log_ice_dmc = false
        log_crcl_dmc = log_hydro_dmc = false
        log_atmos_dmc = log_co2_dmc = log_ocean_dmc = true
        log_qflux_dmc = false
        
        log_clouds_drsp = log_vapor_drsp = false
        log_crcl_drsp = log_hydro_drsp = false
        log_topo_drsp = log_humid_drsp = false
        
        log_ice = log_hdif = log_hadv = false
        log_vdif = log_vadv = log_conv = false
        
        log_rain, log_eva, log_clim = -1, -1, 0
        log_tsurf_ext = log_hwind_ext = log_omega_ext = false
		
	else # custom
		log_clouds_dmc = mean_climate_panel.clouds
        log_vapor_dmc = mean_climate_panel.vapor
        log_ice_dmc = mean_climate_panel.ice
        log_crcl_dmc = mean_climate_panel.crcl
        log_hydro_dmc = mean_climate_panel.hydro
        log_atmos_dmc = mean_climate_panel.atmos
        log_co2_dmc = mean_climate_panel.co2
        log_ocean_dmc = mean_climate_panel.ocean
        log_qflux_dmc = mean_climate_panel.qflux
        
        log_clouds_drsp = co2_response_panel.clouds
        log_vapor_drsp = co2_response_panel.vapor
        log_crcl_drsp = co2_response_panel.crcl
        log_hydro_drsp = co2_response_panel.hydro
        log_topo_drsp = co2_response_panel.topo
        log_humid_drsp = co2_response_panel.humid
        
        log_ice = circulation_panel.ice
        log_hdif = circulation_panel.hdif
        log_hadv = circulation_panel.hadv
        log_vdif = circulation_panel.vdif
        log_vadv = circulation_panel.vadv
        log_conv = circulation_panel.conv
        
        log_rain = hydro_rain_mode
        log_eva = hydro_eva_mode
        log_clim = hydro_clim_dataset
        
        log_tsurf_ext = external_forcing_panel.tsurf
        log_hwind_ext = external_forcing_panel.hwind
        log_omega_ext = external_forcing_panel.omega
	end
end;

# ╔═╡ fee2639d-d8cf-4ee3-b824-73a0b02fef4a
md"""
### Run Duration Configuration (years)
"""

# ╔═╡ 952d5e5d-119a-4ed9-9e22-0b0fe66ae04d
@bind time_flux Slider(0:10, default=0, show_value=true)

# ╔═╡ 584e767e-4dc5-4821-af63-d6a825326d9e
function qflux_correction!(CO2_ctrl, Ts, Ta, q, To, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    for it in 1:(time_flux * ndt_days * ndays_yr)
        timestate.jday = mod((it - 1) ÷ ndt_days, ndays_yr) + 1
        timestate.ityr = mod(it - 1, nstep_yr) + 1
        ityr = timestate.ityr

        tend = tendencies!(CO2_ctrl, Ts, Ta, To, q, ws, timestate, cfg)

        # Views into climatology & correction fields
        Tc   = @view Tclim[:, :, ityr]
        Toc  = @view Toclim[:, :, ityr]
        qc   = @view qclim[:, :, ityr]
        TFc  = @view TF_correct[:, :, ityr]
        ToFc = @view ToF_correct[:, :, ityr]
        qFc  = @view qF_correct[:, :, ityr]

        # ── Surface temperature ──────────────────────────────
        # Uncorrected state (store in workspace buffer)
        @. ws.Ts0_buf = Ts + tend.dT_ocean + Δt * (
            tend.SW + tend.LW_surf - tend.LW_down +
            tend.Q_lat + tend.Q_sens
        ) / cap_surf

        # Correction and corrected state
        @. TFc = (Tc - ws.Ts0_buf) * cap_surf / Δt
        @. ws.Ts0_buf = ws.Ts0_buf + TFc * Δt / cap_surf

        # ── Air temperature ──────────────────────────────────
        @. ws.Ta0_buf = Ta + tend.dTa_crcl + ΔT_AIR_FACTOR * (
            tend.LW_up + tend.LW_down - tend.em * tend.LW_surf +
            tend.Q_lat_air - tend.Q_sens
        )

        # ── Deep ocean ───────────────────────────────────────
        @. ws.To0_buf = To + tend.dTo
        @. ToFc = Toc - ws.To0_buf
        @. ws.To0_buf = ws.To0_buf + ToFc

        # ── Humidity ─────────────────────────────────────────
        @. ws.q0_buf = q + tend.dq_crcl + Δt * (tend.dq_eva + tend.dq_rain)
        @. qFc = qc - ws.q0_buf
        @. ws.q0_buf = ws.q0_buf + qFc

        # Sea ice (updates cap_surf in place)
        seaice!(ws.Ts0_buf, timestate, cfg)

        # Diagnostics
        diagnostics!(it, 0.0, CO2_ctrl, ws.Ts0_buf, ws.Ta0_buf, ws.To0_buf, ws.q0_buf,
                     tend.albedo, tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens, timestate)

        # Advance state (broadcast – same as `.=`)
        @. Ts = ws.Ts0_buf
        @. Ta = ws.Ta0_buf
        @. q  = ws.q0_buf
        @. To = ws.To0_buf
    end
    return nothing
end

# ╔═╡ 92c3bd68-bd07-4381-9c04-e6611650cd1e
function greb_model!(time_flux, time_ctrl, time_scnr, cfg::PhysicsConfig)

	# ── 1. Initialisation ───────────────────────────────────────
	ini = init_model!(cfg)
	Ts_ini   = ini.Ts_ini;	Ta_ini   = ini.Ta_ini
	To_ini   = ini.To_ini;	q_ini    = ini.q_ini
	CO2_ctrl = ini.CO2_ctrl

	# Determine experiment type from cfg
    is_orbital_exp = cfg.experiment in (:obliquity, :eccentricity, :earth_sun_distance)
    is_forced_boundary = cfg.experiment in (:rcp85, :elnino, :lanina)
    is_sst_plus1 = cfg.experiment == :sst_plus1

	# Workspace and accumulator
    ws = CirculationWorkspace()
    acc = MonthlyAccumulator()

	# Initialize time state
    timestate = TimeState(1, 1)
	
	# ── 2. Flux-correction spin-up ──────────────────────────────
	if cfg.log_topo_drsp || cfg.log_qflux_dmc
        if !cfg.log_topo_drsp && cfg.log_qflux_dmc
            println("% loading flux correction fields...")
            load_flux_corrections_jdal2!(jdal2_dir)
        end
        println("% flux correction  CO2 = ", CO2_ctrl)
        qflux_correction!(CO2_ctrl, Ts_ini, Ta_ini, q_ini, To_ini, timestate, cfg, ws)
    else
        println("Flux correction skipped")
    end
	
	# Reset accumulators after spin-up
    reset!(acc)

	# ── 3. Control run ──────────────────────────────────────────
	println("CONTROL RUN: CO2 = ", CO2_ctrl, " time = ", time_ctrl, " yr")
	
	# Initialize state arrays
	Ts = copy(Ts_ini);  Ta = copy(Ta_ini)
	To = copy(To_ini);  q  = copy(q_ini)
	sw_solar_forcing_state[] = 1.0
	mon  = 1;  year = 1970;  irec = 0

	ctrl_output = MonthlyRecord[]
	sizehint!(ctrl_output, time_ctrl * 12)  # Pre-allocate for all months
	timestate = TimeState(1, 1)  # Initialize time state

	for it in 1:(time_ctrl * nstep_yr)
        (mon, irec) = time_loop!(it, year, CO2_ctrl, mon, irec,
                                 Ts, Ta, q, To, ctrl_output, ws, acc, timestate, cfg)
        if mod(it, nstep_yr) == 0
            year += 1
        end
    end

	# ── Build ice climatology from control output ───────────────
    # Compute annual mean ice cover from the stored control monthly means
    # (Assuming control output contains monthly ice cover; adjust if needed)
    ice_forcing = compute_annual_ice_climatology(ctrl_output)   

	# ── 4. Scenario run ─────────────────────────────────────────
	println("SCENARIO: ", cfg.experiment, "  time = ", time_scnr, " yr")
	
	# Reset state to initial conditions
	Ts .= Ts_ini;  Ta .= Ta_ini
	q  .= q_ini;   To .= To_ini
    year = is_orbital_exp ? 1 : 1950
	CO2 = 340.0;  mon = 1;  irec = 0

	sw_solar_forcing_state[] = 1.0
	reset!(acc)  # Use accumulator reset

	scnr_output = MonthlyRecord[]
	if time_scnr > 0
    	sizehint!(scnr_output, time_scnr * 12)
	end

	for it in 1:(time_scnr * nstep_yr)
        # Obtain forcing (CO2 and solar multiplier)
        forcing_result = forcing(it, year, cfg, ice_forcing; nstep_yr=nstep_yr)
        CO2 = forcing_result.CO2
        sw_solar_forcing_state[] = forcing_result.sw_solar_forcing

        # Forced‑boundary experiments: overwrite Ts with climatology
        if is_forced_boundary
            ityr_now = mod(it - 1, nstep_yr) + 1
            Ts .= @view Tclim[:, :, ityr_now]
        end

        # SST+1 K experiment
        if is_sst_plus1
            CO2 = CO2_ctrl
            ityr_now = mod(it - 1, nstep_yr) + 1
            @. Ts = ifelse(z_topo < 0.0, Tclim[:, :, ityr_now] + 1.0, Ts)
        end

        (mon, irec) = time_loop!(it, year, CO2, mon, irec,
                                 Ts, Ta, q, To, scnr_output, ws, acc, timestate, cfg)

        if mod(it, nstep_yr) == 0
            year += 1
        end
    end

	# Post‑processing: anomalies for non‑orbital experiments
    if !is_orbital_exp && !isempty(ctrl_output) && !isempty(scnr_output)
        ctrl_clim = build_monthly_climatology(ctrl_output)
        scnr_output = apply_scenario_anomalies(scnr_output, ctrl_clim)
    end

    return (ctrl = ctrl_output, scnr = scnr_output)
end

# ╔═╡ 82bc0ee6-d0ca-4f47-a43f-67ae6702bb6a
md"""
**Flux correction years:** $(time_flux)
"""

# ╔═╡ 13647b0a-561c-4f65-a64c-b9a4ea76c9f0
@bind time_ctrl Slider(0:100, default=10, show_value=true)

# ╔═╡ 28371449-d04e-4d64-8b4a-cd37ceb25ef5
md"""
**Control run years:** $(time_ctrl)
"""

# ╔═╡ 4d3bedd5-1d4a-4150-828f-c5352c70874b
@bind time_scnr Slider(0:100, default=0, show_value=true)

# ╔═╡ 0ab27970-986e-4359-bed1-8908c5fd109b
md"""
**Scenario run years:** $(time_scnr)
"""

# ╔═╡ 10e20296-cc4f-4c4a-80bf-c39ae6450e85
md"""
### Execute Model

**Current Configuration:**
- Experiment: $(experiment_type)
- Flux correction: $time_flux years
- Control run: $time_ctrl years  
- Scenario run: $time_scnr years
"""

# ╔═╡ 9b7a847e-a0d2-4421-9fb9-ff01b73c4194
@bind run_toggle CheckBox(default=false)

# ╔═╡ f7fb4846-5ee8-4ca5-a76f-1edcd58a0559
md"""
### Accessing Results

After running the model, results are available in the `last_run` variable:

```julia
# Control run monthly records
last_run.ctrl  # Vector{MonthlyRecord}

# Scenario run monthly records  
last_run.scnr  # Vector{MonthlyRecord}

# Each MonthlyRecord contains:
# .Ts      - Surface temperature (96×48)
# .Ta      - Atmospheric temperature (96×48)
# .To      - Ocean temperature (96×48)
# .q       - Atmospheric humidity (96×48)
# .albedo  - Surface albedo (96×48)
# .ice     - Ice coverage (96×48)
# .precip  - Precipitation (96×48)
# .evap    - Evaporation (96×48)
# .qcrcl   - Circulation heat flux (96×48)
# .sw      - Shortwave radiation (96×48)
# .lw      - Longwave radiation (96×48)
# .qlat    - Latent heat flux (96×48)
# .qsens   - Sensible heat flux (96×48)
```

**Example: Plot global mean temperature evolution**
```julia
Ts_global_mean = [mean(rec.Ts) for rec in last_run.ctrl]
plot(Ts_global_mean, xlabel="Month", ylabel="Global Mean Ts [K]")
```
"""

# ╔═╡ 39afec69-f90d-4ff2-99d1-5cfec84d09a1
function current_physics_config()
	return PhysicsConfig(
		# Mean Climate Switches
		log_clouds_dmc = log_clouds_dmc,
		log_vapor_dmc = log_vapor_dmc,
		log_ice_dmc = log_ice_dmc,
		log_crcl_dmc = log_crcl_dmc,
		log_hydro_dmc = log_hydro_dmc,
		log_atmos_dmc = log_atmos_dmc,
		log_co2_dmc = log_co2_dmc,
		log_ocean_dmc = log_ocean_dmc,
		log_qflux_dmc = log_qflux_dmc,
		# CO₂ Response Switches
		log_clouds_drsp = log_clouds_drsp,
		log_vapor_drsp = log_vapor_drsp,
		log_crcl_drsp = log_crcl_drsp,
		log_hydro_drsp = log_hydro_drsp,
		log_topo_drsp = log_topo_drsp,
		log_humid_drsp = log_humid_drsp,
		# Circulation Components
		log_ice = log_ice,
		log_hdif = log_hdif,
		log_hadv = log_hadv,
		log_vdif = log_vdif,
		log_vadv = log_vadv,
		log_conv = log_conv,
		# Hydrology Parameters
		log_rain = log_rain,
		log_eva = log_eva,
		log_clim = log_clim,
		# External Forcing
		log_tsurf_ext = log_tsurf_ext,
		log_hwind_ext = log_hwind_ext,
		log_omega_ext = log_omega_ext,
		experiment = experiment_type,
	)
end

# ╔═╡ cf2a51eb-a661-4c66-9e0a-d1730110e4bc
begin
    # Ensure last_run is always defined globally
    if !@isdefined(last_run)
        global last_run = nothing
    end
    
    # Only run when toggle is ON
    if run_toggle
		cfg = current_physics_config()
		
        # Execute model with current parameters
        result = greb_model!(time_flux, time_ctrl, time_scnr, cfg)
        
        # Store globally for visualization cells (ALWAYS update this)
        global last_run = result
        
        # Show results
        md"""
        **✅ Run complete!**  
        Control months: $(length(result.ctrl))  
        Scenario months: $(length(result.scnr))
        """
    else
        # Show ready message but DON'T change last_run
        md"""**⏸️ Toggle ON to run model**"""
    end
end

# ╔═╡ b0d6a9e3-5c3a-4408-ab87-1067c2dfaeb7
begin
    if @isdefined(last_run) && !isnothing(last_run) && !isempty(last_run.ctrl)
        println("="^50)
        println("GREB MODEL OUTPUT")
        println("="^50)
        
        rec = last_run.ctrl[1]
        
        println("\n🌡️ TEMPERATURE (K):")
        println("   Surface (Ts): min=$(round(minimum(rec.Ts), digits=1)), mean=$(round(mean(rec.Ts), digits=1)), max=$(round(maximum(rec.Ts), digits=1))")
        println("   Air (Ta):     min=$(round(minimum(rec.Ta), digits=1)), mean=$(round(mean(rec.Ta), digits=1)), max=$(round(maximum(rec.Ta), digits=1))")
        println("   Ocean (To):   min=$(round(minimum(rec.To), digits=1)), mean=$(round(mean(rec.To), digits=1)), max=$(round(maximum(rec.To), digits=1))")
        
        println("\n💧 HYDROLOGY (mm/day):")
        println("   Precip: min=$(round(minimum(rec.precip), digits=2)), mean=$(round(mean(rec.precip), digits=2)), max=$(round(maximum(rec.precip), digits=2))")
        println("   Evap:   min=$(round(minimum(rec.evap), digits=2)), mean=$(round(mean(rec.evap), digits=2)), max=$(round(maximum(rec.evap), digits=2))")
        
        println("\n☀️ RADIATION (W/m²):")
        println("   SW: min=$(round(minimum(rec.sw), digits=1)), mean=$(round(mean(rec.sw), digits=1)), max=$(round(maximum(rec.sw), digits=1))")
        println("   LW: min=$(round(minimum(rec.lw), digits=1)), mean=$(round(mean(rec.lw), digits=1)), max=$(round(maximum(rec.lw), digits=1))")
        
        println("\n🪞 ALBEDO:")
        println("   min=$(round(minimum(rec.albedo), digits=2)), mean=$(round(mean(rec.albedo), digits=2)), max=$(round(maximum(rec.albedo), digits=2))")
        
        println("\n❄️ ICE FRACTION:")
        println("   min=$(round(minimum(rec.ice), digits=2)), mean=$(round(mean(rec.ice), digits=2)), max=$(round(maximum(rec.ice), digits=2))")
        
        println("\n✅ All finite: $(all(isfinite, rec.Ts))")
        
    else
        println("❌ No model output found. Run the model first!")
        println("   Set time_ctrl >= 1 and toggle run_toggle ON")
    end
end

# ╔═╡ a9089bf1-3241-4038-81d6-b953a00c3cef
md"""
---
## ⏱️ Benchmarking

Benchmarking components to find optimizations.
"""

# ╔═╡ 9fa7c129-2289-4979-aa1f-8f8112e7875f
md"""### SWradiation:"""

# ╔═╡ caf443c0-5fbc-4fd4-8ef9-3565417d40ab
# ╠═╡ skip_as_script = true
#=╠═╡
# Benchmark SWradiation
@benchmark SWradiation!($bench_data.Ts, $bench_data.timestate, $bench_data.cfg, $bench_data.ws)
  ╠═╡ =#

# ╔═╡ 1d46e4a8-41fc-4b76-bc9d-a115b6b6d9c4
md"""### LWradiation:"""

# ╔═╡ e20657c6-cf03-4f84-a782-2a3e5d6d303c
# ╠═╡ skip_as_script = true
#=╠═╡
# Benchmark LWradiation  
@benchmark LWradiation!($bench_data.Ts, $bench_data.Ta, $bench_data.q, $bench_data.CO2, $bench_data.timestate, $bench_data.cfg, $bench_data.ws)
  ╠═╡ =#

# ╔═╡ e2227570-2fcb-4090-8bbb-ec85af3a519c
md"""### hydro:"""

# ╔═╡ 1b94212a-a5b9-46dc-81ea-896e4c839755
# ╠═╡ skip_as_script = true
#=╠═╡
# Benchmark hydrological cycle
@benchmark hydro!($bench_data.Ts, $bench_data.q, $bench_data.timestate, 
                  $bench_data.cfg, $bench_data.ws)
  ╠═╡ =#

# ╔═╡ 40717ac9-ff3f-4835-ba7e-3cde10b2c3f0
md"""### circulation (temperature):"""

# ╔═╡ fb0719e8-b9fe-4ca3-8c5a-de3586d193a9
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark circulation!($bench_data.Ta, $z_air, $bench_data.ws.dX_crcl, $bench_data.ws, $bench_data.timestate, $bench_data.cfg)
  ╠═╡ =#

# ╔═╡ 0c201c53-0db2-4805-beec-e7705082509d
md"""### circulation (vapor):"""

# ╔═╡ eb10715c-247d-40e4-b11b-aceabb60628c
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark circulation!($bench_data.q, $z_vapor, $bench_data.ws.dX_crcl, $bench_data.ws, $bench_data.timestate, $bench_data.cfg)
  ╠═╡ =#

# ╔═╡ f9da0fed-a0a3-4bfb-8bac-6a839177bc73
md"""### seaice:"""

# ╔═╡ 2e05f3f6-687b-4096-a869-d64564db60da
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark seaice!($bench_data.Ts, $bench_data.timestate, $bench_data.cfg)
  ╠═╡ =#

# ╔═╡ 09cbb86e-64e1-45e5-b0f6-1e8b868506e8
md"""### deep_ocean:"""

# ╔═╡ c186464e-e705-49be-ad5d-4352ce719c1f
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark deep_ocean!($bench_data.Ts, $bench_data.To, $bench_data.timestate, $bench_data.cfg, $bench_data.ws)
  ╠═╡ =#

# ╔═╡ 395242dc-cb73-481f-ad25-71eba93d341f
md"""### diffusion:"""

# ╔═╡ 718e4f4b-af27-4271-a37f-83ea6e30bdbe
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	T1 = rand(Float64, xdim, ydim)
	h_scl = 8400.0  # z_air
	@benchmark diffusion!($T1, $h_scl, $bench_data.ws, $bench_data.timestate)
end
  ╠═╡ =#

# ╔═╡ 59c97f2b-de11-4f42-9fcf-aaf2924d8677
md"""### advection:"""

# ╔═╡ b9b3564e-c59b-4f3a-80d7-4431fd8ce2db
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark advection!($bench_data.T_test, $z_air, $bench_data.ws, $bench_data.timestate, $bench_data.cfg)
  ╠═╡ =#

# ╔═╡ 435bae7c-e0d3-46af-b185-5221f687d6bf
md"""### convergence:"""

# ╔═╡ 1229bf7d-204c-41a2-97dd-9d75e79b18ce
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark convergence!($bench_data.T_test, $omegaclim, $bench_data.timestate, $bench_data.ws)
  ╠═╡ =#

# ╔═╡ 426c1511-f912-4e10-9b8b-f858471019bc
md"""### tendencies:"""

# ╔═╡ 50777c8a-ee6f-449b-b704-bba7e9b36490
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	@benchmark tendencies!($bench_data.CO2, $bench_data.Ts, $bench_data.Ta, $bench_data.To, $bench_data.q, $bench_data.ws, $bench_data.timestate, $bench_data.cfg)
end
  ╠═╡ =#

# ╔═╡ 99b335e7-b28a-44f2-a07d-a84c16366710
# ╠═╡ disabled = true
#=╠═╡
begin
# 1. Profile tendencies! with your benchmark data
@profile for i in 1:10
    tendencies!(bench_data.CO2, bench_data.Ts, bench_data.Ta, bench_data.To, 
                bench_data.q, bench_data.ws, bench_data.timestate, bench_data.cfg)
end

# 2. Flat profile (function-level)
println("\n=== FLAT PROFILE (function-level) ===")
Profile.print(format=:flat, mincount=10, C=false)

# 3. Tree profile (call hierarchy) - MORE USEFUL
println("\n=== TREE PROFILE (call hierarchy) ===")
Profile.print(format=:tree, mincount=10, C=false)

# 4. Allocations profile (find memory hotspots)
println("\n=== ALLOCATION PROFILE ===")
Profile.Allocs.@profile sample_rate=1 begin
    for i in 1:10
        tendencies!(bench_data.CO2, bench_data.Ts, bench_data.Ta, bench_data.To, 
                    bench_data.q, bench_data.ws, bench_data.timestate, bench_data.cfg)
    end
end
Profile.Allocs.print(mincount=100)
end
  ╠═╡ =#

# ╔═╡ d571c448-cea6-49c0-81a4-2376e4fc3f67
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
# profile_greb_model_save
begin
    profile_cfg = PhysicsConfig(experiment=:full_model)
    
    # Run profile and save to file
    Profile.clear()
    @profile greb_model!(0, 1, 0, profile_cfg)
    
    # Save to file
    open("greb_profile.txt", "w") do io
        Profile.print(io, format=:tree, mincount=1, C=false)
    end
    
    println("Profile saved to greb_profile.txt")
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
LoopVectorization = "bdcacae8-1622-11e9-2a5c-532679323890"
NCDatasets = "85f8d34a-cbdd-5861-8df4-14fed0d494ab"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Profile = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
ProfileSVG = "132c30aa-f267-4189-9183-c8a63c7e05e6"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
BenchmarkTools = "~1.8.0"
LoopVectorization = "~0.12.173"
NCDatasets = "~0.14.10"
PlutoUI = "~0.7.65"
ProfileSVG = "~0.2.2"
StaticArrays = "~1.9.13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.6"
manifest_format = "2.0"
project_hash = "04fda499e3cb2e0a3310aec3a5f4f1e0d0da7c3c"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7715e5b2b186c4d9b664d299d2c9e48b9a778c88"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.6.1"

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

    [deps.Adapt.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "3d0cabd25fab32390e3bcb82cd67e700aebd9816"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.25.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceAMDGPUExt = "AMDGPU"
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BenchmarkTools]]
deps = ["Compat", "JSON", "Logging", "PrecompileTools", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "9670d3febc2b6da60a0ae57846ba74670290653f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.8.0"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.Blosc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Lz4_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "535c80f1c0847a4c967ea945fca21becc9de1522"
uuid = "0b7ba130-8d10-5ba8-a3d6-c5182647fed9"
version = "1.21.7+0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CFTime]]
deps = ["Dates", "Printf"]
git-tree-sha1 = "fad2f199d1f1ae0c8e820ab68f81f6dbf62e60b2"
uuid = "179af706-886a-5703-950a-314cd64e0468"
version = "0.2.10"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Preferences", "Static"]
git-tree-sha1 = "f3a21d7fc84ba618a779d1ed2fcca2e682865bab"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.7"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonDataModel]]
deps = ["CFTime", "DataStructures", "Dates", "DiskArrays", "Preferences", "Printf", "Statistics"]
git-tree-sha1 = "bf07704e843daabd2cb2bb1404571656f80bce16"
uuid = "1fbeeb36-5f17-413c-809b-666fb144f157"
version = "0.4.3"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "6fb53a69613a0b2b68a0d12671717d307ab8b24e"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.5"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DiskArrays]]
deps = ["ConstructionBase", "LRUCache", "Mmap", "OffsetArrays"]
git-tree-sha1 = "7821ce71d0b9c2948ab80f86237f1f4212dca861"
uuid = "3c3547ce-8d99-4f5e-a174-61eb10b00ae3"
version = "0.4.21"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "8e9c059d6857607253e837730dbf780b6b151acd"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.19.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.FlameGraphs]]
deps = ["AbstractTrees", "Colors", "FileIO", "FixedPointNumbers", "IndirectArrays", "LeftChildRightSiblingTrees", "Profile"]
git-tree-sha1 = "0166baf81babb91cf78bfcc771d8e87c43d568df"
uuid = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"
version = "1.1.0"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "LibCURL_jll", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "libaec_jll"]
git-tree-sha1 = "e94f84da9af7ce9c6be049e9067e511e17ff89ec"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.14.6+0"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Preferences", "Static"]
git-tree-sha1 = "af9ab7d1f70739a47f03be78771ebda38c3c71bf"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.18"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XML2_jll", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "baaaebd42ed9ee1bd9173cfd56910e55a8622ee1"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.13.0+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LRUCache]]
git-tree-sha1 = "5519b95a490ff5fe629c4a7aa3b3dfc9160498b3"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.6.2"
weakdeps = ["Serialization"]

    [deps.LRUCache.extensions]
    SerializationExt = ["Serialization"]

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "95ba48564903b43b2462318aa243ee79d81135ff"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.2.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "a9fc7883eb9b5f04f46efb9a540833d1fad974b3"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.173"

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    ForwardDiffNNlibExt = ["ForwardDiff", "NNlib"]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.LoopVectorization.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "191686b1ac1ea9c89fc52e996ad15d1d241d1e33"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.10.1+0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "9341048b9f723f2ae2a72a5269ac2f15f80534dc"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.3.2+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "8e98d5d80b87403c311fd51e8455d4546ba7a5f8"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.12"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "675df097f8eeb28998b2cfe3b25655af73d5f7df"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.5.6+0"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bc95bf4149bf535c09602e3acdf950d9b4376227"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+3"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NCDatasets]]
deps = ["CFTime", "CommonDataModel", "DataStructures", "Dates", "DiskArrays", "NetCDF_jll", "NetworkOptions", "Printf"]
git-tree-sha1 = "c82c73e2e0c57a0fe13d3414d7c5a6a821d24016"
uuid = "85f8d34a-cbdd-5861-8df4-14fed0d494ab"
version = "0.14.10"

    [deps.NCDatasets.extensions]
    NCDatasetsMPIExt = "MPI"

    [deps.NCDatasets.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.NetCDF_jll]]
deps = ["Artifacts", "Blosc_jll", "Bzip2_jll", "HDF5_jll", "JLLWrappers", "LazyArtifacts", "LibCURL_jll", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "TOML", "XML2_jll", "Zlib_jll", "Zstd_jll", "libaec_jll", "libzip_jll"]
git-tree-sha1 = "d574803b6055116af212434460adf654ce98e345"
uuid = "7243133f-43d8-5620-bbf4-c2c921802cf3"
version = "401.900.300+0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML", "Zlib_jll"]
git-tree-sha1 = "6d6c0ca4824268c1a7dca1f4721c535ac63d9074"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "5.0.11+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "94ba93778373a53bfd5a0caaf7d809c445292ff4"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.2"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "468dbe2b510c876dc091b2c74ed52c7c34f48b9b"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.5"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3151a0c8061cc3f887019beebf359e6c4b3daa08"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.65"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Profile]]
deps = ["StyledStrings"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.ProfileSVG]]
deps = ["Colors", "FlameGraphs", "Profile", "UUIDs"]
git-tree-sha1 = "95ef58783baaa61cd227342b7b605d8a4dfba4a9"
uuid = "132c30aa-f267-4189-9183-c8a63c7e05e6"
version = "0.2.2"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "72312aa278823c0e99ce31186e22d917d2d11f99"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.46"

[[deps.SciMLPublic]]
git-tree-sha1 = "0ba076dbdce87ba230fff48ca9bca62e1f345c9b"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "bb072715f158b59ad8819ff80da5ffa90cce6ceb"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.4.0"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "SciMLPublic", "Static"]
git-tree-sha1 = "2a635e15d5035c53b345077c947f31ff91744078"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.10.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "7c73336785b21f723f5b143f6e99cf6c43b37dc1"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.6"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "807a234dc5e6132dd6cf4c9317ca0917c4001ab3"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.74"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "58972370b81423fc546c56a60ed1a009450177c3"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.19.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60f4792734488db6f42e2c7699f1d4594780bd03"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.1.7+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libzip_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "OpenSSL_jll", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "86addc139bca85fdf9e7741e10977c45785727b7"
uuid = "337d8026-41b4-5cde-a456-74a10e5b31d1"
version = "1.11.3+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╟─1165eeb4-10d8-4fb5-859d-7fe410189608
# ╟─3f63ce64-effd-49c7-9906-92cda40d59f0
# ╟─0780e492-cdcd-43fb-af85-f9a4bab37450
# ╟─d0ebd34f-a83a-4dc1-9c66-24a3878b681d
# ╟─b303e4e9-49fa-45ad-967e-20f165fdf38c
# ╟─f578f25e-047e-4a7e-8483-d544c7b4bec3
# ╟─2bf0fe8e-5718-4c1e-863b-85db7b3ae7f3
# ╟─8578d6aa-2782-4279-8f6b-78194b8ecc10
# ╟─ebca5877-1305-40f1-ac52-66356dc17661
# ╟─0dbfb663-46e7-4873-ac77-1e8e392fe69d
# ╟─8165ed3f-ede0-45c9-83ac-e4e49262457c
# ╟─07a88c94-93bb-4564-8737-980144c4af43
# ╟─531589ab-c6b5-4048-9ba5-f9ad62ab00a6
# ╟─813d59f7-eeea-402a-96e8-e5cbe2fd4582
# ╟─dc0f1a5e-1e7b-4a9d-816e-731d2747ac4e
# ╟─19f106e4-2b82-47e1-9284-799a105f30cb
# ╟─bb82d1dc-8f08-4351-9371-a8efed1dd9bc
# ╟─54e63724-44ca-42b9-9246-7afb271b8154
# ╟─32b531ab-ee71-4685-af65-a2bae0b868f6
# ╟─d0d74213-96ad-4a45-9a31-37f640a21a45
# ╟─caec8065-9874-4d8c-8e7e-598a97506f06
# ╟─80ac789e-4fe7-4184-9946-b8d7c24b04ea
# ╟─b9f7bde9-0aa4-4075-8fb9-14f84db0b3fa
# ╟─12a89062-bc27-4e36-8b5e-1670b4bd0b28
# ╟─f39520b7-a246-4980-b8b9-0215367d0b46
# ╟─1879fd40-0117-478b-acdd-18112738b81f
# ╟─f8a2c2de-5045-4ab6-a6fa-7bca502afc9b
# ╟─b23bc922-e9bc-4012-9782-1258e3cc8e7b
# ╟─0308940f-16bc-4b03-8321-9ea8c91a21c7
# ╟─898dd5aa-5273-4833-9a4a-0f2b94cc8d38
# ╟─78605b5f-b03d-4851-b339-aea63bad3688
# ╟─75f3b78f-924a-4a65-b2e3-a79c6f2082f9
# ╟─0be9bc61-1a59-4dfe-84f9-bb1a27ca30fc
# ╟─0d200ce7-eadb-41de-8e31-45783a1faab9
# ╟─22d5e751-f095-42d7-9c24-78e546b3ce37
# ╟─047b312f-8d6c-4732-aa0b-bea3de3e99e2
# ╟─d403b01e-7d9f-4778-813d-bbbb0bfdbfb9
# ╟─c94558a5-0ae2-4563-93d7-0609ea3cce41
# ╟─c06deaa7-2e6a-4143-a195-b473cbd84329
# ╟─2d02d36f-efee-411e-befb-52d61b3fb16e
# ╟─a8b5fa01-0526-46f4-9aa0-31b52398a2bd
# ╟─1bdeea50-a73b-47a5-b941-d9d3374f54cf
# ╟─1b38b24d-bb82-4d52-a26d-850d142e25d5
# ╟─0ff4be7c-580c-4db1-a874-175badc8c11a
# ╟─7a06bf0d-a61c-4d28-b144-2725fe90ae62
# ╟─3beff4df-89da-4f25-a0e4-be03048c5f2c
# ╟─d404043f-8080-4262-9ab6-d9bb13eee504
# ╟─fd193b13-336f-42f8-8f93-3ca3bea2e620
# ╟─33fa7b1f-938b-481c-bf25-eca8d7fb33a7
# ╟─241496f6-580b-4e81-acbb-3295696fe40e
# ╟─e493fae7-239a-494c-9a59-728446d70f7a
# ╟─a4ab662c-5876-4861-aecf-2f695d15e15d
# ╟─1df2b91b-be14-427a-87b3-95cdef26ce00
# ╟─2a859d43-5768-4b98-88fa-978b12203513
# ╟─c0c40037-4169-4d38-bebe-2086cebc24f2
# ╟─0febe534-237e-4921-b39c-3828dbae9d19
# ╟─606032a2-b2ca-4fd8-9930-afd83aecec7a
# ╟─a1c04c52-f7c9-430a-8791-7fea15650b2c
# ╟─c6a0d656-8289-4f74-b8d8-f94c236e541d
# ╟─34fb479a-9834-4354-8f36-2680bedec798
# ╟─9bff59c5-4631-4091-8230-989a835788e5
# ╟─ae38f814-fe9c-443f-ae3b-42fa7a7d199a
# ╟─625089e2-ef77-4821-a6d6-d0a0f88207f2
# ╟─a4cd8d40-eedd-4a4b-ac17-51e68334328c
# ╟─db75ea52-9387-4c15-bde5-61777ac9b570
# ╟─cc7db228-0a06-4812-8b32-6541988b2115
# ╟─ba96178d-77d4-4f26-a94f-5ad43c5242db
# ╟─99f0e9ad-438f-4fa5-be9d-36d0fa78d89c
# ╟─2bab06b9-ca98-4142-99cb-d2ad4f1cde93
# ╟─8fdcfefd-0490-434a-a2bb-d171557b6ae7
# ╟─584e767e-4dc5-4821-af63-d6a825326d9e
# ╟─7e514484-9d81-4c5c-83ab-191d68c13043
# ╟─1894ad94-cdf8-4e79-a0e5-b72088db31be
# ╟─7b8e4130-7d8e-492a-a949-bd7fb6808dd6
# ╟─cc0e682c-8767-498e-8178-2de4e796b3a8
# ╟─7281dc60-e0e2-4a34-b4b6-c3f63a97f60f
# ╟─dfdde9f1-b226-4af0-9ac2-36f1b01622fa
# ╟─4f97badf-a501-4c70-a943-d3c86b48f8a1
# ╟─bf5ccbb7-fff9-42bc-8f2f-2e628b448b3a
# ╟─db95cf44-4d37-463c-8fa7-8f084c2a20e8
# ╟─4995d3d8-1f95-41d8-be6c-50663edbce10
# ╟─83e54812-2291-44e6-9b6e-0c0be57865d3
# ╟─867b193e-3390-4b3b-b5fe-5c399a250660
# ╟─92c3bd68-bd07-4381-9c04-e6611650cd1e
# ╟─d17c570f-ad56-4d43-98cb-ccd03a8fcb4a
# ╟─71676720-4685-4216-a55e-5918829d258f
# ╟─b7cef546-f0f4-4f42-a0ed-00ffb92a422e
# ╟─7969e897-121e-4fe0-9b75-36d21931357f
# ╟─dec9de54-eede-45ae-94cd-5bcc477298bf
# ╟─2d2aa153-1b51-4972-89a0-de0f5d803838
# ╟─66527a77-1926-41f7-b86f-84e0d65d9f30
# ╟─393a8bae-c5ee-4a46-85f3-2c6a7cead357
# ╟─506e51a4-d054-4e5d-8741-c4d36a9f2714
# ╟─3ac049a9-50e2-4db6-acb9-aacc20586ca4
# ╟─26b4d172-72ec-4dd6-83c0-0779904634ae
# ╟─2d4f533e-a240-459d-9575-f9982c03184b
# ╟─8ef38e96-42f4-4bb2-9063-88db9842fb47
# ╟─d5331962-6bb1-485d-b538-ff35221e14e6
# ╟─54e19b14-4bc1-4086-9015-ff2cd8f4afbf
# ╟─efb23f98-fb75-42e6-8681-c5147a82a09c
# ╟─d74882fa-3a38-4cc7-b476-e52cb5a7db31
# ╟─fee2639d-d8cf-4ee3-b824-73a0b02fef4a
# ╟─82bc0ee6-d0ca-4f47-a43f-67ae6702bb6a
# ╟─952d5e5d-119a-4ed9-9e22-0b0fe66ae04d
# ╟─28371449-d04e-4d64-8b4a-cd37ceb25ef5
# ╟─13647b0a-561c-4f65-a64c-b9a4ea76c9f0
# ╟─0ab27970-986e-4359-bed1-8908c5fd109b
# ╟─4d3bedd5-1d4a-4150-828f-c5352c70874b
# ╟─10e20296-cc4f-4c4a-80bf-c39ae6450e85
# ╟─9b7a847e-a0d2-4421-9fb9-ff01b73c4194
# ╟─cf2a51eb-a661-4c66-9e0a-d1730110e4bc
# ╟─f7fb4846-5ee8-4ca5-a76f-1edcd58a0559
# ╟─39afec69-f90d-4ff2-99d1-5cfec84d09a1
# ╟─b0d6a9e3-5c3a-4408-ab87-1067c2dfaeb7
# ╟─a9089bf1-3241-4038-81d6-b953a00c3cef
# ╟─698a990b-cd78-43a1-a747-e79102767d62
# ╟─9fa7c129-2289-4979-aa1f-8f8112e7875f
# ╟─caf443c0-5fbc-4fd4-8ef9-3565417d40ab
# ╟─1d46e4a8-41fc-4b76-bc9d-a115b6b6d9c4
# ╟─e20657c6-cf03-4f84-a782-2a3e5d6d303c
# ╟─e2227570-2fcb-4090-8bbb-ec85af3a519c
# ╟─1b94212a-a5b9-46dc-81ea-896e4c839755
# ╟─40717ac9-ff3f-4835-ba7e-3cde10b2c3f0
# ╟─fb0719e8-b9fe-4ca3-8c5a-de3586d193a9
# ╟─0c201c53-0db2-4805-beec-e7705082509d
# ╟─eb10715c-247d-40e4-b11b-aceabb60628c
# ╟─f9da0fed-a0a3-4bfb-8bac-6a839177bc73
# ╟─2e05f3f6-687b-4096-a869-d64564db60da
# ╟─09cbb86e-64e1-45e5-b0f6-1e8b868506e8
# ╟─c186464e-e705-49be-ad5d-4352ce719c1f
# ╟─395242dc-cb73-481f-ad25-71eba93d341f
# ╟─718e4f4b-af27-4271-a37f-83ea6e30bdbe
# ╟─59c97f2b-de11-4f42-9fcf-aaf2924d8677
# ╟─b9b3564e-c59b-4f3a-80d7-4431fd8ce2db
# ╟─435bae7c-e0d3-46af-b185-5221f687d6bf
# ╟─1229bf7d-204c-41a2-97dd-9d75e79b18ce
# ╟─426c1511-f912-4e10-9b8b-f858471019bc
# ╟─50777c8a-ee6f-449b-b704-bba7e9b36490
# ╟─99b335e7-b28a-44f2-a07d-a84c16366710
# ╟─d571c448-cea6-49c0-81a4-2376e4fc3f67
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
