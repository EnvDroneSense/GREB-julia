### A Pluto.jl notebook ###
# v0.20.21

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
	using Plots
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

Functions to read GrADS binary input files (.bin) and populate climate arrays.
"""

# ╔═╡ b303e4e9-49fa-45ad-967e-20f165fdf38c
"""
    read_grads_binary(filepath::String, nx::Int, ny::Int, nt::Int=1)

Read a GrADS binary file (sequential 32-bit floats, Fortran order).

# Arguments
- `filepath`: Full path to the .bin file
- `nx`: Number of longitude points (1 for solar data)
- `ny`: Number of latitude points
- `nt`: Number of time steps (default=1 for static fields)

# Returns
- Array{Float64} of size (nx, ny, nt) for 3D fields
- Array{Float64} of size (ny, nt) for 2D fields (when nx=1, e.g., solar)
- Array{Float64} of size (nx, ny) for static 2D fields (when nt=1)
"""
function read_grads_binary(filepath::String, nx::Int, ny::Int, nt::Int=1)
    # Calculate expected size
    nvals = nx * ny * nt
	expected_bytes = nvals * sizeof(Float32)

	if !isfile(filepath)
		error("Required input file not found: $filepath")
	end

	actual_bytes = filesize(filepath)
	if actual_bytes != expected_bytes
		error("Binary size mismatch for $filepath: expected $expected_bytes bytes, got $actual_bytes bytes")
	end
    
    # Read binary file as Float32 (GrADS default)
    data = open(filepath, "r") do io
        read!(io, Array{Float32}(undef, nvals))
    end
    
    # Reshape: GrADS stores as (lon, lat, time)
    if nt == 1 && nx > 1
        # Static 2D field (topography, glacier)
        return Array{Float64}(reshape(data, nx, ny))
    elseif nx == 1
        # Solar radiation case: (lat, time) only
        return Array{Float64}(reshape(data, ny, nt))
    else
        # Full 3D field (lon, lat, time)
        return Array{Float64}(reshape(data, nx, ny, nt))
    end
end

# ╔═╡ f578f25e-047e-4a7e-8483-d544c7b4bec3
"""
    load_flux_corrections!(input_dir::String)

Load flux correction fields conditionally (called only when needed by qflux_correction!).
Populates: TF_correct, qF_correct, ToF_correct with fallback to zeros if missing.
"""
function load_flux_corrections!(input_dir::String)
    println("  • Loading flux corrections...")
    
    print("    - Surface temperature... ")
    tf_path = joinpath(input_dir, "Tsurf_flux_correction.bin")
    if isfile(tf_path)
        global TF_correct .= read_grads_binary(tf_path, 96, 48, 730)
        println("✓")
    else
        global TF_correct .= 0.0
        @warn "Tsurf flux correction not found, using zeros" file=tf_path
        println("⚠")
    end
    
    print("    - Water vapor... ")
    qf_path = joinpath(input_dir, "vapour_flux_correction.bin")
    if isfile(qf_path)
        global qF_correct .= read_grads_binary(qf_path, 96, 48, 730)
        println("✓")
    else
        global qF_correct .= 0.0
        @warn "Vapour flux correction not found, using zeros" file=qf_path
        println("⚠")
    end
    
    print("    - Deep ocean... ")
    tof_path = joinpath(input_dir, "Tocean_flux_correction.bin")
    if isfile(tof_path)
        global ToF_correct .= read_grads_binary(tof_path, 96, 48, 730)
        println("✓")
    else
        global ToF_correct .= 0.0
        @warn "Tocean flux correction not found, using zeros" file=tof_path
        println("⚠")
    end
    
    return nothing
end

# ╔═╡ 2bf0fe8e-5718-4c1e-863b-85db7b3ae7f3
"""
    load_greb_input_data!(input_dir::String; dataset::Symbol=:ncep)

Load all required GREB input files from the specified directory and populate
the global climate field arrays in-place.

# Arguments
- `input_dir`: Path to the input/ directory containing .bin files
- `dataset`: Climate dataset to use (`:ncep` (default), `:era`, or `:era_ncep`)

# Modifies global arrays
Populates: `z_topo`, `glacier`, `Tclim`, `uclim`, `vclim`, `qclim`, `mldclim`,
`Toclim`, `cldclim`, `swetclim`, `sw_solar`, and MSCM fields (`omegaclim`,
`omegastdclim`, `wsclim`). Flux-correction fields (`TF_correct`, `qF_correct`,
`ToF_correct`) are loaded separately by `load_flux_corrections!()` when needed.
"""
function load_greb_input_data!(input_dir::String; dataset::Symbol=:ncep)
    println("═══════════════════════════════════════════════════════════")
    println("  Loading GREB Climate Data")
    println("  Directory: $input_dir")
    println("  Dataset: $dataset")
    println("═══════════════════════════════════════════════════════════")

	if !isdir(input_dir)
		error("Input directory not found: $input_dir")
	end
    
    # ─── Static 2D fields ───────────────────────────────────────────
    println("\n[1/4] Loading static 2D fields...")
    
    print("  • Topography... ")
    global z_topo .= read_grads_binary(joinpath(input_dir, "global.topography.bin"), 96, 48, 1)
    println("✓")
    
    print("  • Glacier mask... ")
    global glacier .= read_grads_binary(joinpath(input_dir, "greb.glaciers.bin"), 96, 48, 1)
    println("✓")
    
    # ─── 3D climatology fields (96 × 48 × 730) ──────────────────────
    println("\n[2/4] Loading 3D climatology fields (96×48×730)...")
    
    # Select dataset-specific files
    if dataset == :ncep
        tsurf_file = "ncep.tsurf.1948-2007.clim.bin"
        uwind_file = "ncep.zonal_wind.850hpa.clim.bin"
        vwind_file = "ncep.meridional_wind.850hpa.clim.bin"
        humid_file = "ncep.atmospheric_humidity.clim.bin"
        swet_file = "ncep.soil_moisture.clim.bin"
    elseif dataset == :era
        tsurf_file = "erainterim.tsurf.1979-2015.clim.bin"
        uwind_file = "erainterim.zonal_wind.850hpa.clim.bin"
        vwind_file = "erainterim.meridional_wind.850hpa.clim.bin"
        humid_file = "erainterim.atmospheric_humidity.clim.bin"
        swet_file = "ncep.soil_moisture.clim.bin"  # ERA doesn't have soil moisture
    else  # :era_ncep mixed
        tsurf_file = "erainterim.tsurf.1979-2015.clim.bin"
        uwind_file = "ncep.zonal_wind.850hpa.clim.bin"
        vwind_file = "ncep.meridional_wind.850hpa.clim.bin"
        humid_file = "ncep.atmospheric_humidity.clim.bin"
        swet_file = "ncep.soil_moisture.clim.bin"
    end
    
	# Strictly required for a baseline run.
	required_files = [
		"global.topography.bin",
		"greb.glaciers.bin",
		tsurf_file,
		uwind_file,
		vwind_file,
		humid_file,
		"isccp.cloud_cover.clim.bin",
		swet_file,
		"woce.ocean_mixed_layer_depth.clim.bin",
		"Tocean.clim.bin",
		"solar_radiation.clim.bin",
	]

	for req in required_files
		req_path = joinpath(input_dir, req)
		if !isfile(req_path)
			error("Required GREB input file missing: $req_path")
		end
	end

	print("  • Surface temperature ($tsurf_file)... ")
    global Tclim .= read_grads_binary(joinpath(input_dir, tsurf_file), 96, 48, 730)
    println("✓")
    
    print("  • Zonal wind (850 hPa)... ")
    global uclim .= read_grads_binary(joinpath(input_dir, uwind_file), 96, 48, 730)
    println("✓")
    
	print("  • Meridional wind (850 hPa)... ")
	global vclim .= read_grads_binary(joinpath(input_dir, vwind_file), 96, 48, 730)
	println("✓")
    
    print("  • Atmospheric humidity... ")
    global qclim .= read_grads_binary(joinpath(input_dir, humid_file), 96, 48, 730)
    println("✓")
    
    print("  • Cloud cover... ")
    global cldclim .= read_grads_binary(joinpath(input_dir, "isccp.cloud_cover.clim.bin"), 96, 48, 730)
    println("✓")
    
    print("  • Soil moisture... ")
    global swetclim .= read_grads_binary(joinpath(input_dir, swet_file), 96, 48, 730)
    println("✓")
    
    print("  • Ocean mixed-layer depth... ")
    global mldclim .= read_grads_binary(joinpath(input_dir, "woce.ocean_mixed_layer_depth.clim.bin"), 96, 48, 730)
    println("✓")
    
    print("  • Deep ocean temperature... ")
    global Toclim .= read_grads_binary(joinpath(input_dir, "Tocean.clim.bin"), 96, 48, 730)
    println("✓")
    
    # ─── MSCM-specific fields ────────────────────────────────────────
    println("\n[3/4] Loading MSCM-specific fields...")
    
    print("  • Vertical velocity (omega)... ")
    # Use ERA-Interim omega climatology for all datasets
    omega_file = "erainterim.omega.vertmean.clim.bin"
	if isfile(joinpath(input_dir, omega_file))
		global omegaclim .= read_grads_binary(joinpath(input_dir, omega_file), 96, 48, 730)
		println("✓")
	else
		global omegaclim .= 0.0
		@warn "Omega file not found, using zeros" file=joinpath(input_dir, omega_file)
		println("⚠")
	end
    
    print("  • Omega std deviation... ")
    # Use ERA-Interim omega std climatology for all datasets
    omegastd_file = "erainterim.omega_std.vertmean.clim.bin"
	if isfile(joinpath(input_dir, omegastd_file))
		global omegastdclim .= read_grads_binary(joinpath(input_dir, omegastd_file), 96, 48, 730)
		println("✓")
	else
		global omegastdclim .= 0.0
		@warn "Omega std file not found, using zeros" file=joinpath(input_dir, omegastd_file)
		println("⚠")
	end
    
	print("  • Wind speed... ")
	# Load ERA-Interim windspeed file for all datasets
	ws_file = "erainterim.windspeed.850hpa.clim.bin"
	ws_path = joinpath(input_dir, ws_file)
	if isfile(ws_path)
		global wsclim .= read_grads_binary(ws_path, 96, 48, 730)
		println("✓")
	else
		global wsclim .= 0.0
		@warn "Windspeed file not found, using zeros" file=ws_path
		println("⚠")
	end
    
    # ─── Solar radiation ────────────────────────────────────────────────
    println("\n[4/4] Loading solar radiation...")
    # Note: Flux corrections (TF_correct, qF_correct, ToF_correct) are loaded conditionally
    # within greb_model!() when log_exp==1 && cfg.log_qflux_dmc and only as needed.
    
    print("  • Solar radiation (48×730)... ")
    global sw_solar .= read_grads_binary(joinpath(input_dir, "solar_radiation.clim.bin"), 1, 48, 730)
    println("✓")
    
    println("\n═══════════════════════════════════════════════════════════")
    println("  ✓ All climate data loaded successfully!")
    println("═══════════════════════════════════════════════════════════\n")
    
    return nothing
end

# ╔═╡ 8578d6aa-2782-4279-8f6b-78194b8ecc10
"""
    load_solar_forcing(input_dir::String, forcing_type::Symbol, index::Int=0)

Load optional solar forcing scenarios for paleoclimate and orbital forcing experiments.

# Arguments
- `input_dir`: Path to the input/ directory
- `forcing_type`: Type of forcing (`:eccentricity`, `:obliquity`, or `:paleo`)
- `index`: Index for the forcing scenario (e.g., 0-60 for eccentricity, 0-230 for obliquity)

# Returns
- Array{Float64} of size (48, 730) - solar radiation forcing pattern

# Solar Forcing Types
- `:eccentricity` - Eccentricity variations (0-60)
- `:obliquity` - Obliquity variations (0-230 in steps of 5)
- `:paleo` - Paleoclimate 231K Hybers corrected

# Example
```julia
sw_solar_forcing = load_solar_forcing(input_dir, :eccentricity, 10)
```
"""
function load_solar_forcing(input_dir::String, forcing_type::Symbol, index::Int=0)
    solar_dir = joinpath(input_dir, "solar_forcing_scenarios")
    
    if forcing_type == :eccentricity
        filename = "greb.solar.eccentricity.$index.bin"
        println("Loading eccentricity forcing scenario $index...")
    elseif forcing_type == :obliquity
        filename = "greb.solar.obliquity.$index.bin"
        println("Loading obliquity forcing scenario $index...")
    elseif forcing_type == :paleo
        filename = "greb.solar.231K_hybers.corrected.bin"
        println("Loading paleoclimate 231K solar forcing...")
    else
        error("Unknown solar forcing type: $forcing_type. Use :eccentricity, :obliquity, or :paleo")
    end
    
    filepath = joinpath(solar_dir, filename)
    
    if !isfile(filepath)
        error("Solar forcing file not found: $filepath")
    end
    
    # Solar forcing files are (1, 48, 730) but we read as (48, 730)
    solar_forcing = read_grads_binary(filepath, 1, 48, 730)
    println("  ✓ Solar forcing loaded: $(size(solar_forcing))")
    
    return solar_forcing
end

# ╔═╡ ebca5877-1305-40f1-ac52-66356dc17661
md"""
### Solar Forcing Storage

Global variable to store loaded solar forcing scenarios. This is populated on-demand
when running experiments that require orbital or paleoclimate forcing.

**Usage:**
```julia
# Load eccentricity forcing for experiment 36
global sw_solar_forcing_data = load_solar_forcing(input_data_dir, :eccentricity, 36)

# Load obliquity forcing for experiment 35
global sw_solar_forcing_data = load_solar_forcing(input_data_dir, :obliquity, 230)

# Load paleoclimate forcing for experiments 30-32
global sw_solar_forcing_data = load_solar_forcing(input_data_dir, :paleo)
```
"""

# ╔═╡ 0dbfb663-46e7-4873-ac77-1e8e392fe69d
begin
    # ☀️ Solar forcing storage (populated on-demand for orbital/paleo experiments)
    global sw_solar_forcing_data = nothing  # Will hold (48, 730) array when loaded
	global sw_solar_forcing_state = Ref(1.0)  # Runtime solar multiplier used by SWradiation!	
end;

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
- `log_deepocean_dmc`: Deep ocean coupling in mean climate

### 2×CO₂ Response Switches (`*_drsp`) 
- `log_clouds_drsp`: Cloud feedback in CO₂ response
- `log_vapor_drsp`: Water vapor feedback in CO₂ response
- `log_ice_drsp`: Ice-albedo feedback in CO₂ response
- `log_crcl_drsp`: Atmospheric circulation in CO₂ response
- `log_hydro_drsp`: Hydrological cycle in CO₂ response  
- `log_deepocean_drsp`: Deep ocean coupling in CO₂ response

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
	ndt_days = 24 * 3600 ÷ Δt                 # time steps per day
	nstep_yr = ndays_yr * ndt_days            # time steps per year (= 730)

	# 📅 Calendar constants ───────────────────────────────
	jday_mon = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
	jday_mon_cumsum = cumsum(jday_mon)
	
	# 🔢 Physical Constants & Numerical Limits ────────────
	MIN_TEMPERATURE_K = 100.0        # Absolute temperature floor (~-173°C)
	MAX_HUMIDITY_CHANGE = 0.020      # Maximum humidity increment [kg/kg]
	MIN_HUMIDITY_FRACTION = 0.9      # Fraction of humidity that can be removed
	
	# 🔍 Diagnostic print-out location ────────────────────
	ipx = 1          # longitude index for diagnostics
	ipy = 1          # latitude index for diagnostics

end;

# ╔═╡ 813d59f7-eeea-402a-96e8-e5cbe2fd4582
md"""
## 𝚿 `mo_physics` - Physical Parameters & Climate Fields

This module defines all physical constants, model parameters, emissivity coefficients, and declares the climate field arrays.
"""

# ╔═╡ b96052b9-4272-4d78-9e43-927bc872c43a
md"""
### Sensitivity Experiment Control
`log_exp` selects which physical processes are active. It replaces the Fortran `namelist /physics/` input.

| Value | Effect |
|:------|:-------|
| 0 | Full model (default) |
| 1 | Constant topography |
| ≤ 2 | + constant cloud cover |
| ≤ 3 | + constant water vapor |
| ≤ 4 | + no circulation |
| ≤ 5 | + no ice-albedo feedback |
| ≤ 6 | + no hydrology |
| 7 | No water vapor transport |
| 8 | No water vapor advection |
| ≤ 9 | No deep ocean |
| 11 | No deep ocean + linearized emissivity |
| 12, 13 | A1B CO₂ scenario |
| 14–16 | SST+1K experiments |
"""

# ╔═╡ 60f0f581-0bcf-40e7-9a45-86eaceba9040
md"""
### 🔬 Physics Switches

Experiment control (`log_exp`) and advanced MSCM-style physics switches for process deconstruction.
"""

# ╔═╡ 19f106e4-2b82-47e1-9284-799a105f30cb
Base.@kwdef struct PhysicsConfig
	log_clouds_dmc::Bool = true
	log_vapor_dmc::Bool = true
	log_ice_dmc::Bool = true
	log_crcl_dmc::Bool = true
	log_hydro_dmc::Bool = true
	log_deepocean_dmc::Bool = true
	log_atmos_dmc::Bool = true
	log_co2_dmc::Bool = true
	log_ocean_dmc::Bool = true
	log_qflux_dmc::Bool = true
	log_clouds_drsp::Bool = true
	log_vapor_drsp::Bool = true
	log_ice_drsp::Bool = true
	log_crcl_drsp::Bool = true
	log_hydro_drsp::Bool = true
	log_deepocean_drsp::Bool = true
	log_topo_drsp::Bool = true
	log_humid_drsp::Bool = true
	log_ice::Bool = true
	log_hdif::Bool = true
	log_hadv::Bool = true
	log_vdif::Bool = true
	log_vadv::Bool = true
	log_conv::Bool = true
	log_rain::Int = 0
	log_eva::Int = -1
	log_clim::Int = 0
	log_tsurf_ext::Bool = false
	log_hwind_ext::Bool = false
	log_omega_ext::Bool = false
end

# ╔═╡ 54e63724-44ca-42b9-9246-7afb271b8154
md"""
### Hydrology Parameters
Optimized parameter system with lookup table for different parameterization modes and climatology datasets.
"""

# ╔═╡ 32b531ab-ee71-4685-af65-a2bae0b868f6
begin
	# 💧 Optimized Hydrology Parameter Lookup Table ────────────────────────
	const HYDRO_PARAMS = (
		-1 => (1.0, 0.0, 0.0, 0.0),                    # Original GREB
		 1 => (-1.391649, 3.018774, 0.0, 0.0),        # +Relative humidity
		 2 => (0.862162, 0.0, -29.02096, 0.0),        # +Omega convergence
		 3 => (-0.2685845, 1.4591853, -26.9858807, 0.0), # +RH & Omega
		 0 => (-1.88, 2.25, -17.69, 59.07)            # Best GREB (ERA-Interim)
	)
	
	# 🗂️ Pre-allocated Workspace Arrays (reduces dynamic allocation) ─────────
	global temp_2d_1 = zeros(Float64, xdim, ydim)     # General 2D temporary array 1
	global temp_2d_2 = zeros(Float64, xdim, ydim)     # General 2D temporary array 2
	global temp_2d_3 = zeros(Float64, xdim, ydim)     # General 2D temporary array 3
	global temp_2d_4 = zeros(Float64, xdim, ydim)     # General 2D temporary array 4
	
	# 🎯 Cached Weight Arrays (avoid recomputation in diffusion/advection) ───
	global WZ_CACHE = Dict{Float64, Matrix{Float64}}()
	
	# ⚙️ Optimized Parameter Initialization ──────────────────────────────────
	function set_hydrology_parameters!(cfg::PhysicsConfig)
		global c_q, c_rq, c_omega, c_omegastd
		
		# Fast lookup instead of if-else chain
		params = get(HYDRO_PARAMS, cfg.log_rain, (1.0, 0.0, 0.0, 0.0))
		c_q, c_rq, c_omega, c_omegastd = params
		
		# Complete climatology adjustment
		if cfg.log_rain == 0 && cfg.log_clim == 1
			# NCEP parameter set (missing from original Julia version)
			c_q, c_rq, c_omega, c_omegastd = -1.27, 1.99, -16.54, 21.15
		end
		
		@info "⚙️ MSCM hydrology: log_rain=$(cfg.log_rain), log_clim=$(cfg.log_clim) → (c_q=$c_q, c_rq=$c_rq, c_omega=$c_omega, c_omegastd=$c_omegastd)"
	end
end;

# ╔═╡ caec8065-9874-4d8c-8e7e-598a97506f06
md"""
## 🔧 Optimization Structs

Workspace and accumulator structures to eliminate allocations and reduce global variable access.
"""

# ╔═╡ 80ac789e-4fe7-4184-9946-b8d7c24b04ea
begin
	# Circulation Workspace (eliminates allocations in diffusion/advection)
	mutable struct CirculationWorkspace
		T1h_buf::Vector{Float64}      # xdim buffer for polar sub-stepping
		dTxh_buf::Vector{Float64}     # xdim buffer for polar tendencies
		dX_diff::Matrix{Float64}      # xdim × ydim for diffusion output
		dX_adv::Matrix{Float64}       # xdim × ydim for advection output
		dTx::Matrix{Float64}          # xdim × ydim for zonal tendency
		dTy::Matrix{Float64}          # xdim × ydim for meridional tendency
		dX_crcl::Matrix{Float64}      # xdim × ydim for circulation output 
		temp_buf::Matrix{Float64}     # xdim × ydim general workspace
		# Pre-allocated zero buffers 
		zero_buf1::Matrix{Float64}    # xdim × ydim reusable zero buffer 1
		zero_buf2::Matrix{Float64}    # xdim × ydim reusable zero buffer 2
		zero_buf3::Matrix{Float64}    # xdim × ydim reusable zero buffer 3
		zero_buf4::Matrix{Float64}    # xdim × ydim reusable zero buffer 4
		# State buffers (eliminate allocations in time_loop!)
		Ts0_buf::Matrix{Float64}      # Surface temperature output
		Ta0_buf::Matrix{Float64}      # Air temperature output
		To0_buf::Matrix{Float64}      # Ocean temperature output
		q0_buf::Matrix{Float64}       # Humidity output
		a_surf_buf::Matrix{Float64}   # Surface albedo buffer
	end
	
	function CirculationWorkspace()
		CirculationWorkspace(
			zeros(Float64, xdim),
			zeros(Float64, xdim),
			zeros(Float64, xdim, ydim),
			zeros(Float64, xdim, ydim),
			zeros(Float64, xdim, ydim),
			zeros(Float64, xdim, ydim),
			zeros(Float64, xdim, ydim),
			zeros(Float64, xdim, ydim),
			# Zero buffers
			zeros(Float64, xdim, ydim),
			zeros(Float64, xdim, ydim),
			zeros(Float64, xdim, ydim),
			zeros(Float64, xdim, ydim),
			# State buffers
			zeros(Float64, xdim, ydim),
			zeros(Float64, xdim, ydim),
			zeros(Float64, xdim, ydim),
			zeros(Float64, xdim, ydim),
			zeros(Float64, xdim, ydim)
		)
	end
	
	# Global instance (reused across all time steps)
	global circ_workspace = CirculationWorkspace()
end;

# ╔═╡ b9f7bde9-0aa4-4075-8fb9-14f84db0b3fa
begin
	# Monthly Accumulator (eliminates global variable access)
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
	const_pi   = pi                       # π (model precision)
	σ          = 5.6704e-8                # Stefan-Boltzmann constant [W/m²/K⁴]
	ρ_ocean    = 999.1                    # density of water at T=15°C [kg/m³]
	ρ_land     = 2600.0                   # density of solid rock [kg/m³]
	ρ_air      = 1.2                      # density of air at 20°C at sea 
	grav       = 9.80665                  # gravitational acceleration [m/s²]
	cp_ocean   = 4186.0                   # specific heat of water at T=15°C [J/kg/K]
	cp_land    = cp_ocean / 4.5           # specific heat of dry land [J/kg/K]
	cp_air     = 1005.0                   # specific heat of air [J/kg/K]
	ε          = 1.0                      # emissivity for IR
end;

# ╔═╡ b23bc922-e9bc-4012-9782-1258e3cc8e7b
begin
	# ── Column depths ────────────────────────────────────────────────
	d_ocean   = 50.0                          # depth of ocean column [m]
	d_land    = 2.0                           # depth of land column [m]
	d_air     = 5000.0                        # depth of air column [m]

	# ── Heat capacities [J/K/m²] ────────────────────────────────────
	cap_ocean = cp_ocean * ρ_ocean            # heat capacity of 1m ocean [J/K/m²]
	cap_land  = cp_land * ρ_land * d_land     # heat capacity of land column [J/K/m²]
	cap_air   = cp_air * ρ_air * d_air        # heat capacity of air column [J/K/m²]

	# ── Sensible heat ───────────────────────────────────────────────
	ct_sens   = 22.5                          # sensible heat coupling [W/K/m²]

	# ── Albedo parameters ───────────────────────────────────────────
	da_ice    = 0.25                          # albedo increase for ice-covered points
	a_no_ice  = 0.1                           # albedo for non-ice-covered points
	a_cloud   = 0.35                          # cloud albedo

	# ── Ice/snow temperature thresholds [K] ─────────────────────────
	# Land: snow-albedo feedback range
	Tl_ice1   = 273.15 - 10.0                 # lower bound (full ice albedo)
	Tl_ice2   = 273.15                        # upper bound (no ice albedo)
	# Ocean: ice-albedo feedback range
	To_ice1   = 273.15 - 7.0                  # lower bound (full ice)
	To_ice2   = 273.15 - 1.7                  # upper bound (open ocean)

	# ── Deep ocean ──────────────────────────────────────────────────
	co_turb   = 5.0                           # turbulent mixing to deep ocean [W/K/m²]

	# ── Atmospheric transport ───────────────────────────────────────
	κ         = 8e5                           # atmos. diffusion coefficient [m²/s]

	# ── Latent heat / hydrology ─────────────────────────────────────
	ce        = 2e-3                          # latent heat transfer coefficient
	cq_latent = 2.257e6                       # latent heat of evaporation [J/kg]
	cq_rain   = -0.1 / 24.0 / 3600.0          # rain-related vapor decrease [1/s]

	# ── Scaling heights ─────────────────────────────────────────────
	z_air     = 8400.0                        # atmos. scaling height for heat & CO₂ [m]
	z_vapor   = 5000.0                        # scaling height for water vapor diffusion [m]

	# ── Regression factor ───────────────────────────────────────────
	r_qviwv   = 2.6736e3                      # VIWV ↔ q_air regression factor [kg/m³]

	# ── solar factor ────────────────────────────────────────────────
	S0_var = 100.0  						  # variation of solar constant [%], default 100%

	# ── Precomputed transport geometry & coefficients ───────────────
	deg_grid  = 2.0 * const_pi * 6.371e6 / 360.0
	dyy_grid  = dlat * deg_grid
	lat_grid  = [dlat * k - dlat / 2.0 - 90.0 for k in 1:ydim]
	dxlat_grid = [dlon * deg_grid * cos(2.0 * const_pi / 360.0 * lat_grid[k]) for k in 1:ydim]

	ccy_diff  = κ * Δt_crcl / dyy_grid^2
	ccx_diff  = [κ * Δt_crcl / dxlat_grid[k]^2 for k in 1:ydim]
	ccy_adv   = Δt_crcl / dyy_grid / 2.0
	ccx_adv   = [Δt_crcl / dxlat_grid[k] / 2.0 for k in 1:ydim]

	# ── Precomputed periodic longitude neighbour indices ────────────
	lon_jm1 = SVector{96}([j == 1 ? xdim : j - 1 for j in 1:xdim])
	lon_jp1 = SVector{96}([j == xdim ? 1 : j + 1 for j in 1:xdim])
	lon_jm2 = SVector{96}([j <= 2 ? xdim + (j - 2) : j - 2 for j in 1:xdim])
	lon_jp2 = SVector{96}([j >= xdim - 1 ? j + 2 - xdim : j + 2 for j in 1:xdim])
	lon_jm3 = SVector{96}([j <= 3 ? xdim + (j - 3) : j - 3 for j in 1:xdim])
	lon_jp3 = SVector{96}([j >= xdim - 2 ? j + 3 - xdim : j + 3 for j in 1:xdim])
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
	z_topo    = zeros(Float64, xdim, ydim)    # topography [m], <0 = ocean
	glacier   = zeros(Float64, xdim, ydim)    # glacier mask
	z_ocean   = zeros(Float64, xdim, ydim)    # derived ocean depth [m]
	cap_surf  = zeros(Float64, xdim, ydim)    # surface heat capacity [J/K/m²]
	wz_air    = zeros(Float64, xdim, ydim)    # exp(-z_topo / z_air)
	wz_vapor  = zeros(Float64, xdim, ydim)    # exp(-z_topo / z_vapor)
	nothing  # suppress output
end;

# ╔═╡ a291a9c2-b89a-4915-8a19-3f900647fbbb
begin
	# 🔧 3D derived/correction fields (xdim, ydim, nstep_yr) ────────
	TF_correct = zeros(Float64, xdim, ydim, nstep_yr)  # temperature flux correction [W/m²]
	qF_correct = zeros(Float64, xdim, ydim, nstep_yr)  # humidity flux correction
	ToF_correct= zeros(Float64, xdim, ydim, nstep_yr)  # deep-ocean temp flux correction [K/dt]
	dTrad      = zeros(Float64, xdim, ydim, nstep_yr)  # Tatmos-radiation offset
	
	nothing  # suppress output
end;

# ╔═╡ 0d200ce7-eadb-41de-8e31-45783a1faab9
begin
	# ☀️ 2D solar field (ydim, nstep_yr) ─────────────────────────────
	sw_solar   = zeros(Float64, ydim, nstep_yr)  # 24-hr mean solar radiation [W/m²]
	
	nothing  # suppress output
end;

# ╔═╡ 22d5e751-f095-42d7-9c24-78e546b3ce37
md"""
### Time-Step Indices
`jday` and `ityr` track the current calendar day and time-step within the year during integration. In Fortran these were module-level integers; here they start at 1.
"""

# ╔═╡ 047b312f-8d6c-4732-aa0b-bea3de3e99e2
begin
	# 🕐 Time State Struct (eliminates global variables & Pluto conflicts)
	mutable struct TimeState
		jday::Int  # Current calendar day in year [1..365]
		ityr::Int  # Current timestep in year [1..730]
	end

	# 📅 Precomputed Calendar Lookup (avoids mod operations) ────────────────
	const MAX_TIMESTEPS = 200 * nstep_yr  # Support up to 200-year runs
	CALENDAR_LOOKUP = [(
		day = mod((it - 1) ÷ ndt_days, ndays_yr) + 1,
		step = mod(it - 1, nstep_yr) + 1
	) for it in 1:MAX_TIMESTEPS]
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
	# ── Annual-mean diagnostic accumulators (xdim, ydim) ────────────
	Tsmn    = zeros(Float64, xdim, ydim)   # surface temperature
	Tamn    = zeros(Float64, xdim, ydim)   # air temperature
	Tomn    = zeros(Float64, xdim, ydim)   # deep ocean temperature
	qmn     = zeros(Float64, xdim, ydim)   # atmospheric humidity
	amn     = zeros(Float64, xdim, ydim)   # albedo
	
	# MSCM additional diagnostic accumulators
	icmn    = zeros(Float64, xdim, ydim)   # ice cover fraction
	prmn    = zeros(Float64, xdim, ydim)   # precipitation tendency
	evmn    = zeros(Float64, xdim, ydim)   # evaporation tendency
	qcrclmn = zeros(Float64, xdim, ydim)   # humidity circulation tendency
	
	swmn    = zeros(Float64, xdim, ydim)   # shortwave radiation
	lwmn    = zeros(Float64, xdim, ydim)   # longwave radiation
	qlatmn  = zeros(Float64, xdim, ydim)   # latent heat flux
	qsensmn = zeros(Float64, xdim, ydim)   # sensible heat flux
	ftmn    = zeros(Float64, xdim, ydim)   # temperature flux correction
	fqmn    = zeros(Float64, xdim, ydim)   # humidity flux correction
	
	nothing  # suppress output
end;

# ╔═╡ 2d02d36f-efee-411e-befb-52d61b3fb16e
md"""
### Monthly-Mean Output Buffers
These fields accumulate values within each calendar month. At month-end they are divided by the number of time steps in that month and stored (replacing the Fortran binary file writes).

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
	# ── Monthly-mean output buffers (xdim, ydim) ───────────────────
	Tmm     = zeros(Float64, xdim, ydim)   # surface temperature
	Tamm    = zeros(Float64, xdim, ydim)   # air temperature
	Tomm    = zeros(Float64, xdim, ydim)   # deep ocean temperature
	qmm     = zeros(Float64, xdim, ydim)   # atmospheric humidity
	apmm    = zeros(Float64, xdim, ydim)   # albedo
		
	# MSCM additional monthly buffers
	icmm     = zeros(Float64, xdim, ydim)   # ice cover fraction
	prmm     = zeros(Float64, xdim, ydim)   # precipitation tendency
	evmm     = zeros(Float64, xdim, ydim)   # evaporation tendency
	qcrclmm  = zeros(Float64, xdim, ydim)   # humidity circulation tendency
	swmm     = zeros(Float64, xdim, ydim)   # shortwave radiation
	lwmm     = zeros(Float64, xdim, ydim)   # longwave radiation
	qlatmm   = zeros(Float64, xdim, ydim)   # latent heat flux
	qsensmm  = zeros(Float64, xdim, ydim)   # sensible heat flux
	
	nothing  # suppress output
end;

# ╔═╡ 1bdeea50-a73b-47a5-b941-d9d3374f54cf
begin
	# ── Control run monthly storage ──────────────────────
	# Store monthly means from control run for comparison
	Tmn_ctrl    = zeros(Float64, xdim, ydim, 12)   # control surface temperature
	Tamn_ctrl   = zeros(Float64, xdim, ydim, 12)   # control air temperature
	Tomn_ctrl   = zeros(Float64, xdim, ydim, 12)   # control deep ocean temperature
	qmn_ctrl    = zeros(Float64, xdim, ydim, 12)   # control humidity
	icmn_ctrl   = zeros(Float64, xdim, ydim, 12)   # control ice cover
	prmn_ctrl   = zeros(Float64, xdim, ydim, 12)   # control precipitation
	evamn_ctrl  = zeros(Float64, xdim, ydim, 12)   # control evaporation
	qcrclmn_ctrl = zeros(Float64, xdim, ydim, 12)  # control humidity circulation
	
	nothing  # suppress output
end;

# ╔═╡ 1b38b24d-bb82-4d52-a26d-850d142e25d5
md"""
## 🎮 `greb_model` - Main Driver

The driver:

1. **Initialises** derived fields (`dTrad`, `z_ocean`, `cap_surf`) and applies sensitivity-experiment overrides
2. **Runs flux correction** spin-up via `qflux_correction!`
3. **Runs the control integration** at fixed CO₂
4. **Runs the scenario integration** with time-varying CO₂

Fortran binary file I/O (`open(21,...)`, `open(22,...)`) is replaced by in-memory
`Vector` of `NamedTuple` records — one per monthly snapshot.
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
MonthlyRecord = NamedTuple{(:Ts, :Ta, :To, :q, :albedo, :ice, :precip, :evap, :qcrcl, :sw, :lw, :qlat, :qsens), NTuple{13, Matrix{Float64}}}

# ╔═╡ f93e3556-5016-47c0-8ecb-42a9f922d787
md"""
### MSCM Helper Functions
Utility functions for area-weighted diagnostics and ice fraction computation.
"""

# ╔═╡ 5842dab4-b354-40a4-9057-a3752d525c10
begin
	# ── Ice fraction computation from temperature thresholds ────────
	function ice_fraction_field(Ts, ws::CirculationWorkspace=circ_workspace)
		ice = ws.temp_buf
		fill!(ice, 0.0)  # Clear workspace
		
		@inbounds for j in 1:ydim, i in 1:xdim
			if glacier[i, j] > 0.5
				ice[i, j] = 1.0
			elseif z_topo[i, j] >= 0.0  # land
				if Ts[i, j] <= Tl_ice1
					ice[i, j] = 1.0
				elseif Ts[i, j] >= Tl_ice2
					ice[i, j] = 0.0
				else
					ice[i, j] = 1.0 - (Ts[i, j] - Tl_ice1) / (Tl_ice2 - Tl_ice1)
				end
			else  # ocean
				if Ts[i, j] <= To_ice1
					ice[i, j] = 1.0
				elseif Ts[i, j] >= To_ice2
					ice[i, j] = 0.0
				else
					ice[i, j] = 1.0 - (Ts[i, j] - To_ice1) / (To_ice2 - To_ice1)
				end
			end
		end
		return ice
	end
end;

# ╔═╡ 3beff4df-89da-4f25-a0e4-be03048c5f2c
md"""
### `init_model!` - Initialisation & Sensitivity Overrides

Mutates the global climate arrays in-place:
- Initializes MSCM hydrology parameters via `set_hydrology_parameters!()`
- Computes `dTrad` (radiation offset)
- Derives `z_ocean` (3× max mixed-layer depth)
- Applies sensitivity experiment overrides based on `log_exp`
- Sets `cap_surf` (surface heat capacity) from topography
- Extracts initial conditions from the last time step of climatology
- Returns `(Ts_ini, Ta_ini, To_ini, q_ini, CO2_ctrl)` as a `NamedTuple`
"""

# ╔═╡ fd193b13-336f-42f8-8f93-3ca3bea2e620
md"""
## 🔄 `time_loop` - Single Time-Step Integration

Each call advances the model by one time step `Δt`:

1. Compute calendar indices `jday` and `ityr`
2. Call `tendencies!` to get all physics tendencies
3. Update **surface temperature** `Ts0` from SW, LW, latent/sensible heat, flux correction
4. Update **air temperature** `Ta0` from LW balance, latent heat, circulation
5. Update **deep ocean temperature** `To0` from turbulent mixing + flux correction
6. Update **humidity** `q0` from evaporation, rain, circulation + flux correction (clamped ≥ 0)
7. Call `seaice!` to adjust heat capacity for sea ice
8. Call `output!` and `diagnostics!`
"""

# ╔═╡ 241496f6-580b-4e81-acbb-3295696fe40e
md"""
## ⚔️ `tendencies` - Physics Aggregation

Aggregates all physics tendencies for one time step.
"""

# ╔═╡ a4ab662c-5876-4861-aecf-2f695d15e15d
md"""
## ☀️ `SWradiation` - Short-Wave Radiation

Computes surface albedo from ice-cover thresholds (separate for land and ocean), atmospheric albedo from cloud climatology, and the net downward short-wave flux.

**Key logic:**
- Land surface albedo ramps linearly between `Tl_ice1` (full ice) and `Tl_ice2` (ice-free)
- Ocean surface albedo ramps between `To_ice1` and `To_ice2`
- Glaciers always get full ice albedo
- Sensitivity experiments (`log_exp ≤ 5`) override surface albedo to `a_no_ice`
- Combined albedo: `α = α_surf + α_atmos − α_surf x α_atmos`
- SW flux broadcast: each longitude row gets `SW_solar[:,ityr] x (1 − α[i,:])`
"""

# ╔═╡ 1df2b91b-be14-427a-87b3-95cdef26ce00
function SWradiation!(Ts1, timestate::TimeState, log_exp::Int, cfg::PhysicsConfig,                         ws::CirculationWorkspace)

	# ── Atmospheric albedo from cloud climatology ──────────────
	cld = @view cldclim[:, :, timestate.ityr]
	a_atmos = @. cld * a_cloud

	# ── Surface albedo ─────────────────────────────────────────
	a_surf = ws.a_surf_buf

	# Land: albedo = f(Ts) with ice thresholds Tl_ice1 / Tl_ice2
	# Use static scheduling for deterministic results
	Threads.@threads :static for j in 1:ydim
		for i in 1:xdim
			if z_topo[i,j] >= 0.0                         # land
				if Ts1[i,j] <= Tl_ice1
					a_surf[i,j] = a_no_ice + da_ice             # full ice
				elseif Ts1[i,j] >= Tl_ice2
					a_surf[i,j] = a_no_ice                      # ice-free
				else
					a_surf[i,j] = a_no_ice + da_ice * (1.0 - (Ts1[i,j] - Tl_ice1) / (Tl_ice2 - Tl_ice1))
				end
			else                                              # ocean
				if Ts1[i,j] <= To_ice1
					a_surf[i,j] = a_no_ice + da_ice             # full ice
				elseif Ts1[i,j] >= To_ice2
					a_surf[i,j] = a_no_ice                      # ice-free
				else
					a_surf[i,j] = a_no_ice + da_ice * (1.0 - (Ts1[i,j] - To_ice1) / (To_ice2 - To_ice1))
				end
			end
		end
	end

	# Glacier cells always get full ice albedo
	@. a_surf = ifelse(glacier > 0.5, a_no_ice + da_ice, a_surf)

	# Sensitivity experiments: no ice-albedo feedback
	if log_exp <= 5
		a_surf .= a_no_ice
	end

	# SWradiation applies log_ice as the direct ice-albedo gate.
	if !cfg.log_ice
		a_surf .= a_no_ice
	end

	# ── Combined albedo & SW flux ──────────────────────────────
	albedo = @. a_surf + a_atmos - a_surf * a_atmos

	global temp_2d_4
	sw = temp_2d_4
	@simd for i in 1:xdim
		# Apply solar forcing for experiments (e.g., 27-28: solar constant changes)
		solar_flux = sw_solar[:, timestate.ityr]
		solar_flux = solar_flux .* sw_solar_forcing_state[] .* (0.01 * S0_var)
		@views sw[i, :] .= solar_flux .* (1.0 .- albedo[i, :])
	end

	return (SW = sw, albedo = albedo)
end

# ╔═╡ 2a859d43-5768-4b98-88fa-978b12203513
md"""
## 🌙 `LWradiation` - Long-Wave Radiation

Computes emissivity from CO₂, water-vapour column and cloud cover, then derives the upward/downward long-wave fluxes.

**Key logic:**
- Effective CO₂ and vapour columns are pressure-scaled via `exp(-z_topo / z_air)`
- Emissivity `em` is a log-regression of CO₂ + vapour + cloud adjustment  
  (10-parameter vector `p_emi`)
- Sensitivity experiment `log_exp == 11`: vapour from climatology + linear correction
- `LW_surf = -σ x Ts⁴`
- `LWair_down = -em x σ x (Ta + dTrad)⁴` and `LWair_up = LWair_down`
"""

# ╔═╡ c0c40037-4169-4d38-bebe-2086cebc24f2
function LWradiation!(Ts1, Ta1, q1, CO2, timestate::TimeState, log_exp::Int, cfg::PhysicsConfig, ws::CirculationWorkspace)
	
	# ── Effective column amounts (pressure-scaled) ─────────────
	e_co2   = zeros(Float64, xdim, ydim)       # CO₂ (spatial)
	e_vapor = @. wz_air * r_qviwv * q1         # water vapour
	e_cloud = @view cldclim[:, :, timestate.ityr]        # clouds

	# Apply spatial CO₂ masking for regional experiments
	# Spatial threading (grid-independent)
	Threads.@threads :static for j in 1:ydim
		for i in 1:xdim
			e_co2[i, j] = wz_air[i, j] * CO2 * co2_part[i, j]
		end
	end

	# ── Emissivity (log-regression with 10 parameters) ─────────
	em = @. p_emi[4] * log(p_emi[1] * e_co2 + p_emi[2] * e_vapor + p_emi[3]) +
		    p_emi[7] +
		    p_emi[5] * log(p_emi[1] * e_co2   + p_emi[3]) +
		    p_emi[6] * log(p_emi[2] * e_vapor + p_emi[3])

	# Cloud adjustment (in-place)
	@. em = (p_emi[8] - e_cloud) / p_emi[9] * (em - p_emi[10]) + p_emi[10]

	# ── Fluxes ──────────────────────────────────────────────────
	LW_surf    = @. -σ * Ts1^4
	dTr        = @view dTrad[:, :, timestate.ityr]
	LWair_down = @. -em * σ * (Ta1 + dTr)^4
	LWair_up   = LWair_down

	# Mean-state atmosphere deconstruction parity with Fortran
	if !cfg.log_atmos_dmc
		LWair_down .= 0.0
	end

	return (LW_surf = LW_surf, LWair_up = LWair_up, LWair_down = LWair_down, em = em)
end

# ╔═╡ 0febe534-237e-4921-b39c-3828dbae9d19
md"""
## 💧 `hydro` - Hydrological Cycle

Computes latent heat fluxes and water-vapour tendencies (evaporation and rain).

**Key logic:**
- Early return with zeros when `log_exp ≤ 6`, `== 13`, or `== 15` (no hydrology feedback)
- Absolute wind speed from `(u,v)` climatology + gustiness (2 m/s land, 3 m/s ocean)
- Saturated humidity: Clausius–Clapeyron with topo-pressure scaling
- Latent heat: `Q_lat = (q − q_s) x |wind| x cq_latent x ρ_air x ce x swet`
- Evaporation: `dq_eva = -Q_lat / (cq_latent x r_qviwv)`
- Rain: `dq_rain = cq_rain x q`
- Atmospheric latent heat: `Q_lat_air = -dq_rain x cq_latent x r_qviwv`
"""

# ╔═╡ 606032a2-b2ca-4fd8-9930-afd83aecec7a
function hydro!(Ts1, q1, timestate::TimeState, log_exp::Int, cfg::PhysicsConfig, ws::CirculationWorkspace=circ_workspace)

	# Use pre-allocated zero buffers 
	fill!(ws.zero_buf1, 0.0)
	fill!(ws.zero_buf2, 0.0)
	fill!(ws.zero_buf3, 0.0)
	fill!(ws.zero_buf4, 0.0)

	# Process control switches
	!cfg.log_atmos_dmc && return (Q_lat = ws.zero_buf1, Q_lat_air = ws.zero_buf2, dq_eva = ws.zero_buf3, dq_rain = ws.zero_buf4)
	!cfg.log_hydro_dmc && return (Q_lat = ws.zero_buf1, Q_lat_air = ws.zero_buf2, dq_eva = ws.zero_buf3, dq_rain = ws.zero_buf4)
	!cfg.log_hydro_drsp && return (Q_lat = ws.zero_buf1, Q_lat_air = ws.zero_buf2, dq_eva = ws.zero_buf3, dq_rain = ws.zero_buf4)
	
	# Sensitivity experiments with no hydrological feedback
	if log_exp <= 6 || log_exp == 13 || log_exp == 15
		return (Q_lat = ws.zero_buf1, Q_lat_air = ws.zero_buf2, dq_eva = ws.zero_buf3, dq_rain = ws.zero_buf4)
	end

	# ── Get climatology fields ──────────────────────────
	u = @view uclim[:, :, timestate.ityr]
	v = @view vclim[:, :, timestate.ityr]
	swet = @view swetclim[:, :, timestate.ityr]
	omega = @view omegaclim[:, :, timestate.ityr]
	omegastd = @view omegastdclim[:, :, timestate.ityr]

	# ── Saturated humidity (Clausius–Clapeyron + topo-pressure) ─
	qs = @. 3.75e-3 * exp(17.08085 * (Ts1 - 273.15) / (Ts1 - 273.15 + 234.175)) *
		wz_air

	# ── Enhanced Evaporation Modes ──────────────────────────
	if cfg.log_eva == -1  # Original GREB model
		# +2 m/s land gustiness, +3 m/s ocean gustiness, uniform ce.
		wind_eff = @. sqrt(u^2 + v^2)
		@. wind_eff = ifelse(z_topo > 0.0,
			sqrt(wind_eff^2 + 2.0^2),
			sqrt(wind_eff^2 + 3.0^2))
		Q_lat = @. (q1 - qs) * wind_eff * cq_latent * ρ_air * ce * swet
		
	elseif cfg.log_eva == 0  # Skin temperature version
		# Enhanced with skin temperature offset and wsclim integration
		Tskin_offset = @. (z_topo > 0.5) ? 5.0 : 1.0  # Land +5K, Ocean +1K   
		Tskin = @. max(Ts1 + Tskin_offset, 200.0)  # Prevent unrealistic temps
		
		# Recalculate saturation with skin temperature
		qs_skin = @. 3.75e-3 * exp(17.08085 * (Tskin - 273.15) / (Tskin - 273.15 + 234.175)) * wz_air
		
		# Use wsclim if available, otherwise computed wind
		ws_view = @view wsclim[:, :, timestate.ityr]
		if maximum(abs, ws_view) > 1e-10  # Check if real data
			ws_base = ws_view
		else
			ws_base = @. sqrt(u^2 + v^2)
		end
		gustiness_land = 11.5  
		gustiness_ocean = 5.4  
		wind_eff = @. (z_topo > 0.5) ? sqrt(ws_base^2 + gustiness_land^2) : sqrt(ws_base^2 + gustiness_ocean^2)

		# Mode 0 coefficients 
		cE_land = 0.25 * ce  
		cE_ocean = 0.58 * ce 
		cE = @. (z_topo > 0.5) ? cE_land : cE_ocean 
		Q_lat = @. cE * wind_eff * ρ_air * cq_latent * (qs_skin - q1) * swet
		
	elseif cfg.log_eva == 1  # Enhanced land/ocean coefficients 
		# Improved parameterization with different land/ocean physics
		gustiness_land = 144.0   
		gustiness_ocean = 7.1    
		wind_eff = @. (z_topo > 0.5) ? sqrt(u^2 + v^2 + gustiness_land^2) : sqrt(u^2 + v^2 + gustiness_ocean^2)
		
		# Enhanced coefficients from MSCM calibration
		cE_land = 0.04 * ce   
		cE_ocean = 0.73 * ce  
		cE = @. (z_topo > 0.5) ? cE_land : cE_ocean
		Q_lat = @. cE * wind_eff * ρ_air * cq_latent * (qs - q1) * swet
		
	elseif cfg.log_eva == 2  # Wind climatology version  
		# Always uses wsclim instead of computed wind components
		ws_view = @view wsclim[:, :, timestate.ityr]
		if maximum(abs, ws_view) < 1e-10
			@warn "Mode 2 requires wsclim data. Falling back to computed wind."
			ws_base = @. sqrt(u^2 + v^2)
		else
			ws_base = ws_view  # Use wind speed climatology
		end
		
		# gustiness values: land=9.0, ocean=4.0
		gustiness_land = 9.0    
		gustiness_ocean = 4.0    
		wind_eff = @. (z_topo > 0.5) ? sqrt(ws_base^2 + gustiness_land^2) : sqrt(ws_base^2 + gustiness_ocean^2)
		
		# Mode 2 coefficients (wsclim-optimized)
		cE_land = 0.56 * ce   
		cE_ocean = 0.79 * ce  
		cE = @. (z_topo > 0.5) ? cE_land : cE_ocean 
		Q_lat = @. cE * wind_eff * ρ_air * cq_latent * (qs - q1) * swet
	else
		@warn "Unknown log_eva mode: $(cfg.log_eva). Using original GREB."
		# Fallback to original GREB
		gustiness = @. (z_topo > 0.5) ? 3.0 : 10.0
		wind_eff = @. sqrt(u^2 + v^2 + gustiness^2)
		cE = @. (z_topo > 0.5) ? 0.65e-3 : 1.15e-3
		Q_lat = @. cE * wind_eff * ρ_air * cq_latent * (qs - q1) * swet
	end

	# ── Precipitation system (Eq. 11 in Stassen et al 2019) ────────────────────────
	if cfg.log_rain != -1
		rq = @. q1 / qs
		dq_rain = @. (c_q + c_rq * rq + c_omega * omega + c_omegastd * omegastd) * cq_rain * q1

		if cfg.log_rain == 1
			limit_value = @. -0.0015 / (wz_vapor * r_qviwv * 86400.0)
			# Fortran parity: where dq_rain >= limit, set to limit.
			dq_rain = @. min(dq_rain, limit_value)
		end
	else
		# Original GREB precipitation (log_rain == -1)
		dq_rain = @. cq_rain * q1
	end

	# ── Water-vapour tendencies ─────────────────────────────
	dq_eva  = @. -Q_lat / cq_latent / r_qviwv        # evaporation
	
	# Clamp precipitation to avoid negative humidity
	dq_rain = @. max(dq_rain, -0.9 * q1 / Δt)

	# ── Atmospheric latent heat ────────────────────────────
	Q_lat_air = @. -dq_rain * cq_latent * r_qviwv

	return (Q_lat = Q_lat, Q_lat_air = Q_lat_air, dq_eva = dq_eva, dq_rain = dq_rain)
end;

# ╔═╡ a1c04c52-f7c9-430a-8791-7fea15650b2c
md"""
## 🌫️ `convergence` - Moisture Flux Convergence

Calculates divergence (convergence) of humidity field using vertical velocity (omega).
Implements Eq. 18 from Stassen et al 2019 for MSCM moisture transport.
"""

# ╔═╡ c6a0d656-8289-4f74-b8d8-f94c236e541d
function convergence!(T1, omegaclim, ityr)
	"""Calculate moisture flux convergence using omega vertical velocity.
	
	Based on Fortran subroutine in greb.model.mscm.f90 lines 1238-1262.
	Implements Eq. 18 from Stassen et al 2019.
	
	Args:
		T1: Input field (typically specific humidity) [kg/kg]
		omegaclim: Vertical velocity climatology [Pa/s] 
		ityr: Time step index within year
		
	Returns:
		div: Convergence field [kg/kg/s]
	"""
	div = zeros(Float64, xdim, ydim)
	
	for j in 1:ydim
		for i in 1:xdim
			# Vertical velocity omega (Pa/s) to m/s
			w = -omegaclim[i,j,ityr] / (ρ_air * grav)
			# Convergence calculation
			div[i,j] = T1[i,j] * w * Δt_crcl / z_vapor * 2.5
		end
	end
	
	return div
end

# ╔═╡ 34fb479a-9834-4354-8f36-2680bedec798
md"""
## ❄️ `seaice` - Sea-Ice Heat Capacity

Adjusts the surface heat capacity `cap_surf` over ocean grid cells depending on sea-ice state.

**Key logic (mutates `cap_surf` in-place):**
- Ocean + fully frozen (`Ts ≤ To_ice1`): `cap_surf = cap_land`
- Ocean + ice-free (`Ts ≥ To_ice2`): `cap_surf = cap_ocean x mld`
- Ocean + partial ice: linear ramp between the two
- Sensitivity experiments (`log_exp ≤ 5`): no ice feedback — land gets `cap_land`, ocean gets `cap_ocean x mld`
- Glacier cells always get `cap_land`
"""

# ╔═╡ ae38f814-fe9c-443f-ae3b-42fa7a7d199a
md"""
## 🌊 `deep_ocean` - Deep Ocean Coupling

Computes mixed-layer ↔ deep-ocean heat exchange via entrainment/detrainment and turbulent mixing.

**Key logic:**
- Early return with zeros when `log_exp ≤ 9`, `== 11`, or `14–16`
- `dmld`: change in mixed-layer depth between consecutive time steps
- Entrainment (`dmld < 0`): deep ocean warms as mixed layer shallows
- Detrainment (`dmld > 0`): surface cools as mixed layer deepens
- Both scaled by `c_effmix = 0.5`
- Turbulent mixing: `co_turb x (T_x − To)` with `T_x = max(To_ice2, Ts)`
"""

# ╔═╡ 625089e2-ef77-4821-a6d6-d0a0f88207f2
function deep_ocean!(Ts1, To1, timestate::TimeState, log_exp::Int, cfg::PhysicsConfig, ws::CirculationWorkspace=circ_workspace)
	
	# Use pre-allocated zero buffers 
	fill!(ws.zero_buf1, 0.0)
	fill!(ws.zero_buf2, 0.0)
	dT_ocean = ws.zero_buf1
	dTo = ws.zero_buf2
	
	# Process control switches
	!cfg.log_deepocean_dmc && return (dT_ocean = dT_ocean, dTo = dTo)
	!cfg.log_deepocean_drsp && return (dT_ocean = dT_ocean, dTo = dTo)

	# Sensitivity experiments with no deep-ocean coupling
	if log_exp <= 9 || log_exp == 11
		return (dT_ocean = dT_ocean, dTo = dTo)
	end
	if 14 <= log_exp && log_exp <= 16
		return (dT_ocean = dT_ocean, dTo = dTo)
	end

	# ── Change in mixed-layer depth ─────────────────────────
	mld_now  = @view mldclim[:, :, timestate.ityr]
	mld_prev = timestate.ityr > 1 ? (@view mldclim[:, :, timestate.ityr - 1]) : (@view mldclim[:, :, nstep_yr])
	dmld = @. mld_now - mld_prev

	# ── Entrainment & detrainment ─────────────────────────
	# Thread over grid points (independent)
	Threads.@threads :static for j in 1:ydim
		for i in 1:xdim
			z_topo[i,j] >= 0.0 && continue                   # skip land
			Ts1[i,j] < To_ice2 && continue                   # skip ice-covered
			if dmld[i,j] < 0.0
				# Entrainment: mixed layer shallows → deep ocean warms
				dTo[i,j] = -dmld[i,j] / (z_ocean[i,j] - mld_now[i,j]) * (Ts1[i,j] - To1[i,j])
			elseif dmld[i,j] > 0.0
				# Detrainment: mixed layer deepens → surface cools
				dT_ocean[i,j] = dmld[i,j] / mld_now[i,j] * (To1[i,j] - Ts1[i,j])
			end
		end
	end

	# Efficiency factor
	c_effmix = 0.5
	dTo      .*= c_effmix
	dT_ocean .*= c_effmix

	# ── Turbulent mixing ─────────────────────────────────
	Threads.@threads :static for j in 1:ydim
		for i in 1:xdim
			z_topo[i,j] >= 0.0 && continue
			Tx = max(To_ice2, Ts1[i,j])
			dTo[i,j] += Δt * co_turb * (Tx - To1[i,j]) / (cap_ocean * (z_ocean[i,j] - mld_now[i,j]))
			dT_ocean[i,j] += Δt * co_turb * (To1[i,j] - Tx) / (cap_ocean * mld_now[i,j])
		end
	end

	return (dT_ocean = dT_ocean, dTo = dTo)
end

# ╔═╡ a4cd8d40-eedd-4a4b-ac17-51e68334328c
md"""
## 🌀 `circulation` - Atmospheric Circulation

Runs the diffusion + advection sub-stepping loop at the finer circulation time step `Δt_crcl`.

**Key logic:**
- Early return with zeros for `log_exp ≤ 4`, or when transporting water vapour
  (`h_scl == z_vapor`) under `log_exp == 7` or `== 16`
- Number of sub-steps: `time = max(1, round(Δt / Δt_crcl))`
- Each sub-step: `X₁ += diffusion!(X₁, h_scl) + advection!(X₁, h_scl)`
- Returns total change `dX_crcl = X₁ − X_in`
"""

# ╔═╡ cc7db228-0a06-4812-8b32-6541988b2115
md"""
## 🌀 `diffusion` - 3rd-Order Diffusion

Applies topography-weighted diffusion on a lat–lon grid with periodic longitude boundaries.

**Key logic:**
- Latitude diffusion: 2nd-order centred differences, zero-flux at poles
- Longitude at moderate latitudes (`dxlat > 250 km`): 3rd-order weighted stencil
  using ±1, ±2, ±3 neighbours with weights 10/4/1 (normalised by 20)
- Longitude near poles (`dxlat ≤ 250 km`): sub-time-stepping with the same stencil
  + stability clamp (`dTxh ≥ -0.9 T1h`)
- Final result scaled by `exp(-z_topo / h_scl)` (pressure weighting)
"""

# ╔═╡ ba96178d-77d4-4f26-a94f-5ad43c5242db
function diffusion!(T1::Matrix{Float64}, h_scl::Float64, ws::CirculationWorkspace, timestate::TimeState)::Matrix{Float64}
	# Use workspace buffers instead of temp_2d_*
	dX_diffuse = ws.dX_diff
	fill!(dX_diffuse, 0.0)  # Clear workspace

	# ── Precomputed geometry/coefficients ───────────────────
	dxlat = dxlat_grid
	ccy  = ccy_diff
	ccx  = ccx_diff
	# Use cached weight arrays
	wz = if h_scl == z_air
		wz_air
	elseif h_scl == z_vapor
		wz_vapor
	else
		get!(WZ_CACHE, h_scl) do
			@. exp(-z_topo / h_scl)
		end
	end

	dTx = ws.dTx
	dTy = ws.dTy
	fill!(dTx, 0.0)
	fill!(dTy, 0.0)

	@inbounds for k in 1:ydim
		# ── Latitudinal diffusion (2nd-order centred) ────────
		if k >= 2 && k <= ydim - 1
			@inbounds for j in 1:xdim
				dTy[j,k] = ccy * (wz[j,k-1] * (T1[j,k-1] - T1[j,k]) +
				                   wz[j,k+1] * (T1[j,k+1] - T1[j,k]))
			end
		elseif k == 1
			@inbounds for j in 1:xdim
				dTy[j,k] = ccy * wz[j,k+1] * (-T1[j,k] + T1[j,k+1])
			end
		else  # k == ydim
			@inbounds for j in 1:xdim
				dTy[j,k] = ccy * wz[j,k-1] * (T1[j,k-1] - T1[j,k])
			end
		end

		# ── Longitudinal diffusion ───────────────────────────
		if dxlat[k] > 2.5e5   # moderate latitudes: single-step 3rd-order
			@inbounds for j in 1:xdim
				jm1 = lon_jm1[j]; jp1 = lon_jp1[j]
				jm2 = lon_jm2[j]; jp2 = lon_jp2[j]
				jm3 = lon_jm3[j]; jp3 = lon_jp3[j]

				dTx[j,k] = ccx[k] * (
					10 * (wz[jm1,k] * (T1[jm1,k] - T1[j,k]) + wz[jp1,k] * (T1[jp1,k] - T1[j,k])) +
					 4 * (wz[jm2,k] * (T1[jm2,k] - T1[jm1,k]) + wz[jm1,k] * (T1[j,k] - T1[jm1,k])) +
					 4 * (wz[jp1,k] * (T1[j,k] - T1[jp1,k]) + wz[jp2,k] * (T1[jp2,k] - T1[jp1,k])) +
					 1 * (wz[jm3,k] * (T1[jm3,k] - T1[jm2,k]) + wz[jm2,k] * (T1[jm1,k] - T1[jm2,k])) +
					 1 * (wz[jp2,k] * (T1[jp1,k] - T1[jp2,k]) + wz[jp3,k] * (T1[jp3,k] - T1[jp2,k]))
				) / 20.0
			end
		else  # near-polar latitudes: sub-time-stepping
			dd     = max(1, round(Int, Δt_crcl / (1.0 * dxlat[k]^2 / κ)))
			dtdff2 = Δt_crcl / dd
			time2  = max(1, round(Int, Δt_crcl / dtdff2))
			ccx2   = κ * dtdff2 / dxlat[k]^2

			# Use workspace buffers (no allocation)
			copyto!(ws.T1h_buf, view(T1, :, k))  # Copy latitude row
			fill!(ws.dTxh_buf, 0.0)

			for _tt2 in 1:time2
				for j in 1:xdim
					jm1 = lon_jm1[j]; jp1 = lon_jp1[j]
					jm2 = lon_jm2[j]; jp2 = lon_jp2[j]
					jm3 = lon_jm3[j]; jp3 = lon_jp3[j]

					ws.dTxh_buf[j] = ccx2 * (
						10 * (wz[jm1,k] * (ws.T1h_buf[jm1] - ws.T1h_buf[j]) + wz[jp1,k] * (ws.T1h_buf[jp1] - ws.T1h_buf[j])) +
						 4 * (wz[jm2,k] * (ws.T1h_buf[jm2] - ws.T1h_buf[jm1]) + wz[jm1,k] * (ws.T1h_buf[j] - ws.T1h_buf[jm1])) +
						 4 * (wz[jp1,k] * (ws.T1h_buf[j] - ws.T1h_buf[jp1]) + wz[jp2,k] * (ws.T1h_buf[jp2] - ws.T1h_buf[jp1])) +
						 1 * (wz[jm3,k] * (ws.T1h_buf[jm3] - ws.T1h_buf[jm2]) + wz[jm2,k] * (ws.T1h_buf[jm1] - ws.T1h_buf[jm2])) +
						 1 * (wz[jp2,k] * (ws.T1h_buf[jp1] - ws.T1h_buf[jp2]) + wz[jp3,k] * (ws.T1h_buf[jp3] - ws.T1h_buf[jp2]))
					) / 20.0
				end
				# Stability clamp: no negative q
				@. ws.dTxh_buf = ifelse(ws.dTxh_buf <= -ws.T1h_buf, -0.9 * ws.T1h_buf, ws.dTxh_buf)
				ws.T1h_buf .+= ws.dTxh_buf
			end
			dTx[:, k] .= ws.T1h_buf .- T1[:, k]
		end
	end

	# Final result: weighted by pressure (reuse wz = exp(-z_topo/h_scl))
	@. dX_diffuse = wz * (dTx + dTy)
	return dX_diffuse
end

# ╔═╡ 99f0e9ad-438f-4fa5-be9d-36d0fa78d89c
md"""
## ➡️ `advection` - Upwind Advection

Upwind scheme for meridional and zonal transport driven by climatological wind fields `uclim`, `vclim`.

**Key logic:**
- Meridional (v): 2nd-order upwind with ±1/±2 neighbours; zero at pole boundaries
- Zonal at moderate latitudes (`dxlat > 250 km`): 2nd-order upwind (±1/±2)
- Zonal near poles: sub-time-stepping with 3rd-order upwind (weights 10/4/1, /20)
  + stability clamp
- Sensitivity experiment `log_exp == 8` with `h_scl == z_vapor`: early return
"""

# ╔═╡ 2bab06b9-ca98-4142-99cb-d2ad4f1cde93
function advection!(T1::Matrix{Float64}, h_scl::Float64, ws::CirculationWorkspace, timestate::TimeState, log_exp::Int)::Matrix{Float64}
	# Use workspace buffers
	dX_advec = ws.dX_adv
	fill!(dX_advec, 0.0)  # Clear workspace

	# Sensitivity experiment 8: no vapour advection
	(log_exp == 8 && h_scl == z_vapor) && return dX_advec

	# ── Precomputed geometry/coefficients ───────────────────
	dxlat = dxlat_grid
	ccy   = ccy_adv
	ccx   = ccx_adv
	# Use cached weight arrays
	wz = if h_scl == z_air
		wz_air
	elseif h_scl == z_vapor
		wz_vapor
	else
		get!(WZ_CACHE, h_scl) do
			@. exp(-z_topo / h_scl)
		end
	end

	dTx = ws.dTx
	dTy = ws.dTy
	fill!(dTx, 0.0)
	fill!(dTy, 0.0)

	@inbounds for k in 1:ydim
		# ── Meridional (v) advection ───────────────────────
		@inbounds for j in 1:xdim
			v = vclim[j, k, timestate.ityr]
			if v >= 0.0
				if k >= 3
					dTy[j,k] = -v * ccy * (wz[j,k-1] * (T1[j,k] - T1[j,k-1]) +
					                        wz[j,k-2] * (T1[j,k] - T1[j,k-2])) / 3.0
				elseif k == 2
					dTy[j,k] = -v * ccy * wz[j,k-1] * (T1[j,k] - T1[j,k-1])
				else  # k == 1
					dTy[j,k] = 0.0
				end
			else
				if k <= ydim - 2
					dTy[j,k] = v * ccy * (wz[j,k+1] * (T1[j,k] - T1[j,k+1]) +
					                       wz[j,k+2] * (T1[j,k] - T1[j,k+2])) / 3.0
				elseif k == ydim - 1
					dTy[j,k] = v * ccy * wz[j,k+1] * (-T1[j,k] + T1[j,k+1])
				else  # k == ydim
					dTy[j,k] = 0.0
				end
			end
		end

		# ── Zonal (u) advection ───────────────────────────
		if dxlat[k] > 2.5e5   # moderate latitudes
			@inbounds for j in 1:xdim
				jm1 = lon_jm1[j]; jp1 = lon_jp1[j]
				jm2 = lon_jm2[j]; jp2 = lon_jp2[j]

				u = uclim[j, k, timestate.ityr]
				if u >= 0.0
					dTx[j,k] = -u * ccx[k] * (wz[jm1,k] * (T1[j,k] - T1[jm1,k]) +
					                          wz[jm2,k] * (T1[j,k] - T1[jm2,k])) / 3.0
				else
					dTx[j,k] =  u * ccx[k] * (wz[jp1,k] * (T1[j,k] - T1[jp1,k]) +
					                          wz[jp2,k] * (T1[j,k] - T1[jp2,k])) / 3.0
				end
			end
		else  # near-polar: sub-time-stepping
			dd     = max(1, round(Int, Δt_crcl / (dxlat[k] / 10.0)))
			dtdff2 = Δt_crcl / dd
			time2  = max(1, round(Int, Δt_crcl / dtdff2))
			ccx2   = dtdff2 / dxlat[k] / 2.0

			# Use workspace buffers (no allocation)
			copyto!(ws.T1h_buf, view(T1, :, k))
			fill!(ws.dTxh_buf, 0.0)

			for _tt2 in 1:time2
				for j in 1:xdim
					jm1 = lon_jm1[j]; jp1 = lon_jp1[j]
					jm2 = lon_jm2[j]; jp2 = lon_jp2[j]
					jm3 = lon_jm3[j]; jp3 = lon_jp3[j]

					u = uclim[j, k, timestate.ityr]
					if u >= 0.0
						ws.dTxh_buf[j] = -u * ccx2 * (
							10 * wz[jm1,k] * (ws.T1h_buf[j] - ws.T1h_buf[jm1]) +
							 4 * wz[jm2,k] * (ws.T1h_buf[jm1] - ws.T1h_buf[jm2]) +
							 1 * wz[jm3,k] * (ws.T1h_buf[jm2] - ws.T1h_buf[jm3])
						) / 20.0
					else
						ws.dTxh_buf[j] =  u * ccx2 * (
							10 * wz[jp1,k] * (ws.T1h_buf[j] - ws.T1h_buf[jp1]) +
							 4 * wz[jp2,k] * (ws.T1h_buf[jp1] - ws.T1h_buf[jp2]) +
							 1 * wz[jp3,k] * (ws.T1h_buf[jp2] - ws.T1h_buf[jp3])
						) / 20.0
					end
				end
				# Stability clamp
				@. ws.dTxh_buf = ifelse(ws.dTxh_buf <= -ws.T1h_buf, -0.9 * ws.T1h_buf, ws.dTxh_buf)
				ws.T1h_buf .+= ws.dTxh_buf
			end
			dTx[:, k] .= ws.T1h_buf .- T1[:, k]
		end
	end

	@. dX_advec = dTx + dTy
	return dX_advec
end

# ╔═╡ 8fdcfefd-0490-434a-a2bb-d171557b6ae7
md"""
## 🔄 `qflux_correction` - Flux Correction

Runs a `time_flux`-year spin-up loop that computes the correction fields `TF_correct`, `ToF_correct`, and `qF_correct` so the model stays close to observed climatology.

**Key logic (mutates global arrays in-place):**
1. Loop over `time_flux × ndt_days × ndays_yr` steps
2. Compute tendencies, then surface/air/ocean temperature and humidity *without* correction
3. Derive correction = climatology − uncorrected state, store into `TF_correct[:,:,ityr]` etc.
4. Re-compute state *with* correction, call `seaice!` and `diagnostics!`
"""

# ╔═╡ 7e514484-9d81-4c5c-83ab-191d68c13043
md"""
## 🧪 `forcing` - Enhanced Forcing System

Comprehensive forcing dispatcher supporting 50+ experiments. Replaces simple `co2_level()`.

**CO₂ Scaling Experiments**:
- 20: 2×CO₂ (680 ppm), 21: 4×CO₂ (1360 ppm), 22: 10×CO₂ (3400 ppm)
- 23: 0.5×CO₂ (170 ppm), 24: 0×CO₂ (0 ppm)

**Time-Varying CO₂ Experiments**:
- 25: CO₂ sine wave (30-year period), 26: CO₂ step function (2×→1× at 1980)

**Solar Experiments**: 
- 27: +27 W/m² solar constant, 28: 11-year solar cycle

**Paleoclimate Experiments**:
- 30: Paleo solar + CO₂=200ppm, 31: Paleo solar + modern CO₂, 32: Modern solar + paleo CO₂

**Orbital Forcing**:
- 35: Obliquity changes, 36: Eccentricity changes, 37: Earth-Sun distance effects

**Scenario Experiments**:
- 12/13: Legacy A1B trajectory (1950–2100)
- 95: Enhanced A1B scenario mode
- 96-99: IPCC RCP scenarios (requires external files)
- 100: Custom CO₂ trajectory (requires file)

**Advanced Forcing**:
- 40-47: Regional/partial CO₂ experiments (NH/SH/tropics/ocean/land/seasonal)
- 230, 240-241: Forced boundary conditions (requires external data)
"""

# ╔═╡ 1894ad94-cdf8-4e79-a0e5-b72088db31be
function forcing(it, year, log_exp, icmn_ctrl=zeros(xdim,ydim,12); nstep_yr=nstep_yr)
	# Default CO₂ concentration
	CO2 = 680.0
	sw_solar_forcing = 1.0  # Solar forcing multiplier (ratio to reference)
	
	# 📜 Legacy experiments (preserve exact backward compatibility) ───────────
	if log_exp == 1
		CO2 = 550.0  # 550 ppm CO₂ steady state
	elseif log_exp == 12 || log_exp == 13
		CO2_1950 = 310.0;  CO2_2000 = 370.0;  CO2_2050 = 520.0
		if year <= 2000
			CO2 = CO2_1950 + 60.0 / 50.0 * (year - 1950)
		elseif year <= 2050
			CO2 = CO2_2000 + 150.0 / 50.0 * (year - 2000)
		elseif year <= 2100
			CO2 = CO2_2050 + 180.0 / 50.0 * (year - 2050)
		end
	
	# 💨 CO₂ scaling experiments ──────────────────────────────────────────────
	elseif log_exp == 20
		CO2 = 2 * 340.0  # 2×CO₂
	elseif log_exp == 21  
		CO2 = 4 * 340.0  # 4×CO₂
	elseif log_exp == 22
		CO2 = 10 * 340.0  # 10×CO₂ 
	elseif log_exp == 23
		CO2 = 0.5 * 340.0  # 0.5×CO₂
	elseif log_exp == 24
		CO2 = 0.0  # 0×CO₂ (no greenhouse effect)
	
	# ☀️ Solar forcing experiments ───────────────────────────────────────────
	elseif log_exp == 27
		CO2 = 340.0  # Baseline CO₂
		sw_solar_forcing = (1365.0 + 27.0) / 1365.0  # +27 W/m² solar constant
	elseif log_exp == 28
		CO2 = 340.0  # Baseline CO₂
		# 11-year solar cycle: ±1 W/m² amplitude
		sw_solar_forcing = (1365.0 + 1.0 * sin(2π * year / 11.0)) / 1365.0
	
	# 📈 Enhanced A1B scenario ──────────────────────────────────────────────
	elseif log_exp == 95
		CO2_1950 = 310.0;  CO2_2000 = 370.0;  CO2_2050 = 520.0
		if year <= 2000
			CO2 = CO2_1950 + 60.0 / 50.0 * (year - 1950)
		elseif year <= 2050
			CO2 = CO2_2000 + 150.0 / 50.0 * (year - 2000)  
		elseif year <= 2100
			CO2 = CO2_2050 + 180.0 / 50.0 * (year - 2050)
		end

		# ── Time-varying CO₂ experiments ────────────
	elseif log_exp == 25
		# CO₂ sine wave (30-year period)
		CO2 = 340.0 + 170.0 + 170.0 * cos(2π * (year - 13.0) / 30.0)
	elseif log_exp == 26
		# CO₂ step function: 2×CO₂ until 1980, then 1×CO₂
		CO2 = year >= 1980 ? 340.0 : 2 * 340.0
	
	# ── Paleoclimate experiments ────────────────────   
	elseif log_exp == 30
		# Paleo: 231kyr BP solar + CO₂=200ppm
		CO2 = 200.0
	elseif log_exp == 31  
		# Paleo: 231kyr BP solar + modern CO₂
		CO2 = 340.0
	elseif log_exp == 32
		# Paleo: modern solar + CO₂=200ppm (231kyr BP)
		CO2 = 200.0
	
	# ── Orbital forcing experiments ─────────────────
	elseif log_exp == 35
		# Solar obliquity changes
		CO2 = 340.0  # Baseline CO₂
	elseif log_exp == 36
		# Solar eccentricity changes  
		CO2 = 340.0  # Baseline CO₂
	elseif log_exp == 37
		# Solar constant as function of Earth-Sun distance
		CO2 = 340.0  # Baseline CO₂
	
	# 📂 File I/O dependent experiments (placeholders) ───────────────────────
	elseif 96 <= log_exp && log_exp <= 99
		# IPCC RCP scenarios - requires external CO₂ data files
		error("Experiment $log_exp (RCP scenario) requires external data file. Not yet implemented.")
	elseif log_exp == 100
		# Custom CO₂ scenario - requires external trajectory file  
		error("Experiment $log_exp (custom scenario) requires external data file. Not yet implemented.")

	# 🌍 Regional/partial CO₂ experiments ────────────────────────────────────
	elseif 40 <= log_exp && log_exp <= 47  
		# Reset co2_part to full CO₂ first
		co2_part .= 1.0

		if log_exp == 40
			# 2×CO₂ Northern Hemisphere only
			CO2 = 2 * 340.0
			co2_part[:, 1:24] .= 0.5
		elseif log_exp == 41
			# 2×CO₂ Southern Hemisphere only
			CO2 = 2 * 340.0
			co2_part[:, 25:48] .= 0.5
		elseif log_exp == 42
			# 2×CO₂ Tropics only
			CO2 = 2 * 340.0
			co2_part[:, 1:15] .= 0.5
			co2_part[:, 33:48] .= 0.5
			for i in 4:4:96
				co2_part[i, 33] = 1.0
				co2_part[i, 15] = 1.0
			end
		elseif log_exp == 43
			# 2×CO₂ Extratropics only
			CO2 = 2 * 340.0
			co2_part[:, 16:32] .= 0.5
			for i in 4:4:96
				co2_part[i, 32] = 1.0
				co2_part[i, 16] = 1.0
			end
		elseif log_exp == 44
			# 2×CO₂ Ocean only (land + perennial ice reduced)
			CO2 = 2 * 340.0
			for j in 1:ydim, i in 1:xdim
				if z_topo[i, j] > 0.0
					co2_part[i, j] = 0.5
				end
			end
			# Caller passes annual-mean ice in icmn_ctrl[:, :, 1] (replicated across 12 months).
			icmn_ctrl1 = @view icmn_ctrl[:, :, 1]
			for j in 1:ydim, i in 1:xdim
				if icmn_ctrl1[i, j] >= 0.5
					co2_part[i, j] = 0.5
				end
			end
		elseif log_exp == 45
			# 2×CO₂ Land/Ice only (ocean reduced, perennial ice restored to full)
			CO2 = 2 * 340.0
			for j in 1:ydim, i in 1:xdim
				if z_topo[i, j] <= 0.0
					co2_part[i, j] = 0.5
				end
			end
			# Caller passes annual-mean ice in icmn_ctrl[:, :, 1] (replicated across 12 months).
			icmn_ctrl1 = @view icmn_ctrl[:, :, 1]
			for j in 1:ydim, i in 1:xdim
				if icmn_ctrl1[i, j] >= 0.5
					co2_part[i, j] = 1.0
				end
			end
		elseif log_exp == 46
			# 2×CO₂ Boreal Winter only (seasonal)
			ityr_step = mod(it - 1, nstep_yr) + 1
			CO2 = (ityr_step <= 181 || ityr_step >= 547) ? 2 * 340.0 : 340.0
		elseif log_exp == 47
			# 2×CO₂ Boreal Summer only (seasonal)
			ityr_step = mod(it - 1, nstep_yr) + 1
			CO2 = (ityr_step <= 181 || ityr_step >= 547) ? 340.0 : 2 * 340.0
		end
	
	# 📊 Forced boundary condition experiments (handled in scenario loop) ───── 
	elseif log_exp == 230 || log_exp == 240 || log_exp == 241
		CO2 = 340.0
	end
	
	return (CO2 = CO2, sw_solar_forcing = sw_solar_forcing)
end

# ╔═╡ 7b8e4130-7d8e-492a-a949-bd7fb6808dd6
md"""
## 📊 `diagnostics` - Annual Diagnostics

Accumulates annual-mean fields and prints a summary line at the end of each year.

**Key logic (mutates global accumulators in-place):**
- Accumulates `Tsmn`, `Tamn`, `Tomn`, `qmn`, `amn`, `swmn`, `lwmn`, `qlatmn`, `qsensmn`, `ftmn`, `fqmn`
- At `ityr == nstep_yr`: divides by `nstep_yr`, prints global-mean / sample-point temperatures, then resets
"""

# ╔═╡ cc0e682c-8767-498e-8178-2de4e796b3a8
function diagnostics!(it, year, CO2, Ts0, Ta0, To0, q0, albedo, sw, lw_surf, q_lat, q_sens, timestate::TimeState)
	# Accumulate
	Tsmn    .+= Ts0;   Tamn    .+= Ta0;   Tomn    .+= To0
	qmn     .+= q0;    amn     .+= albedo
	swmn    .+= sw;    lwmn    .+= lw_surf
	qlatmn  .+= q_lat; qsensmn .+= q_sens
	ftmn    .+= @view TF_correct[:, :, timestate.ityr]
	fqmn    .+= @view qF_correct[:, :, timestate.ityr]

	if timestate.ityr == nstep_yr
		# Compute annual means
		Tsmn    ./= nstep_yr;  Tamn    ./= nstep_yr;  Tomn    ./= nstep_yr
		qmn     ./= nstep_yr;  amn     ./= nstep_yr
		swmn    ./= nstep_yr;  lwmn    ./= nstep_yr
		qlatmn  ./= nstep_yr;  qsensmn ./= nstep_yr
		ftmn    ./= nstep_yr;  fqmn    ./= nstep_yr

		# Print summary: global mean, sample points (Fortran indices 48,27 and 16,38)
		println(year, "  ", sum(Tsmn) / (xdim * ydim) - 273.15,
		        "  ", Tsmn[48, 24 + 3] - 273.15,
		        "  ", Tsmn[16, 24 + 14] - 273.15)

		# Reset accumulators
		Tsmn .= 0.0;  Tamn .= 0.0;  Tomn .= 0.0;  qmn .= 0.0;  amn .= 0.0
		swmn .= 0.0;  lwmn .= 0.0;  qlatmn .= 0.0;  qsensmn .= 0.0
		ftmn .= 0.0;  fqmn .= 0.0
	end
	return nothing
end

# ╔═╡ 7281dc60-e0e2-4a34-b4b6-c3f63a97f60f
md"""
## 💾 `output` - Monthly Output

Accumulates monthly-mean fields and pushes a `MonthlyRecord` to the output buffer at month boundaries.

**Key logic:**
- Accumulates `Tmm`, `Tamm`, `Tomm`, `qmm`, `apmm` every time step
- At the last time step of each calendar month: divides by `ndm = jday_mon[mon] × ndt_days`,
  pushes a `MonthlyRecord` NamedTuple, resets accumulators, advances `mon`
- Returns `(irec = irec, mon = mon)` for caller bookkeeping
"""

# ╔═╡ dfdde9f1-b226-4af0-9ac2-36f1b01622fa
function output!(it::Int, irec::Int, mon::Int, 
                Ts0, Ta0, To0, q0, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, 					qsens, output_buf::Vector{MonthlyRecord}, 
                acc::MonthlyAccumulator, timestate::TimeState)::NamedTuple{(:irec, :mon),Tuple{Int, Int}}
	# Use MonthlyAccumulator struct
	accumulate!(acc, Ts0, Ta0, To0, q0, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens)

	# Check: last time step of current calendar month?
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
		# Reset accumulator
		reset!(acc)
		mon += 1
		if mon == 13
			mon = 1
		end
	end
	return (irec = irec, mon = mon)
end

# ╔═╡ 4f97badf-a501-4c70-a943-d3c86b48f8a1
function build_monthly_climatology(records::Vector{MonthlyRecord})::Vector{MonthlyRecord}
	isempty(records) && return MonthlyRecord[]

	fields = propertynames(records[1])
	counts = zeros(Int, 12)
	clim_acc = [Dict{Symbol, Matrix{Float64}}() for _ in 1:12]

	for idx in eachindex(records)
		mon = mod(idx - 1, 12) + 1
		rec = records[idx]
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

	clim = MonthlyRecord[]
	for mon in 1:12
		if counts[mon] == 0
			push!(clim, records[1])
		else
			push!(clim,
				NamedTuple{fields}(Tuple(clim_acc[mon][fld] ./ counts[mon] for fld in fields))
			)
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

	for idx in eachindex(scnr_records)
		mon = mod(idx - 1, 12) + 1
		rec = scnr_records[idx]
		ref = ctrl_clim[mon]
		push!(anom,
			NamedTuple{fields}(Tuple(getfield(rec, fld) .- getfield(ref, fld) for fld in fields))
		)
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
	input_data_dir = joinpath(@__DIR__, "greb-official-official", "greb-official-official", "input")
	
	# Load the data (choose dataset: :ncep, :era, or :era_ncep)
	load_greb_input_data!(input_data_dir, dataset=:ncep)
end

# ╔═╡ 83e54812-2291-44e6-9b6e-0c0be57865d3
md"""
---
### `greb_model` — Main Integration Driver

Orchestrates the three phases of the GREB model:
1. Flux-correction spin-up (`time_flux` years)
2. Control run at constant CO₂ (`time_ctrl` years)
3. Scenario run with time-varying CO₂ (`time_scnr` years)

Returns `(ctrl = Vector{MonthlyRecord}, scnr = Vector{MonthlyRecord})`.
"""

# ╔═╡ d17c570f-ad56-4d43-98cb-ccd03a8fcb4a
md"""
---
## 🎛️ Interactive Model Control Panel

Configure and run GREB model experiments interactively using the controls below.

**Usage:**
1. Select experiment type from dropdown
2. Configure model settings (years, CO₂ concentration)
3. Click "Run Model" to execute
4. View results in plots below

**⚠️ Troubleshooting:**

**If Precipitation = 0 everywhere:**
- Check that **Hydrological Cycle** is enabled in Mean Climate settings
- Verify `log_exp` is not ≤ 6 (which disables hydrology)
- Ensure climate input data is loaded properly

**If Albedo appears constant/uniform:**
- Check that `log_exp` is not ≤ 5 (which disables ice-albedo feedback)
- Verify surface temperature varies enough to trigger albedo changes
- Check that **Ice-Albedo Feedback** is enabled in Mean Climate settings
"""

# ╔═╡ 71676720-4685-4216-a55e-5918829d258f
md"""
### Experiment Selection
"""

# ╔═╡ b7cef546-f0f4-4f42-a0ed-00ffb92a422e
begin
@bind log_exp Select(
	[
		0 => "Full Model (default)",
		1 => "Constant Topography",
		20 => "2×CO₂",
		21 => "4×CO₂",
		27 => "Solar +27W/m²",
		240 => "El Niño",
		241 => "La Niña",
		30 => "Paleo 231kyr",
		230 => "RCP8.5 Climate Change"
	],
	default=0
)
end

# ╔═╡ 9bff59c5-4631-4091-8230-989a835788e5
function seaice!(Ts0, timestate::TimeState, cfg::PhysicsConfig)
	mld = @view mldclim[:, :, timestate.ityr]

	for j in 1:ydim, i in 1:xdim
		z_topo[i,j] >= 0.0 && continue          # skip land
		if Ts0[i,j] <= To_ice1
			cap_surf[i,j] = cap_land                             # full sea ice
		elseif Ts0[i,j] >= To_ice2
			cap_surf[i,j] = cap_ocean * mld[i,j]                 # open ocean
		else
			cap_surf[i,j] = cap_land + (cap_ocean * mld[i,j] - cap_land) /
				(To_ice2 - To_ice1) * (Ts0[i,j] - To_ice1)       # partial ice
		end
	end

	# Process control switches
	if !cfg.log_ice_dmc || !cfg.log_ice_drsp
		# No ice feedback: use simple land/ocean distinction
		for j in 1:ydim, i in 1:xdim
			if z_topo[i,j] >= 0.0
				cap_surf[i,j] = cap_land
			else
				cap_surf[i,j] = cap_ocean * mld[i,j]
			end
		end
		return
	end

	# Sensitivity experiments: no ice feedback
	if log_exp <= 5
		for j in 1:ydim, i in 1:xdim
			if z_topo[i,j] >= 0.0
				cap_surf[i,j] = cap_land
			else
				cap_surf[i,j] = cap_ocean * mld[i,j]
			end
		end
	end

	# Glacier cells always land-ice capacity
	@. cap_surf = ifelse(glacier > 0.5, cap_land, cap_surf)

	return nothing
end

# ╔═╡ 7969e897-121e-4fe0-9b75-36d21931357f
md"""
#### 🎛️ Configuration Preset
"""

# ╔═╡ dec9de54-eede-45ae-94cd-5bcc477298bf
@bind config_preset Select([
	"full" => "🌍 Full Physics (All processes active)",
	"no_feedbacks" => "🔒 No Feedbacks (Fixed albedo, clouds, vapor)",
	"mscm" => "🌦️ MSCM Original (Mean State Climate Model)",
	"sensitivity" => "🧪 Sensitivity Test (Minimal processes)",
	"custom" => "⚙️ Custom Configuration"
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
	$(Child("ice", CheckBox(default=true))) Ice Feedback  
	$(Child("crcl", CheckBox(default=true))) Circulation  
	$(Child("hydro", CheckBox(default=true))) Hydrology  
	$(Child("deepocean", CheckBox(default=true))) Deep Ocean  
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
	$(Child("ice", CheckBox(default=true))) Ice Feedback  
	$(Child("crcl", CheckBox(default=true))) Circulation  
	$(Child("hydro", CheckBox(default=true))) Hydrology  
	$(Child("deepocean", CheckBox(default=true))) Deep Ocean  
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
		# 🌍 Full Physics: All processes active
		log_clouds_dmc, log_vapor_dmc, log_ice_dmc = true, true, true
		log_crcl_dmc, log_hydro_dmc, log_deepocean_dmc = true, true, true
		log_atmos_dmc, log_co2_dmc, log_ocean_dmc, log_qflux_dmc = true, true, true, true
		
		log_clouds_drsp, log_vapor_drsp, log_ice_drsp = true, true, true
		log_crcl_drsp, log_hydro_drsp, log_deepocean_drsp = true, true, true
		log_topo_drsp, log_humid_drsp = true, true
		
		log_ice, log_hdif, log_hadv = true, true, true
		log_vdif, log_vadv, log_conv = true, true, true
		
		log_rain, log_eva, log_clim = 0, -1, 0
		log_tsurf_ext, log_hwind_ext, log_omega_ext = false, false, false
		
	elseif config_preset == "no_feedbacks"
		# 🔒 No Feedbacks: Fixed albedo, clouds, vapor
		log_clouds_dmc, log_vapor_dmc, log_ice_dmc = false, false, false
		log_crcl_dmc, log_hydro_dmc, log_deepocean_dmc = true, true, true
		log_atmos_dmc, log_co2_dmc, log_ocean_dmc, log_qflux_dmc = true, true, true, true
		
		log_clouds_drsp, log_vapor_drsp, log_ice_drsp = false, false, false
		log_crcl_drsp, log_hydro_drsp, log_deepocean_drsp = true, true, true
		log_topo_drsp, log_humid_drsp = true, true
		
		log_ice, log_hdif, log_hadv = false, true, true
		log_vdif, log_vadv, log_conv = true, true, true
		
		log_rain, log_eva, log_clim = 0, -1, 0
		log_tsurf_ext, log_hwind_ext, log_omega_ext = false, false, false
		
	elseif config_preset == "mscm"
		# 🌦️ MSCM Original: Mean State Climate Model configuration
		log_clouds_dmc, log_vapor_dmc, log_ice_dmc = true, true, true
		log_crcl_dmc, log_hydro_dmc, log_deepocean_dmc = true, true, true
		log_atmos_dmc, log_co2_dmc, log_ocean_dmc, log_qflux_dmc = true, true, true, true
		
		log_clouds_drsp, log_vapor_drsp, log_ice_drsp = true, true, true
		log_crcl_drsp, log_hydro_drsp, log_deepocean_drsp = true, true, true
		log_topo_drsp, log_humid_drsp = true, true
		
		log_ice, log_hdif, log_hadv = true, true, true
		log_vdif, log_vadv, log_conv = true, true, true
		
		log_rain, log_eva, log_clim = 0, -1, 0
		log_tsurf_ext, log_hwind_ext, log_omega_ext = false, false, false
		
	elseif config_preset == "sensitivity"
		# 🧪 Sensitivity Test: Minimal processes for testing
		log_clouds_dmc, log_vapor_dmc, log_ice_dmc = false, false, false
		log_crcl_dmc, log_hydro_dmc, log_deepocean_dmc = false, false, false
		log_atmos_dmc, log_co2_dmc, log_ocean_dmc, log_qflux_dmc = true, true, true, false
		
		log_clouds_drsp, log_vapor_drsp, log_ice_drsp = false, false, false
		log_crcl_drsp, log_hydro_drsp, log_deepocean_drsp = false, false, false
		log_topo_drsp, log_humid_drsp = false, false
		
		log_ice, log_hdif, log_hadv = false, false, false
		log_vdif, log_vadv, log_conv = false, false, false
		
		log_rain, log_eva, log_clim = -1, -1, 0
		log_tsurf_ext, log_hwind_ext, log_omega_ext = false, false, false
		
	else  # custom
		# ⚙️ Custom: Read from interactive panels
		log_clouds_dmc = mean_climate_panel.clouds
		log_vapor_dmc = mean_climate_panel.vapor
		log_ice_dmc = mean_climate_panel.ice
		log_crcl_dmc = mean_climate_panel.crcl
		log_hydro_dmc = mean_climate_panel.hydro
		log_deepocean_dmc = mean_climate_panel.deepocean
		log_atmos_dmc = mean_climate_panel.atmos
		log_co2_dmc = mean_climate_panel.co2
		log_ocean_dmc = mean_climate_panel.ocean
		log_qflux_dmc = mean_climate_panel.qflux
		
		log_clouds_drsp = co2_response_panel.clouds
		log_vapor_drsp = co2_response_panel.vapor
		log_ice_drsp = co2_response_panel.ice
		log_crcl_drsp = co2_response_panel.crcl
		log_hydro_drsp = co2_response_panel.hydro
		log_deepocean_drsp = co2_response_panel.deepocean
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

# ╔═╡ 0be9bc61-1a59-4dfe-84f9-bb1a27ca30fc
begin
	# 🌡️ 3D climate fields (xdim, ydim, nstep_yr) ───────────────────
	# Input fields (to be loaded from data)
	Tclim     = zeros(Float64, xdim, ydim, nstep_yr)   # surface temperature [K]
	uclim     = zeros(Float64, xdim, ydim, nstep_yr)   # zonal wind [m/s]
	vclim     = zeros(Float64, xdim, ydim, nstep_yr)   # meridional wind [m/s]
	qclim     = zeros(Float64, xdim, ydim, nstep_yr)   # atmospheric humidity [kg/kg]
	mldclim   = zeros(Float64, xdim, ydim, nstep_yr)   # mixed-layer depth [m]
	
	# MSCM additional climatology fields
	omegaclim    = zeros(Float64, xdim, ydim, nstep_yr) # vertical velocity [Pa/s]
	omegastdclim = zeros(Float64, xdim, ydim, nstep_yr) # omega std deviation [Pa/s]
	wsclim       = zeros(Float64, xdim, ydim, nstep_yr) # wind speed [m/s]

	# 📊 Anomaly Fields for ENSO/Climate Change Experiments ────────────────
	# ENSO anomaly fields (log_exp 240-241)
	Tclim_anom_enso  = zeros(Float64, xdim, ydim, nstep_yr) # surface temperature [K]
	uclim_anom_enso  = zeros(Float64, xdim, ydim, nstep_yr) # zonal wind [m/s]
	vclim_anom_enso  = zeros(Float64, xdim, ydim, nstep_yr) # meridional wind [m/s]
	omegaclim_anom_enso = zeros(Float64, xdim, ydim, nstep_yr) # vertical velocity [Pa/s]
	wsclim_anom_enso = zeros(Float64, xdim, ydim, nstep_yr) # wind speed [m/s]
	
	# Climate change anomaly fields (log_exp 230)
	Tclim_anom_cc    = zeros(Float64, xdim, ydim, nstep_yr) # surface temperature [K]
	uclim_anom_cc    = zeros(Float64, xdim, ydim, nstep_yr) # zonal wind [m/s] 
	vclim_anom_cc    = zeros(Float64, xdim, ydim, nstep_yr) # meridional wind [m/s]
	omegaclim_anom_cc = zeros(Float64, xdim, ydim, nstep_yr) # vertical velocity [Pa/s]
	wsclim_anom_cc   = zeros(Float64, xdim, ydim, nstep_yr) # wind speed [m/s]
	
	# Apply mean climate deconstruction switches
	if !log_hydro_dmc
		qclim .= 0.0  # zero out humidity climatology
	end
	
	# 🌬️ Precomputed wind sign splits (MSCM optimization) ───────────────────
	uclim_m   = zeros(Float64, xdim, ydim, nstep_yr)   # negative u components
	uclim_p   = zeros(Float64, xdim, ydim, nstep_yr)   # positive u components  
	vclim_m   = zeros(Float64, xdim, ydim, nstep_yr)   # negative v components
	vclim_p   = zeros(Float64, xdim, ydim, nstep_yr)   # positive v components
	
	# Initialize wind component separation (CRITICAL: affects advection)
	@. uclim_m = ifelse(uclim >= 0.0, uclim, 0.0)  # positive winds only
	@. uclim_p = ifelse(uclim < 0.0, uclim, 0.0)   # negative winds only
	@. vclim_m = ifelse(vclim >= 0.0, vclim, 0.0)  # positive winds only
	@. vclim_p = ifelse(vclim < 0.0, vclim, 0.0)   # negative winds only
	Toclim    = zeros(Float64, xdim, ydim, nstep_yr)   # deep ocean temperature [K]
	cldclim   = zeros(Float64, xdim, ydim, nstep_yr)   # cloud cover fraction
	swetclim  = zeros(Float64, xdim, ydim, nstep_yr)   # soil wetness [0-1]
	
	nothing  # suppress output
end;

# ╔═╡ d404043f-8080-4262-9ab6-d9bb13eee504
function init_model!(log_exp, cfg::PhysicsConfig)

	# ── MSCM Hydrology Parameter Initialization ────────────
	set_hydrology_parameters!(cfg)
	
	# ── dTrad: offset between Tatmos and radiation temperature ────
	@. dTrad = -0.16 * Tclim - 5.0

	# ── z_ocean: 3× maximum mixed-layer depth over the year ──────
	z_ocean .= 3.0 .* dropdims(maximum(mldclim; dims=3); dims=3) # vectorized

	# ── Sensitivity experiment overrides ─────────────────────────
	if log_exp == 1
		@. z_topo = min(z_topo, 1.0)       # constant topography
	end
	
	if log_exp <= 2
		cldclim .= 0.7                     # constant cloud cover
	end

	# Apply cloud deconstruction switch
	if !cfg.log_clouds_dmc
		cldclim .= 0.0  # zero cloud climatology
	end
	
	# Apply flux correction conditional zeroing (MSCM feature)
	if log_exp == 1 && !cfg.log_qflux_dmc
		TF_correct .= 0.0
		qF_correct .= 0.0
		ToF_correct .= 0.0
	end

	if !cfg.log_hydro_dmc
		qclim .= 0.0  # zero out humidity climatology
	end
	
	if log_exp <= 3
		qclim .= 0.0052                    # constant water vapor
	end
	
	if log_exp <= 9 || log_exp == 11
		mldclim .= d_ocean                 # no deep ocean
	end

	# ── MSCM Experiment Handler ──────────────────────────────────────────
	# Apply advanced experiment forcing (230=climate change, 240/241=ENSO)
	if log_exp == 230
		@info "Applying CMIP5 RCP8.5 climate change forcing (log_exp=230)"
		Tclim      .+= Tclim_anom_cc
		uclim      .+= uclim_anom_cc
		vclim      .+= vclim_anom_cc
		omegaclim  .+= omegaclim_anom_cc
		wsclim     .+= wsclim_anom_cc
	elseif log_exp == 240 || log_exp == 241
		sign = (log_exp == 240) ? 1.0 : -1.0  # El Niño vs La Niña
		type_str = (log_exp == 240) ? "El Niño" : "La Niña"
		@info "Applying ERA-Interim $type_str forcing (log_exp=$log_exp)"
		
		Tclim      .+= sign .* Tclim_anom_enso
		uclim      .+= sign .* uclim_anom_enso
		vclim      .+= sign .* vclim_anom_enso 
		omegaclim  .+= sign .* omegaclim_anom_enso
		wsclim     .+= sign .* wsclim_anom_enso
	end
	
	# ── Topography pressure weights (used in multiple kernels) ─────
	@. wz_air = exp(-z_topo / z_air)
	@. wz_vapor = exp(-z_topo / z_vapor)

	# ── Surface heat capacity ────────────────────────────────────
	for j in 1:ydim, i in 1:xdim
		if z_topo[i, j] > 0.0
			cap_surf[i, j] = cap_land
		else
			# Apply ocean deconstruction switch (MSCM feature)
			if cfg.log_ocean_dmc
				cap_surf[i, j] = cap_ocean * mldclim[i, j, 1]  # normal ocean
			else
				cap_surf[i, j] = cap_land  
			end
		end
	end

	# ── Initial conditions from last time step of climatology ────
	Ts_ini  = Tclim[:, :, nstep_yr]       |> copy   # surface temperature
	Ta_ini  = copy(Ts_ini)                          # air temperature = Tsurf
	To_ini  = Toclim[:, :, nstep_yr]      |> copy   # deep ocean temperature
	q_ini   = qclim[:, :, nstep_yr]       |> copy   # atmospheric water vapor

	# ── Control CO₂ level ───────────────────────────────────────
	CO2_ctrl = 340.0
	if !cfg.log_co2_dmc
		CO2_ctrl = 0.0  # Zero CO2 for deconstruction experiments
	end
	if log_exp == 12 || log_exp == 13
		CO2_ctrl = 298.0                             # A1B scenario
	end
	if 95 <= log_exp && log_exp <= 100
		CO2_ctrl = 280.0  # IPCC scenarios baseline
	end

	return (Ts_ini = Ts_ini, Ta_ini = Ta_ini, To_ini = To_ini,
	        q_ini  = q_ini,  CO2_ctrl = CO2_ctrl)
end

# ╔═╡ db75ea52-9387-4c15-bde5-61777ac9b570
begin
	function circulation!(X_in, h_scl, ws::CirculationWorkspace, timestate::TimeState, log_exp::Int, cfg::PhysicsConfig)
		# Use pre-allocated workspace buffer
		dX_crcl = ws.dX_crcl
		fill!(dX_crcl, 0.0)

		# Process control switches
		!cfg.log_atmos_dmc && return dX_crcl  
		!cfg.log_crcl_dmc && return dX_crcl
		!cfg.log_crcl_drsp && return dX_crcl

		# Sensitivity experiments: no circulation
		log_exp <= 4 && return dX_crcl
		(log_exp == 7  && h_scl == z_vapor) && return dX_crcl
		(log_exp == 16 && h_scl == z_vapor) && return dX_crcl

		# Number of sub-steps within one main time step
		ntime = max(1, round(Int, Δt / Δt_crcl))

		X1 = copy(X_in)
		for _tt in 1:ntime
    		fill!(ws.dX_diff, 0.0)
    		fill!(ws.dX_adv, 0.0)
    		fill!(ws.zero_buf1, 0.0)
    		dx_diff = ws.dX_diff
    		dx_adv = ws.dX_adv
    		dx_conv = ws.zero_buf1
			
			# Process-specific switches (pass workspace to diffusion!/advection!)
			if cfg.log_vdif && h_scl == z_vapor
				dx_diff = diffusion!(X1, h_scl, ws, timestate)
			end
			if cfg.log_vadv && h_scl == z_vapor
				dx_adv = advection!(X1, h_scl, ws, timestate, log_exp)
			end
			# CONVERGENCE INTEGRATION (log_conv == 0 enables convergence)
			if !cfg.log_conv && h_scl == z_vapor
				dx_conv = convergence!(X1, omegaclim, timestate.ityr)
			end
			if cfg.log_hdif && h_scl == z_air
				dx_diff = diffusion!(X1, h_scl, ws, timestate)
			end
			if cfg.log_hadv && h_scl == z_air
				dx_adv = advection!(X1, h_scl, ws, timestate, log_exp)
			end
			
			# Update field with all tendencies
			@. X1 = X1 + dx_diff + dx_adv + dx_conv
		end

		@. dX_crcl = X1 - X_in
		return dX_crcl
	end
end

# ╔═╡ e493fae7-239a-494c-9a59-728446d70f7a
function tendencies!(CO2, Ts1, Ta1, To1, q1, ws::CirculationWorkspace,
					 timestate::TimeState, log_exp::Int, cfg::PhysicsConfig)
	
	# Short-wave radiation → albedo, SW flux
	sw_out   = SWradiation!(Ts1, timestate, log_exp, cfg, ws)
	
	# Long-wave radiation → LW_surf, LWair_up, LWair_down, emissivity
	lw_out   = LWradiation!(Ts1, Ta1, q1, CO2, timestate, log_exp, cfg, ws)
	
	# Sensible heat flux
	Q_sens   = @. ct_sens * (Ta1 - Ts1)
	if !cfg.log_atmos_dmc
		Q_sens .= 0.0
	end
	
	# Hydrological cycle → latent heat + evaporation/rain tendencies
	hy_out   = hydro!(Ts1, q1, timestate, log_exp, cfg, ws)
	
	# Atmospheric circulation — temperature diffusion/advection (pass workspace)
	dTa_crcl = circulation!(Ta1, z_air, ws, timestate, log_exp, cfg)
	
	# Atmospheric circulation — water-vapour diffusion/advection (pass workspace)
	dq_crcl  = circulation!(q1, z_vapor, ws, timestate, log_exp, cfg)
	
	# Deep ocean coupling
	do_out   = deep_ocean!(Ts1, To1, timestate, log_exp, cfg, ws)

		return (albedo     = sw_out.albedo,
	        SW         = sw_out.SW,
	        LW_surf    = lw_out.LW_surf,
	        Q_lat      = hy_out.Q_lat,
	        Q_sens     = Q_sens,
	        Q_lat_air  = hy_out.Q_lat_air,
	        dq_eva     = hy_out.dq_eva,
	        dq_rain    = hy_out.dq_rain,
	        dq_crcl    = dq_crcl,
	        dTa_crcl   = dTa_crcl,
	        dT_ocean   = do_out.dT_ocean,
	        dTo        = do_out.dTo,
	        LWair_down = lw_out.LWair_down,
	        LWair_up   = lw_out.LWair_up,
	        em         = lw_out.em)
	end

# ╔═╡ 33fa7b1f-938b-481c-bf25-eca8d7fb33a7
function time_loop!(it::Int, year::Int, CO2::Float64, mon::Int, irec::Int,
				   Ts1::Matrix{Float64}, Ta1::Matrix{Float64},
				   q1::Matrix{Float64}, To1::Matrix{Float64},
				   output_buf::Vector{MonthlyRecord},
				   ws::CirculationWorkspace, acc::MonthlyAccumulator,
				   timestate::TimeState, log_exp::Int, cfg::PhysicsConfig)
	
	# ── Calendar indices (use precomputed lookup) ───────────────
	cal = it <= MAX_TIMESTEPS ? CALENDAR_LOOKUP[it] : (
		day = mod((it - 1) ÷ ndt_days, ndays_yr) + 1,
		step = mod(it - 1, nstep_yr) + 1
	)
	# Update time state
	timestate.jday = cal.day
	timestate.ityr = cal.step

	# ── Compute all tendencies (pass workspace) ─────────────────
	tend = tendencies!(CO2, Ts1, Ta1, To1, q1, ws, timestate, log_exp, cfg)

	# ── Surface temperature ─────────────────────────────────────
	Ts0 = ws.Ts0_buf
	TFc = @view TF_correct[:, :, timestate.ityr]  # Phase 4: Use view
	@inbounds @turbo for j in 1:ydim
		for i in 1:xdim
			Ts0[i,j] = Ts1[i,j] + tend.dT_ocean[i,j] + Δt * (
				tend.SW[i,j] + tend.LW_surf[i,j] - tend.LWair_down[i,j] +
				tend.Q_lat[i,j] + tend.Q_sens[i,j] + TFc[i,j]
			) / cap_surf[i,j]
		end
	end

	# ── Air temperature ─────────────────────────────────────────
	Ta0 = ws.Ta0_buf
	@inbounds @turbo for j in 1:ydim
		for i in 1:xdim
			Ta0[i,j] = Ta1[i,j] + tend.dTa_crcl[i,j] + Δt * (
				tend.LWair_up[i,j] + tend.LWair_down[i,j] - tend.em[i,j] * tend.LW_surf[i,j] +
				tend.Q_lat_air[i,j] - tend.Q_sens[i,j]
			) / cap_air
		end
	end

	# ── Numerical stability clamps ──────────────────────────────
	@. Ts0 = max(Ts0, MIN_TEMPERATURE_K)  # Temperature floor
	@. Ta0 = max(Ta0, MIN_TEMPERATURE_K)  # Atmosphere temperature floor

	# ── Deep ocean temperature ──────────────────────────────────
	To0 = ws.To0_buf
	@. To0 = To1 + tend.dTo + ToF_correct[:, :, timestate.ityr]

	# ── Apply process control switches to tendencies ────────────
	fill!(ws.zero_buf2, 0.0)
	fill!(ws.zero_buf3, 0.0)
	dq_eva_final = cfg.log_hydro_dmc ? tend.dq_eva : ws.zero_buf2
	dq_rain_final = cfg.log_hydro_dmc ? tend.dq_rain : ws.zero_buf2
	dq_crcl_final = cfg.log_crcl_dmc ? tend.dq_crcl : ws.zero_buf3
	
	# ── Humidity ────────────────────────────────────────────────
	q0 = ws.q0_buf
	qFc = @view qF_correct[:, :, timestate.ityr]  # Phase 4: Use view
	dq = @. Δt * (dq_eva_final + dq_rain_final) + dq_crcl_final + qFc
	# Enhanced humidity bounds
	@. dq = clamp(dq, -MIN_HUMIDITY_FRACTION * q1, MAX_HUMIDITY_CHANGE)  # Prevent negative q, limit positive changes

	# Update humidity state
	@. q0 = q1 + dq

	# ── Sea-ice heat capacity adjustment ─────────────────────────
	seaice!(Ts0, timestate, cfg)

	# ── Compute ice fraction ─────────────────────────────────────
	ice = ice_fraction_field(Ts0, ws)

	# Analysis-friendly units (mm/day)
		# Precipitation: -dq_rain (negative tendency) → positive mm/day
		# Evaporation: dq_eva (positive tendency) → positive mm/day
		# Conversion: dq [kg/kg/timestep] × r_qviwv [kg/m³] × wz_vapor × 86400 [s/day]
	precip_out = @. -dq_rain_final * wz_vapor * r_qviwv * 86400.0 
	evap_out = @. dq_eva_final * wz_vapor * r_qviwv * 86400.0
	qcrcl_out = @. dq_crcl_final

	# ── Output + diagnostics ─────────────────────────────────────
	res_out = output!(it, irec, mon, Ts0, Ta0, To0, q0, tend.albedo, 
	                  ice, precip_out, evap_out, qcrcl_out,
	                  tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens,
	                  output_buf, acc, timestate)
	diagnostics!(it, year, CO2, Ts0, Ta0, To0, q0, tend.albedo,
	             tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens, timestate)

	return (Ts0 = Ts0, Ta0 = Ta0, q0 = q0, To0 = To0,
	        albedo = tend.albedo, mon = res_out.mon, irec = res_out.irec)
end

# ╔═╡ fee2639d-d8cf-4ee3-b824-73a0b02fef4a
md"""
### Run Duration Configuration (years)
"""

# ╔═╡ 952d5e5d-119a-4ed9-9e22-0b0fe66ae04d
@bind time_flux Slider(0:10, default=0, show_value=true)

# ╔═╡ 584e767e-4dc5-4821-af63-d6a825326d9e
function qflux_correction!(CO2_ctrl, Ts1, Ta1, q1, To1, log_exp::Int, cfg::PhysicsConfig)
	
	# Create workspace for circulation
	ws = CirculationWorkspace()
	# Create local time state
	timestate = TimeState(1, 1)
	
	for it in 1:(time_flux * ndt_days * ndays_yr)
		# Update time state
		timestate.jday = mod((it - 1) ÷ ndt_days, ndays_yr) + 1
		timestate.ityr = mod(it - 1, nstep_yr) + 1
		
		# All physics tendencies (pass workspace)
		tend = tendencies!(CO2_ctrl, Ts1, Ta1, To1, q1, ws, timestate, log_exp, cfg)

		# Views into climatology & correction fields (avoid copies)
		Tc   = @view Tclim[:, :, timestate.ityr]
		Toc  = @view Toclim[:, :, timestate.ityr]
		qc   = @view qclim[:, :, timestate.ityr]
		TFc  = @view TF_correct[:, :, timestate.ityr]
		ToFc = @view ToF_correct[:, :, timestate.ityr]
		qFc  = @view qF_correct[:, :, timestate.ityr]

		# ── Uncorrected updates ───────────────────────────────
		# Surface temperature (no flux correction)
		dTs = @. Δt * (tend.SW + tend.LW_surf - tend.LWair_down +
			           tend.Q_lat + tend.Q_sens) / cap_surf
		Ts0 = @. Ts1 + dTs + tend.dT_ocean

		# Air temperature
		dTa = @. Δt * (tend.LWair_up + tend.LWair_down -
			           tend.em * tend.LW_surf + tend.Q_lat_air - tend.Q_sens) / cap_air
		Ta0 = @. Ta1 + dTa + tend.dTa_crcl

		# Deep ocean (no correction)
		To0 = @. To1 + tend.dTo

		# Humidity (no correction)
		dq = @. Δt * (tend.dq_eva + tend.dq_rain)
		q0 = @. q1 + dq + tend.dq_crcl

		# ── Compute correction fields ────────────────────────
		# Surface temperature correction
		@. TFc = (Tc - Ts0) * cap_surf / Δt

		# Re-compute Ts0 with correction
		Ts0 = @. Ts1 + dTs + tend.dT_ocean + TFc * Δt / cap_surf

		# Deep ocean correction
		@. ToFc = Toc - To0
		To0 = @. To1 + tend.dTo + ToFc

		# Humidity correction
		@. qFc = qc - q0
		q0 = @. q1 + dq + tend.dq_crcl + qFc

		# Sea ice
		seaice!(Ts0, timestate, cfg)

		# Diagnostics
		diagnostics!(it, 0.0, CO2_ctrl, Ts0, Ta0, To0, q0, tend.albedo,
			         tend.SW, tend.LW_surf, tend.Q_lat, tend.Q_sens, timestate)

		# Advance state
		Ts1 .= Ts0;  Ta1 .= Ta0
		q1  .= q0;   To1 .= To0
	end
	return nothing
end

# ╔═╡ 92c3bd68-bd07-4381-9c04-e6611650cd1e
function greb_model!(log_exp, time_flux, time_ctrl, time_scnr, cfg::PhysicsConfig)

	# ── 1. Initialisation ───────────────────────────────────────
	ini = init_model!(log_exp, cfg)
	Ts_ini   = ini.Ts_ini
	Ta_ini   = ini.Ta_ini
	To_ini   = ini.To_ini
	q_ini    = ini.q_ini
	CO2_ctrl = ini.CO2_ctrl

	# Optimization: Create workspace and accumulator
	ws = CirculationWorkspace()
	acc = MonthlyAccumulator()
	
	# ── 2. Flux-correction spin-up ──────────────────────────────
	if log_exp != 1 || cfg.log_qflux_dmc
		# Load flux corrections when needed (log_exp==1 && cfg.log_qflux_dmc)
		if log_exp == 1 && cfg.log_qflux_dmc
			println("% loading flux correction fields...")
			load_flux_corrections!(input_data_dir)
		end
		println("% flux correction  CO2 = ", CO2_ctrl)
		qflux_correction!(CO2_ctrl, Ts_ini, Ta_ini, q_ini, To_ini, log_exp, cfg)
	else
		println("% flux correction skipped (log_exp=1 and log_qflux_dmc=false)")
	end

	# ── 3. Control run ──────────────────────────────────────────
	println("% CONTROL RUN  CO2 = ", CO2_ctrl, "  time = ", time_ctrl, " yr")
	Ts1 = copy(Ts_ini);  Ta1 = copy(Ta_ini)
	To1 = copy(To_ini);  q1  = copy(q_ini)
	sw_solar_forcing_state[] = 1.0
	mon  = 1;  year = 1970;  irec = 0
	reset!(acc)  # Use accumulator reset

	ctrl_output = MonthlyRecord[]
	timestate = TimeState(1, 1)  # Initialize time state

	for it in 1:(time_ctrl * nstep_yr)
		res = time_loop!(it, year, CO2_ctrl, mon, irec, Ts1, Ta1, q1, To1, ctrl_output, ws, acc, timestate, log_exp, cfg)
		Ts1 .= res.Ts0;  Ta1 .= res.Ta0
		q1  .= res.q0;   To1 .= res.To0
		mon  = res.mon;  irec = res.irec
	end

	# ── 4. Scenario run ─────────────────────────────────────────
	println("% SCENARIO EXP: ", log_exp, "  time = ", time_scnr, " yr")
	Ts1 .= Ts_ini;  Ta1 .= Ta_ini
	q1  .= q_ini;   To1 .= To_ini
	year = 1950;  CO2 = 340.0;  mon = 1;  irec = 0
	if 35 <= log_exp && log_exp <= 37
		year = 1
	end
	sw_solar_forcing_state[] = 1.0
	reset!(acc)  # Use accumulator reset

	scnr_output = MonthlyRecord[]
	timestate = TimeState(1, 1)  # Initialize time state
	# Precompute annual-mean control ice once and replicate
	# across month dimension to match forcing(icmn_ctrl) interface.
	icmn_ctrl_annual = dropdims(sum(icmn_ctrl; dims=3) ./ 12.0; dims=3)
	icmn_ctrl_forcing = repeat(reshape(icmn_ctrl_annual, xdim, ydim, 1), 1, 1, 12)

	for it in 1:(time_scnr * nstep_yr)
		# Get forcing from enhanced dispatcher
		forcing_result = forcing(it, year, log_exp, icmn_ctrl_forcing; nstep_yr=nstep_yr)
		CO2 = forcing_result.CO2
		sw_solar_forcing_state[] = forcing_result.sw_solar_forcing

		# Forced-boundary experiments: keep Tsurf on external climatological boundary
		if log_exp == 230 || log_exp == 240 || log_exp == 241
			ityr_now = mod(it - 1, nstep_yr) + 1
			Ts1 .= @view Tclim[:, :, ityr_now]
		end

		# Sensitivity experiment: SST+1 K over ocean
		if 14 <= log_exp && log_exp <= 16
			CO2 = CO2_ctrl
			ityr_now = mod(it - 1, nstep_yr) + 1
			for j in 1:ydim, i in 1:xdim
				if z_topo[i, j] < 0.0
					Ts1[i, j] = Tclim[i, j, ityr_now] + 1.0
				end
			end
		end

		res = time_loop!(it, year, CO2, mon, irec, Ts1, Ta1, q1, To1, scnr_output, ws, acc, timestate, log_exp, cfg)
		Ts1 .= res.Ts0;  Ta1 .= res.Ta0
		q1  .= res.q0;   To1 .= res.To0
		mon  = res.mon;  irec = res.irec

		# Advance year at year boundary
		if mod(it, nstep_yr) == 0
			year += 1
		end
	end 

	if !(35 <= log_exp && log_exp <= 37) && !isempty(ctrl_output) && !isempty(scnr_output)
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
- Experiment: log_exp = $log_exp
- Flux correction: $time_flux years
- Control run: $time_ctrl years  
- Scenario run: $time_scnr years

⚠️ **Warning:** Runs with total time > 50 years may take 10+ minutes.
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

# ╔═╡ 9d3eea3d-1ba6-4650-9187-85fa5aeca210
md"""
---
### 💾 Export Data

Export model results to NetCDF format for further analysis.
"""

# ╔═╡ 3d93d6ce-6f7a-4083-ade9-89ba0058bae0
md"""
**Export Controls:**
"""

# ╔═╡ 47d3626b-463c-4090-b12a-5ed5f5e8ce62
@bind export_controls PlutoUI.combine() do Child
	md"""
	**Dataset:** $(Child("dataset", Select(["ctrl" => "Control Run", "scnr" => "Scenario Run"], default="ctrl")))

	**Export Mode:** $(Child("mode", Select(["single" => "Single Month", "timeseries" => "Full Time Series"], default="single")))
	
	**Export Month:** $(Child("month", Slider(1:12, default=1, show_value=true)))
	
	**Filename:** $(Child("filename", TextField(default="greb_output")))
	"""
end

# ╔═╡ 93913dcc-aae0-40fe-9837-92e9a251aadf
md"""
**Export Summary:**

Click the button below to export either a single month or the full time series. Files will be saved in the workspace directory.
"""

# ╔═╡ 3f7349aa-da11-4d27-a637-5286950d3244
@bind do_export CheckBox()

# ╔═╡ 39afec69-f90d-4ff2-99d1-5cfec84d09a1
function current_physics_config()
	return PhysicsConfig(
		log_clouds_dmc = log_clouds_dmc,
		log_vapor_dmc = log_vapor_dmc,
		log_ice_dmc = log_ice_dmc,
		log_crcl_dmc = log_crcl_dmc,
		log_hydro_dmc = log_hydro_dmc,
		log_deepocean_dmc = log_deepocean_dmc,
		log_atmos_dmc = log_atmos_dmc,
		log_co2_dmc = log_co2_dmc,
		log_ocean_dmc = log_ocean_dmc,
		log_qflux_dmc = log_qflux_dmc,
		log_clouds_drsp = log_clouds_drsp,
		log_vapor_drsp = log_vapor_drsp,
		log_ice_drsp = log_ice_drsp,
		log_crcl_drsp = log_crcl_drsp,
		log_hydro_drsp = log_hydro_drsp,
		log_deepocean_drsp = log_deepocean_drsp,
		log_topo_drsp = log_topo_drsp,
		log_humid_drsp = log_humid_drsp,
		log_ice = log_ice,
		log_hdif = log_hdif,
		log_hadv = log_hadv,
		log_vdif = log_vdif,
		log_vadv = log_vadv,
		log_conv = log_conv,
		log_rain = log_rain,
		log_eva = log_eva,
		log_clim = log_clim,
		log_tsurf_ext = log_tsurf_ext,
		log_hwind_ext = log_hwind_ext,
		log_omega_ext = log_omega_ext,
	)
end

# ╔═╡ cf2a51eb-a661-4c66-9e0a-d1730110e4bc
begin
    # Ensure last_run is always defined globally
    if !@isdefined(last_run)
        global last_run = nothing
    end
    
    # Initialize result variable
    local result = nothing
    
    # Only run when toggle is ON
    if run_toggle
		cfg = current_physics_config()
		
        # Execute model with current parameters
        result = greb_model!(log_exp, time_flux, time_ctrl, time_scnr, cfg)
        
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

# ╔═╡ 3f7349aa-da11-4d27-a637-5286950d3245
begin
    # Define lon/lat centers for NetCDF export
    lon_centers = range(0, length=xdim, step=dlon)
    lat_centers = range(-90 + dlat/2, length=ydim, step=dlat)
    
    const EXPORT_FIELDS = (:Ts, :Ta, :To, :q, :albedo, :ice, :precip, :evap, :qcrcl, 						   :sw, :lw, :qlat, :qsens)
    const EXPORT_UNITS = Dict(
        :Ts => "K", :Ta => "K", :To => "K", :q => "kg/kg", :albedo => "1", 
		:ice => "1", :precip => "mm/day", :evap => "mm/day", :qcrcl => "mm/day",
        :sw => "W m-2", :lw => "W m-2", :qlat => "W m-2", :qsens => "W m-2",
    )

	function export_to_netcdf(record::MonthlyRecord, filename::String; dataset_name::String="", month_index::Int=0)
        base = endswith(lowercase(filename), ".nc") ? filename[1:end-3] : filename
        nc_file = abspath(base * ".nc")
        out_dir = dirname(nc_file)
        isdir(out_dir) || mkpath(out_dir)

        NCDataset(nc_file, "c") do ds
            defDim(ds, "lon", xdim)
            defDim(ds, "lat", ydim)

            vlon = defVar(ds, "lon", Float64, ("lon",))
            vlat = defVar(ds, "lat", Float64, ("lat",))
            vlon[:] = lon_centers
            vlat[:] = lat_centers
            vlon.attrib["units"] = "degrees_east"
            vlat.attrib["units"] = "degrees_north"

            for fld in EXPORT_FIELDS
                var_data = getfield(record, fld)
                size(var_data) == (xdim, ydim) || error("Unexpected data shape for $fld: $(size(var_data)), expected ($(xdim), $(ydim))")
                v = defVar(ds, String(fld), Float64, ("lon", "lat"))
                v[:, :] = var_data
                v.attrib["units"] = EXPORT_UNITS[fld]
				v.attrib["long_name"] = "$(fld) monthly field"
            end

				ds.attrib["title"] = "GREB monthly export"
				ds.attrib["dataset"] = dataset_name
            ds.attrib["month_index"] = string(month_index)
            ds.attrib["export_mode"] = "single"
        end

        return nc_file
    end

	function export_timeseries_to_netcdf(records::Vector{MonthlyRecord}, filename::String;
		dataset_name::String="")
		isempty(records) && error("Cannot export empty time series")

		base = endswith(lowercase(filename), ".nc") ? filename[1:end-3] : filename
		nc_file = abspath(base * ".nc")
		out_dir = dirname(nc_file)
		isdir(out_dir) || mkpath(out_dir)
		nmonths = length(records)

		NCDataset(nc_file, "c") do ds
			defDim(ds, "lon", xdim)
			defDim(ds, "lat", ydim)
			defDim(ds, "time", nmonths)

			vlon = defVar(ds, "lon", Float64, ("lon",))
			vlat = defVar(ds, "lat", Float64, ("lat",))
			vtime = defVar(ds, "time", Int32, ("time",))
			vlon[:] = lon_centers
			vlat[:] = lat_centers
			vtime[:] = collect(Int32(1):Int32(nmonths))
			vlon.attrib["units"] = "degrees_east"
			vlat.attrib["units"] = "degrees_north"
			vtime.attrib["units"] = "month_index"
			vtime.attrib["long_name"] = "monthly record index"

			for fld in EXPORT_FIELDS
				v = defVar(ds, String(fld), Float64, ("lon", "lat", "time"))
				v.attrib["units"] = EXPORT_UNITS[fld]
				v.attrib["long_name"] = "$(fld) monthly time series"
				for t in 1:nmonths
					var_data = getfield(records[t], fld)
					size(var_data) == (xdim, ydim) || error("Unexpected data shape for $fld at time $t: $(size(var_data)), expected ($(xdim), $(ydim))")
					v[:, :, t] = var_data
				end
			end

			ds.attrib["title"] = "GREB monthly export"
			ds.attrib["dataset"] = dataset_name
			ds.attrib["export_mode"] = "timeseries"
			ds.attrib["n_months"] = string(nmonths)
		end

		return nc_file
	end
end;

# ╔═╡ 3f7349aa-da11-4d27-a637-5286950d3246
begin
	if do_export > 0
		if isnothing(last_run)
			md"**Export failed:** no model run available. Run the model first."
		else
			dataset = Symbol(export_controls.dataset)
			export_mode = String(export_controls.mode)
			records = getfield(last_run, dataset)
			if isempty(records)
				md"**Export failed:** selected dataset is empty."
			else
				month_idx = clamp(export_controls.month, 1, length(records))
				filename = String(strip(export_controls.filename))
				if isempty(filename)
					md"**Export failed:** please provide a filename."
				else
					try
						nc_file = if export_mode == "single"
							export_to_netcdf(records[month_idx], filename;
								dataset_name=string(dataset), month_index=month_idx)
						elseif export_mode == "timeseries"
							export_timeseries_to_netcdf(records, filename;
								dataset_name=string(dataset))
						else
							error("Unknown export mode: $export_mode")
						end

						md"""
						**Export complete:**
						- Dataset: `$(dataset)`
						- Mode: `$(export_mode)`
						- Record index: `$(export_mode == "single" ? month_idx : "all")`
						- Format: `NetCDF (.nc)`
						- File: `$(nc_file)`
						"""
					catch err
						md"**Export failed:** $(sprint(showerror, err))"
					end
				end
			end
		end
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LoopVectorization = "bdcacae8-1622-11e9-2a5c-532679323890"
NCDatasets = "85f8d34a-cbdd-5861-8df4-14fed0d494ab"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
LoopVectorization = "~0.12.173"
NCDatasets = "~0.14.10"
Plots = "~1.41.1"
PlutoUI = "~0.7.65"
StaticArrays = "~1.9.13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.0"
manifest_format = "2.0"
project_hash = "20d3a42f85a62697022a1cd21243a9a61d75403e"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "d81ae5489e13bc03567d4fbbb06c546a5e53c857"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.22.0"

    [deps.ArrayInterface.extensions]
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

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

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
git-tree-sha1 = "0836c647014903bedccf23ba72b5ebb8c89a7db8"
uuid = "179af706-886a-5703-950a-314cd64e0468"
version = "0.2.5"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Preferences", "Static"]
git-tree-sha1 = "f3a21d7fc84ba618a779d1ed2fcca2e682865bab"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.7"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonDataModel]]
deps = ["CFTime", "DataStructures", "Dates", "DiskArrays", "Preferences", "Printf", "Statistics"]
git-tree-sha1 = "cd10f8b38725a6458dd971464daa5a751a67e6b0"
uuid = "1fbeeb36-5f17-413c-809b-666fb144f157"
version = "0.4.2"

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

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

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

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "6c72198e6a101cccdd4c9731d3985e904ba26037"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiskArrays]]
deps = ["ConstructionBase", "LRUCache", "Mmap", "OffsetArrays"]
git-tree-sha1 = "e8c9406f3164633756a0db93334ef23102933eef"
uuid = "3c3547ce-8d99-4f5e-a174-61eb10b00ae3"
version = "0.4.18"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7bb1361afdb33c7f2b085aa49ea8fe1b0fb14e58"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.1+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "83dc665d0312b41367b7263e8a4d172eac1897f4"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.4"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3a948313e7a41eb1db7a1e733e6335f17b4ab3c4"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "7.1.1+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "1828eb7275491981fa5f1752a5e126e8f26f8741"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.17"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "27299071cc29e409488ada41ec7643e0ab19091f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.17+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "50c11ffab2a3d50192a228c313f05b5b5dc5acb2"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.0+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "LibCURL_jll", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "libaec_jll"]
git-tree-sha1 = "e94f84da9af7ce9c6be049e9067e511e17ff89ec"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.14.6+0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Preferences", "Static"]
git-tree-sha1 = "af9ab7d1f70739a47f03be78771ebda38c3c71bf"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.18"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XML2_jll", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "3d468106a05408f9f7b6f161d9e7715159af247b"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.12.2+0"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LRUCache]]
git-tree-sha1 = "5519b95a490ff5fe629c4a7aa3b3dfc9160498b3"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.6.2"
weakdeps = ["Serialization"]

    [deps.LRUCache.extensions]
    SerializationExt = ["Serialization"]

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.11.1+1"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

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
git-tree-sha1 = "c105fe467859e7f6e9a852cb15cb4301126fac07"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.11"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "e214f2a20bdd64c04cd3e4ff62d3c9be7e969a59"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.5.4+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3cce3511ca2c6f87b19c34ffc623417ed2798cbd"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.10+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bc95bf4149bf535c09602e3acdf950d9b4376227"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.NCDatasets]]
deps = ["CFTime", "CommonDataModel", "DataStructures", "Dates", "DiskArrays", "NetCDF_jll", "NetworkOptions", "Printf"]
git-tree-sha1 = "c82c73e2e0c57a0fe13d3414d7c5a6a821d24016"
uuid = "85f8d34a-cbdd-5861-8df4-14fed0d494ab"
version = "0.14.10"

    [deps.NCDatasets.extensions]
    NCDatasetsMPIExt = "MPI"

    [deps.NCDatasets.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

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

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML", "Zlib_jll"]
git-tree-sha1 = "ab6596a9d8236041dcd59b5b69316f28a8753592"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "5.0.9+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "f1a7e086c677df53e064e0fdd2c9d0b0833e3f6e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.5.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.1+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1f7f9bbd5f7a2e5a9f7d96e51c9754454ea7f60b"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.4+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "12ce661880f8e309569074a61d3767e5756a199f"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.1"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

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
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "34f7e5d2861083ec7596af8b8c092531facf2192"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.8.2+2"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "da7adf145cce0d44e892626e647f9dcbe9cb3e10"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.8.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "9eca9fc3fe515d619ce004c83c31ffd3f85c7ccf"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.8.2+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "8f528b0851b5b7025032818eb5abbeb8a736f853"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.8.2+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

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
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SciMLPublic]]
git-tree-sha1 = "ed647f161e8b3f2973f24979ec074e8d084f1bee"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "49440414711eddc7227724ae6e570c7d5559a086"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.3.1"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
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
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2c962245732371acd51700dbb268af311bddd719"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.6"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "d969183d3d244b6c33796b5ed01ab97328f2db85"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.5"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

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

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "d1d9a935a26c475ebffd54e9c7ad11627c43ea85"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.72"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1aa23f01927b2dac46db77a56b31088feee0a491"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.1.4+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.13.1+1"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libzip_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "OpenSSL_jll", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "86addc139bca85fdf9e7741e10977c45785727b7"
uuid = "337d8026-41b4-5cde-a456-74a10e5b31d1"
version = "1.11.3+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.5.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "fbf139bce07a534df0e699dbb5f5cc9346f95cc1"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.9.2+0"
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
# ╟─b96052b9-4272-4d78-9e43-927bc872c43a
# ╟─60f0f581-0bcf-40e7-9a45-86eaceba9040
# ╟─dc0f1a5e-1e7b-4a9d-816e-731d2747ac4e
# ╟─19f106e4-2b82-47e1-9284-799a105f30cb
# ╟─54e63724-44ca-42b9-9246-7afb271b8154
# ╟─32b531ab-ee71-4685-af65-a2bae0b868f6
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
# ╟─a291a9c2-b89a-4915-8a19-3f900647fbbb
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
# ╟─f93e3556-5016-47c0-8ecb-42a9f922d787
# ╟─5842dab4-b354-40a4-9057-a3752d525c10
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
# ╠═606032a2-b2ca-4fd8-9930-afd83aecec7a
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
# ╟─9d3eea3d-1ba6-4650-9187-85fa5aeca210
# ╟─3d93d6ce-6f7a-4083-ade9-89ba0058bae0
# ╟─47d3626b-463c-4090-b12a-5ed5f5e8ce62
# ╟─93913dcc-aae0-40fe-9837-92e9a251aadf
# ╟─3f7349aa-da11-4d27-a637-5286950d3244
# ╟─39afec69-f90d-4ff2-99d1-5cfec84d09a1
# ╟─3f7349aa-da11-4d27-a637-5286950d3245
# ╟─3f7349aa-da11-4d27-a637-5286950d3246
# ╟─b0d6a9e3-5c3a-4408-ab87-1067c2dfaeb7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
