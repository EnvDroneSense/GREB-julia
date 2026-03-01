# GREB Model Input Data Guide

This document describes the input data required to run the GREB climate model.

## 📁 Data Directory Structure

Create the following directory structure in your project:

```
ClimaModel/
└── Data/
    └── input/
        ├── [climate data files - see below]
        └── solar_forcing_scenarios/
            └── [optional solar forcing files]
```

## 📋 Required Input Files

### Static 2D Fields

| File | Size | Description |
|------|------|-------------|
| `global.topography.bin` | 96×48 | Global topography (m), ocean points < 0 |
| `greb.glaciers.bin` | 96×48 | Glacier mask (0 or 1) |

### 3D Climatology Fields (96×48×730)

All climatology files contain 730 time steps (12-hour intervals over one year).

#### NCEP Dataset Files
| File | Description |
|------|-------------|
| `ncep.tsurf.1948-2007.clim.bin` | Surface temperature climatology (K) |
| `ncep.zonal_wind.850hpa.clim.bin` | Zonal wind at 850 hPa (m/s) |
| `ncep.meridional_wind.850hpa.clim.bin` | Meridional wind at 850 hPa (m/s) |
| `ncep.atmospheric_humidity.clim.bin` | Atmospheric specific humidity (kg/kg) |
| `ncep.soil_moisture.clim.bin` | Soil moisture fraction (0-1) |

#### ERA-Interim Dataset Files (Alternative)
| File | Description |
|------|-------------|
| `erainterim.tsurf.1979-2015.clim.bin` | Surface temperature climatology (K) |
| `erainterim.zonal_wind.850hpa.clim.bin` | Zonal wind at 850 hPa (m/s) |
| `erainterim.meridional_wind.850hpa.clim.bin` | Meridional wind at 850 hPa (m/s) |
| `erainterim.atmospheric_humidity.clim.bin` | Atmospheric specific humidity (kg/kg) |
| `erainterim.omega.vertmean.clim.bin` | Vertical velocity (Pa/s) |
| `erainterim.omega_std.vertmean.clim.bin` | Vertical velocity std deviation |
| `erainterim.windspeed.850hpa.clim.bin` | Wind speed at 850 hPa (m/s) |
| `erainterim.evaporation.clim.bin` | Evaporation rate |

#### Common Files (Required for Both Datasets)
| File | Description |
|------|-------------|
| `isccp.cloud_cover.clim.bin` | Cloud cover fraction (0-1) from ISCCP |
| `woce.ocean_mixed_layer_depth.clim.bin` | Ocean mixed layer depth (m) |
| `Tocean.clim.bin` | Deep ocean temperature (K) |

### MSCM-Specific Fields (96×48×730)

| File | Description |
|------|-------------|
| `cmip5.omega.rcp85.ensmean.forcing.new.bin` | Vertical velocity from CMIP5 |
| `cmip5.omegastd.rcp85.ensmean.forcing.new.bin` | Vertical velocity std dev |

### Flux Correction Fields (96×48×730)

| File | Description |
|------|-------------|
| `Tsurf_flux_correction.bin` | Surface temperature flux correction (W/m²) |
| `vapour_flux_correction.bin` | Water vapor flux correction (kg/m²/s) |
| `Tocean_flux_correction.bin` | Deep ocean flux correction (W/m²) |

### Solar Radiation (48×730)

| File | Description |
|------|-------------|
| `solar_radiation.clim.bin` | 24-hour mean solar radiation at top of atmosphere (W/m²) |

### Scenario Forcing Files

Text files containing time series of forcing values:

| File | Description |
|------|-------------|
| `ipcc.scenario.rcp26.forcing.txt` | RCP 2.6 scenario (low emissions) |
| `ipcc.scenario.rcp45.forcing.txt` | RCP 4.5 scenario (moderate emissions) |
| `ipcc.scenario.rcp6.forcing.txt` | RCP 6.0 scenario |
| `ipcc.scenario.rcp85.forcing.txt` | RCP 8.5 scenario (high emissions) |

### Optional: Solar Forcing Scenarios

Located in `solar_forcing_scenarios/` subdirectory:

- Eccentricity variations: `greb.solar.eccentricity.{0-60}.bin`
- Obliquity variations: `greb.solar.obliquity.{0-230}.bin` (steps of 5)
- Paleoclimate: `greb.solar.231K_hybers.corrected.bin`

## 🔧 Data Format Specifications

### GrADS Binary Format

All `.bin` files are in GrADS binary format:
- **Data type**: 32-bit floating point (Float32)
- **Byte order**: Little-endian (standard)
- **Storage order**: Column-major (Fortran-style): longitude, latitude, time

### Dimensions

- **Longitude**: 96 points (3.75° resolution)
- **Latitude**: 48 points (3.75° resolution)
- **Time**: 730 steps (12-hour intervals, one year)

### Control Files (.ctl)

Some datasets include `.ctl` files that describe the binary file structure for GrADS. These are not required by the Julia code but can be useful for verification.

## ✅ Verification

To verify your data files are correctly formatted, the model will print diagnostic information when loading:

```julia
load_greb_input_data!("Data/input"; dataset=:ncep)
```

Expected output:
```
═══════════════════════════════════════════════════════════
  Loading GREB Climate Data
  Directory: Data/input
  Dataset: ncep
═══════════════════════════════════════════════════════════

[1/4] Loading static 2D fields...
  • Topography... ✓
  • Glacier mask... ✓

[2/4] Loading 3D climatology fields (96×48×730)...
  • Surface temperature... ✓
  ...
```

## 🐛 Troubleshooting

### File Not Found Errors

- Ensure all required files are in `Data/input/`
- Check file names match exactly (case-sensitive on Linux/Mac)
- Verify file extensions are `.bin` (not `.BIN` or `.binary`)

### Size Mismatch Errors

- Verify files have correct dimensions (96×48×730 for 3D fields)
- Check that files are in Float32 format
- Ensure no header bytes in binary files

### NaN or Unrealistic Values

- Check that ocean/land masks are correctly applied
- Verify flux correction files are available
- Ensure topography file has negative values for ocean points

---

**Last Updated**: March 2026
