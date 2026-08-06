# ── notebook cell b303e4e9-49fa-45ad-967e-20f165fdf38c  (orig lines 73-112) ──
"""
    read_jld2(filepath::String)

Read a `.jld2` field file written by `scripts/convert_greb_to_jld2.jl`.

# Returns
- named tuple `(data, dim_names, coords, ctl)` where:
  - `data`: Array{Float32} with shape as stored
  - `dim_names`: Vector{String} of dimension names (e.g., ["lon", "lat", "time"])
  - `coords`: `Dict{Int,Vector{Float64}}` of physical coordinate values per
    dimension index, or `nothing` if the file has none
  - `ctl`: raw GrADS `.ctl` metadata text, or `nothing` if the file has none
"""
function read_jld2(filepath::String)
    jldopen(filepath, "r") do file
        return (
            data=file["data"],
            dim_names=file["dim_names"],
            coords=haskey(file, "coords") ? file["coords"] : nothing,
            ctl=haskey(file, "ctl") ? file["ctl"] : nothing,
        )
    end
end

# ── notebook cell 8578d6aa-2782-4279-8f6b-78194b8ecc10  (orig lines 113-140) ──
function load_solar_forcing_jld2(jld2_dir::String, forcing_type::Symbol, index::Int=0)

    if forcing_type == :paleo
        filepath = joinpath(jld2_dir, "solar_scenarios", "solar_paleo.jld2")
        result = read_jld2(filepath)
        return result.data

    elseif forcing_type == :eccentricity
        filepath = joinpath(jld2_dir, "solar_scenarios", "solar_eccentricity.jld2")
        result = read_jld2(filepath)
        values = Int.(result.coords[1])
        pos = findfirst(==(index), values)
        @assert pos !== nothing "Eccentricity index $index not found in $(values)"
        return result.data[pos, :, :]

    elseif forcing_type == :obliquity
        filepath = joinpath(jld2_dir, "solar_scenarios", "solar_obliquity.jld2")
        result = read_jld2(filepath)
        values = Int.(result.coords[1])
        pos = findfirst(==(index), values)
        @assert pos !== nothing "Obliquity index $index not found in $(values)"
        return result.data[pos, :, :]

    else
        error("Unknown forcing type: $forcing_type. Use :paleo, :eccentricity, or :obliquity")
    end
end

# ── notebook cell f578f25e-047e-4a7e-8483-d544c7b4bec3  (orig lines 750-772) ──
"""Load flux corrections from JLD2 files into `fields` (zeros if missing)"""
function load_flux_corrections_jld2!(jld2_dir::String, fields::ClimateFields)
    correction_files = Dict(
        "Tsurf_flux_correction.jld2" => fields.TF_correct,
        "vapour_flux_correction.jld2" => fields.qF_correct,
        "Tocean_flux_correction.jld2" => fields.ToF_correct
    )

    for (filename, array) in correction_files
        filepath = joinpath(jld2_dir, "climatology", filename)
        if isfile(filepath)
            result = read_jld2(filepath)
            array .= result.data
            println("✅ Loaded $filename")
        else
            fill!(array, 0.0)
            @warn "$filename not found, using zeros"
        end
    end
end

# ── notebook cell 2bf0fe8e-5718-4c1e-863b-85db7b3ae7f3  (orig lines 877-977) ──
"""Load all GREB input data from JLD2 formatted files, returning a fresh [`ClimateFields`](@ref)"""
function load_greb_jld2!(jld2_dir::String; dataset::Symbol=:ncep)
    if !isdir(jld2_dir)
        error("JLD2 directory not found: $jld2_dir")
    end

    fields = ClimateFields()

    # Static 2D files
    println("📂 Loading static fields...")
    topo_result = read_jld2(joinpath(jld2_dir, "static", "global.topography.jld2"))
    fields.z_topo .= topo_result.data

    glacier_result = read_jld2(joinpath(jld2_dir, "static", "greb.glaciers.jld2"))
    fields.glacier .= glacier_result.data

    # Dataset-specific file mapping
    file_map = Dict(
        :ncep => Dict(
            "Tclim" => "ncep.tsurf.1948-2007.clim.jld2",
            "uclim" => "ncep.zonal_wind.850hpa.clim.jld2",
            "vclim" => "ncep.meridional_wind.850hpa.clim.jld2",
            "qclim" => "ncep.atmospheric_humidity.clim.jld2",
            "swetclim" => "ncep.soil_moisture.clim.jld2"
        ),
        :era => Dict(
            "Tclim" => "erainterim.tsurf.1979-2015.clim.jld2",
            "uclim" => "erainterim.zonal_wind.850hpa.clim.jld2",
            "vclim" => "erainterim.meridional_wind.850hpa.clim.jld2",
            "qclim" => "erainterim.atmospheric_humidity.clim.jld2",
            "swetclim" => "ncep.soil_moisture.clim.jld2"
        )
    )

    # Use mixed dataset as fallback
    files = get(file_map, dataset, file_map[:ncep])

    println("📂 Loading 3D climatology ($dataset dataset)...")
    climatology_dir = joinpath(jld2_dir, "climatology")

    # Load each variable individually with unique result names
    tsurf_result = read_jld2(joinpath(climatology_dir, files["Tclim"]))
    fields.Tclim .= tsurf_result.data

    uwind_result = read_jld2(joinpath(climatology_dir, files["uclim"]))
    fields.uclim .= uwind_result.data

    vwind_result = read_jld2(joinpath(climatology_dir, files["vclim"]))
    fields.vclim .= vwind_result.data

    humid_result = read_jld2(joinpath(climatology_dir, files["qclim"]))
    fields.qclim .= humid_result.data

    swet_result = read_jld2(joinpath(climatology_dir, files["swetclim"]))
    fields.swetclim .= swet_result.data

    # Common climatology files
    println("📂 Loading common climatology fields...")

    cld_result = read_jld2(joinpath(climatology_dir, "isccp.cloud_cover.clim.jld2"))
    fields.cldclim .= cld_result.data

    mld_result = read_jld2(joinpath(climatology_dir, "woce.ocean_mixed_layer_depth.clim.jld2"))
    fields.mldclim .= mld_result.data

    tocean_result = read_jld2(joinpath(climatology_dir, "Tocean.clim.jld2"))
    fields.Toclim .= tocean_result.data

    omega_result = read_jld2(joinpath(climatology_dir, "erainterim.omega.vertmean.clim.jld2"))
    fields.omegaclim .= omega_result.data

    omegastd_result = read_jld2(joinpath(climatology_dir, "erainterim.omega_std.vertmean.clim.jld2"))
    fields.omegastdclim .= omegastd_result.data

    ws_result = read_jld2(joinpath(climatology_dir, "erainterim.windspeed.850hpa.clim.jld2"))
    fields.wsclim .= ws_result.data

    # Solar radiation (special: lat × time)
    println("📂 Loading solar radiation...")
    solar_path = joinpath(jld2_dir, "solar", "solar_radiation.clim.jld2")
    if isfile(solar_path)
        solar_result = read_jld2(solar_path)
        @assert size(solar_result.data) == (ydim, nstep_yr) "Wrong solar dimensions"
        fields.sw_solar .= solar_result.data
    else
        error("Solar radiation file not found: $solar_path")
    end

    # Optional: Load flux corrections
    println("📂 Loading flux corrections...")
    load_flux_corrections_jld2!(jld2_dir, fields)

    # Update wind sign splits
    @. fields.uclim_m = ifelse(fields.uclim >= 0.0, fields.uclim, 0.0)
    @. fields.uclim_p = ifelse(fields.uclim < 0.0, fields.uclim, 0.0)
    @. fields.vclim_m = ifelse(fields.vclim >= 0.0, fields.vclim, 0.0)
    @. fields.vclim_p = ifelse(fields.vclim < 0.0, fields.vclim, 0.0)

    println("✅ All GREB data loaded successfully from JLD2")
    return fields
end
