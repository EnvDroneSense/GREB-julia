### Convert GREB binary files to JLD2 (JuliaIO/JLD2.jl) files ###
#
# Converts the raw GREB `.bin` input files (see DATA_README.md for the
# expected layout, normally under `Data/input/`) into the `.jld2` dataset
# that `load_greb_jld2!`/`GREB.read_jld2` read (see README.md's "Input Data"
# section for the resulting directory structure).
#
# Each field is written as its own `.jld2` file storing plain Julia values
# under the keys "data" (Array{Float32}), "dim_names", and optionally
# "coords" (physical coordinate values per dimension, e.g. the actual
# eccentricity/obliquity value behind a stacked scenario slice) and "ctl"
# (the original GrADS .ctl metadata text, if present) — no custom binary
# format, no header to hand-parse.
#
# Solar forcing scenario families (eccentricity, obliquity) are discovered
# by globbing `greb.solar.<prefix>.<N>.bin` and sorting numerically, NOT by
# a hardcoded stride — an earlier version of this pipeline (JDAL2) hardcoded
# `0:25:230` for obliquity, which only matched 10 of the 47 files that
# actually exist (steps of 5, not 25), silently discarding ~79% of the
# obliquity scenario data. Globbing means this can't happen again even if
# the input file set changes.
#
# Usage:
#   julia --project=. scripts/convert_greb_to_jld2.jl [input_dir] [output_dir]
#   # defaults: input_dir  = Data/input        (see DATA_README.md)
#   #           output_dir = greb_dataset_jld2  (see README.md)

using JLD2

# ============================================================================
# CONSTANTS
# ============================================================================

const LON, LAT, TIME = 96, 48, 730
const SIZE_2D = LON * LAT * 4
const SIZE_3D = LON * LAT * TIME * 4
const SIZE_SOLAR = LAT * TIME * 4

# ============================================================================
# CLASSIFICATION
# ============================================================================

function classify_by_size(size_bytes::Int)
    if size_bytes == SIZE_2D
        return (layout="lon_lat", dims=(LON, LAT), dim_names=["lon", "lat"])
    elseif size_bytes == SIZE_3D
        return (layout="lon_lat_time", dims=(LON, LAT, TIME), dim_names=["lon", "lat", "time"])
    elseif size_bytes == SIZE_SOLAR
        return (layout="lat_time", dims=(LAT, TIME), dim_names=["lat", "time"])
    else
        error("Unknown file size: $size_bytes")
    end
end

# ============================================================================
# RAW BINARY READER
# ============================================================================

function read_bin(bin_path::String, dims::Tuple)
    data = Vector{Float32}(undef, prod(dims))
    open(bin_path, "r") do io; read!(io, data); end
    return reshape(data, dims)
end

# ============================================================================
# JLD2 FIELD WRITER / READER
# ============================================================================

function write_jld2_field(filepath::String, data, dim_names::Vector{String};
                           coords::Union{Nothing,Dict{Int,Vector{Float64}}}=nothing,
                           ctl::Union{Nothing,String}=nothing)
    mkpath(dirname(filepath))
    jldopen(filepath, "w") do file
        file["data"] = data
        file["dim_names"] = dim_names
        coords !== nothing && !isempty(coords) && (file["coords"] = coords)
        ctl !== nothing && (file["ctl"] = ctl)
    end
end

function read_jld2_field(filepath::String)
    jldopen(filepath, "r") do file
        return (
            data=file["data"],
            dim_names=file["dim_names"],
            coords=haskey(file, "coords") ? file["coords"] : nothing,
            ctl=haskey(file, "ctl") ? file["ctl"] : nothing,
        )
    end
end

function read_ctl_text(bin_path::String)
    ctl_path = splitext(bin_path)[1] * ".ctl"
    return isfile(ctl_path) ? read(ctl_path, String) : nothing
end

# ============================================================================
# BATCH CONVERSION (main dataset)
# ============================================================================

function convert_all(input_path::String, output_dir::String)
    println("🔄 CONVERTING TO JLD2")
    println("="^50)
    println("Input:  $input_path")
    println("Output: $output_dir\n")

    bin_files = filter(endswith(".bin"), readdir(input_path; join=true))
    println("Found $(length(bin_files)) .bin files\n")

    success, failed, total_bytes = 0, 0, 0

    for bin_path in sort(bin_files)
        name = splitext(basename(bin_path))[1]
        size_bytes = stat(bin_path).size
        class = classify_by_size(size_bytes)

        if occursin(r"topography|glacier", name)
            out_path = joinpath(output_dir, "static", "$name.jld2")
        elseif occursin(r"solar_radiation", name)
            out_path = joinpath(output_dir, "solar", "$name.jld2")
        else
            out_path = joinpath(output_dir, "climatology", "$name.jld2")
        end

        try
            arr = read_bin(bin_path, class.dims)
            write_jld2_field(out_path, arr, class.dim_names; ctl=read_ctl_text(bin_path))
            total_bytes += size_bytes
            println("  ✓ $name → $(size_bytes ÷ 1024) KB [$(join(class.dim_names, "×"))]")
            success += 1
        catch e
            println("  ❌ $name: $e")
            failed += 1
        end
    end

    println("\n" * "="^50)
    println("✅ CONVERSION COMPLETE")
    println("  Successful: $success files")
    println("  Failed: $failed files")
    println("  Total data: $(round(total_bytes / 1e6, digits=1)) MB")
    println("  Output: $output_dir")
end

# ============================================================================
# SOLAR SCENARIOS (GROUPED) - discovered by glob, not a hardcoded stride
# ============================================================================

function read_solar_bin(filepath::String)
    return read_bin(filepath, (LAT, TIME))
end

"""
Find all files matching `greb.solar.<prefix>.<number>.bin` in `dir`, and
return them sorted by the numeric value (not lexicographically), along
with the parsed numeric values. Errors if no files match, since a missing
prefix almost certainly means the input layout changed.
"""
function discover_indexed_scenario_files(dir::String, prefix::String)
    pattern = Regex("^greb\\.solar\\.$(prefix)\\.(\\d+)\\.bin\$")
    matches = Tuple{Float64,String}[]
    for fname in readdir(dir)
        m = match(pattern, fname)
        m === nothing && continue
        push!(matches, (parse(Float64, m.captures[1]), joinpath(dir, fname)))
    end
    isempty(matches) && error("No files found matching greb.solar.$(prefix).<N>.bin in $dir")
    sort!(matches, by=first)
    values = first.(matches)
    paths = last.(matches)
    return values, paths
end

function convert_scenario_family(dir::String, prefix::String, output_path::String, index_name::String)
    values, paths = discover_indexed_scenario_files(dir, prefix)

    stack = Vector{Any}(undef, length(paths))
    missing_idx = Int[]
    for (i, path) in enumerate(paths)
        if isfile(path)
            stack[i] = read_solar_bin(path)
        else
            push!(missing_idx, i)
        end
    end

    if !isempty(missing_idx)
        error("$(length(missing_idx))/$(length(paths)) expected $prefix files went missing between listing and reading: " *
              join(paths[missing_idx], ", "))
    end

    arr_3d = permutedims(cat(stack..., dims=3), [3, 1, 2])
    write_jld2_field(output_path, arr_3d, [index_name, "lat", "time"];
                      coords=Dict(1 => values))
    return length(paths)
end

function convert_solar_scenarios(input_path::String, output_dir::String)
    println("🔄 CONVERTING SOLAR SCENARIOS (GROUPED)")
    println("="^60)

    solar_dir = joinpath(input_path, "solar_forcing_scenarios")
    if !isdir(solar_dir)
        println("  ⚠ $solar_dir not found — skipping solar scenarios (optional, see DATA_README.md)")
        return
    end
    mkpath(output_dir)

    # Paleo (single file, not an indexed family)
    paleo_path = joinpath(solar_dir, "greb.solar.231K_hybers.corrected.bin")
    if isfile(paleo_path)
        try
            write_jld2_field(joinpath(output_dir, "solar_paleo.jld2"),
                              read_solar_bin(paleo_path), ["lat", "time"];
                              ctl=read_ctl_text(paleo_path))
            println("  ✓ Paleo (48×730)")
        catch e
            println("  ❌ Paleo: $e")
        end
    else
        println("  ⚠ Paleo file not found, skipping: $paleo_path")
    end

    # Eccentricity (all greb.solar.eccentricity.<N>.bin files, whatever N range exists)
    try
        n = convert_scenario_family(solar_dir, "eccentricity",
                                     joinpath(output_dir, "solar_eccentricity.jld2"), "ecc_index")
        println("  ✓ Eccentricity ($(n)×48×730)")
    catch e
        println("  ❌ Eccentricity: $e")
    end

    # Obliquity (all greb.solar.obliquity.<N>.bin files)
    try
        n = convert_scenario_family(solar_dir, "obliquity",
                                     joinpath(output_dir, "solar_obliquity.jld2"), "obl_index")
        println("  ✓ Obliquity ($(n)×48×730)")
    catch e
        println("  ❌ Obliquity: $e")
    end

    println("\n✅ Solar scenarios converted to 3 grouped files")
end

# ============================================================================
# VERIFICATION
# ============================================================================

function verify(output_dir::String)
    println("\n" * "="^60)
    println("🔍 VERIFICATION")
    println("="^60)

    for (name, path) in [
        ("Static", joinpath(output_dir, "static", "global.topography.jld2")),
        ("Solar", joinpath(output_dir, "solar", "solar_radiation.clim.jld2")),
        ("3D", joinpath(output_dir, "climatology", "ncep.tsurf.1948-2007.clim.jld2")),
        ("Paleo", joinpath(output_dir, "solar_scenarios", "solar_paleo.jld2")),
        ("Eccentricity", joinpath(output_dir, "solar_scenarios", "solar_eccentricity.jld2")),
        ("Obliquity", joinpath(output_dir, "solar_scenarios", "solar_obliquity.jld2")),
    ]
        if isfile(path)
            result = read_jld2_field(path)
            println("✅ $name: $(size(result.data)) → $(join(result.dim_names, "×"))")
        end
    end

    println("\n✅ Done!")
end

# ============================================================================
# ENTRY POINT
# ============================================================================

function main(input_path::String, output_dir::String)
    if !isdir(input_path)
        error("Input directory not found: $input_path (see DATA_README.md for the expected layout)")
    end
    convert_all(input_path, output_dir)
    convert_solar_scenarios(input_path, joinpath(output_dir, "solar_scenarios"))
    verify(output_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    input_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "Data", "input")
    output_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "..", "greb_dataset_jld2")
    main(input_path, output_dir)
end
