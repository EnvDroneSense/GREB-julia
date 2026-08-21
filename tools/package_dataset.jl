# =============================================================================
# package_dataset.jl - build the distributable dataset archive and its SHA256.
#
# MAINTAINER TOOL. Produces the `.tar.gz` that `src/data.jl`'s DataDep points
# at, and prints the hash to paste into `DATA_SHA256`.
#
#   julia --project=. tools/package_dataset.jl [dataset_dir] [output_path]
#   # defaults: dataset_dir = greb_input_data
#   #           output_path = greb_input_data-v1.tar.gz (in the current dir)
#
# The build is reproducible: entries are sorted, owner/group are zeroed, and
# gzip runs with -n so no timestamp is embedded. Rebuilding from an identical
# tree gives a byte-identical archive, so the recorded SHA256 stays valid.
#
# Requires `tar`, `gzip` and `sha256sum` (Git Bash / WSL / any Unix).
# =============================================================================

using SHA

const REPO = normpath(joinpath(@__DIR__, ".."))

# ── validate the tree before packaging ───────────────────────────────────────
"""
Check `dir` against the converter's allowlist: every field the model reads must
be present, and nothing else may be. Returns the file count and total bytes.
"""
function validate_dataset(dir::AbstractString)
    isdir(dir) || error("dataset directory not found: $dir")

    # MODEL_FIELD_NAMES lives in the converter; read it rather than duplicating.
    conv = read(joinpath(REPO, "tools", "convert_greb_to_jld2.jl"), String)
    m = match(r"const MODEL_FIELD_NAMES = Set\{String\}\(\[(.*?)\n\]\)"s, conv)
    m === nothing && error("could not parse MODEL_FIELD_NAMES from the converter")
    body = m.captures[1]
    cut = findfirst("(\"erainterim.", body)
    cut === nothing || (body = body[1:first(cut)-1])

    expected = Set{String}()
    for lit in eachmatch(r"\"([^\"]+)\"", body)
        s = lit.captures[1]
        occursin('$', s) && continue
        occursin('.', s) && !occursin(' ', s) && push!(expected, s)
    end
    for f in ("tsurf", "zonal.wind", "meridional.wind", "windspeed", "omega"),
        s in ("elnino", "lanina")
        push!(expected, "erainterim.$f.$s.forcing")
    end

    # Combined multi-field files the converter writes but that are not per-field
    # allowlist entries.
    combined = Set(["flux_corrections", "ipcc_scenarios", "solar_paleo",
                    "solar_eccentricity", "solar_obliquity",
                    "historical_emissions_population"])

    present, nbytes, nfiles = Set{String}(), 0, 0
    for (root, _, files) in walkdir(dir), f in files
        endswith(f, ".jld2") || continue
        push!(present, first(splitext(f)))
        nbytes += stat(joinpath(root, f)).size
        nfiles += 1
    end

    extra = setdiff(present, union(expected, combined))
    absent = setdiff(expected, present)
    isempty(absent) || error("dataset is missing $(length(absent)) field(s) the model reads:\n  " *
                             join(sort(collect(absent)), "\n  "))
    if !isempty(extra)
        error("""
              dataset contains $(length(extra)) file(s) the model never reads:
                $(join(sort(collect(extra)), "\n  "))
              Regenerate without --all, or add them to MODEL_FIELD_NAMES if they
              are now genuinely used. See .claude/notes/data-distribution.md.
              """)
    end
    return nfiles, nbytes
end

function build_archive(dir::AbstractString, out::AbstractString)
    # Reproducible tar: sorted names, no owner info, gzip without a timestamp.
    cmd = pipeline(`tar --sort=name --owner=0 --group=0 --numeric-owner -cf - -C $dir .`,
                   `gzip -n -6`)
    open(out, "w") do io
        run(pipeline(cmd; stdout = io))
    end
    return out
end

# Read the expected tag out of src/data.jl rather than hardcoding it twice.
function data_release_tag()
    s = read(joinpath(REPO, "src", "data.jl"), String)
    m = match(r"const DATA_RELEASE_TAG = \"([^\"]+)\"", s)
    m === nothing ? "<tag>" : m.captures[1]
end

function main(dir::AbstractString, out::AbstractString)
    println("Validating $dir ...")
    nfiles, nbytes = validate_dataset(dir)
    println("  ", nfiles, " .jld2 files, ", round(nbytes / 1048576; digits = 1), " MB unpacked")

    println("Building $out (reproducible tar.gz) ...")
    build_archive(dir, out)

    sz = stat(out).size
    hash = bytes2hex(open(sha256, out))
    println()
    println("archive : ", out)
    println("size    : ", sz, " bytes (", round(sz / 1048576; digits = 1), " MB)")
    println("sha256  : ", hash)
    println()
    println("Next steps:")
    println("  1. Update DATA_SHA256 in src/data.jl to the hash above.")
    println("  2. Attach the archive to the '", data_release_tag(), "' GitHub release")
    println("     as '", basename(out), "' (the name must match DATA_ARCHIVE_NAME).")
    println("  3. Verify end-to-end on a clean machine, or by clearing the")
    println("     datadeps cache and calling greb_data_dir() with no local")
    println("     greb_input_data/ and no GREB_DATA set.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    dataset_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(REPO, "greb_input_data")
    output_path = length(ARGS) >= 2 ? ARGS[2] : "greb_input_data-v1.tar.gz"
    main(dataset_dir, output_path)
end
