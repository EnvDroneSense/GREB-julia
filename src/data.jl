# =============================================================================
# data.jl - locating the JLD2 input dataset.
#
# The dataset is ~439 MB and is not shipped with the package. `greb_data_dir()`
# resolves where it lives, falling back to a DataDeps.jl download.
#
# =============================================================================

"""
Name of the DataDep, and the folder name the dataset is cached under inside
DataDeps' storage (`~/.julia/scratchspaces/<uuid>/datadeps/` by default).
"""
const DATA_DEP_NAME = "GREB-input-data"

"""
Release tag and asset name for the published dataset bundle.

The bundle is built from a validated `greb_input_data/` tree with
`tools/package_dataset.jl`, which also prints the SHA256 below. 
"""
const DATA_RELEASE_TAG = "data-v1"
const DATA_ARCHIVE_NAME = "greb_input_data-v1.tar.gz"

const DATA_URL = "https://github.com/EnvDroneSense/GREBClimate.jl/releases/download/" *
                 "$DATA_RELEASE_TAG/$DATA_ARCHIVE_NAME"

"""
SHA256 of `DATA_ARCHIVE_NAME`, 370563584 bytes.

Reproducible: `tools/package_dataset.jl` builds the archive with sorted entries,
zeroed owner/group, **pinned entry timestamps** and `gzip -n`, so the archive
depends only on the dataset's contents. Regenerating the `.jld2` tree from the
raw `.bin` files and repackaging reproduces this exact hash.

(The timestamp pin matters: without it tar stored each file's mtime, so a
regenerated-but-byte-identical dataset produced a different archive.)
"""
const DATA_SHA256 = "a9799ecb2e50d6f01517c69e2a7c6a8f646887123271d83eb3538990597f65c5"

const _DATA_MESSAGE = """
The GREB input dataset (~353 MB download, ~439 MB unpacked).

Derived from the input files of the original GREB model by Dietmar Dommenget
and colleagues (Monash University), converted to JLD2. It combines fields from
several upstream sources, each with its own terms of use and acknowledgement
requirements:

  NCEP/NCAR reanalysis      surface temperature, winds, humidity, soil moisture
  ERA-Interim (ECMWF)       surface temperature, winds, humidity, omega
  ISCCP                     cloud cover
  WOCE                      ocean mixed-layer depth
  CMIP5 (RCP8.5 ens. mean)  climate-change anomaly forcing
  IPCC                      RCP/SSP CO2 scenario tables

"""

"""
    greb_data_dir(path = nothing; allow_download = true) -> String or nothing

Return the directory holding the JLD2 input dataset, resolving in this order:

1. `path`, if given and non-empty.
2. `ENV["GREB_DATA"]`, if set.
3. `greb_input_data/` next to the package, if it exists.
4. The `$DATA_DEP_NAME` DataDep - downloading it on first use, after asking.

Pass `allow_download = false` to stop after step 3 and return `nothing` when no
local dataset is found. Test suites and benchmarks use this so that running them
can never pull 353 MB over the network as a side effect.

The result is a plain path, suitable for [`load_greb_jld2!`](@ref) and
`greb_model!`'s `jld2_dir`:

```julia
dir    = greb_data_dir()
fields = load_greb_jld2!(dir; dataset = :ncep)
result = greb_model!(RunSpec(), cfg; jld2_dir = dir, fields = fields)
```
"""
function greb_data_dir(path::Union{Nothing,AbstractString} = nothing;
                       allow_download::Bool = true)
    if path !== nothing && !isempty(path)
        isdir(path) || error("GREB dataset directory not found: $path")
        return String(path)
    end

    env = get(ENV, "GREB_DATA", "")
    if !isempty(env)
        isdir(env) || error("GREB_DATA points at a directory that does not exist: $env")
        return env
    end

    local_dir = normpath(joinpath(@__DIR__, "..", "greb_input_data"))
    isdir(local_dir) && return local_dir

    allow_download || return nothing

    # Resolved lazily: this is the only branch that can trigger a download, and
    # it must never run at precompile time
    return @datadep_str DATA_DEP_NAME
end

"""
    register_greb_datadep()

Register the dataset with DataDeps. Called from `GREBClimate.__init__`;
registration itself is cheap and downloads nothing.
"""
function register_greb_datadep()
    register(DataDep(DATA_DEP_NAME, _DATA_MESSAGE, DATA_URL, DATA_SHA256;
                     post_fetch_method = unpack))
    return nothing
end
