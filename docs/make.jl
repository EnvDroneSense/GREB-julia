# =============================================================================
# docs/make.jl — builds the GREBClimate.jl Documenter site.
#
# Local build (once, to link the docs env to the local package source):
#   julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#   julia --project=docs docs/make.jl
#
# CI (see .github/workflows/docs.yml) uses julia-actions/julia-docdeploy@v1,
# which runs the equivalent `Pkg.develop` step automatically before this file.
# =============================================================================

using Documenter
using GREBClimate

makedocs(
    sitename = "GREBClimate.jl",
    modules = [GREBClimate],
    authors = "Thomas Struys",
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Physics Switches" => "switches.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/EnvDroneSense/GREBClimate.jl.git",
    devbranch = "main",
)
