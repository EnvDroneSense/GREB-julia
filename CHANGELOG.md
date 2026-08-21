# Changelog

Notable changes to GREBClimate.jl, in roughly chronological order. Loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) conventions;
since the package hasn't had its first registered release yet, entries are
grouped by development phase rather than version number. For the detailed,
pass-by-pass discovery narrative behind any entry below — which bug, how it
was found, how it was verified against the Fortran reference — see
[`claude/AUDIT_LOG.md`](claude/AUDIT_LOG.md).

## [Unreleased] — toward v1.0.0 / first registration

### Renamed
- Package `GREB` → `GREBClimate`; repository `GREB-julia` → `GREBClimate.jl`
  (the General registry's AutoMerge requires package names of at least 5
  characters).

### Added
- Documenter.jl site — tutorial, API reference, and a physics-switches guide
  — deployed to GitHub Pages.
- CI: matrix over Julia `1.10`/`1`, split into `light`/`heavy` test shards.
- A benchmark suite (`benchmark/run_benchmarks.jl`) with `year`/`stages`/
  `threads`/`alloc` modes.
- `ClimateFields`/`ModelState`/`SurfaceState` structs, replacing ~40 mutable
  module-level globals with explicit, passed-in state.
- `RunSpec` for flux-correction/control/scenario run lengths.
- All four IPCC RCP scenarios and five SSP scenarios, a historical CO₂
  (1850–2017) hindcast, and a user-supplied CO₂ trajectory (`:custom_co2`).
- "Deconstruct" experiment presets (`:decon_mean_climate`, `:decon_2xco2`)
  toggling individual feedback processes.
- 3-way threaded physics: the temperature and humidity `circulation!` calls
  run concurrently via `Threads.@spawn`.
- `Float32` compute path throughout, matching the JLD2 input data's native
  precision.

### Fixed
18 correctness bugs found by direct comparison against the Fortran
reference (`greb.model.mscm.f90`); the ones with the widest-reaching impact:

- `circulation!` left a stale moisture-convergence term in every
  temperature sub-step of every default run.
- `set_hydrology_parameters!` wrote to disconnected module globals instead
  of the config struct — the `log_rain` switch had zero effect on any run.
- `hydro!`'s `log_eva` modes `1`/`2` silently duplicated mode `-1`'s formula
  instead of using their own coefficients.
- `deep_ocean!` cut off ocean-atmosphere heat exchange under sea ice and in
  high-latitude winters.
- 19 of 36 exported functions had docstrings silently detached from their
  definitions by a comment-placement pattern that breaks Julia's `@doc`
  binding — found while wiring up the Documenter.jl site.

The remaining 13 fixes, plus one investigated-and-reverted finding, are
detailed in `claude/AUDIT_LOG.md`.

### Performance
- ~2.31× faster per simulated year (2.7s → 1.17s) from 3-way threading of
  `tendencies!` combined with 4 `@turbo`-rewritten hot-loop kernels.
- `Float32` throughout: a further ~1.6× on top of threading, with output
  validated against the previous `Float64` path to well under 0.01 K.
- Flux-correction files (`Tsurf`/`vapour`/`Tocean` corrections) merged into
  a single `flux_corrections.jld2` — ~35% faster to load, no size penalty.

### Changed
- `src/GREB.jl`, originally a single 2,245-line file, split into topical
  files (`constants`, `config`, `state`, `io`, `physics/`, `circulation`,
  `tendencies`, `output`, `postprocess`, `model`).
- Test suite split into `light`/`heavy` shards (`GREB_TEST_SHARD` env var),
  matching the CI matrix.

## [0.1.0] - 2026-08-06
Initial extraction from the interactive Pluto notebook into a standard Julia
package layout (`Project.toml`, `src/`, `test/`).
