# Changelog

Notable changes to GREBClimate.jl, in roughly chronological order. Loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) conventions;
since the package hasn't had its first registered release yet, entries are
grouped by development phase rather than version number. For the detailed,
pass-by-pass discovery narrative behind any entry below - which bug, how it
was found, how it was verified against the Fortran reference - see
[`.claude/notes/audit-history.md`](.claude/notes/audit-history.md).

## [Unreleased] - toward v1.0.0 / first registration

### Renamed
- Package `GREB` → `GREBClimate`; repository `GREB-julia` → `GREBClimate.jl`
  (the General registry's AutoMerge requires package names of at least 5
  characters).

### Added
- **Automatic dataset download.** The JLD2 input dataset is now fetched and
  cached on first use via [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl).
  The new exported `greb_data_dir()` resolves an explicit path, then
  `$GREB_DATA`, then a local `greb_input_data/`, and only then downloads - so
  existing setups keep working untouched and offline use is unaffected. Pass
  `allow_download=false` to stop before the network; the test suite and the
  benchmarks do exactly that, so neither can pull 353 MB as a side effect.
  Set `DATADEPS_ALWAYS_ACCEPT=true` for non-interactive sessions.
- `tools/package_dataset.jl` builds the published dataset archive
  reproducibly and prints its SHA256, validating the tree against the
  converter's allowlist first so a stray field cannot enlarge the download.
- `greb_model!` now refuses an uninitialized `ClimateFields` (all-zero
  climatology) rather than silently simulating a physically meaningless
  ~-40 °C world. `ClimateFields` carries a `loaded` flag set by
  `load_greb_jld2!`; genuinely data-free runs - package precompilation and
  config/scenario-plumbing tests - opt in with `allow_uninitialized=true`.
- `tools/convert_greb_to_jld2.jl` filters on an explicit `MODEL_FIELD_NAMES`
  allowlist, so it no longer emits 11 `.jld2` files (~148 MB) that nothing
  reads. `--all` restores the previous convert-everything behaviour. A test
  asserts the allowlist and `src/io.jl`'s loads agree in both directions.
- Documenter.jl site - tutorial, API reference, and a physics-switches guide
  - deployed to GitHub Pages.
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
- **The multithreaded circulation path was never tested.** `tendencies!` runs
  `circulation!(Ta)`/`circulation!(q)` concurrently only when
  `Threads.nthreads() > 1`, and the test process was single-threaded, so that
  branch had only ever been validated by benchmarking. A heavy-shard test now
  spawns `-t 1` and `-t 2` subprocesses and asserts bit-identical monthly
  means, and CI sets `JULIA_NUM_THREADS=2`.
- **README quick-start produced wrong results.** The documented first run
  called `load_greb_jld2!` and discarded its return value, so the model ran on
  a zero climatology and reported a global-mean surface temperature of
  233 K (-40 °C) instead of 277 K (14.8 °C) - while printing
  `✅ All GREB data loaded successfully` and passing `all(isfinite, Ts)`. The
  snippet now threads `fields` through, and the API refuses the mistake.
- A second instance of the same pattern in README §Loading Data.
- `tools/convert_greb_to_jld2.jl` defaulted its input directory to
  `Data/input`, which does not exist - the layout is flat `Data/` plus
  `Data/solar_forcing_scenarios/` - so the documented no-argument invocation
  always failed. Corrected, with an actionable error when the directory is
  missing.
- `docs/src/index.md` claimed Julia 1.9; `Project.toml` requires 1.10.
- `docs/src/tutorial.md` referenced a `claude/BENCHMARKS.md` that never existed.
18 correctness bugs found by direct comparison against the Fortran
reference (`greb.model.mscm.f90`); the ones with the widest-reaching impact:

- `circulation!` left a stale moisture-convergence term in every
  temperature sub-step of every default run.
- `set_hydrology_parameters!` wrote to disconnected module globals instead
  of the config struct - the `log_rain` switch had zero effect on any run.
- `hydro!`'s `log_eva` modes `1`/`2` silently duplicated mode `-1`'s formula
  instead of using their own coefficients.
- `deep_ocean!` cut off ocean-atmosphere heat exchange under sea ice and in
  high-latitude winters.
- 19 of 36 exported functions had docstrings silently detached from their
  definitions by a comment-placement pattern that breaks Julia's `@doc`
  binding - found while wiring up the Documenter.jl site.

The remaining 13 fixes, plus one investigated-and-reverted finding, are
detailed in `.claude/notes/audit-history.md`.

### Performance
- ~2.31× faster per simulated year (2.7s → 1.17s) from 3-way threading of
  `tendencies!` combined with 4 `@turbo`-rewritten hot-loop kernels.
- `Float32` throughout: a further ~1.6× on top of threading, with output
  validated against the previous `Float64` path to well under 0.01 K.
- Flux-correction files (`Tsurf`/`vapour`/`Tocean` corrections) merged into
  a single `flux_corrections.jld2` - ~35% faster to load, no size penalty.

### Changed
- The dataset shrank from 580 MB / 49 files to **439 MB / 39 files**: 11 files
  that no code reads were removed, mostly CMIP5 `.new` variants of fields the
  model reads in their non-`.new` form. The official MSCM Fortran GREB opens
  the non-`.new` names, so this changes no results; see
  `.claude/notes/data-distribution.md`.
- User-facing data documentation now describes obtaining the prepared `.jld2`
  bundle. `DATA_README.md` and the `.bin` converter are labelled as maintainer
  tooling, which is what they are - the raw inputs are collated from several
  upstream sources and are not redistributed. Distribution via DataDeps.jl is
  planned; see `.claude/notes/data-distribution.md`.
- Development notes moved from `claude/` to `.claude/notes/`, split by topic
  with an explicit status (Done / Investigated / Open / Planned) per file and
  an `INDEX.md`. This also removes the `claude/` vs `.claude/` name collision.
- `src/GREB.jl`, originally a single 2,245-line file, split into topical
  files (`constants`, `config`, `state`, `io`, `physics/`, `circulation`,
  `tendencies`, `output`, `postprocess`, `model`).
- Test suite split into `light`/`heavy` shards (`GREB_TEST_SHARD` env var),
  matching the CI matrix.

## [0.1.0] - 2026-08-06
Initial extraction from the interactive Pluto notebook into a standard Julia
package layout (`Project.toml`, `src/`, `test/`).
