# DataDeps distribution

**Status:** Done — implemented 2026-08-21. **One prerequisite outstanding:** the release asset is not uploaded, so the URL 404s on a machine with no local dataset.

How the `.jld2` dataset is distributed and located. Sections 1-5 are the
context that informed the design, measured or quoted rather than assumed; the
sections after them record what was decided and built. Pairs with
[`data-distribution.md`](data-distribution.md), which holds the dataset audit.

---

## 1. What we are shipping

After removing the 11 unreferenced files: **39 files, 439 MB**.

| Group | Files | Size | Needed for |
|:------|------:|-----:|:-----------|
| static + solar + NCEP + common climatology + flux corrections | 15 | **180 MB** | any `:full_model` run with `dataset=:ncep` |
| ERA-Interim climatology | 4 | 52 MB | `dataset=:era` only |
| CMIP5 RCP8.5 anomalies | 5 | 65 MB | `:rcp85` / climate-change experiments |
| ENSO anomalies (elnino + lanina) | 10 | 129 MB | `:elnino` / `:lanina` only |
| solar scenarios (paleo, eccentricity, obliquity) | 3 | 15 MB | orbital / paleo experiments |
| CO2 scenario tables | 2 | <1 MB | scenario runs |

The headline number: **a default NCEP run needs 180 MB of the 439 MB.** That is
the strongest argument for splitting the bundle rather than shipping one blob —
it more than halves the download for the common case.

`climatology/` alone is 424 MB (97% of the payload), so any split has to divide
*that*, not the top-level directories.

## 2. Archive format — use `.tar.gz`

Measured on the real 439 MB tree:

| Format | Size | Notes |
|:-------|-----:|:------|
| `tar.gz` (gzip -6) | **353 MB** | 17s to build |
| `tar.xz` (xz -2 -T0) | 327 MB | only 7% better |

Prior measurement in [`data-organization.md`](data-organization.md) reached the
same conclusion from the other direction: compress **only** for distribution,
never the copy `load_greb_jld2!` reads, since JLD2-internal compression costs
8–20× on read time for a ~17% size win.

**Use `.tar.gz`, not `.tar.xz.`** Not for the size — for correctness. DataDeps'
`unpack` special-cases only `.tar`-secondary-extension, `.tgz` and `.tbz`; a
`.tar.xz` falls through to the generic 7z path, which strips the `.xz` layer and
leaves an un-extracted `.tar` behind. The 26 MB is not worth a two-pass
`post_fetch_method`.

## 3. DataDeps API (quoted from source)

Constructor, from `src/types.jl`:

```julia
DataDep(
    name::String,
    message::String,
    remotepath,
    hash=nothing;
    fetch_method=fetch,
    post_fetch_method=identity
)
```

- `name` — also the folder name the data lands in.
- `message` — shown when asking the user to consent to the download. This is
  where upstream attribution and licence text belongs.
- `hash` — SHA256 hex string is the normal form. Multiple files' checksums are
  XORed. `Any` suppresses the checksum warning (don't).
- `post_fetch_method` — `unpack` extracts and **deletes the archive** unless
  `keep_originals` is set.

Registration goes in the module's `__init__()`. `datadep"Name"` resolves lazily
at *evaluation* time, so nothing downloads at `using` time — only when the macro
is actually reached.

Environment variables that matter:

| Variable | Effect |
|:---------|:-------|
| `DATADEPS_ALWAYS_ACCEPT` | Skips the consent prompt. **Must** be set in CI or the process hangs waiting on stdin. |
| `DATADEPS_DISABLE_DOWNLOAD` | Any triggered download throws instead. Useful for CI jobs that should use a cache. |
| `DATADEPS_LOAD_PATH` | Extra colon-separated search paths, checked before the default. |

Default cache location: `~/.julia/scratchspaces/<DataDeps-UUID>/datadeps`.

## 4. Integration surface — smaller than expected

`load_greb_jld2!(jld2_dir)` already takes a plain path, and every downstream
loader (`load_flux_corrections_jld2!`, `load_cc_anomaly_jld2!`,
`load_enso_anomaly_jld2!`, `load_solar_forcing_jld2`, `load_co2_scenario_jld2`)
derives its path from that one argument. **No loader signature needs to change.**

What needed to change was only how the *default* path is resolved. Five places
hardcoded or half-resolved it:

| Site | Behaviour before this change |
|:-----|:-----------------------------|
| `examples/run_greb.jl` | `ARGS[1]` → `GREB_DATA` → `../greb_input_data` |
| `benchmark/run_benchmarks.jl` | `ARGS[2]` → `GREB_DATA` → `../greb_input_data` |
| `test/runtests.jl:9` | `DATA_DIR = ../greb_input_data`, hardcoded |
| `notebooks/GREB_julia.jl:159` | `../greb_input_data`, hardcoded |
| `README.md` / docs | tell the user to pass the path |

Shape adopted: one exported `greb_data_dir()` resolving in priority order —
explicit argument, then `GREB_DATA`, then a local `greb_input_data/` if present,
then the DataDep. All five sites call it. That keeps every current workflow
working (including a manually-unpacked bundle and an offline machine) and makes
DataDeps the fallback rather than a new requirement.

Interaction with precompilation: the `@compile_workload` block runs with
`jld2_dir=""` and `allow_uninitialized=true`, so it needs no data. Registration
in `__init__` is safe — but `datadep"..."` must never be evaluated at
precompile time, or building the package would trigger a 353 MB download.

`Project.toml` gains DataDeps in `[deps]` and `[compat]` (`"1"`, resolved to
v1.0.0). It is a light dependency, and making it a weak dependency/extension
would defeat the purpose — the whole point is that resolution works by default.

## 5. CI implications

`.github/workflows/ci.yml` shards tests via `GREB_TEST_SHARD` and has no
dataset, so every data-dependent test currently `@test_skip`s — including the
golden regression.

DataDeps makes it *possible* to run those in CI, which is a genuine gain in
coverage and a 353 MB download per job. Recommended: keep CI skipping by
default, and add one opt-in job with `DATADEPS_ALWAYS_ACCEPT=true` plus a cache
on the datadeps directory. Do not turn it on for the whole matrix.

## Decisions taken

- **Hosting: GitHub Releases.** Tag `data-v1`, asset `greb_input_data-v1.tar.gz`.
  No DOI, so citation stays informal — the DataDep `message` carries the
  attribution text instead.
- **Single bundle**, with the split kept as a future option (see below).

## What was implemented

- `src/data.jl` — `DATA_DEP_NAME`/`DATA_RELEASE_TAG`/`DATA_ARCHIVE_NAME`/
  `DATA_URL`/`DATA_SHA256`, the consent `message` carrying upstream
  attribution, `register_greb_datadep()`, and the exported `greb_data_dir()`.
- `GREBClimate.__init__` registers the DataDep. Registration is cheap and
  downloads nothing; `datadep"..."` is only reached inside `greb_data_dir`'s
  last branch. Verified that `Pkg.precompile` does not trigger a download.
- `tools/package_dataset.jl` — validates the tree against the converter's
  allowlist, builds a reproducible `tar.gz`, prints size + SHA256 and the
  release steps.
- Call sites wired: `examples/run_greb.jl` (download allowed),
  `benchmark/run_benchmarks.jl` and `test/runtests.jl`
  (`allow_download=false`, so measuring or testing can never pull 353 MB),
  `notebooks/GREB_julia.jl`.
- `DATADEPS_ALWAYS_ACCEPT=false` + `DATADEPS_DISABLE_DOWNLOAD=true` in both CI
  workflows, as a guard against a stray resolution hanging on stdin.
- Tests: resolution order and priority, error cases, and coherence of the
  URL/tag/asset/SHA constants (a typo there 404s only on a clean machine).

### The archive

```
greb_input_data-v1.tar.gz   370563615 bytes (353.4 MB)
sha256 370732f4166af8f9a5dfb2fd160c035b0881c79d1acfc53a0e59bc9bdb040315
```

Reproducible — `--sort=name --owner=0 --group=0 --numeric-owner` plus `gzip -n`.
Built twice from the same tree, byte-identical both times, and the Julia
packager independently reproduced the hash the shell build produced.

## Still outstanding

1. **Upload the asset.** The `data-v1` release does not exist yet, so
   `greb_data_dir()` currently 404s on a machine with no local dataset. Every
   other resolution path works. End-to-end verification needs the release, or a
   local file:// URL substituted temporarily.
2. ~~Redistribution terms~~ — **reviewed and cleared 2026-08-21** by the
   maintainer: the upstream data is Creative Commons and this is non-commercial
   research use. The DataDep `message` names all six sources and asks for
   acknowledgement, which stays the right thing to show at the consent prompt.

   One thing to revisit if the bundle is ever used commercially or
   redistributed onward: if any upstream component is CC **NC**, that
   restriction passes to downstream users of the bundle, and the `message`
   currently does not say so. Not a problem for the current use.

## Future option: splitting the bundle

Kept deliberately easy. A default NCEP run needs 180 MB of the 439 MB, so a
split would cut the common-case download by 59%. It would mean several
`DataDep` registrations and a rule mapping `dataset=`/`experiment` to the deps
each requires — plus a decent error when a group is absent. `greb_data_dir()`
returning a single directory is the interface that would have to change (or
grow a per-group variant), which is why it is worth doing only once there is
evidence users care about the download size.
