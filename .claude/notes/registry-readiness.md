# Julia General Registry publication readiness

**Status:** Open — checklist partially complete.

<!-- Split out of the former monolithic claude/IMPROVEMENTS.md on 2026-08-21.
     Section numbers in cross-references (§N.M) refer to that original file;
     see INDEX.md for where each section now lives. -->

---

## 9. Julia General Registry publication readiness

### 9.1 Data distribution (Artifacts.jl) — deferred by decision

The model needs ~1.16 GB of climatology/scenario data (`Data/` 581 MB +
`greb_dataset_jld2/` 580 MB) to actually run, but neither directory is
git-tracked — both are gitignored, so the git tree at any tag (what `Pkg`'s
package tarball is built from) already excludes them. That's correct: a
General-registry package must stay small, and nothing in the registry
process requires or expects large data to live inside the package tree
itself.

What's missing is the other half — an automated way for a `Pkg.add("GREB")`
user to actually *get* the data. Today it's "ask the maintainer, or
regenerate from raw `.bin` files via `scripts/convert_greb_to_jld2.jl`"
(README's Input Data section). There is no `Artifacts.toml` and no
`artifact"..."` usage anywhere in the repo — every loader in `src/io.jl`
(`load_greb_jld2!`, `load_cc_anomaly_jld2!`, `load_solar_forcing_jld2`,
`load_co2_scenario_jld2`, etc.) takes a plain required `jld2_dir::String`
positional argument, joined with hardcoded subdirectory names
(`climatology/`, `static/`, `solar/`, `scenario/`).

**Deferred plan, not implemented in this pass:**
1. Build the `greb_dataset_jld2/` tree, tar it, and attach it as an asset on
   a GitHub Release (the repo already publishes releases — see `v0.1.0`,
   `v1.0.0`).
2. Generate `Artifacts.toml` at the repo root via
   `Pkg.Artifacts.create_artifact`/`bind_artifact!`, pointing at that release
   asset URL with its SHA256 + git-tree-sha1.
3. Change `src/io.jl`'s loaders to resolve `artifact"greb_dataset"` instead
   of requiring a manually-supplied `jld2_dir`, falling back to an explicit
   path override for anyone who wants to point at a local/custom dataset
   (e.g. during development, or a different reanalysis product).
4. This also closes a real coverage gap noted during the registry audit:
   `.github/workflows/ci.yml` never fetches `greb_dataset_jld2/`, so the
   "heavy" test shard's data-gated tests — including the golden-regression
   check — currently always hit the `@test_skip "greb_dataset_jld2/ not
   present"` branch in `test/runtests.jl` and silently never run in CI.
   Once artifact-based loading exists, CI can resolve the artifact like any
   other user and actually exercise that coverage.

### 9.2 Repo rename (`GREB-julia` → `GREB.jl`) — deferred by decision

Registry convention expects the GitHub repo backing a package named `GREB`
to be called `GREB.jl`; this repo is `GREB-julia`. Not an AutoMerge
blocker, but commonly expected by reviewers. Rename is deferred for now —
revisit before actually submitting the registration PR
(`docs/make.jl`'s `deploydocs(repo = "github.com/EnvDroneSense/GREB-julia.git", ...)`
and any other hardcoded repo URLs would need updating alongside the
GitHub-side rename).

**Resolved 2026-08-21:** the package name itself also needed to change —
`GREB` is 4 characters, below AutoMerge's 5-character minimum (not
identified as a blocker in the 9.3 audit below; found afterward). Renamed
the package to `GREBClimate` and the repo to `GREBClimate.jl`, updating
`docs/make.jl`'s `deploydocs` call, `Project.toml`, the module name, and
every `using GREB`/`GREB.`-qualified reference across `src/`, `test/`,
`examples/`, `benchmark/`, and `docs/`.

### 9.3 Registry audit summary

A full AutoMerge-requirements pass (2026-08-12) otherwise found the package
in good shape: `Project.toml` has compat bounds on every dependency
including `julia` itself; `LICENSE` (MIT) exists; `test/runtests.jl` runs
cleanly with zero external data required (data-gated tests skip gracefully,
per §9.1 above); and a live check against `JuliaRegistries/General`
confirmed the name `GREB` is not yet taken. Fixed as part of this pass:
version reset to `1.0.0` for a clean first registration (the package had
briefly been bumped to `1.1.0`), a new `CITATION.cff`, a
`Pkg.add("GREB")` forward-reference in the README, and a git-history
rewrite (`git filter-repo`) to strip ~26 MB of now-deleted `Data/input/*.bin`
files that were still permanently baked into old commits.
