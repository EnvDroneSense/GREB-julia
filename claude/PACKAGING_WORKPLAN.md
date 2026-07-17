# Workplan: From "works as a package" to a polished Julia package

**Status:** Draft for review.

**Context:** The notebook→package extraction (see [WORKPLAN.md](WORKPLAN.md)) is done:
`Project.toml` has a name/UUID/version, the model lives in `src/GREB.jl`
(`module GREB`), smoke tests run via `Pkg.test()`, and there are `examples/`,
`benchmark/`, and `notebooks/` directories. This plan covers what remains to make
GREB.jl a *proper* Julia package: installable by others, CI-tested, documented,
and (optionally) registerable in the General registry.

**Guiding principle:** each phase leaves the repo in a working, committable state.
Phases are ordered by impact; later phases are optional depending on ambition.

---

## Phase 0 — Housekeeping (~30 min)

Small cleanups so the repo history is tidy before real changes.

- [x] Commit the pending `.gitignore` change and the untracked `benchmark/`
      directory (`benchmarks.jl` + `benchmark/Project.toml`; its `Manifest.toml`
      is already gitignored).
- [x] Fix the README clone URL (still points at `EnvDroneSense/GREB-julia`) and
      align install instructions with the package layout
      (`Pkg.add(url="...")` / `Pkg.develop` instead of "activate the environment").
- [x] Decide what to do with `claude/` planning docs — **kept in `claude/`**;
      they document the port's history. Revisit before a public release.
- [x] Consider renaming `notebooks/PultoUI.jl` → `PlutoUI.jl` or deleting it
      (it is an unrelated PlutoUI demo notebook, not GREB code). **Deleted** —
      recoverable from git history if ever needed.

**Done when:** `git status` is clean and the README describes the actual repo.

---

## Phase 1 — Project.toml compat bounds (~30 min) ⚠ registration blocker

The General registry's AutoMerge requires `[compat]` entries for **every**
dependency. This is also just good practice: it protects users from breaking
upstream releases.

- [ ] Add compat bounds for all four deps, e.g.:

  ```toml
  [compat]
  LoopVectorization = "0.12"
  NCDatasets = "0.14"
  StaticArrays = "1"
  Statistics = "1"
  julia = "1.9"
  ```

  Check the currently-resolved versions in `Manifest.toml` and bound to those
  (semver-compatible caret bounds are the default).
- [ ] Consider raising `julia = "1.10"` (current LTS) — LoopVectorization support
      on future Julia versions is a known risk; verify it precompiles on the
      target version.
- [ ] Decide: keep `LoopVectorization` or replace `@turbo` with
      `@inbounds @simd`? LV is a heavy dependency with an uncertain maintenance
      future. **Recommendation:** keep for now, benchmark the difference in
      Phase 5, decide then.
- [ ] Re-resolve and run tests: `julia --project=. -e 'using Pkg; Pkg.update(); Pkg.test()'`.

**Done when:** every `[deps]` entry has a `[compat]` bound and tests pass.

---

## Phase 2 — Code hygiene: constants and exports (~2–4 h)

The notebook extraction left module-level `begin ... end` blocks defining **bare
(non-`const`) globals**: `xdim`, `ydim`, `dlon`, `dlat`, `Δt`, `nstep_yr`,
`jday_mon`, physical constants (`σ`, `min_T_K`, …), grid lookup tables, and more.
Only the `lon_j*` index vectors got `const`. Untyped globals defeat type
inference in every function that reads them — a real performance and correctness
hazard for a SIMD-optimized model.

- [ ] Sweep `src/GREB.jl` and mark every module-level assignment `const`.
      Scalars and arrays alike (`const` on an array fixes the *binding*, contents
      stay mutable — fine for read-only tables).
- [ ] Run the benchmark suite before/after (`benchmark/benchmarks.jl`) and record
      the speedup in the PR description.
- [ ] Trim the export list: stop exporting bare grid constants
      (`xdim`, `ydim`, `nstep_yr`) — users should write `GREB.xdim`. Keep
      exporting types and model functions.
- [ ] Optional: split the 2 250-line `src/GREB.jl` into
      `src/{constants,types,io,physics,circulation,model}.jl` with `include`s
      in `GREB.jl`. Pure file moves, no logic changes. Improves navigability;
      defer if it would disrupt in-flight work.

**Done when:** no non-`const` module-level globals remain; tests + benchmarks pass;
exports are types & functions only.

---

## Phase 3 — Continuous integration (~1 h)

No `.github/workflows/` exists yet.

- [ ] `CI.yml`: `julia-actions/setup-julia` + `julia-actions/julia-runtest`,
      matrix over `['1.10', '1']` (LTS + latest) on `ubuntu-latest`
      (add `macos-latest` if cheap). Trigger on push to `main` + PRs.
- [ ] `CompatHelper.yml` (weekly cron) — keeps compat bounds current.
- [ ] `TagBot.yml` — only useful once registered; add now, it's inert until then.
- [ ] Add coverage upload (Codecov) + badge, replacing/joining the current
      Julia/Pluto badges in the README.
- [ ] Tip: `PkgTemplates.jl` can generate all of these into an existing repo;
      or copy from any modern registered package.

**Done when:** a PR to `main` runs tests automatically and the README shows a
passing CI badge.

---

## Phase 4 — Input data as artifacts (~3–6 h, biggest UX win)

The model needs the external `greb_dataset_jdal2/` directory; today users must
obtain it manually and pass a path. Idiomatic Julia packages ship or fetch data
automatically.

- [ ] Host the JDAL2 dataset as a versioned tarball (GitHub release asset of
      this repo works well; check size limits — 2 GB/asset).
- [ ] Wire it up with **Artifacts.jl** (`Artifacts.toml`, lazy download) or
      **DataDeps.jl** (runtime download with user consent). Artifacts is the
      standard choice for fixed, versioned data.
- [ ] Change `load_greb_jdal2!()` to default to the artifact path, keeping the
      explicit-path method for custom data.
- [ ] Create a **tiny synthetic test fixture** (a few-KB JDAL2 file with the
      right header + small dims) checked into `test/fixtures/`, so CI can
      exercise `read_jdal2` and ideally a short model run on real code paths —
      today's tests are smoke-level only (types, constants, NaN-run).
- [ ] Add an integration test gated on the full artifact
      (`if haskey(ENV, "GREB_FULL_DATA")` or lazy-artifact download in a
      separate CI job): one control year, assert output shapes and plausible
      temperature ranges (regression guard against the Fortran reference).

**Done when:** `using GREB; load_greb_jdal2!()` works on a fresh machine with no
manual downloads, and CI tests real I/O paths.

---

## Phase 5 — Documentation (~2–4 h, optional)

Many docstrings already exist; `.gitignore` already anticipates `docs/build/`.

- [ ] Set up Documenter.jl: `docs/Project.toml`, `docs/make.jl`,
      `docs/src/index.md` (adapted from README) + an API-reference page
      (`@autodocs`).
- [ ] Add a `docs.yml` workflow deploying to GitHub Pages (`gh-pages` branch).
- [ ] Port the best explanatory markdown from the Pluto notebook into a
      "Model description" docs page (energy balance equations, experiment
      catalogue, scenario options).
- [ ] Add doctests to key docstrings where feasible.

**Done when:** docs build in CI and are browsable at
`https://<user>.github.io/GREB.jl`.

---

## Phase 6 — Registration decision (~1 h + review latency, optional)

Only if the goal is `Pkg.add("GREB")` for the world; a GitHub-installable
package (`Pkg.add(url=...)`) needs none of this.

- [ ] **Name check:** General-registry AutoMerge wants names ≥ 5 characters;
      "GREB" is 4. Options: (a) request a manual-merge exception (justifiable —
      GREB is the model's established name), or (b) rename to e.g. `GREBModel`.
      Decide before tagging v0.1.0.
- [ ] Verify AutoMerge guidelines: compat bounds (Phase 1) ✓, license ✓ (MIT),
      repo name ends in `.jl` ✓, sensible sequence of versions.
- [ ] Register via JuliaRegistrator (`@JuliaRegistrator register` comment on the
      release commit); TagBot then creates the GitHub release.
- [ ] Post-registration: subsequent releases are just version bump + re-register.

**Done when:** `Pkg.add("GREB")` (or the renamed package) works from a clean
Julia install.

---

## Suggested order & effort summary

| Phase | What | Effort | Blocking for |
|:--|:--|:--|:--|
| 0 | Housekeeping | 30 min | — |
| 1 | Compat bounds | 30 min | Registration |
| 2 | `const` globals + exports | 2–4 h | Performance credibility |
| 3 | CI | 1 h | Everything after it |
| 4 | Data artifacts + real tests | 3–6 h | Usability by others |
| 5 | Documenter docs | 2–4 h | — |
| 6 | Registration | 1 h + review | `Pkg.add("GREB")` |

Phases 0–3 (~1 day) give a solid, CI-tested package installable from GitHub.
Phase 4 is the biggest usability win. Phases 5–6 are polish for a public release.
