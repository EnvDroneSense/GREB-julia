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

- [x] Add compat bounds for all four deps:

  ```toml
  [compat]
  LoopVectorization = "0.12"
  NCDatasets = "0.14"
  StaticArrays = "1"
  Statistics = "1"
  julia = "1.10"
  ```

  Bounds derived from the resolved versions in `Manifest.toml`
  (LV 0.12.174, NCDatasets 0.14.15, StaticArrays 1.9.17, Statistics 1.11.1).
- [x] Consider raising `julia = "1.10"` (current LTS) — **done**; precompiles
      and tests pass on Julia 1.12.6.
- [x] Decide: keep `LoopVectorization` or replace `@turbo` with
      `@inbounds @simd`? **Kept for now** — it precompiles and works on 1.12;
      benchmark the difference in Phase 5, revisit then.
- [x] Re-resolve and run tests — `Pkg.update()` + `Pkg.test()`: 31/31 pass.

**Done when:** every `[deps]` entry has a `[compat]` bound and tests pass.

---

## Phase 2 — Code hygiene: constants, parameters, and exports (~4–8 h)

The notebook extraction left module-level `begin ... end` blocks defining **bare
(non-`const`) globals**: `xdim`, `ydim`, `dlon`, `dlat`, `Δt`, `nstep_yr`,
`jday_mon`, physical constants (`σ`, `min_T_K`, …), grid lookup tables, and more.
Only the `lon_j*` index vectors got `const`. Untyped globals defeat type
inference in every function that reads them — a real performance and correctness
hazard for a SIMD-optimized model.

### Design: three kinds of "constants", two destinations

Not everything should become `const`, and nothing should become a `Dict` or a
pile of kwargs (hash lookups in hot loops, no typo checking, no documented
defaults, and threading dozens of kwargs through the kernel call chain is
error-prone). Sort every module-level scalar/table into one of:

1. **Natural constants** — `σ`, `grav`, `cp_*`, densities, `const_pi`, `p_emi`.
   Nobody tunes these → plain `const` module globals.
2. **Structural constants** — `xdim`, `ydim`, `Δt`, `Δt_crcl`, `nstep_yr`,
   calendar tables, neighbour-index tables. Array sizes and precomputed lookups
   depend on them; changing them means changing the code → plain `const`.
3. **Tunable process parameters** (~20: transport `κ`, exchange coefficients
   `ct_sens`/`ce`/`co_turb`/`c_effmix`, albedo & ice parameters `da_ice`/
   `a_no_ice`/`a_cloud`/ice thresholds, scale heights `z_air`/`z_vapor`,
   hydrology coefficients, `S0_var`, clamp limits) → an **immutable
   `Base.@kwdef struct ModelParams`** with the current values as documented
   defaults:

   ```julia
   Base.@kwdef struct ModelParams
       κ::Float64       = 8e5     # atmospheric diffusion coefficient [m²/s]
       ct_sens::Float64 = 22.5    # sensible heat coupling [W/K/m²]
       ce::Float64      = 2e-3    # latent heat transfer coefficient
       da_ice::Float64  = 0.25    # albedo increase for ice cover
       # ...
   end
   p = ModelParams(κ = 1e6)       # sensitivity run: tweak one, keep the rest
   ```

   Field access on a concrete immutable struct compiles to a constant — zero
   overhead, `@turbo`-compatible (keep the existing pattern of copying fields
   into locals before hot loops, as `hydro!` does with `cfg.c_q`).

**Derived parameters follow their bases.** `turb_coeff`, `ccx_diff`/`ccy_diff`,
`ccx_adv`/`ccy_adv`, `const_factor`, `cap_ocean`/`cap_land`/`cap_air`,
`inv_*_ice_range` are precomputed at module load from tunables — today, tuning
`κ` would silently leave `ccx_diff` stale. Compute all derived quantities in
the `ModelParams` inner constructor (or a `derive(p)` builder) so base and
derived values can never desynchronize. This closes a latent footgun, not just
a style issue.

**Plumbing:** since `cfg::PhysicsConfig` already reaches every kernel, nest the
parameters as a field — `params::ModelParams = ModelParams()` inside
`PhysicsConfig`; kernels read `cfg.params.κ`. No signature changes anywhere.
(A separate argument is conceptually cleaner — switches vs. numbers — but
touches every call site; split later only if it earns it.) End state for
sensitivity studies:
`greb_model!(...; cfg = PhysicsConfig(params = ModelParams(ce = 2.5e-3)))`
instead of editing source.

### Tasks (in order — the model stays runnable after each step)

- [ ] Sweep `src/GREB.jl` and mark every module-level assignment `const`
      (categories 1–3 alike, as a mechanical first step; `const` on an array
      fixes the *binding*, contents stay mutable — fine for read-only tables).
      Benchmark before/after with `benchmark/benchmarks.jl`, record the speedup.
- [ ] Introduce `ModelParams` (with derived quantities in the constructor) and
      nest it in `PhysicsConfig`. Migrate kernels one at a time — each migrated
      parameter's `const` global is deleted in the same commit; run tests +
      benchmarks per kernel.
- [ ] Trim the export list: stop exporting bare grid constants
      (`xdim`, `ydim`, `nstep_yr`) — users should write `GREB.xdim`. Keep
      exporting types and model functions. Export `ModelParams`.
- [ ] Optional: split the 2 250-line `src/GREB.jl` into
      `src/{constants,types,io,physics,circulation,model}.jl` with `include`s
      in `GREB.jl`. Pure file moves, no logic changes. Improves navigability;
      defer if it would disrupt in-flight work.

**Done when:** no non-`const` module-level globals remain; tunable parameters
live in `ModelParams` with derived values computed from bases; tests +
benchmarks pass; exports are types & functions only.

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
| 2 | `const` globals + `ModelParams` + exports | 4–8 h | Performance credibility |
| 3 | CI | 1 h | Everything after it |
| 4 | Data artifacts + real tests | 3–6 h | Usability by others |
| 5 | Documenter docs | 2–4 h | — |
| 6 | Registration | 1 h + review | `Pkg.add("GREB")` |

Phases 0–3 (~1 day) give a solid, CI-tested package installable from GitHub.
Phase 4 is the biggest usability win. Phases 5–6 are polish for a public release.
