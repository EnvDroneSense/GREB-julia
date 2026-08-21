# Fortran switch/experiment validation

**Status:** Done — switches validated against the original Fortran GREB.

<!-- Split out of the former monolithic claude/IMPROVEMENTS.md on 2026-08-21.
     Section numbers in cross-references (§N.M) refer to that original file;
     see INDEX.md for where each section now lives. -->

---

## 7. Fortran switch/experiment validation

Cross-referenced every `PhysicsConfig` field and `cfg.experiment` branch
against two `greb.shell.mscm.f90` variants and four `run.greb.*.csh`
scripts. 7.1–7.3 were initially report-only by explicit user decision, then
actually implemented in a follow-up pass (2026-08-12); 7.4–7.6 remain
documentation-only findings.

### 7.1 `:rcp26`/`:rcp45`/`:rcp60` ✅ implemented
Fortran `log_exp=96/97/98` load their CO2 files directly; Julia's
`forcing()` used to error unconditionally for these three. Fixed: they now
share the existing `:ssp*`/`:historical_co2` `cfg.co2_scenario[year]` lookup
dispatch in `forcing()` (`tendencies.jl`), loaded at scenario start via
`_CO2_SCENARIO_SYMBOLS`/`_CO2_SCENARIO_KEY` in `model.jl` (`:rcp60` maps to
the jld2 key `"rcp6"`, `:rcp26`/`:rcp45` map 1:1) from the already-working
`ipcc_scenarios.jld2`. `create_experiment_config` gained a matching preset.
Test: the former "not yet implemented" placeholder loop in
`test/runtests.jl` was replaced with a real synthetic-`ipcc_scenarios.jld2`
run for each of the three symbols.

### 7.2 `:custom_co2` (EXP=100) ✅ implemented
Fortran's `log_exp=100` lets a user point at an arbitrary `year CO2`
text file. Fixed: `PhysicsConfig` gained a `custom_co2_path::String` field;
`src/io.jl` gained `load_custom_co2_scenario(path)`, a plain-text parser
matching Fortran's format (skips blank/`#`-comment lines); `model.jl` loads
it into `cfg.co2_scenario` at scenario start (erroring clearly if the path
is unset), reusing the same `forcing()` dispatch as §7.1/`:ssp*`.
`create_experiment_config(:custom_co2; co2_path=...)` sets the path.

### 7.3 `_dmc`/`_drsp` switch families — combined presets added ✅ implemented
Every individual switch was already a real, wired `PhysicsConfig` field;
`create_experiment_config` had no preset that named a whole family at once
the way the Fortran run scripts do. Fixed (pure discoverability/convenience
addition, no dispatch-logic changes — every switch already worked):
`create_experiment_config(:decon_mean_climate; ...)` mirrors
`run.greb.decon_mean_climate.csh` (`log_exp=1`) — the 6 `_dmc` switches plus
the 5 circulation switches (`log_ice`/`log_hdif`/`log_hadv`/`log_vdif`/
`log_vadv`), each keyword-overridable, defaulting to `true` exactly as the
Fortran script's namelist defaults. `create_experiment_config(:decon_2xco2;
...)` mirrors `run.greb.decon2xco2.csh` (`log_exp=10`) — the 5 `_drsp`
switches (`log_crcl_drsp` excluded, no Fortran counterpart in this family)
plus the same 5 circulation switches, with `co2_concentration` fixed at
680.0.

### 7.4 `:sst_plus1` has no Fortran `log_exp` counterpart — confirmed intentional ✅
`model.jl:273,374-378` implements it correctly as a Julia-only addition
(logic lives entirely in `model.jl`, `forcing()` needs no special case).
Confirmed deliberate by explicit user decision, not a mistranslated
experiment number — no action taken.

### 7.5 `:a1b_scenario`/`:a1b_enhanced` — hardcoded formula, not a data file
No `ipcc.scenario.a1b*.txt` exists anywhere; Julia hardcodes a
piecewise-linear 1950→2000→2050→2100 (310→370→520 ppm) approximation
(`tendencies.jl:80-90,118-128`) — reasonable given no source table was
supplied, but an approximation, not a validated match.

### 7.6 Confirmed correct / already-resolved — no action
`log_clim`'s Fortran-convention mismatch (already resolved, §4.2);
obliquity/eccentricity range (Fortran comment is stale, glob-based file
discovery already handles the real file set correctly);
`earth_sun_distance_pct`'s `dradius` translation (confirmed correct);
regional CO2 experiments' dynamic-vs-static mask asymmetry between the two
halves of the family (both paths produce correct masks, architectural note
only).

---
