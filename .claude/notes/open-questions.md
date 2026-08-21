# Open questions needing domain input

**Status:** Open — investigated and documented, deliberately not changed.

These need a climate-science or design decision rather than more code reading.
Nothing here is a bug with a known fix.

<!-- Split out of the former monolithic claude/IMPROVEMENTS.md on 2026-08-21.
     Section numbers in cross-references (§N.M) refer to that original file;
     see INDEX.md for where each section now lives. -->

---

## 4. Findings needing domain/design input — documented, not fixed

### 4.1 Seven config switches wired to nothing ✅ resolved
`log_vapor_drsp`/`log_ice_drsp`/`log_ice_dmc`/`solar_multiplier` had zero
Fortran basis — **removed** from `PhysicsConfig`. `log_tsurf_ext`/
`log_hwind_ext`/`log_omega_ext` are declared in Fortran but genuinely dead
there too — **kept**, with a comment, for structural fidelity.

### 4.2 `log_clim` doesn't select the ERA/NCEP dataset its label implies ✅ resolved
`cfg.log_clim` (hydrology regression coefficients) and `load_greb_jld2!`'s
`dataset` kwarg (ERA vs. NCEP file selection) are fully disconnected. Per
user decision, kept orthogonal (a caller can combine `log_clim=1` with
`dataset=:era` for a sensitivity experiment) and documented the relationship
with cross-referencing comments on both (`src/config.jl`, `src/io.jl`).

### 4.6 More `const`-ify candidates — resolved by §1.1
The original target (module-level globals never rebound) no longer exists;
nothing to do.

### 4.8 `forcing()`'s dispatch style is inconsistent with the codebase — deferred by decision
Its ~180-line `if/elseif` chain over `cfg.experiment` is the odd one out
next to the Dict-based dispatch used elsewhere. **User decision: leave
as-is** — zero known bugs after two Fortran-audit passes, and refactoring
risks breaking the "missing branch = valid no-op" semantics §0.6 protects,
for a style-only win.

### 4.9 Testing gaps — ✅ closed, see §3

### 4.10 CI has no lint/format or doc-build check — deferred by decision
**User decision: skip for now** — `JuliaFormatter` isn't present anywhere in
the repo; adding one means picking a style and reformatting the whole
codebase, bigger scope than a mechanical addition.

### 4.11 `output.jl`'s `global_mean` was an unweighted grid mean ✅ fixed
Fortran's `gmean()` (`greb.model.mscm.f90:1497-1513`) area-weights by
`cos(lat)`; Julia's plain `sum/count` didn't. Fixed by reusing the existing
`dxlat_grid` (already `cos(lat)`-proportional) as the weight. Diagnostic-only
(`println`, never stored in `MonthlyRecord`) — golden regression test
unaffected.

---
