# Code structure

**Status:** Done — the refactor described here has landed.

How the model went from a single 2245-line file with ~40 mutable globals to
topical files under `src/` with explicit `ClimateFields`/`ModelState` structs.

<!-- Split out of the former monolithic claude/IMPROVEMENTS.md on 2026-08-21.
     Section numbers in cross-references (§N.M) refer to that original file;
     see INDEX.md for where each section now lives. -->

---

## 1. Structure

### 1.1 Module-global mutable state → state structs ✅ done
Replaced by `ClimateFields`/`ModelState` (`src/state.jl`), threaded
explicitly through every physics/circulation/tendencies/output/model
function. `greb_model!` defaults to a fresh `ClimateFields()` per call,
making cross-run leaks (§0.5) structurally impossible for that path.
Validated to ~7e-12 absolute agreement against pre-refactor output; the same
1-year run went from 49s/17.1 GiB to 7.0s/100 MiB (unplanned bonus that also
fixed §2.2 and unblocked §2.4).

### 1.2 Module split into files ✅ done
`src/GREB.jl` is now a shell with `include()`s in dependency order:
`constants.jl` → `config.jl` → `state.jl` → `io.jl` →
`physics/{radiation,hydrology,ocean}.jl` → `circulation.jl` →
`tendencies.jl` → `output.jl` → `postprocess.jl` → `model.jl`.

### 1.3 Positional argument lists → structs ✅ done
`SurfaceState` (`src/state.jl`) wraps `Ts0`/`Ta0`/`To0`/`q0` by reference;
the rest of `output!`/`diagnostics!`'s former 9–13 positional args already
lived on the existing `tend` `NamedTuple` or `ws::CirculationWorkspace`, so
those are passed through instead of a second struct. Zero new allocations.
Bonus: fixed a real ordering footgun between `time_loop!`'s and
`diagnostics!`/`output!`'s old positional argument order (`q`/`To` swapped).

### 1.4 Run durations → `RunSpec` ✅ done
`greb_model!(run::RunSpec, cfg; ...)`, `RunSpec(; flux=0, ctrl=1, scnr=1)`
(`src/config.jl`), replacing three easily-swapped positional ints. All call
sites updated.

### 1.5 Reunite notebook-only helpers with the package ✅ done
`examples/run_greb.jl` provides a non-Pluto path.

### 1.6 Remove hidden/dead globals/fields ✅ done
Removed: `WZ_CACHE`, `circ_workspace`, two shadowed monthly-mean buffer
blocks (21 arrays), 4 dead accumulator entries, `co2_part_scn`,
`sw_solar_forcing_data`, `CirculationWorkspace`'s `dX_crcl`/`ws.eva`/
`ws.rain`. All confirmed dead by grep before removal.

---
