# GREB.jl — Potential Improvements (Structure & Performance)

Observations and concrete suggestions after extracting the Pluto notebook into a
package. Nothing here has been applied — it is a menu for future work. Line
references are to `src/GREB.jl` as of the extraction (2253 lines, single file).

The single theme underneath most of this: **the model runs on ~40 mutable
module-level globals** (grid arrays, climatology fields, work buffers, flux
corrections). That one design choice drives most of the structural *and*
performance issues below — it's why `time_flux`/`jdal2_dir` leaked as globals,
why the module eats ~750 MB at load, and why the code is hard to thread or test.

---

## 1. Structure

### 1.1 Replace module-global mutable state with state structs  ⭐ highest impact
Today the model's entire state lives as module globals: climatology fields
(`Tclim`, `uclim`, `qclim`, …), work buffers (`Tsmn`, `Tmm`, …), masks
(`co2_part`), and flux corrections (`TF_correct`, …). Functions read/write them
by name (`load_greb_jdal2!` fills them; `tendencies!`/`forcing` read them).

Problems this causes:
- **No reentrancy / concurrency** — you cannot run two experiments at once, or
  thread over ensembles; they'd clobber shared globals.
- **Hidden coupling** — exactly what produced the `time_flux`/`jdal2_dir` bugs.
- **Type instability** (see §2.1) and **eager 750 MB allocation** (see §2.2).
- **Hard to test** — you can't construct a small, isolated state.

Suggested shape:
```julia
struct ClimateFields          # the loaded climatology (immutable container, mutable arrays)
    Tclim::Array{Float64,3}; uclim::Array{Float64,3}; ...
end
mutable struct ModelState     # evolving state + reusable work buffers
    Ts::Matrix{Float64}; Ta::Matrix{Float64}; To::Matrix{Float64}; q::Matrix{Float64}
    ws::CirculationWorkspace; acc::MonthlyAccumulator; ...
end
greb_model(fields::ClimateFields, cfg::PhysicsConfig, run::RunSpec) -> Result
```
`CirculationWorkspace` and `MonthlyAccumulator` already do this well for their
slices — extend the same pattern to the rest. This is a large but mechanical
refactor and unlocks almost everything else.

### 1.2 Split the 2253-line module into files
`include()` topical files from `GREB.jl`:
```
src/GREB.jl            # module, usings, includes, exports
src/io.jl              # read_jdal2, load_*_jdal2!
src/config.jl          # PhysicsConfig, create_experiment_config, set_hydrology_parameters!
src/physics/radiation.jl   # SWradiation!, LWradiation!
src/physics/hydrology.jl   # hydro!
src/physics/ocean.jl       # seaice!, deep_ocean!
src/circulation.jl     # diffusion!, advection!, convergence!, circulation!
src/tendencies.jl      # tendencies!, forcing
src/output.jl          # diagnostics!, output!, time_loop!, accumulators
src/model.jl           # init_model!, qflux_correction!, greb_model!
src/postprocess.jl     # build_monthly_climatology, apply_scenario_anomalies, ...
```
No logic changes; purely navigability.

### 1.3 Collapse giant positional argument lists into structs
Several kernels thread ~10–20 positional arrays:
- `output!(it, irec, mon, Ts0, Ta0, To0, q0, albedo, ice, precip, evap, qcrcl, sw, lw, qlat, qsens, …)` (line 1899)
- `diagnostics!(it, year, CO2, Ts0, Ta0, To0, q0, albedo, sw, lw_surf, q_lat, q_sens, timestate)` (line 1844)
- `time_loop!(…)` (line 1932)

These are error-prone (positional mixups) and hard to read. Passing a
`ModelState`/`Diagnostics` struct (from §1.1) removes most of them.

### 1.4 Turn run durations into a `RunSpec`, not positional ints
`greb_model!(time_flux, time_ctrl, time_scnr, cfg)` — three bare ints whose
order is easy to swap. A `RunSpec(; flux=0, ctrl=1, scnr=1)` keyword struct is
self-documenting and pairs with §1.1.

### 1.5 Reunite the two notebook-only helpers with the package
`current_physics_config` (reads `@bind` widgets) and `setup_benchmark` were left
in the notebook. Once §1.4 exists, provide a small non-Pluto constructor
(`PhysicsConfig` from a preset + overrides) so the notebook and scripts share one
path instead of the widget-coupled function.

### 1.6 Remove remaining hidden globals
- `WZ_CACHE = Dict{Float64,Matrix{Float64}}()` (line ~304) is a global memoization
  cache — move into `ModelState` or a proper memoized function to keep it from
  leaking across runs.
- `sw_solar_forcing_state = Ref(1.0)` — same; belongs to run state.

---

## 2. Performance

### 2.1 `const`-ify the constants  ⭐ biggest easy win
Grid dimensions and physical constants are **non-`const` module globals**:
`xdim = 96`, `ydim = 48`, `nstep_yr` (lines 104–114); `σ`, `ρ_ocean`, `grav`,
`cp_air`, … (lines ~513–523). Every function reading them incurs a dynamic,
type-unstable global lookup, and the compiler cannot constant-fold loop bounds
like `for i in 1:xdim`.

Fix: mark the true constants `const` (or, better, make grid size a type
parameter / pass via a `Grid` struct). Verify the win with `@code_warntype` on a
hot kernel (`tendencies!`, `advection!`) — expect fewer `Any`/`Box` and faster
inner loops. Low effort, high payoff, and it composes with `@turbo`.

### 2.2 Stop allocating ~750 MB at module load
There are **28 top-level `zeros(Float64, xdim, ydim, nstep_yr)` allocations** that
run when the module loads — ~27 MB each ⇒ **~750 MB reserved even for
`using GREB` with no run** (this is why `using GREB` is heavy and precompile is
slow). Move these into the state structs from §1.1 so they're allocated only when
a model is actually set up. Also shrinks the test/precompile footprint.

### 2.3 Consider `Float32` for climatology fields
The JDAL2 files are `Float32` and are up-converted to `Float64` on load. Keeping
the large 3D fields `Float32` (or making element type a parameter) halves memory
(~375 MB) and roughly doubles SIMD lane throughput under `@turbo`. Keep
accumulators/tendencies in `Float64` if numerics require it — benchmark both.

### 2.4 Thread the spatial operators
Time-stepping is inherently sequential, but per-step spatial kernels
(`diffusion!`, `advection!`, `SWradiation!`, `LWradiation!`, `hydro!`) are
independent across grid columns. Once state is passed explicitly (§1.1, no shared
globals) these become safe to `@threads` / `@batch` (Polyester) over latitude.
Combined with existing `@turbo`, this is the main path to multi-core speedup that
the README's "parallelisation" future-plan asks for.

### 2.5 Cut allocations in hot paths
Audit with `@time`/`--track-allocation` or `BenchmarkTools`:
- `forcing(it, year, cfg, icmn_ctrl = zeros(xdim, ydim, 12); …)` (line 1658)
  **allocates a fresh `zeros` default every call** — hoist it into state.
- Watch for temporaries in broadcast/`@.` expressions in `tendencies!`/`hydro!`;
  prefer in-place `@turbo` loops writing to preallocated buffers.

### 2.6 Reduce time-to-first-run (TTFX)
Add a `PrecompileTools.@compile_workload` running a tiny 1-step model on synthetic
fields, so the heavy kernels compile at precompile time rather than on the user's
first call.

---

## 3. Testing, tooling & reproducibility

- **Reference/regression tests**: once the JDAL2 dataset is available, snapshot a
  short control run and assert against it (tolerance-based). Today's tests are
  smoke-only (they run on zeroed fields → `NaN`, so they can't catch numerical
  regressions).
- **Per-kernel unit tests** with small synthetic inputs (e.g. `diffusion!`
  conserves the field mean; `seaice!` clamps to physical bounds).
- **CI**: GitHub Actions matrix (Julia 1.9 → latest) running `Pkg.test()`.
- **`[compat]` bounds** for `NCDatasets`, `StaticArrays`, `LoopVectorization`
  (currently only `julia = "1.9"` is bounded) — needed before any registration.
- **Data acquisition**: a `scripts/get_data.jl` (or documented URL) so
  `examples/run_greb.jl` is runnable end-to-end; note the JDAL2 data is currently
  "available on request" only.
- **Docs**: a `Documenter.jl` site would suit the physics-heavy API; the notebook
  can become a rendered tutorial page.
- **`NCDatasets` dependency**: confirm it's actually used by the kept functions;
  if NetCDF I/O isn't reached, drop it to slim the dependency tree.

---

## 4. Suggested order (low-risk → high-value)

1. `const`-ify constants (§2.1) — tiny diff, immediate type-stability win.
2. Fix per-call allocations in `forcing`/hot loops (§2.5).
3. Split into files (§1.2) — mechanical, improves everything after.
4. Add `[compat]`, CI, per-kernel unit tests (§3).
5. **Introduce state structs** (§1.1, §1.3, §1.4) — the big one; enables §2.2,
   §2.4, and removes the global-coupling class of bugs for good.
6. Lazy allocation (§2.2), threading (§2.4), `Float32`/TTFX (§2.3, §2.6).

> ⚠️ Every performance change must be validated against a reference run
> (§3) — this is a numerical model, so "faster" only counts if the output is
> unchanged within tolerance.
