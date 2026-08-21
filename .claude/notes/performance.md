# Performance

**Status:** Done — optimizations landed; benchmark methodology and numbers recorded.

The largest note by far: every optimization pass, what was measured, and what
was rejected. Run the suite via the `benchmark` skill
(`.claude/skills/benchmark/SKILL.md`) rather than reading numbers off this page,
since they are machine-specific.

<!-- Split out of the former monolithic claude/IMPROVEMENTS.md on 2026-08-21.
     Section numbers in cross-references (§N.M) refer to that original file;
     see INDEX.md for where each section now lives. -->

---

## 2. Performance

### 2.1 Constants ✅ done
Grid dimensions and physical constants in `src/constants.jl` are `const`.

### 2.2 Stop allocating ~750 MB at module load ✅ done
Direct consequence of §1.1 — climatology arrays are `ClimateFields` fields,
allocated only when constructed, not at `using GREB` time.

### 2.3 `Float32` for climatology/working fields — ✅ done (2026-08-13)
Retyping just the static (JLD2-native-`Float32`) climatology is safe, but
`init_model!`'s `copy(Tclim[:,:,end])`-style inits would silently downgrade
the *working* simulation state (`Ts`/`Ta`/`To`/`q`) too without an explicit
`Float64.(...)` guardrail. Measured the actual arithmetic payoff of a full
switch (standalone clone of the circulation kernels, both eltypes):

| | `circulation!(Ta)` | `circulation!(q)` | total |
|---|---|---|---|
| `Float64` | 1301–1385 µs | 1360–1455 µs | 2661–2840 µs |
| `Float32` | 1210–1218 µs | 1337–1370 µs | 2554–2580 µs |

Only **~5–10%** faster, not 2× — the zonal terms gather through permuted
neighbor-index arrays rather than reading contiguous memory, and the 96×48
grid already fits in cache at either precision. `CirculationWorkspace` is
also a concrete non-parametric `Float64` struct, so a real switch needs a
`CirculationWorkspace{T}` refactor first. **Not worth it for this kernel
shape** given §2.11 already gets 2.0× more cheaply; canceled by user
request. Could matter more for a GPU/ensemble-batched use case (§2.11),
untested here. Full methodology: `AUDIT_LOG.md`.

**Follow-up, benchmarked (`CircularArrays.jl` for the gather pattern itself) — tried, measured, rejected.**
An outside suggestion proposed `CircularArrays.jl` as a direct fix for the
gather pattern identified above: wrap `T1`/`wz` so periodic-longitude
neighbor reads (`T1[j-1,k]`/`T1[j+1,k]`) become ordinary offset accesses
instead of indexing through `lon_jm1..jp3`. Built a standalone clone of
`diffusion!`'s zonal stencil (`src/circulation.jl:73-94`) two ways —
today's `Matrix{Float64}` + index-array approach vs. a `CircularVector`
wrapping each row (`CircularArrays.jl` circularizes *every* dimension of a
`CircularArray`, so a whole-2D-array wrap would incorrectly wrap the
non-periodic latitude axis too — a per-row `CircularVector` is the only
structurally correct option) — and measured both, `@turbo`-annotated
identically:

| | time/call | vs. baseline |
|---|---|---|
| Baseline (`Matrix` + index arrays, `@turbo`) | 16–30 µs | 1× |
| `CircularVector`-per-row (`@turbo`) | 340–426 µs | **0.03–0.09× (10–20× slower)** |

Output was correct (diffs `1.6e-16`, float-reassociation only, 0 bytes/call
either way), but `@turbo` **silently fails to compile** over the
`CircularVector`-backed loop (`LoopVectorization.check_args` rejects it,
falling back to a plain `@inbounds @fastmath` loop with a runtime warning)
— confirmed directly in `code_native`: the baseline's inner loop is real
packed-SIMD code (`vgatherqpd`×12, `vmovupd`/`vmovapd`×26, `...pd` = packed
double), while the circular version's inner loop is **100% scalar**
(`vmovsd`/`vmulsd`/`vsubsd`/`vaddsd`, `...sd` = scalar double, zero packed
instructions of any kind). The fix removes the gather instructions
entirely, exactly as claimed — but takes the SIMD vectorization down with
it, which costs far more than the gather ever did. **Rejected**: this is
worse than doing nothing on every axis (correctness aside), confirming the
`@turbo`-compatibility risk flagged before running the experiment. Script:
`benchmark/circular_arrays_experiment.jl`.

**Re-review (2026-08-13) — different verdict for non-gather kernels, still
not implemented.** The original analysis above only ever tested
`circulation!` (gather-bound via the periodic-longitude index arrays). Redid
the Float32-vs-Float64 comparison fresh on this machine/Julia version,
*and* extended it to `hydro!`/`LWradiation!`-style kernels, which are pure
elementwise (no neighbor gather at all — every read is at `[i,j]`, nothing
offset):

| Kernel shape | `Float64` | `Float32` | Speedup |
|---|---|---|---|
| Gather-bound (`circulation!`-style zonal diffusion) | 210.5 ns | 220.9 ns | **0.95× (no gain, slightly slower)** |
| Pure elementwise, no gather (`hydro!`-style) | 8733 ns | 2889 ns | **3.02×** |
| Transcendental-heavy elementwise (`LWradiation!`-style, `log()` calls) | 29800 ns | 11300 ns | **2.64×** |

So the original "not worth it" conclusion still holds for `circulation!`
itself — on this machine Float32 doesn't even keep its earlier claimed
5-10%, confirming it's not the lever for the dominant 65-75%-of-runtime
kernel. But `hydro!`/`LWradiation!`/`ocean!`'s elementwise kernels (no
gather to fight) show a genuine **2.6-3× speedup**, never tested before —
plausibly double SIMD lane width plus faster `exp`/`log`/`sqrt` at single
precision. These three together are a meaningful chunk of the ~25-35% of
per-timestep cost that isn't circulation (§2.10), so a 2.6-3× speedup there
could plausibly be a double-digit-percent full-model win — bigger than
anything in §2.15 and in the same league as §2.11's threading win.
Precision check: accumulating 730 `Float32` additions of O(250) values (a
proxy for `diagnostics!`/`output!`'s annual-mean accumulators) against the
same in `Float64` gives a relative error of ~2e-7 — negligible if the
*accumulator* arrays stay `Float64` and only per-timestep working buffers
go to `Float32`.
**Not implemented, and not a quick follow-up**: this would need a real
design decision, not a drop-in change — `ClimateFields`/`CirculationWorkspace`
are concrete non-parametric `Float64` structs (same blocker as before), a
mixed-precision split (elementwise physics in `Float32`, circulation and
accumulators in `Float64`) means crossing precision boundaries at several
call sites, and any change here needs validation against the Fortran
reference over full multi-decade runs, not just a standalone kernel diff.
Worth a dedicated pass if the user wants to pursue it — flagging it here as
reopened, not reclosed.

**Implemented (2026-08-13) — full conversion, plus the fix that made the
gather kernel benefit too.** Investigating the "not worth it for
`circulation!`" result further (per user question) turned up the actual
cause: the gather-kernel comparison above changed the data to `Float32`
but left the periodic-longitude index arrays (`lon_jm1`/`lon_jp1`/`lon_jm2`/
`lon_jp2`/`lon_jm3`/`lon_jp3`) as `Int64` — a width mismatch.
`LoopVectorization` needs matching-width gather addresses; `Float32` data
gathered through `Int64` indices nearly *doubled* the gather-instruction
count (15→28) instead of halving it. Narrowing the index arrays to `Int32`
alongside the data dropped the count back to 15 and made the kernel
genuinely faster than the `Float64`/`Int64` baseline (~245ns → ~175ns,
confirmed reproducible across 3 runs) — reversing the "not worth it"
verdict for `circulation!` itself, not just the elementwise kernels.

Given that, converted the whole model to `Float32` end-to-end, per explicit
user decision on scope (not a mixed-precision split):
- `ClimateFields`/`CirculationWorkspace`/`SurfaceState`/`MonthlyAccumulator`/
  `ModelState`/`MonthlyRecord` (`state.jl`): every field retyped, including
  the public output (`MonthlyRecord` is now `Matrix{Float32}`, not promoted
  back to `Float64`).
- `PhysicsConfig`'s numeric fields (`co2_concentration`, `c_q`/`c_rq`/
  `c_omega`/`c_omegastd`, `earth_sun_distance_pct`, `co2_scenario::Dict{Int,
  Float32}`) — also converted, not just cast locally, so a stray `Float64`
  scalar can't re-promote an elementwise kernel back to double precision.
- `lon_jm1`/`lon_jp1`/`lon_jm2`/`lon_jp2`/`lon_jm3`/`lon_jp3` → `Int32` (the
  fix above), plus every other `constants.jl` physical/derived constant.
- Every bare `Float64` literal inside a `@turbo` loop across
  `circulation.jl`, `physics/hydrology.jl`, `physics/radiation.jl`,
  `physics/ocean.jl`, `output.jl`, `model.jl`, `tendencies.jl` got an `f0`
  suffix — the highest-risk, most mechanical part of the change, since one
  missed literal silently re-promotes that loop to `Float64` and reproduces
  the original false-negative. One simplification found along the way:
  plain integer literals (`3`, `4`, `20`) do **not** need suffixing —
  `Float32 op Int` promotes to `Float32` in Julia, not `Float64`; only
  literals with a decimal point/exponent do.
- `load_greb_jld2!` no longer upconverts at all — JLD2 data is already
  `Float32`, so loading is now a same-type copy instead of a promotion (a
  small free win on top).

**Verification**: full 438-test suite passes (6 tests needed tolerance
loosening — `atol=1e-9`/`rtol=1e-10` style assertions to `~1e-5`–`1e-4`,
an expected consequence of the precision switch, not a workaround; two
`==` comparisons against `Float64` literal tuples/dicts needed `isapprox`
since a round-tripped `Float32` literal doesn't compare exactly equal to
the original `Float64` literal). Numeric sanity check (1-year control run,
`Float32` vs. a saved `Float64` reference from before the change): max
absolute difference 0.0036 K on `Ts`, 0.00086 K on `Ta`, 0.003 K on `To`,
2e-7 kg/kg on `q` — physically negligible, and the printed annual
diagnostic line matched the `Float64` baseline to 2 decimal places.
`@code_native` on the real (not standalone) `diffusion!` confirmed no
`cvtss2sd`/`cvtsd2ss`-style conversion instructions anywhere in the
compiled loop — no literal was missed.

**Real full-model result**: `benchmark/run_benchmarks.jl` (real NCEP data,
`-t 3`) — **0.68-0.72s/year, mean 0.70s**, vs. the ~1.1-1.2s/year `Float64`
baseline measured earlier this session. **~1.6× faster**, matching the
weighted estimate from the kernel-level numbers above and landing in the
same league as §2.11's threading win — the single biggest lever found this
session.

### 2.4 Thread the spatial operators — tried, measured, reverted
`Threads.@threads :static` on `diffusion!`/`advection!`'s outer row loop,
validated bit-identical output, but **3–6× slower** end-to-end at `-t auto`/14
threads (`tendencies!`: 4.2ms→23.7ms) — the grid (96×48, ~1µs/row) is too
small for the scheduling overhead to pay off at any thread count tried.
Reverted to plain `for k in 1:ydim` loops. **Takeaway**: don't retry plain
`Threads.@threads` at this grid/kernel shape without a reason to expect a
different result; `Polyester.@batch` or a much larger grid might fare
better.

**Follow-up, benchmarked (`Polyester.jl` retry) — real but narrow, fragile win; not implemented.**
Retried the reverted idea with `Polyester.@batch` instead of
`Threads.@threads`, on a standalone clone of `diffusion!`'s *full* outer-`k`
loop (mid-latitude **and** polar branches — the polar branches reuse a
single shared `ws.T1h`/`ws.dTxh` scratch buffer across an inner sub-step
loop in the real code, `src/state.jl:10-11`, which would race under
concurrent `k`; the clone gives each lane a private `(xdim, nlanes)`
scratch matrix indexed by `Threads.threadid()` instead). Verified
bit-identical output vs. sequential across every thread count tested,
including 5 repeated trials per count to rule out an intermittent race —
none found. Timed at every thread count from 1 to the machine's full 14
logical cores:

| `-t` | time/call | speedup |
|---|---|---|
| 1 | 61.6 µs | 1.00× |
| 2 | 65.1 µs | 0.92× |
| 3 | 54.3 µs | 1.17× |
| **4** | **35.4 µs** | **1.79×** |
| 5 | 110.5 µs | 0.55× |
| 6 | 158.7 µs | 0.37× |
| 7 | 166.6 µs | 0.35× |
| 8 | 137.6 µs | 0.43× |
| 14 (all logical cores) | 16,900–18,000 µs | **0.003–0.004×** (270–340× *slower*) |

Real result, not a wash: **`-t 4` gives a genuine 1.79× speedup** with
correct output — a materially better outcome than `Threads.@threads`'s flat
3–6× slowdown at every thread count tried in the original attempt. But the
curve is sharply non-monotonic and **not safe to run at `-t auto`** — every
count above 4 degrades, and running at the machine's full logical-core
count (14) doesn't just lose the speedup, it collapses catastrophically
(reproduced twice, both ~300× slower than sequential), consistent with
Polyester's static per-lane scheduling overhead dominating once lane count
exceeds what a 48-row loop can actually use, compounded by likely
oversubscription past physical-core count on a hyperthreaded machine. **Not
implemented**: a real win exists, but only if the thread count is
explicitly pinned (e.g. `-t 4`, not `-t auto`) rather than left to the
caller's environment — that's a behavioral/deployment decision, not a pure
code change, and is left for a follow-up once a specific thread-count
policy is decided. Script: `benchmark/polyester_experiment.jl`.

### 2.5 Cut per-call allocations in `forcing` ✅ done
`icmn_ctrl` is a required argument computed once, not fresh zeros every call.

### 2.6 Reduce time-to-first-run (TTFX) ✅ done
`PrecompileTools.@compile_workload` compiles the heavy kernels at package
build time instead of on every user's first call.

### 2.7 Cut `SWradiation!`'s per-call allocation ✅ done
`@views`-wrapped the loop body — was materializing a fresh vector per
iteration (40704 bytes/48 allocs), the *entire* allocation footprint of
`tendencies!`. Every kernel now benchmarks at 0 bytes.

### 2.8 Hot-loop field/config localization + type stability — audited, clean
`fields.xyz`/`cfg.xyz` are already hoisted to locals before every hot loop;
no struct field is untyped/abstract. Nothing to change.

### 2.9 Precomputation/caching — audited, clean
Every candidate for hoisting is already hoisted (`regional_co2_*` masks
where possible — the other 2 genuinely can't be, they depend on
`icmn_ctrl`) or already a precomputed lookup table/`const` array. Nothing to
do.

**Flagged, not fixed, while auditing**: `output.jl`'s `global_mean` was an
unweighted grid mean — resolved, see §4.11.

### 2.10 Where per-timestep time goes vs. Fortran — profiled, documented
A 1-simulated-year `:full_model` control run takes ~2.7s wall-clock
(post-precompile). Profiling (`Profile.@profile` + manual `@elapsed`) over
5000 warm `time_loop!` calls attributes it as:

| Stage | Time | Share |
|---|---|---|
| `circulation!` (`Ta`, `z_air`, 24 substeps) | 1370 µs | ~34% |
| `circulation!` (`q`, `z_vapor`, 24 substeps) | 1498 µs | ~37% |
| everything else (radiation, hydrology, ocean, output) | ~1100–1800 µs | ~26–29% |

**~65–75% of every timestep is `circulation!`**, run twice per main step.
Ruled out as causes: allocations (96 bytes/step total) and the polar
sub-stepping (only the 2 outermost rows need >1 sub-step). The remaining gap
to native Fortran is likely structural — 24×2 separate `@turbo`-swept
function calls per step rarely fuse as tightly as one hand-written Fortran
subroutine.

**The Fortran-comparison assumptions were checked directly, not assumed:**
the reference build (`gfortran -Ofast -ffast-math -funroll-loops -fopenmp`,
found in `run.greb.*.csh` next to the Fortran source) isn't apples-to-apples
with a single-threaded Julia run. `greb.model.mscm.f90`'s `tendencies`
subroutine (lines 512–535) wraps `SW`/`LW`+sensible/`hydro`/`circulation`(Ta)/
`circulation`(q)/`deep_ocean` each in an `!$omp section` — the **same six
stages** §2.11 found and benchmarked independently — so the folklore "<1s/yr"
figure already assumes multi-threaded execution; the same folder's own
`readme.md` states unoptimized/unparallelized Fortran takes **~30s/simulated
year**, slower than this Julia port's current ~2.7s/year single-threaded.
Checked the other flags too: `julia -O3` vs. default `-O2` showed no
measurable difference (2.13–2.22s/year, within noise — `@turbo`'s codegen is
largely independent of the global `-O` level); `-march=native` is already
Julia's default (`Base.JLOptions().cpu_target == "native"`, confirmed on
this machine). `-ffast-math` has no Julia equivalent applied globally, but
see §2.12 for its per-kernel impact where it'd actually matter.

### 2.12 Missing `@turbo` spots — ✅ done (4 of 4 candidates)
Every hot kernel in `circulation.jl`/`hydrology.jl`/`ocean.jl` uses `@turbo`.
The 4 candidates below (all rewritten as explicit `@turbo` loops) were found
by checking whether `LWradiation!`'s fix generalizes, prioritizing by call
frequency (per-substep code compounds: `circulation!` runs 24 substeps × 2
fields = 48×/main-timestep). Verified against the full test suite (incl.
golden regression); real 1-year `greb_model!` run: 2.7s → 2.212s.

| Spot | Frequency | Baseline | `@turbo` | Speedup | ~Saved/year |
|---|---|---|---|---|---|
| `diffusion!` mid-lat branch (`circulation.jl:63-66`) | ~46 rows × 48/step | 9.2 µs/call | 2.3 µs/call | 75% | 0.22s |
| `circulation!` substep accumulate (`circulation.jl:326`) | 48/step | 208.8 µs/step | 45.6 µs/step | 78% | 0.12s |
| `LWradiation!` (`radiation.jl:82`, 3 `log()`/cell) | 1/step | 145.7 µs | 36.3 µs | 75% | 0.08s |
| `time_loop!` state update (`output.jl:137-169`) | 1/step | 63.3 µs | 12.3 µs | 80% | 0.04s |

All four verified numerically identical to baseline (diffs `1e-14`–`1e-18`,
pure float reassociation) — the first correctness pass hit a spurious `NaN`
traced to the *test harness* forgetting to call `init_model!`/`seaice!` first
(`cap_surf` still all-zero → division by zero), not a real bug. **Total
~0.46s/simulated year** across these four — composes with (doesn't
double-count) §2.11's ~1.0s/year, since §2.11 parallelizes *across* the two
`circulation!` calls while these speed up code *inside* each call.

**One negative result, useful as a boundary case**: `circulation!`'s final
subtraction (`circulation.jl:330`, a single elementwise op, only 2 calls/step
— not per-substep) was **17.6% slower** under `@turbo` — the macro's own
overhead isn't worth paying below some size/frequency threshold. Confirms
this isn't "add `@turbo` everywhere," it's "add it where the call is large
or frequent enough" — `LWradiation!`/mid-lat-diffusion/substep-accumulate all
clear that bar; this one line doesn't. `@fastmath` (tested only on
`LWradiation!`, 45% faster there) is the softer middle option when a spot is
too awkward to rewrite as an explicit `@turbo` loop.

**Follow-up sweep (2026-08-12), 4 more confirmed wins + 1 more negative
result** — a broader pass (2 parallel Explore agents across all of `src/`)
resolved the two candidates above and found 3 more of the same shape:

| Spot | Baseline | `@turbo` | Speedup |
|---|---|---|---|
| `diagnostics!` accumulate (`output.jl:11-21`, 11 ops, 1×/step) | 22.65 µs | 9.75 µs | 57% |
| `state.jl accumulate!` (13 ops, 1×/step) | 24.2 µs | 11.85 µs | 51% |
| `SWradiation!`'s `a_atmos`/albedo-combine/`sw`-flux-loop (1×/step) | 5.5 µs | 3.4 µs | 38% |
| `qflux_correction!`'s per-timestep update (12 ops, spin-up only) | 57.0 µs | 18.0 µs | 68% |

All four **implemented** (`output.jl`, `state.jl`, `radiation.jl`,
`model.jl`), verified numerically identical (diffs `1e-16`ish) and against
the full test suite incl. golden regression. Combined contribution to a
real `:full_model` `flux=0` run is small — the first 3 apply there (~20ms/
simulated year estimated from the per-call numbers × 730 steps); the 4th
only runs during the flux-correction spin-up (`run.flux > 0`, off by
default) so it contributes nothing to that particular run. This ~20ms/year
is genuine but small enough to be swamped by this machine's own run-to-run
wall-clock noise (~300ms swings observed across repeated 1-year timings) —
unlike §2.11/§2.12's fixes, which were large enough to show up cleanly
end-to-end, this batch's win is only visible in the isolated per-call
microbenchmarks above, not in an aggregate before/after timing. Documented
honestly rather than reporting a noisy or cherry-picked aggregate number.

**One more confirmed negative result**: `seaice!`'s trailing glacier-override
line (`ocean.jl:46`, single op) was flat-to-slower under `@turbo` (−2.5% to
−47% across repeated runs, i.e. noise-dominated but never a real win) —
**not implemented**, a second instance of the same "too small/infrequent"
pattern as `circulation.jl:330`. `tendencies.jl:39`'s `Q_sens = ct_sens *
(Ta - Ts)` (single op, 1×/step) was also consistently slower under `@turbo`
(−13% to −16%) — **not implemented**, third instance of the pattern. Both
below the profitability threshold: a single elementwise op over 96×48
doesn't carry enough work to amortize `@turbo`'s own dispatch/setup cost.

### 2.13 `forcing()`'s regional-CO₂ mask recomputed every scenario timestep — ✅ done
`:regional_co2_ocean`/`:regional_co2_land_ice` (`tendencies.jl`) recomputed
their `co2_part` mask from scratch on every scenario timestep, even though
it depends only on `z_topo` and `icmn_ctrl` — both fixed for the whole
scenario run. Not a `@turbo` candidate (branchy, not a clean elementwise
op) — fixed algorithmically by guarding the mask computation with `it ==
1`; `fields.co2_part` correctly holds the static mask for every later
timestep since nothing else touches it for these two experiments. Measured
in isolation (`forcing()` called directly, `:regional_co2_ocean`): 7.65 µs
baseline → 8.35 µs at `it==1` (pays the real cost once) → **0.1 µs at
`it>1`** (98.7% faster) — ~6ms/simulated year saved for these two
experiments (small in absolute terms, since the mask itself is cheap; the
point is skipping it 729/730 times per year is free). Verified via the full
test suite's `"greb_model! reaches every :experiment symbol..."` test,
which exercises both experiments directly.

### 2.15 Four small-scale fixes (2026-08-13) — ✅ done
All verified against the full 438-test suite and a full-model benchmark run
(`benchmark/run_benchmarks.jl`, real NCEP data, `-t 3`) before/after —
output unchanged, no regression. Each fix was also benchmarked **in
isolation** against its pre-fix version (old code pulled from git `HEAD`,
run side-by-side with the fixed version in the same process, real NCEP
data where applicable):

| Fix | Old | New | Speedup | Est. saving/year |
|---|---|---|---|---|
| `output!` copy elimination (per field) | 23.6 µs | 3.2 µs | **7.4×** | ~3.2 ms (×13 fields, ×12 months) |
| `LWradiation!` fused `LW_up` write | 36.6 µs | 34.3 µs | 1.07× | ~1.7 ms (×730 calls) |
| `hydro!` loop fusion | 16.8 µs | 11.5 µs | **1.46×** | ~3.9 ms (×730 calls) |
| polar sub-step constants (`diffusion!`) | 36.9 ns | 2.2 ns | **16.8×** | ~22 ms (×700,800 calls) |
| polar sub-step constants (`advection!`) | 33.3 ns | 2.4 ns | **17.7×** | ~22 ms (×700,800 calls) |

Combined estimated saving: **~53 ms/simulated year**, out of a ~1.1-1.2
s/year baseline (~4-5%) — real, but still inside the full-model benchmark's
own run-to-run noise band (0.94-1.22s/year measured on this machine), so it
doesn't show up as a clean before/after delta at the whole-run level; the
per-fix isolation above is what actually demonstrates each one's effect.
Detail per fix:
- **`output!`'s redundant double allocation** (`output.jl`): each of the 13
  monthly-record fields was built as `copy(acc.Tmm ./ ndm)` — `./` already
  returns a fresh, non-aliased array, so `copy(...)` allocated a second
  array and immediately discarded the first. Removed the redundant `copy`.
  Biggest *relative* win of the four (7.4×) but only runs at month
  boundaries (12×/year), so modest in absolute terms.
- **Extra full-grid pass in `LWradiation!`** (`radiation.jl`): `LW_up .=
  LW_down` copied the whole grid in a second pass right after the `@turbo`
  loop that computed `LW_down`. Merged the write into the same loop
  (`LW_up[i,j] = LW_down_val` alongside `LW_down[i,j] = LW_down_val`). Only
  a 1.07× win in isolation — copying a cache-resident 96×48 array is
  already fast (basically a 36 KB memcpy), so removing the second pass
  saves little; the value here is code cleanliness, not raw speed.
- **`hydro!`'s 6 separate `@turbo` loops over the same grid** (default
  config runs 5 of them; `log_rain==1` adds the 6th): saturation humidity,
  relative humidity, evaporation (one of 4 `log_eva` branches), precipitation,
  the optional rain-limit clamp, and water-vapor tendencies were each a
  separate full-grid pass, despite being purely elementwise (no
  cross-grid-point coupling) — nothing here needs neighbor values the way
  `diffusion!`/`advection!` do. Fused each `log_eva` branch into a single
  `@turbo` loop that walks the whole per-point pipeline (including the
  rain-limit clamp, now an inline `ifelse` gated on `cfg.log_rain==1` rather
  than a separate conditional pass) in one grid traversal. Output matches
  the pre-fusion version to floating-point roundoff (~1e-13, from
  reassociated `exp`/`sqrt` calls). `hydro!` runs once per main timestep
  (not per circulation sub-step), so this is ~730×/year.
- **Polar sub-step constants recomputed every call** (`circulation.jl`'s
  `diffusion!`/`advection!`): `dd`/`dtdff2`/`time2`/`ccx2` in each polar-row
  branch depend only on `k` (via `dxlat_grid[k]`) and fixed module
  constants — never on the per-call `T1` data — yet were recomputed from
  scratch on every call. There are **20 of 48 rows** flagged `IS_POLAR`
  (not a handful, as the raw formula count might suggest) — 20 rows × 2
  fields (Ta, q) × 24 sub-steps × 730 steps/year = **700,800 calls/year**
  for diffusion alone, and the same count again for advection (different
  formula, same waste). Precomputed as `POLAR_DIFF_TIME2`/`POLAR_DIFF_CCX2`/
  `POLAR_ADV_TIME2`/`POLAR_ADV_CCX2` (`constants.jl`, same pattern as the
  existing `ccx_diff`/`ccx_adv`/`dxlat_grid` lookup tables), indexed by `k`.
  Individually this was only a handful of scalar ops (~35 ns), but at
  700,800× the call count this turned out to be the **largest absolute
  saving of the four** (~44 ms/year combined) — the "small-scale fix" label
  undersold this one; the call-count multiplier is what makes it matter.

### 2.14 Other sweep findings — investigated, not implemented
A wider pass (2026-08-12) also looked beyond `@turbo` candidates, at
architecture/I-O patterns:
- **CMIP5/ENSO anomaly reload** (`model.jl`'s `greb_model!`,
  `load_cc_anomaly_jld2!`/`load_enso_anomaly_jld2!`): always re-reads from
  disk even when the same `fields` is reused across multiple `greb_model!`
  calls with the same experiment — real for a repeated-run/ensemble use
  case, but a safe fix needs a cache-key/invalidation design (what if
  `jld2_dir`/`cfg` differs between calls?) that's a genuine design decision,
  not a quick win — left for whoever actually hits this in an ensemble
  workflow.
- **`postprocess.jl`'s runtime-`Symbol` `getfield` loop**
  (`build_monthly_climatology`/`apply_scenario_anomalies`): a real dynamic
  dispatch per field per record, but dwarfed by the ~4608-element array op
  it wraps at current `MonthlyRecord` counts (12×`time_ctrl`/`time_scnr`) —
  not worth an `@generated`/`Val`-based rewrite at this scale.
- **`CirculationWorkspace()`'s ~40 separate allocations**: only happens 3×
  per `greb_model!` call (`ws`/`ws_a`/`ws_q`), not per-timestep — negligible
  (~1.44 MiB, low-µs) next to a multi-second run; consolidating into fewer
  larger arrays would save allocation count but not meaningfully change
  wall-clock time for the current usage pattern, and would hurt readability
  for no real gain absent a hypothetical high-frequency-construction use
  case that doesn't currently exist in this codebase.
- **Rust FFI for hot kernels** (2026-08-13): benchmarked a Rust `cdylib`
  (`opt-level=3`, `lto`, both default and `target-cpu=native`) against the
  existing `@turbo` mid-latitude zonal-diffusion kernel
  (`circulation.jl:72-94`), calling it via `ccall` with real NCEP climatology
  data and byte-identical output. Result: Rust was **slower**, not
  faster — 374-413 ns/call vs. 233 ns/call for the current `@turbo` version.
  The `ccall` round-trip itself is negligible (~2 ns, measured with a no-op
  control), so FFI overhead isn't the cause: this kernel is a gather-heavy
  stencil (5 pairs of non-contiguous periodic-neighbour lookups per
  element), and `LoopVectorization.@turbo` already emits SIMD-gather code
  tuned for exactly that pattern, which plain `rustc -O3` autovectorization
  doesn't match without hand-written SIMD intrinsics. At this grid size
  (96×48) and kernel granularity, adding a second-language build (Rust
  toolchain in CI, `ccall` bindings, binary-artifact distribution for a
  registered package) would add real complexity for a net slowdown — not
  pursued further.
- **Loop iteration order vs. column-major memory layout** (2026-08-13):
  `physics/ocean.jl:20,40,81` and `physics/radiation.jl:20,39` write
  `@turbo for i in 1:xdim, j in 1:ydim` (`j`, not the memory-contiguous
  `i`, is the fastest-varying/innermost index) while other kernels nest the
  other way (`i` innermost). Benchmarked both orderings standalone on two
  kernels: `seaice!`'s first loop body (cheap, few branches) and
  `deep_ocean!`'s loop body (heavier — divisions, more branches). For
  **both**, with `@turbo` present, source order made no reliable
  difference — medians swapped which order "won" between repeated runs
  (e.g. `deep_ocean!`-style: 9193 vs. 8100 ns one run, 4243 vs. 4329 ns the
  next), and `@code_native` showed packed SIMD loads (`vmovup*`), **no
  gather instructions**, for either order in both kernels.
  `LoopVectorization.@turbo` picks the true unit-stride axis for
  vectorization itself, independent of source loop nesting order — so this
  mismatch costs nothing today in the actual `@turbo`'d code.
  - However, the same `deep_ocean!`-style body **without** `@turbo` (plain
    `@inbounds` loop) is genuinely order-sensitive: correct
    (memory-contiguous) order measured 8400-9800 ns, mismatched order
    measured 18300-21200 ns — a reproducible **~2.2× slowdown** — and
    `@code_native` confirmed why: the mismatched plain loop compiles to
    `vgather`/scalar loads instead of packed `vmovup*`. So the general
    principle ("wrong order costs nothing") only holds *because* `@turbo`
    is present; it would be a real bug in any loop that lacks it.
  - Checked every double-nested spatial loop in `src/` for this: every
    site using the "wrong" `i`-outer/`j`-inner order already has `@turbo`
    (the two files above); every plain/`@inbounds`-only double loop
    (`model.jl:119` in one-time `init_model!` setup, `tendencies.jl:206-233`
    in rare `it==1` static-mask branches for `:regional_co2_ocean`/
    `:regional_co2_land_ice`) already uses the correct `j`-outer/`i`-inner
    order. No live instance of "no `@turbo` + wrong order" exists in the
    codebase today. Not changed; the `i`/`j` nesting inconsistency in
    `ocean.jl`/`radiation.jl` remains a readability nit only — but the
    ~2.2× plain-loop finding is worth keeping in mind for any *future*
    spatial loop that can't use `@turbo` (e.g. calls a function
    `LoopVectorization` doesn't support): get the nesting order right for
    those, since there's no compiler safety net without `@turbo`.

### 2.11 Parallelize `tendencies!`'s independent stages — ✅ done
All six stages (`SW`/`LW`/`hydro`/`deep_ocean`/`circulation!(Ta)`/
`circulation!(q)`) read only pre-timestep state and write disjoint buffers,
except the two `circulation!` calls which need separate workspace instances.

| Grouping | Time (2000-call avg) | Speedup |
|---|---|---|
| Sequential (all 6 stages) | 2747.5–2887.5 µs | 1× |
| 2-way: `circulation!(Ta)` \| `circulation!(q)` | 1540.0 µs | 1.88× |
| 6-way: every stage individually spawned | 1801.5 µs | 1.56× (worse) |
| **3-way**: `circulation!(Ta)` \| `circulation!(q)` \| rest bundled | **1376.0 µs** | **2.0×** |

Bare `Threads.@spawn`+`wait` overhead (~19.4 µs) exceeds several stages'
own cost, so spawning all 6 individually is a net loss. The 3-way split
captures the full ceiling at 2.0× (~1.0s/simulated year), verified
bit-identical vs. sequential.

**Implemented** in `tendencies!` (`src/tendencies.jl`): `ws_a`/`ws_q`
keyword args (default `=ws`, i.e. sequential — strictly opt-in), gated on
`Threads.nthreads() > 1 && ws_a !== ws_q`. `time_loop!`/`qflux_correction!`
forward the same kwargs; `greb_model!` is the only call site that allocates
distinct `ws_a`/`ws_q` and thus actually parallelizes. Verified: full suite
(319/319, incl. golden regression) passes unchanged at `-t 1`/`-t 2`/`-t 3`.

**Needs `-t 3`, not `-t 2`, to realize the full ceiling** — this is a true
3-lane split (`circulation!(Ta)` \| `circulation!(q)` \| rest), and with
only 2 OS threads the two spawned circulation tasks queue onto the single
non-main worker thread and run *sequentially* there, overlapping with "rest"
on the main thread but not with each other. Real 1-year `greb_model!`
timing (§2.12's turbo fixes already included in every number; `deepcopy`
of the large climatology arrays kept strictly outside every timed region —
an earlier pass of this benchmark hid this exact 2-vs-3-thread distinction
by leaving `deepcopy` inside the timed block, diluting the real
per-thread-count differences with a large constant):

| Threads | Time | vs. `-t 1` |
|---|---|---|
| `-t 1` | 1.912s | 1× |
| `-t 2` | 1.363s | 1.40× |
| **`-t 3`** | **1.170s** | **1.63×** |
| `-t 4` | 1.294s | no further gain (within noise) |

`-t 3` is the minimum thread count that gives the full benefit — `-t 4`+
adds nothing further since only 3 concurrent pieces of work exist per
timestep. **~2.31× total** vs. the original ~2.7s/year baseline at `-t 3`,
consistent with (slightly under) the ~2.8× the isolated integration test
below measured, the remaining gap being `output!`/`diagnostics!`/
bookkeeping outside `tendencies!` that this change doesn't touch (Amdahl's
law, given `tendencies!`'s ~70–75% share, §2.10).

**Re-review (2026-08-13, post-`Float32`) — `-t 2` is now the reliable
choice; `-t 3`'s advantage is gone.** The `-t 3` recommendation above was
measured when circulation was ~65-75% of per-timestep cost and "the rest"
(SW/LW radiation, `hydro!`, `deep_ocean!`) was a real ~25-35% — a genuine
third lane worth its own thread. §2.3's `Float32` conversion changed that
balance: the new `stages` mode in `benchmark/run_benchmarks.jl` shows
circulation(Ta)+circulation(q) at **~98%** of the pipeline and "the rest"
at **~2%** (`Float32` gave the small elementwise kernels a bigger relative
speedup — 2.6-3× — than circulation's ~1.4×, so its tiny absolute cost
shrank further while circulation's share of the *total* actually grew).

With "the rest" now nearly free, `-t 2` already captures almost all the
real parallelism: the main task finishes its ~30µs of synchronous work,
blocks on `wait()`, and Julia's work-stealing scheduler picks up whichever
of the two spawned circulation tasks hasn't started yet — so the two big
tasks still run genuinely concurrently on 2 threads, just staggered by that
~30µs instead of starting simultaneously. The third thread that `-t 3` adds
now sits idle for ~98% of each timestep while still costing task-spawn/
sync overhead 730×/year, and having more active threads exposes more
scheduling variance to background system load (`benchmark/run_benchmarks.jl`
`threads` mode spawns real OS processes per thread count, so it does see
this).

Measured across 5 independent `threads`-mode sweeps (real NCEP dataset,
`-t 3` reps=3-4 each; machine had the usual OneDrive-sync/SearchIndexer
background load — see the benchmark skill's noise-source list, not
suppressed here since the point is `-t 3`'s *sensitivity* to it):

| Threads | Speedup vs. `-t 1` across 5 sweeps | Consistency |
|---|---|---|
| `-t 2` | 1.24×, 1.48×, 1.50×, 1.58×, 1.62× (mean ~1.48×) | Low variance — always a clear win |
| `-t 3` | 1.02×, 1.14×, 1.38×, 1.52×, 1.65× (mean ~1.34×) | High variance — sometimes ties `-t 2`, sometimes barely beats serial |
| `-t 4` | 0.79×, 1.12×, 1.14×, 1.15×, 1.53× (mean ~1.15×) | High variance, lowest mean |

`-t 3` is *not* reliably worse than `-t 2` in every single run (it won
outright in 2 of the 5 sweeps) — the honest finding is that `-t 2` is now
the **consistent, low-risk default**, while `-t 3`/`-t 4` no longer have a
dependable edge to justify recommending them over `-t 2` by default. This
reverses the earlier "`-t 3`, not `-t 2`" guidance for this codebase's
current (`Float32`) state — a direct, load-bearing consequence of §2.3,
not a new independent optimization. Updated: `benchmark/run_benchmarks.jl`'s
`threads` mode/skill guidance, `README.md`'s threading blurb.

Other identified, unimplemented opportunities (not benchmarked further):
ensemble/parameter-sweep parallelism via `Threads`/`Distributed`/GPU arrays
(likely the bigger practical payoff given how the model is actually used —
sensitivity sweeps, CO₂ scenarios); low-cost `Float32`/GPU-array swaps via
multiple dispatch (the physics kernels have no hardcoded `Float64`, only
`CirculationWorkspace`'s fields do, per §2.3); AD-based calibration via
`Enzyme.jl`/`ForwardDiff.jl` for gradient-based parameter tuning.

**Pre-implementation integration test** (kept for methodology reference):
combining the 3-way split with §2.12's turbo fixes in an isolated
`tendencies!`+`time_loop!` harness measured ~2.8× (2.83–2.88s →
1.00–1.04s) — higher than the ~1.78× (2.7s → 1.513s) the real
`greb_model!` shows post-implementation, since the isolated harness doesn't
carry `output!`/`diagnostics!`/bookkeeping overhead. One measurement gotcha
caught along the way, since fixed in the real implementation too: the
first, un-warmed-up timing showed the combo 2.5× *slower*, purely because
its new functions hadn't been JIT-compiled yet (unlike `greb_model!`'s
`PrecompileTools`-baked ones) — resolved with a warm-up call before timing.

---
