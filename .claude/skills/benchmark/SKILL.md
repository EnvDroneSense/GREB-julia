---
name: benchmark
description: Run GREB.jl's timing/allocation benchmarks (benchmark/run_benchmarks.jl) and report results honestly, accounting for machine noise. Use when the user asks to benchmark, time, profile, or check performance/allocations of the model, or to compare thread counts.
---

# GREB benchmark

`benchmark/run_benchmarks.jl` is a dependency-free harness with four modes,
all against the real dataset, post-JIT-warmup:

| Mode | What it does |
|---|---|
| `year` (default) | Wall-clock time for a real 1-simulated-year `:full_model` control run. |
| `stages` | Per-timestep breakdown: times each of `tendencies!`'s stages (`circulation!` for `Ta`/`q`, `SWradiation!`, `LWradiation!`, `hydro!`, `deep_ocean!`) individually and reports each one's share of the total. |
| `threads` | Runs `year` in separate `-t N` subprocesses (thread count is fixed at Julia startup) and reports relative speedup vs. `-t 1`. |
| `alloc` | Reports bytes allocated by one `tendencies!` call - a regression guard for the zero-allocation hot path. |

The model runs natively in `Float32` throughout (climatology, workspace
buffers, model state - see `claude/IMPROVEMENTS.md` §2.3); current baseline
on this machine is **~0.65-0.75s/simulated year at `-t 2`** (down from the
pre-`Float32` ~1.1-1.2s/year). Treat any reading far outside that band as
worth double-checking against the noise sources below before reporting it
as a real regression.

**Default to `-t 2`, not `-t 3`.** This reverses older guidance in this
repo (and in git history of this file) - `-t 3` used to be the recommended
count. `Float32` shifted the per-timestep cost balance enough that the
third thread `-t 3` adds no longer has reliable work to do (see the
"Threading recommendation" note below) - `-t 2` is now the consistent,
low-variance choice.

## Steps

1. Pick the mode that answers the actual question:
   - Plain "how fast is the model" / comparing before-vs-after a change →
     `year` (the default).
   - "Where does the time actually go" / did a kernel-specific change move
     the needle → `stages`.
   - "Is the 3-way threading split still working" → `threads` (this spawns
     its own subprocesses at each thread count, so run it with any `-t`;
     the outer process's thread count doesn't matter for this mode).
   - "Did this change add allocations to the hot path" → `alloc`.

   ```bash
   julia --project=. -t 2 benchmark/run_benchmarks.jl year
   julia --project=. -t 1 benchmark/run_benchmarks.jl stages
   julia --project=. -t 1 benchmark/run_benchmarks.jl threads
   julia --project=. -t 1 benchmark/run_benchmarks.jl alloc
   ```

   Pass a JLD2 data directory as the next positional argument if it isn't
   at the default `../greb_dataset_jld2` location, or set `GREB_DATA`. The
   legacy single-argument call form (`... run_benchmarks.jl <dir>`, no mode)
   still works and defaults to `year`.

2. **Before trusting a slow or surprising `year`/`threads` reading**, rule
   out known sources of noise on this machine - don't report a single
   reading as fact:
   - **Stale precompile cache**: if source files changed since the last
     `using GREB`, the first run in a fresh process pays a one-time recompile.
     Check with `julia --project -e 'using Pkg; Pkg.precompile()'` - if it
     recompiles `GREB` (not just prints the manifest-resolution warning),
     re-run the benchmark afterward.
   - **Post-restart background load**: right after a reboot, `OneDrive.exe`,
     `OneDrive.Sync.Service.exe`, and `SearchIndexer.exe` (this repo lives in
     a OneDrive-synced folder) contend for disk/CPU for a couple of minutes.
     Check with `tasklist | grep -i -E "onedrive|searchindexer"` and, if
     present and the machine was recently restarted, wait ~1-2 min and re-run.
   - General wall-clock noise on this shared/dev machine is real even without
     either cause above - a ~50-150ms swing across repeated `year` runs is
     normal at the current ~0.7s baseline (proportionally similar to the
     ~300-600ms swing seen at the older ~1.1-1.2s baseline).
   - `stages`/`alloc` are much less noise-prone (µs/byte-scale, no dataset
     I/O in the timed region) - trust a single reading there more readily
     than a single `year`/`threads` reading.

3. If comparing thread counts, prefer `threads` over running `year` manually
   at each count - it does exactly that, back-to-back, and prints the
   speedup table directly. **Threading recommendation (re-reviewed
   2026-08-13, post-`Float32` — see `claude/IMPROVEMENTS.md` §2.11): `-t 2`
   is the reliable default.** Circulation is now ~98% of per-timestep cost
   (`Float32` sped up the other stages proportionally more than
   circulation), so the third lane `-t 3` used to justify is nearly free -
   `-t 2` already captures almost all the real parallelism via work-stealing
   once the tiny synchronous "rest" finishes. Measured across 5 independent
   `threads` sweeps: `-t 2` gave a consistent ~1.5× (range 1.24-1.62×,
   *low* variance - always a clear win); `-t 3`/`-t 4` ranged anywhere from
   ~1.0× (barely better than serial) to ~1.65× (occasionally best) - *high*
   variance, no longer a dependable edge over `-t 2`. Don't take one
   `threads` run as settling a `-t 2` vs `-t 3` question either way -
   `-t 3`/`-t 4` are the ones sensitive enough to background load
   (OneDrive/SearchIndexer, see below) to flip either direction; run it a
   few times before concluding a real regression or improvement at those
   counts specifically. `-t 4`+ has no reliable reason to help further on
   this grid size regardless. Compare *relative* speedup, not just absolute
   numbers - absolute baselines drift with machine load, but the relative
   improvement from a real code change should hold up.

4. Report mean/min/max (for `year`/`threads`) or the per-stage table (for
   `stages`) or the byte count (for `alloc`), and call out plainly if noise
   (any of the causes above, or high run-to-run variance) makes a reading
   untrustworthy rather than presenting it as a clean result.

## Other ways to check performance (not scripted, use when they fit better)

- **`@code_native`/gather-instruction count**: when a change touches a
  `@turbo` loop with neighbor-index lookups (`circulation.jl`'s
  `diffusion!`/`advection!`), a wall-clock number alone can hide a width
  mismatch (e.g. `Float32` data gathered through `Int64` indices - see
  `claude/IMPROVEMENTS.md` §2.3's `lon_jm1`/etc. `Int32` fix). Check with
  `InteractiveUtils.code_native` and grep the output for `"gather"` and for
  `cvtss2sd`/`cvtsd2ss`-style conversions - a rising gather count or any
  conversion instruction in a loop that shouldn't have one is a real signal
  wall-clock timing alone can miss or mis-attribute to noise.
- **Isolated kernel microbenchmarks**: for a change scoped to one function
  (not the whole pipeline), a standalone `@benchmark` (BenchmarkTools, not
  a dependency of this package but fine to add temporarily to a scratch
  script) on just that function - with realistic, not synthetic-zero, input
  data - gives a cleaner signal than `year`'s whole-model number, which
  mixes in dataset loading and every other stage's noise. This is how every
  number in `claude/IMPROVEMENTS.md` §2.3/§2.15 was actually produced before
  being confirmed against the real `year`/`stages` numbers here.
- **Numeric regression alongside any timing change**: a faster wrong answer
  is not a win. Pair any timing comparison with a correctness check -
  `test/runtests.jl`'s golden-regression test, or a direct diff against a
  saved reference run - especially for anything touching numerical
  precision (`Float32` conversions, reassociated arithmetic).
