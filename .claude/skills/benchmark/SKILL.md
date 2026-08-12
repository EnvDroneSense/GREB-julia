---
name: benchmark
description: Run the GREB.jl 1-simulated-year timing benchmark (benchmark/run_benchmarks.jl) and report results honestly, accounting for machine noise. Use when the user asks to benchmark, time, or check performance of the model, or to compare thread counts.
---

# GREB benchmark

Runs `benchmark/run_benchmarks.jl`, a simple dependency-free script that times
a real 1-simulated-year `:full_model` control run (post-JIT-warmup) a few
times and reports mean/min/max.

## Steps

1. Run the benchmark at the requested thread count (default `-t 3`, the
   recommended count since `tendencies!`'s 3-way split needs
   `Threads.nthreads() >= 3` to realize full concurrency):

   ```bash
   julia --project -t 3 benchmark/run_benchmarks.jl
   ```

   Pass a JLD2 data directory as an argument if it isn't at the default
   `../greb_dataset_jld2` location, or set `GREB_DATA`.

2. **Before trusting a slow or surprising reading**, rule out known sources of
   noise on this machine — don't report a single reading as fact:
   - **Stale precompile cache**: if source files changed since the last
     `using GREB`, the first run in a fresh process pays a one-time recompile.
     Check with `julia --project -e 'using Pkg; Pkg.precompile()'` — if it
     recompiles `GREB` (not just prints the manifest-resolution warning),
     re-run the benchmark afterward.
   - **Post-restart background load**: right after a reboot, `OneDrive.exe`,
     `OneDrive.Sync.Service.exe`, and `SearchIndexer.exe` (this repo lives in
     a OneDrive-synced folder) contend for disk/CPU for a couple of minutes.
     Check with `tasklist | grep -i -E "onedrive|searchindexer"` and, if
     present and the machine was recently restarted, wait ~1-2 min and re-run.
   - General wall-clock noise on this shared/dev machine is real even without
     either cause above — a ~300-600ms swing across repeated runs is normal.

3. If comparing thread counts (e.g. verifying the 3-way threading speedup is
   still intact), run both back-to-back and compare *relative* speedup, not
   just absolute numbers — absolute baselines drift with machine load, but
   the relative improvement from a real code change should hold up:

   ```bash
   julia --project -t 1 benchmark/run_benchmarks.jl
   julia --project -t 3 benchmark/run_benchmarks.jl
   ```

4. Report mean/min/max for each run, and call out plainly if noise (any of
   the causes above, or high run-to-run variance) makes a reading
   untrustworthy rather than presenting it as a clean result.
