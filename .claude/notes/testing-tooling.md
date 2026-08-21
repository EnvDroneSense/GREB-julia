# Testing, tooling & reproducibility

**Status:** Done — CI, tests and tooling in place.

<!-- Split out of the former monolithic claude/IMPROVEMENTS.md on 2026-08-21.
     Section numbers in cross-references (§N.M) refer to that original file;
     see INDEX.md for where each section now lives. -->

---

## 3. Testing, tooling & reproducibility

- **Golden/snapshot regression test** ✅ done: real 1yr control + 1yr scenario
  `:full_model` run against `greb_dataset_jld2/`, checked against reference
  monthly global-mean `Ts`/`Ta`/`q` within `atol=1e-6`; `@test_skip`s when the
  dataset isn't present. `RUN_GOLDEN` env var (default on) additionally
  skips it locally even when the dataset is present.
- **Branch-coverage tests**: `log_eva`/`log_rain` sweep uses a zipped 5-run
  sweep instead of a 20-run cross product (independent branches, same
  coverage, half the runtime); a direct test confirms `log_eva` values
  produce genuinely different output.
- **CI** ✅, **`[compat]` bounds** ✅, **dead dependencies removed**
  (`NCDatasets`/`StaticArrays`/`Statistics`) ✅, **misplaced docstrings**
  fixed (§0.19) ✅, **benchmark suite** (`benchmark/run_benchmarks.jl`, driven by
  `.claude/skills/benchmark/SKILL.md`; results recorded in
  `.claude/notes/performance.md`) ✅, **`Documenter.jl` site** ✅ (`docs/`, deployed
  to GitHub Pages on push to `main` via `.github/workflows/docs.yml`).
- **Testing gaps** ✅ closed: `qflux_correction!`'s asymmetric `Ta` handling
  confirmed intentional (matches Fortran, `greb.model.mscm.f90:540-597`),
  regression test added; the ~20 `:experiment` symbols unreachable via
  `create_experiment_config` now covered by a direct-construction sweep;
  `build_monthly_climatology`/`apply_scenario_anomalies` have direct unit
  tests; both JLD2 loaders' "file present" branches covered via synthetic
  `mktempdir`+`jldopen` datasets.

### 3.0 Parallelizing the test suite across processes — ✅ done
Splitting `runtests.jl`'s 24 testsets into 2 process-level shards (not
`Threads.@threads` — `Test.jl`'s testset tracking isn't thread-safe) along
an existing boundary:

| | Time |
|---|---|
| Sequential (one process) | 106.5s |
| Split A alone (17 cheap testsets) | 34.9s |
| Split B alone (7 heavy testsets, incl. golden regression) | 65.1s |
| **Parallel (2 concurrent `julia` processes)** | **68.1s** |

**1.56×** on 14 available cores — short of 2× because the split is
unbalanced (B takes ~2× A). A rebalanced/finer split would likely get
closer to the ceiling.

**Implemented**: `test/runtests.jl`'s 24 testsets are now split into
`run_light_tests()` (17) and `run_heavy_tests()` (7, incl. golden
regression) at this exact boundary, gated by `GREB_TEST_SHARD`
(`"all"`/`"light"`/`"heavy"`, default `"all"` — plain `Pkg.test()`/`]test`
is unaffected). `.github/workflows/ci.yml` gained a `shard: [light, heavy]`
matrix axis. Verified locally: 142 + 177 = 319, matching the unsharded
total exactly; both shards pass independently.

---

---

## 2026-08-21: test process was single-threaded

`Pkg.test()` passed no `julia_args`, so the test subprocess ran with one
thread. Since `tendencies!` gates its `Threads.@spawn` branch on
`Threads.nthreads() > 1`, the multithreaded circulation path was never
executed by any test or CI job — it had only ever been validated by
benchmarking.

Fixed by a subprocess-based equivalence test plus `JULIA_NUM_THREADS=2` in CI.
Details and the reasoning in [`performance.md`](performance.md).

General lesson for this repo: anything gated on a startup-fixed setting
(`nthreads`, `--check-bounds`, `JULIA_*`) cannot be covered by an in-process
test. It needs a spawned subprocess, or it is not covered at all.
