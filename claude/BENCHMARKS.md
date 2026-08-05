# GREB.jl — Benchmark Log

Running log of `BenchmarkTools.@benchmark` results for GREB's hot per-timestep kernels, appended by `benchmark/run_benchmarks.jl` each time it's run. Compare entries across runs to see the effect of each optimization pass.

Reproduce: `julia --project=benchmark benchmark/run_benchmarks.jl "<tag>"`

---

## 2026-08-05 23:41 — baseline — before safe-mechanical optimization pass

- branch: `optimize/safe-mechanical-pass`  commit: `609d48f` (dirty working tree)
- data: real `greb_dataset_jld2/` (NCEP)
- julia: `1.12.6`

| kernel | min | median | memory | allocs |
|---|---|---|---|---|
| `convergence!` | 1.160 μs | 1.190 μs | 0 bytes | 0 |
| `diffusion!` | 216.700 μs | 231.800 μs | 140.00 KiB | 2698 |
| `advection!` | 259.900 μs | 284.500 μs | 114.20 KiB | 3416 |
| `SWradiation!` | 25.600 μs | 27.100 μs | 48.91 KiB | 388 |
| `LWradiation!` | 161.200 μs | 172.200 μs | 1.22 KiB | 28 |
| `hydro!` | 30.500 μs | 32.800 μs | 4.17 KiB | 88 |
| `seaice!` | 7.380 μs | 7.700 μs | 1.44 KiB | 38 |
| `deep_ocean!` | 10.700 μs | 11.200 μs | 1.77 KiB | 34 |
| `circulation!` | 16.064 ms | 18.973 ms | 5.96 MiB | 146736 |
| `tendencies!` | 26.881 ms | 33.108 ms | 11.97 MiB | 294010 |

