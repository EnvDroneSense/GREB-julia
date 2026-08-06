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

---

## 2026-08-06 07:39 — after second bug-hunting pass (do_conv/co2_part/log_rain fixes, dead-code cleanup)

- branch: `improvements`  commit: `56421bd`
- data: real `greb_dataset_jld2/` (NCEP)
- julia: `1.12.6`

| kernel | min | median | memory | allocs |
|---|---|---|---|---|
| `convergence!` | 1.240 μs | 1.310 μs | 0 bytes | 0 |
| `diffusion!` | 256.700 μs | 279.500 μs | 140.00 KiB | 2698 |
| `advection!` | 324.600 μs | 673.500 μs | 114.20 KiB | 3416 |
| `SWradiation!` | 27.200 μs | 62.600 μs | 48.91 KiB | 388 |
| `LWradiation!` | 167.300 μs | 276.400 μs | 1.22 KiB | 28 |
| `hydro!` | 33.400 μs | 108.200 μs | 4.17 KiB | 88 |
| `seaice!` | 20.525 μs | 22.200 μs | 1.44 KiB | 38 |
| `deep_ocean!` | 38.500 μs | 40.300 μs | 1.77 KiB | 34 |
| `circulation!` | 28.009 ms | 34.539 ms | 5.96 MiB | 146736 |
| `tendencies!` | 55.108 ms | 64.352 ms | 11.97 MiB | 294010 |

---

## 2026-08-06 10:46 — post state-struct refactor (ClimateFields/ModelState) + Fortran-audit fixes 0.7-0.10

- branch: `improvements`  commit: `4a2f30b` (dirty working tree)
- data: real `greb_dataset_jld2/` (NCEP)
- julia: `1.12.6`

| kernel | min | median | memory | allocs |
|---|---|---|---|---|
| `convergence!` | 1.130 μs | 2.250 μs | 0 bytes | 0 |
| `diffusion!` | 45.400 μs | 116.700 μs | 0 bytes | 0 |
| `advection!` | 22.600 μs | 74.200 μs | 0 bytes | 0 |
| `SWradiation!` | 16.400 μs | 18.600 μs | 38.25 KiB | 48 |
| `LWradiation!` | 165.900 μs | 176.200 μs | 0 bytes | 0 |
| `hydro!` | 25.000 μs | 25.600 μs | 0 bytes | 0 |
| `seaice!` | 5.086 μs | 5.514 μs | 0 bytes | 0 |
| `deep_ocean!` | 7.775 μs | 8.150 μs | 0 bytes | 0 |
| `circulation!` | 1.791 ms | 1.922 ms | 0 bytes | 0 |
| `tendencies!` | 3.689 ms | 4.247 ms | 38.25 KiB | 48 |

---

## 2026-08-06 11:42 — threading (diffusion!/advection!) -t 1, pre-existing baseline for comparison

- branch: `improvements`  commit: `d169451` (dirty working tree)
- data: real `greb_dataset_jld2/` (NCEP)
- julia: `1.12.6`

| kernel | min | median | memory | allocs |
|---|---|---|---|---|
| `convergence!` | 1.240 μs | 1.310 μs | 0 bytes | 0 |
| `diffusion!` | 47.800 μs | 54.400 μs | 512 bytes | 7 |
| `advection!` | 21.700 μs | 68.000 μs | 672 bytes | 7 |
| `SWradiation!` | 16.100 μs | 21.400 μs | 38.25 KiB | 48 |
| `LWradiation!` | 188.600 μs | 205.400 μs | 0 bytes | 0 |
| `hydro!` | 25.600 μs | 26.900 μs | 0 bytes | 0 |
| `seaice!` | 4.867 μs | 5.050 μs | 0 bytes | 0 |
| `deep_ocean!` | 7.250 μs | 7.575 μs | 0 bytes | 0 |
| `circulation!` | 1.753 ms | 1.893 ms | 27.75 KiB | 336 |
| `tendencies!` | 3.788 ms | 4.223 ms | 93.75 KiB | 720 |

---

## 2026-08-06 11:43 — threading (diffusion!/advection!) -t auto (14 logical CPUs)

- branch: `improvements`  commit: `d169451` (dirty working tree)
- data: real `greb_dataset_jld2/` (NCEP)
- julia: `1.12.6`

| kernel | min | median | memory | allocs |
|---|---|---|---|---|
| `convergence!` | 1.170 μs | 1.240 μs | 0 bytes | 0 |
| `diffusion!` | 140.600 μs | 220.000 μs | 6.30 KiB | 72 |
| `advection!` | 132.100 μs | 205.900 μs | 8.48 KiB | 72 |
| `SWradiation!` | 16.000 μs | 18.400 μs | 38.25 KiB | 48 |
| `LWradiation!` | 158.300 μs | 174.000 μs | 0 bytes | 0 |
| `hydro!` | 25.600 μs | 26.800 μs | 0 bytes | 0 |
| `seaice!` | 4.886 μs | 5.000 μs | 0 bytes | 0 |
| `deep_ocean!` | 7.800 μs | 8.150 μs | 0 bytes | 0 |
| `circulation!` | 8.443 ms | 10.871 ms | 354.75 KiB | 3456 |
| `tendencies!` | 17.791 ms | 23.703 ms | 747.75 KiB | 6960 |

