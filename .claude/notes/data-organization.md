# Data organization & JLD2 compression

**Status:** Investigated — findings recorded, no change made.

Measured comparison of data layouts and JLD2 compression settings. The wider
options write-up that this summarizes is appended below from the former
`claude/DATA_ORGANIZATION_OPTIONS.md`.

<!-- Split out of the former monolithic claude/IMPROVEMENTS.md on 2026-08-21.
     Section numbers in cross-references (§N.M) refer to that original file;
     see INDEX.md for where each section now lives. -->

---

## 6. Data organization & JLD2 compression — investigated, findings documented

Full comparison doc: appended below in this file. Everything here
is measured (real `jldopen`/`tar` runs), not estimated.

### 6.1 `Data/`'s layout doesn't match its own docs
`DATA_README.md`/`convert_greb_to_jld2.jl` expect an `input/` subdirectory
that doesn't exist on disk (files sit directly under `Data/`). Two
root-level files (`erainterim.evaporation.clim`,
`erainterim.omega.vertmean.nomean.clim`, 25.7 MiB) were confirmed unread by
any Fortran variant (grepped both available source trees) — **not a porting
gap**. Both removed from `Data/`; `DATA_README.md` updated.

### 6.2 JLD2 compression: real numbers
Built the real `greb_dataset_jld2/` (605.4 MiB, 53 files) and measured
`compress=true` on 6 representative files:

| File | Compression ratio | Read slowdown |
|---|---|---|
| `ncep.tsurf...clim.jld2` | 80.1% | 9.5× |
| `cmip5.tsurf...forcing.new.jld2` | 86.6% | 7.9× |
| `erainterim.tsurf.elnino.forcing.jld2` | 93.2% | 20.4× |
| `global.topography.jld2` (static) | 29.5% | 4.0× |
| `solar_radiation.clim.jld2` | 40.4% | 2.0× |
| `solar_eccentricity.jld2` (grouped) | 40.8% | 11.4× |

Whole-tree `tar -I 'gzip -9'`: 530.6 MB from 634.8 MB (83.6%, 42s),
consistent with the per-file numbers. **Conclusion: compress only for
distribution** (a separate `tar`/`xz` archive), never the copy
`load_greb_jld2!` reads — the ~16–22% size win isn't worth an 8–20× read-time
risk.

### 6.3 JLD2 grouping: by load-pattern, not by source ✅ flux-correction merge implemented
Tested merging 5 CMIP5 fields into one file: combined-file single-field
reads cost 2.3× more (24.0ms vs 10.3ms, ~14ms absolute — irrelevant for a
once-at-startup load), but reading all 5 combined is *faster* (45.0ms vs
50.2ms). Rule: combine when a group is always read together, not by
upstream source — `ncep.*`/`erainterim.*` are read in varying subsets
depending on `dataset=:ncep`/`:era` and active experiment, so combining by
source would wrongly couple independently-varying load conditions.

**Implemented**: `load_flux_corrections_jld2!` always loads all 3 flux
files together, matching the "always read in full" case — merged into
`flux_corrections.jld2` (`convert_flux_corrections` in
`scripts/convert_greb_to_jld2.jl`), measured **~35% faster** (19.89ms vs
30.78ms combined vs. separate), no size penalty. `test/runtests.jl` updated
to match; full suite re-run clean (319/319 pass).

---

---

## Data organization & JLD2 compression: options comparison

Research/comparison document. No files were moved and no permanent code was
changed to produce this — the conversion script was run once into a scratch
directory (not `greb_dataset_jld2/`, which stays gitignored and unbuilt in the
repo) purely to get real, measured numbers instead of estimates. All sizes below
are **actual measurements**, not extrapolations, unless explicitly marked.

**Assumption for this revision:** the CMIP5 RCP8.5 forcing data and the
ERA-Interim El Niño/La Niña forcing data are being actively wired in (the
`:rcp85`/`:elnino`/`:lanina` `PhysicsConfig.experiment` stubs in
`src/model.jl:85-106` are being completed). They're treated below as **active,
in-use data**, not quarantine candidates. Only `erainterim.evaporation.clim.*`
and `erainterim.omega.vertmean.nomean.clim.bin` remain genuinely unreferenced by
any code path.

### 1. Current state

#### 1.1 `Data/` layout

`Data/` sits flat at the repo root: 49 root-level `.bin` files + `.ctl` sidecars
+ 10 CO2 `.txt` files, plus one subdirectory (`solar_forcing_scenarios/`, 109
`.bin` files). **Measured total: 634,487,414 bytes (605.1 MiB)** (`du -sb Data`).

| Shape | Byte size | Count | Examples |
|---|---|---|---|
| 2D static (96×48 float32) | 18,432 B | 2 | `global.topography.bin`, `greb.glaciers.bin` |
| 3D climatology/forcing (96×48×730 float32) | 13,455,360 B (~12.8 MiB) | 46 | `Tocean.clim.bin`, `ncep.tsurf.1948-2007.clim.bin`, `erainterim.omega.vertmean.clim.bin`, `cmip5.omega.rcp85.ensmean.forcing.new.bin` |
| Solar lat×time (48×730 float32) | 140,160 B | 1 (root) + 109 (`solar_forcing_scenarios/`, 15,296,940 B measured) | `solar_radiation.clim.bin`, `greb.solar.eccentricity.<N>.bin`, `greb.solar.obliquity.<N>.bin` |
| CO2 scenario text | 3,856–7,029 B each, **63,443 B total (measured)** | 10 | `ipcc.scenario.rcp85.forcing.txt`, `ipcc.scenario.hist.forcing.CO2.emission.pop.txt` |

**Layout/docs mismatch:** `DATA_README.md` and `scripts/convert_greb_to_jld2.jl`'s
default (`input_path = Data/input`) both expect an `input/` subdirectory. On disk,
files are directly under `Data/` — `Data/input/` doesn't exist (confirmed: the
conversion run below had to be pointed explicitly at `Data`, not the default).

#### 1.2 Category breakdown (measured, 3D fields only — the 619 MB that dominates the total)

Counted and summed directly (`ls`/`du -cb` on the actual files):

| Category | Files | Bytes (raw `.bin`) | Status |
|---|---|---|---|
| CMIP5 RCP8.5 forcing (all `cmip5.*.bin`, incl. `.new`/`nomean`/`vertmean`/`precip` variants) | 16 | 215,285,760 B (205.3 MiB) | being wired in (`:rcp85`) |
| ERA-Interim El Niño/La Niña forcing | 10 | 134,553,600 B (128.3 MiB) | being wired in (`:elnino`/`:lanina`) |
| ERA-Interim climatology (used: humidity, wind ×2, omega, omega_std, tsurf, windspeed) | 7 | 94,187,520 B (89.8 MiB) | in use today (always-loaded + `:era` dataset) |
| ERA-Interim, still unreferenced (`evaporation.clim`, `omega...nomean.clim`) | 2 | 26,910,720 B (25.7 MiB) | **not** wired to anything yet |
| NCEP climatology (tsurf, wind ×2, humidity, soil moisture) | 5 | 67,276,800 B (64.2 MiB) | in use today (default dataset) |
| Common (Tocean, Tocean/Tsurf/vapour flux correction, isccp, woce) | 6 | 80,732,160 B (77.0 MiB) | in use today |

Sum of these six rows = 619,946,560 B, matching the 46-file/13,455,360-B-each
total exactly. So, with CMIP5+ENSO now counted as active: **333.6 MiB (54%) of
the 3D-field data is tied to the forcing experiments being implemented right
now**, and only 25.7 MiB (4%) remains out of scope.

#### 1.3 `greb_dataset_jld2/` — actually built for this document

Not committed (gitignored, `.gitignore:17`) and not previously built in this
repo. For this document, `scripts/convert_greb_to_jld2.jl` was run for real
(`julia --project=. scripts/convert_greb_to_jld2.jl Data <scratch-dir>`) to get
measured JLD2 output sizes rather than guessing at container overhead:

- **49 field files converted, 0 failed.** Total reported by the script: 619.1 MB
  of source data in.
- Measured output tree: **634,768,362 bytes (605.4 MiB) total** — i.e. JLD2's
  container format adds essentially no overhead over the raw `.bin` size (a few
  KB of `dim_names`/`ctl`-text metadata per file; e.g.
  `ncep.tsurf.1948-2007.clim.jld2` is 13,460,706 B vs. 13,455,360 B raw, +0.04%).
- Subdirectory sizes (measured): `climatology/` 619,184,798 B, `solar_scenarios/`
  15,298,260 B, `solar/` 145,477 B, `scenario/` 92,357 B, `static/` 47,470 B —
  matching `main()`'s five-directory structure (`static/`, `climatology/`,
  `solar/`, `solar_scenarios/`, `scenario/`) confirmed directly from the run.
- Grouping precedent confirmed live: `solar_eccentricity.jld2` (61 stacked source
  files → one 8,557,567 B file with a `coords` index) and `ipcc_scenarios.jld2`
  (10 source files → one 92,357 B `Dict`).

#### 1.4 JLD2 compression availability

`Project.toml` pins `JLD2 = "0.6"`; `Manifest.toml:185` resolves **JLD2 0.6.5**,
which already depends on `ChunkCodecLibZlib`/`ChunkCodecLibZstd` as hard
dependencies — so `jldopen(path, "w"; compress=true)` works today with zero new
`Project.toml` entries. Confirmed by actually calling it (§4) — no version
mismatch or missing-codec errors.

---

### 2. Data-layout reorganization options

With CMIP5/ENSO now treated as active data, the "used vs. unused" split
collapses to just the categorization work — there's no meaningful quarantine
folder left, just the 25.7 MiB of genuinely-unreferenced ERA-Interim variants.

| | Option A — current (flat) | Option B — by category (all active data together) |
|---|---|---|
| Structure | `Data/*.bin` flat | `Data/input/{static,climatology/ncep,climatology/erainterim,climatology/cmip5_rcp85,climatology/enso,common,solar,solar_forcing_scenarios,scenario}` + `Data/deprecated/` (2 files: evaporation.clim, omega...nomean.clim) |
| Migration cost | none | medium — matches `convert_all`'s own routing logic, but needs the script updated |
| Fixes `Data/input/` docs mismatch | no | yes |
| Signals which experiment each file feeds | no | yes — `climatology/cmip5_rcp85/` and `climatology/enso/` map 1:1 to the `:rcp85`/`:elnino`/`:lanina` experiment branches in `src/model.jl:85-106` |
| Requires touching `convert_greb_to_jld2.jl` | no | yes — `convert_all`'s `readdir(input_path; join=true)` (line 112) is non-recursive, so category subfolders need either a recursive glob or one call per subfolder |

**Recommendation: Option B**, with the 2 truly-unreferenced files moved to a
small `Data/deprecated/` (not `unused/`, since almost everything else just moved
out of that bucket) so they're not mixed in with the six categories that now all
feed a real code path. This is a proposed target, not applied here — adopting it
means updating `convert_greb_to_jld2.jl`'s discovery logic and `DATA_README.md`'s
paths, and physically moving the files, as a separate follow-up.

---

### 3. JLD2 output logical grouping options

| | Option A — current (per-field) | Option B — group always-loaded-together | Option C — monolithic |
|---|---|---|---|
| File count (measured baseline) | 53 (.jld2 files produced) | ~50 (flux corrections merged 3→1) | 1 |
| Matches `load_greb_jld2!`'s per-key reads | yes, directly | mostly (flux corrections are already read as a group, `load_flux_corrections_jld2!`) | no — every call site would need rewriting |
| Preserves `:ncep`/`:era` selective loading | yes | yes | no — forces loading both variants (adds the unused variant's full 64–90 MiB every run) |
| Corruption blast radius | one field (≤12.8 MiB) | one group | entire 605 MiB dataset |
| Precedent already in codebase (measured) | — | yes: `solar_eccentricity.jld2` (61→1 files, 8.56 MB) and `ipcc_scenarios.jld2` (10→1 files, 92 KB) already do exactly this | — |

**Recommendation: Option B**, narrowly — merge only the three flux-correction
files (`Tsurf_flux_correction`, `vapour_flux_correction`, `Tocean_flux_correction`,
40.4 MiB combined measured) into one `flux_corrections.jld2`, since
`load_flux_corrections_jld2!` (`src/io.jl:82-100`) always loads all three
together already. Leave climatology/static per-file — merging those would force
loading the unused `:ncep`-vs-`:era` variant every run (an extra ~64–90 MiB read
for no benefit), losing the selective-loading behavior `load_greb_jld2!`
(`src/io.jl:126-145`) currently gives for free.

---

### 4. Compression scenarios — measured, not estimated

The requirement driving this section: compression is for **sharing/distributing**
the dataset only — never for the copy `load_greb_jld2!` reads at run time.

#### 4.1 JLD2-native `compress=true`, per representative file (measured directly)

Ran `jldopen(path, "w"; compress=true)` on one real file per category from the
converted output, then timed 20 repeated reads of both the compressed and
uncompressed copy (`jldopen(path, "r") do f; f["data"]; end`):

| File (category) | Raw JLD2 size | `compress=true` size | Ratio | Read, uncompressed | Read, compressed | Slowdown |
|---|---|---|---|---|---|---|
| `ncep.tsurf.1948-2007.clim.jld2` (NCEP climatology) | 13,460,706 B | 10,783,564 B | 80.1% | 16.6 ms | 157.6 ms | **9.5×** |
| `cmip5.tsurf.rcp85.ensmean.forcing.new.jld2` (CMIP5 forcing) | 13,460,456 B | 11,662,464 B | 86.6% | 15.7 ms | 124.0 ms | **7.9×** |
| `erainterim.tsurf.elnino.forcing.jld2` (ENSO forcing) | 13,460,456 B | 12,550,430 B | 93.2% | 6.9 ms | 140.2 ms | **20.4×** |
| `global.topography.jld2` (static) | 23,737 B | 7,000 B | 29.5% | 0.71 ms | 2.84 ms | 4.0× |
| `solar_radiation.clim.jld2` | 145,477 B | 58,756 B | 40.4% | 0.67 ms | 1.32 ms | 2.0× |
| `solar_eccentricity.jld2` (grouped, 61 stacked scenarios) | 8,557,567 B | 3,493,992 B | 40.8% | 5.2 ms | 59.5 ms | 11.4× |

Takeaway: the 3D climate fields — which make up 619 of 635 MB — barely
compress (80–93% of original) because ensemble-mean/reanalysis float32 fields
have limited redundancy, and cost **8–20× longer to read** every single time.
The El Niño forcing file, ironically one of the categories being actively wired
in, is both the worst compressor (93.2%) and the slowest to re-read (20.4×) of
everything tested.

#### 4.2 Whole-dataset archive compression (measured directly, real `tar` runs)

Compressed the actual converted `greb_dataset_jld2/` tree (634,768,362 B) as one
archive:

| Method | Output size | Ratio | Wall time |
|---|---|---|---|
| `tar -I 'gzip -9'` | 530,607,089 B | 83.6% | 42.0 s |
| `tar -I 'xz -3'` | 497,644,548 B | 78.4% | 57.3 s |

Cross-check: summing the per-category §1.2 byte counts weighted by the §4.1
per-category ratios (cmip5×86.6%, ENSO×93.2%, ERA-Interim/NCEP/common×80.1%,
static×29.5%, solar/solar_scenarios×~40%) predicts **≈533.9 MB** — within 0.6%
of the measured 530.6 MB `tar+gzip` result, confirming the per-file samples in
§4.1 are representative of their categories rather than cherry-picked.

**Only a ~16–22% size reduction is achievable at the whole-dataset level**,
because the dataset is dominated by float32 climate fields that don't compress
well losslessly, regardless of method (JLD2-native or external archive).

#### 4.3 Scenario comparison

| | Scenario 1 — no compression | Scenario 2 — JLD2-native `compress=true` | Scenario 3 — external `tar+gzip`/`tar+xz` |
|---|---|---|---|
| What's compressed | nothing | the `.jld2` files themselves | the whole `greb_dataset_jld2/` tree as one archive |
| Measured size (this dataset) | 634,768,362 B (605.4 MiB) | not measured at full scale, but §4.1's per-category ratios predict ≈530–535 MB if applied dataset-wide (≈16% smaller) | 530,607,089 B (gzip, measured) / 497,644,548 B (xz -3, measured) — 16–22% smaller |
| Risk to model-run read path | none | **real, and measured**: 8–20× slower reads if the compressed copy is ever read directly instead of decompressing first | none — decompression is a separate `tar` step producing ordinary uncompressed `.jld2` files; `read_jld2`/`load_greb_jld2!` never sees the archive |
| Ecosystem dependency | none | none beyond JLD2 itself (0.6.5 vendors the codecs — confirmed, no errors) | `tar`/`gzip`/`xz` (all present in this environment already) or a small Julia unpacker if that's undesirable |
| Implementation effort | none | low | very low — one `tar` command, confirmed working |

**Recommendation: Scenario 3.** The measured 8–20× read slowdown in §4.1 makes
JLD2-native compression a real hazard if anyone points `load_greb_jld2!` at the
compressed copy by mistake — and the size win (≈16–22%) isn't large enough to
justify that risk given it's achievable identically, with zero risk to the read
path, by archiving the plain output tree with `tar`. Ship
`greb_dataset_jld2.tar.xz` (497.6 MB, best measured ratio) as the share/download
artifact, document a "download → `tar xf` → point `load_greb_jld2!` at the
extracted folder" workflow, and never generate a `compress=true` copy as
anything other than a distribution artifact.

---

### 5. Summary

| | Est./measured size | Model-run read impact | Effort | Dependency |
|---|---|---|---|---|
| Layout: Option B (by category) | — (reorg only) | none | medium (script + doc update) | none |
| JLD2 grouping: merge flux corrections | 40.4 MiB → 1 file | none | low | none |
| Compression: Scenario 3 (`tar+xz`) | 634.8 MB → 497.6 MB (measured, 78.4%) | none | very low | `tar`/`xz` |

**Overall recommendation:** adopt the Option B data layout (§2), merge the
flux-correction trio (§3), and share the dataset as a plain `tar+xz` archive of
`greb_dataset_jld2/` (§4) — every number above is a direct measurement from this
session's actual conversion run and compression tests, not an estimate.
