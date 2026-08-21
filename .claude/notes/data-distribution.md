# Data distribution & JLD2 dataset audit

**Status:** Mixed — audit **done**; `.new` question **resolved** (2026-08-21); DataDeps distribution **planned**.

How users get the `.jld2` dataset, and what is actually in it. For the *layout*
and compression trade-offs of the dataset itself, see
[`data-organization.md`](data-organization.md).

---

## Planned: distribution via DataDeps.jl

**Status: planned, not started.**

Today a user must request the ~580 MB `.jld2` bundle from the maintainers,
unpack it, and pass the path in (`load_greb_jld2!(dir)` / `GREB_DATA`). The raw
`.bin` files are collated from several upstream sources and are not
redistributed, so regenerating from source is a maintainer path only — see
`DATA_README.md`.

The intent is to serve the prepared `.jld2` bundle through
[DataDeps.jl](https://github.com/oxinabox/DataDeps.jl), so first use fetches
and caches it and the manual step disappears.

Not yet decided / needed before this can land:

- A hosting location with a stable URL and a SHA256 for the bundle.
- Whether to ship one bundle or split it (e.g. core climatology vs. optional
  solar/scenario families) so a minimal run doesn't pull 580 MB.
- Licence/attribution wording for the redistributed derivative, given the
  mixed upstream provenance (NCEP, ERA-Interim, ISCCP, WOCE, CMIP5, IPCC).

Docs are deliberately written against today's manual path, so nothing in
`README.md` promises a mechanism that does not exist yet.

## Done: dataset audit (2026-08-21)

Audited all 49 `.jld2` files in `greb_input_data/` (580 MB).

**No byte-identical duplicates.** All 49 files md5-hash uniquely.

**11 files (148 MB, ~25% of the dataset) are read by nothing** — not `src/`,
`test/`, `benchmark/`, `examples/`, `scripts/`, or the Pluto notebook:

```
cmip5.evaporation.rcp85.ensmean.forcing.jld2          cmip5.omegastd.rcp85.ensmean.forcing.new.jld2
cmip5.evaporation.rcp85.ensmean.forcing.new.jld2      cmip5.omegastd.rcp85.ensmean.vertmean.forcing.jld2
cmip5.meridional.wind.rcp85.ensmean.forcing.new.jld2  cmip5.precip.rcp85.ensmean.forcing.new.jld2
cmip5.omega.rcp85.ensmean.forcing.new.jld2            cmip5.tsurf.rcp85.ensmean.forcing.new.jld2
cmip5.omega.rcp85.ensmean.forcing.nomean.new.jld2     cmip5.windspeed.rcp85.ensmean.forcing.new.jld2
                                                      cmip5.zonal.wind.rcp85.ensmean.forcing.new.jld2
```

**Root cause:** `scripts/convert_greb_to_jld2.jl` converts indiscriminately —
`filter(endswith(".bin"), readdir(input_path))`. Every `.bin` present in the
input directory becomes a `.jld2`, whether the model reads it or not.

Note the ENSO anomaly files (`erainterim.*.{elnino,lanina}.forcing.jld2`) are
**not** in that list despite never appearing as literals in the source: they
are built by string interpolation in `src/io.jl` (`"erainterim.tsurf.$suffix.forcing.jld2"`).
Any future dead-file check has to resolve those, or it will report ten false
positives.

## Resolved: the non-`.new` CMIP5 anomalies are correct

**Status: closed 2026-08-21. The Julia port was already reading the right files.**

`src/io.jl` reads the **non-`.new`** CMIP5 anomaly files. The concern was that
in upstream GREB a `.new` suffix often marks a *corrected* field, which would
have made this a silent scientific error rather than housekeeping. (A previous
pass had already flagged the pair as undecided — see `audit-history.md`,
2026-08-06 entry, where it was a reason not to bake CMIP5 into a shared blob.)

The two sets are genuinely different data, not duplicate copies. Measured over
the full `(96, 48, 730)` fields:

| field | variant | mean | std | corr(old,new) |
|:------|:--------|-----:|----:|--------------:|
| tsurf [K] | non-`.new` (used) | 3.578 | 2.434 | 0.961 |
| | `.new` | 4.206 | 2.471 | |
| zonal wind [m/s] | non-`.new` (used) | 0.1597 | 0.7638 | 0.899 |
| | `.new` | 0.1673 | 0.8387 | |
| omega [Pa/s] | non-`.new` (used) | -5.0e-5 | 0.00197 | 0.556 |
| | `.new` | 1.5e-4 | 0.00282 | |

Both are plausible RCP8.5 anomaly fields, but the choice is not cosmetic: the
`.new` surface-temperature anomaly is **0.63 K warmer in the global mean**, and
the two omega fields correlate at only 0.56.

**Resolution.** The official MSCM Fortran GREB
([christianstassen/greb-official](https://github.com/christianstassen/greb-official),
branch `official`) hardcodes the non-`.new` names in `greb.shell.mscm.f90`:

```fortran
open(31,file='../input/cmip5.tsurf.rcp85.ensmean.forcing.bin'...
open(32,file='../input/cmip5.zonal.wind.rcp85.ensmean.forcing.bin'...
open(33,file='../input/cmip5.meridional.wind.rcp85.ensmean.forcing.bin'...
open(34,file='../input/cmip5.omega.rcp85.ensmean.forcing.bin'...
open(35,file='../input/cmip5.windspeed.rcp85.ensmean.forcing.bin'...
```

That is exactly the set `src/io.jl` loads, so the port matches the reference
implementation and **no code change is needed**. The same file confirms the
ENSO forcing names the Julia code builds by interpolation
(`erainterim.{tsurf,zonal.wind,meridional.wind,omega,windspeed}.{elnino,lanina}.forcing.bin`).

Caveat worth keeping: this establishes what the *published MSCM* GREB uses. If
the `.new` files came from a later unpublished revision, that revision is not
what this port targets — and switching to it would be a deliberate scientific
change, not a fix.

## Done: converter allowlist + drift test

Both follow-ups were blocked on the `.new` question and landed once it closed.

- `scripts/convert_greb_to_jld2.jl` now filters on an explicit
  `MODEL_FIELD_NAMES` allowlist (33 entries) instead of globbing every `.bin`.
  Against the real `Data/` directory it keeps 33 and skips exactly the 11 dead
  files, and warns if an allowlisted field is absent from the input. `--all`
  restores the old convert-everything behaviour for inspecting an extra field.
- `test/runtests.jl` asserts the allowlist and `src/io.jl`'s actual loads agree
  **in both directions**, so neither can drift. Both sides are parsed as text:
  including the converter would run its top-level code, and io.jl's loads are
  spread across several functions with no single introspectable list. The test
  expands the ENSO `$suffix` interpolation on both sides.
- `.claude/skills/docs-check/check_docs.jl` reports unreferenced `.jld2` files
  as a NOTE, so the situation stays visible without failing CI.

## Not done: deleting the 11 dead files

Left alone deliberately — deleting from a 580 MB local dataset is the
maintainer's call, not a cleanup to do silently. They cost 148 MB and nothing
reads them; regenerating the bundle with the new allowlist simply won't emit
them. Keeping the `.new` files somewhere out-of-band has some value as a record
of the variant that was evaluated and rejected.
