---
name: docs-check
description: Verify GREBClimate.jl's documentation matches the code — run every Julia snippet in README.md and docs/src/ against the real package, and check file references, exported names and version claims. Use before claiming docs are updated, after renaming or changing any exported function's signature, and when reviewing changes to README.md, DATA_README.md or docs/.
---

# Docs check

Documentation in this repo has produced a silent scientific error, not just
confusion. The README's headline quick-start once ran the model on an all-zero
climatology and reported a −40 °C world as success, because the snippet
discarded `load_greb_jld2!`'s return value. It printed a ✅ and finished with
`all(isfinite, Ts) == true`. Nobody following the README could have known.

Reading a snippet is not checking it. **Run it.**

Full account of that bug and the rest of the 2026-08-21 sweep:
`.claude/notes/onboarding.md`.

## When to use this

- Before telling the user the docs are correct or updated.
- After changing the signature, name, or return value of anything exported
  from `GREBClimate` — snippets in four files may reference it.
- When reviewing edits to `README.md`, `DATA_README.md`, or `docs/src/*.md`.
- After a rename. The `GREB` → `GREBClimate` rename left a stale
  `docs/Manifest.toml` that broke the documented docs build for weeks.

## The checks

### 1. Run the snippets

Every runnable Julia block in `README.md` and `docs/src/*.md` must execute
against the real package. There are four places that document the same flow —
`README.md`'s quick-start, `README.md` §Running the Model, `docs/src/tutorial.md`,
and `examples/run_greb.jl` — and they drift independently.

```bash
julia --project=. examples/run_greb.jl greb_input_data
```

That covers the canonical path. For a documented snippet that
`examples/run_greb.jl` does not cover, paste it verbatim into a scratch file
and run it. Verbatim matters: the bug above only reproduces if you keep the
missing `fields=`.

**Sanity-check the numbers, not just the exit code.** Two independent tells,
on a `:full_model` control run against the real dataset:

| | healthy | ran on a zero climatology |
|---|---|---|
| the model's own printed line | `1970  14.77  29.82  3.69` | `1970  -40.0  -40.0  -40.0` |
| `mean(result.ctrl[1].Ts)` | 276.64 K | 233.15 K |

(The printed 14.77 °C and the 276.64 K array mean are different quantities —
the printout is the model's own global/land/ocean summary, the other is an
unweighted mean over grid cells. Don't try to reconcile them; just compare each
against its own column.)

Both cases exit 0, print `✅ All GREB data loaded successfully`, and satisfy
`all(isfinite, Ts)`.

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

### 2-4. Static checks (scripted)

```bash
julia .claude/skills/docs-check/check_docs.jl
```

Stdlib-only, does not load the package, runs in seconds, exits nonzero on
failure. It checks:

- **Local markdown links resolve**, relative to the linking file. (`claude/BENCHMARKS.md`
  was referenced from the tutorial and had never existed.)
- **Julia version claims agree** across `README.md`, `docs/src/index.md` and
  `Project.toml`'s compat bound. They disagreed — 1.10 / 1.9 / 1.10 — until
  2026-08-21.
- **`@ref`'d names are actually exported** from `GREBClimate`.
- **No dead `.jld2` files.** Reported as a NOTE, not a failure: 11 files
  (148 MB) are currently unreferenced, pending the open `.new` question in
  `.claude/notes/data-distribution.md`. This check is interpolation-aware —
  `src/io.jl` builds ENSO filenames as `"erainterim.tsurf.$suffix.forcing.jld2"`,
  so a naive grep reports ten false positives. If you write your own version of
  this check, expand `$suffix` first.

Do not replace this with an ad-hoc `grep` over filenames. The obvious one
resolves every path from the repo root regardless of which file the link is in,
and matches bare package names, so it emits ~70 false positives and gets
ignored.

### 5. Documenter cross-references and docstrings

The script cannot validate `@ref` targets or missing docstrings — only a real
build can:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

`checkdocs = :exports` is set, so a missing docstring for an exported name
fails here. A clean build prints no warnings other than the expected
"could not auto-detect the building environment. Skipping deployment."

If `Pkg.develop` errors with "Refusing to add package ... with the same UUID
already exists", `docs/Manifest.toml` is stale — delete it, it is gitignored
and regenerates. That stale manifest, left by the `GREB` → `GREBClimate`
rename, broke the documented docs build for weeks.

## Reporting

Say which checks you ran and what they output. "Docs look right" is not a
result; "ran the quick-start, global-mean Ts 14.77 °C, docs build clean" is.

If a check fails, fix the docs rather than the check — unless the code is what
is wrong, which happens: the converter's documented default input directory
did not exist, and the fix belonged in the script.

## Related

- `.claude/notes/onboarding.md` — what this skill was built from, including
  what was checked and found *accurate* (don't redo those).
- `.claude/skills/dev-notes/SKILL.md` — where to record what you find.
