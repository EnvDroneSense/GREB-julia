# Dev notes index

Working notes for GREBClimate.jl — investigations, measurements, and decisions
that the code and git history do not record on their own. Not user
documentation: users want `README.md` and the
[Documenter site](https://EnvDroneSense.github.io/GREBClimate.jl/).

One file per topic. Every file opens with a **Status:** line, so this index
stays a map rather than a second copy of the content. See
[`../skills/dev-notes/SKILL.md`](../skills/dev-notes/SKILL.md) for the
conventions when adding to these.

## Done — landed, recorded for the reasoning

| Note | What it covers |
|:-----|:---------------|
| [`bugs-fixed.md`](bugs-fixed.md) | Every bug found and fixed, across two sweeps. The outcomes. |
| [`audit-history.md`](audit-history.md) | Pass-by-pass forensic narrative of *how* those bugs were found. |
| [`structure.md`](structure.md) | The split from one 2245-line file with ~40 globals into `src/` topical files with explicit state structs. |
| [`performance.md`](performance.md) | Every optimization pass: threading, `Float32`, `@turbo` — what was measured, what was rejected. The largest note. |
| [`testing-tooling.md`](testing-tooling.md) | Test suite, CI, reproducibility tooling. |
| [`fortran-validation.md`](fortran-validation.md) | Switch/experiment behaviour validated against the original Fortran GREB. |
| [`onboarding.md`](onboarding.md) | 2026-08-21 usability review: the silent −40 °C README bug and the doc-accuracy fixes it triggered. |

## Investigated — findings recorded, no change made

| Note | What it covers |
|:-----|:---------------|
| [`data-organization.md`](data-organization.md) | Dataset layout and JLD2 compression options, measured rather than estimated. |
| [`ensembles-and-autodiff.md`](ensembles-and-autodiff.md) | Ensemble/multi-run costs, and what specifically blocks differentiability. |

## Open — needs a decision, not more reading

| Note | What it covers |
|:-----|:---------------|
| [`open-questions.md`](open-questions.md) | Findings needing climate-science or design input. |
| [`registry-readiness.md`](registry-readiness.md) | Julia General Registry checklist, partially complete. |

## Planned — not started

| Note | What it covers |
|:-----|:---------------|
| [`planned.md`](planned.md) | Next steps not yet begun. |
| [`data-distribution.md`](data-distribution.md) | DataDeps distribution (planned). Also holds the JLD2 dataset audit — 11 dead files, 148 MB — and the resolved `.new` CMIP5 question. |

---

## History

Until 2026-08-21 these lived in a top-level `claude/` directory as four files,
two of which exceeded 750 lines and interleaved fixed work, open questions and
future plans in the same document. They were split by topic and moved here,
which also removed the `claude/` vs `.claude/` name collision. Section
references of the form §N.M inside these notes point at the old
`claude/IMPROVEMENTS.md` numbering:

| Old location | Now |
|:-------------|:----|
| `claude/IMPROVEMENTS.md` §0, §8 | `bugs-fixed.md` |
| `claude/IMPROVEMENTS.md` §1 | `structure.md` |
| `claude/IMPROVEMENTS.md` §2 | `performance.md` |
| `claude/IMPROVEMENTS.md` §3 | `testing-tooling.md` |
| `claude/IMPROVEMENTS.md` §4 | `open-questions.md` |
| `claude/IMPROVEMENTS.md` §5 | `planned.md` |
| `claude/IMPROVEMENTS.md` §6 + `claude/DATA_ORGANIZATION_OPTIONS.md` | `data-organization.md` |
| `claude/IMPROVEMENTS.md` §7 | `fortran-validation.md` |
| `claude/IMPROVEMENTS.md` §9 | `registry-readiness.md` |
| `claude/AUDIT_LOG.md` | `audit-history.md` |
| `claude/EXPANSION.md` | `ensembles-and-autodiff.md` |
