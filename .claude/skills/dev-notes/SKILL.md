---
name: dev-notes
description: Record an investigation, measurement, or decision in GREBClimate.jl's .claude/notes/ — what was tried, what the result was, and what is planned. Use after finishing any investigation, benchmark pass, audit, or design decision worth keeping, and when asked what has been tried before or what is planned next.
---

# Dev notes

`.claude/notes/` holds what the code and git history cannot say on their own:
what was tried, what the measurement was, what was rejected and why, and what
is still open. `.claude/notes/INDEX.md` is the map.

**These notes are gitignored** — local working material, not part of the
package. A fresh clone will not have them, so never point user-facing
documentation (`README.md`, `docs/src/`, `CHANGELOG.md`) at a path under
`.claude/notes/`: it would be a dead link for everyone but the author. Refer to
"the maintainers' working notes" in prose instead.

Read the relevant note **before** starting an investigation. Several things in
here were already measured and deliberately rejected; re-deriving them wastes a
session and risks contradicting a recorded decision.

## The structure

One file per topic. Every file opens with:

```markdown
# Title

**Status:** <one of the four below>

<one or two sentences on what this note is for>

---
```

Four statuses, and they are the point of the layout — the previous structure
interleaved fixed work, open questions and future plans inside single
750+ line documents, so you could not tell state from a heading:

| Status | Means |
|:-------|:------|
| **Done** | Landed. Kept for the reasoning, not the to-do. |
| **Investigated** | Measured, findings recorded, deliberately no change made. |
| **Open** | Needs a decision — domain input or a design call, not more reading. |
| **Planned** | Not started. |

`INDEX.md` groups the notes under those same four headings. A note can be
Mixed (see `data-distribution.md`) when one topic genuinely spans states.

## Adding to the notes

**Prefer extending an existing note over adding a file.** Twelve topics is
about right; twenty single-purpose files is the sprawl this replaced.

1. Check `INDEX.md` for a note that already owns the topic.
2. Append your finding there, with the date and what you actually ran.
3. If the note's status changed (Open → Done, Planned → Investigated), update
   its **Status:** line *and* move its row to the right section of `INDEX.md`.
4. Only create a new file for a genuinely new topic — then add it to `INDEX.md`
   in the same edit. A note not in the index is a note nobody will read.

When something in `planned.md` gets done, move it out to the note for its
topic and record the *result* there. `planned.md` should only ever contain
things that have not started.

## What makes a note worth keeping

Record the reasoning and the measurement — the things that are expensive to
recover:

- **Numbers you actually observed**, with the command that produced them.
  Not "faster" — `0.65 s/simulated year at -t 2, down from 1.1 s`.
- **What was rejected, and why.** More valuable than what was adopted. The
  reason `-t 2` beats `-t 3` is recorded; without it someone re-runs that
  experiment every six months.
- **Blockers by name.** "Differentiability is hard" is useless;
  "`@turbo` kernels are opaque to `ForwardDiff`" is actionable.
- **Constraints that are not visible in the code.** That the precompile
  workload must run data-free is why `allow_uninitialized` exists; nothing in
  `src/` says so.

Do **not** record what the repo already tells you: file layout, what a
function does, that a test passes, or a summary of a diff. Those are in the
code and the git log, and duplicating them means two things to keep in sync.

## Cross-references

Section numbers of the form §N.M inside these notes refer to the *old*
`claude/IMPROVEMENTS.md` numbering, from before the 2026-08-21 split. The
mapping table is at the bottom of `INDEX.md`. Leave those references alone —
rewriting them would break the audit trail in `audit-history.md`.

When adding a cross-reference between notes, use a plain relative link
(`[performance.md](performance.md)`) so it resolves in a local editor. These
files are not published anywhere, so nothing needs to resolve on GitHub.

## Related

- `.claude/skills/docs-check/SKILL.md` — verifying *user-facing* docs. Notes are
  for maintainers; the README and Documenter site are for users, and the two
  should not be confused for one another.
- `.claude/skills/benchmark/SKILL.md` — produces the numbers that belong in
  `performance.md`.
