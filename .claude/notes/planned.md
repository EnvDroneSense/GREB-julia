# Planned next steps

**Status:** Planned — not started.

Forward-looking only. When something here is done, move it to the note for its
topic and record the result there.

<!-- Split out of the former monolithic claude/IMPROVEMENTS.md on 2026-08-21.
     Section numbers in cross-references (§N.M) refer to that original file;
     see INDEX.md for where each section now lives. -->

---

## 5. Next up

Nothing blocking. Remaining open items:
- §4.8 (`forcing()` dispatch style) and §4.10 (CI lint/format) — both
  explicitly deferred by user decision, not forgotten.
- Ensemble/GPU parallelism and AD-based calibration (§2.11) remain
  documented ideas only, not benchmarked.
- §2.14's 3 investigated-not-implemented findings (CMIP5/ENSO anomaly
  reload caching, `postprocess.jl`'s Symbol dispatch, `CirculationWorkspace()`
  allocation count) — each needs a design decision or isn't worth it at
  current scale, not quick safe wins.
- §8's 5 newly found bugs remain documented, not fixed (explicit user
  decision to keep that pass report-only).

Landed across the last two passes: §2.12's original 4 `@turbo` fixes plus
its 4 follow-up fixes, §2.13's `forcing()` mask-recomputation fix, §3.0's
test-suite sharding, and §2.11's 3-way `tendencies!` thread split (originally
needed `-t 3`, not `-t 2`, to realize its full ceiling — since reversed by
the `Float32` conversion, §2.3; `-t 2` is now the recommended default — see
§2.11's post-`Float32` re-review) — see each section above and
`AUDIT_LOG.md` for the real (not just benchmarked) numbers.

> ⚠️ Every performance change must be validated against a reference run —
> "faster" only counts if output is unchanged within tolerance. Every
> physics change must be validated against the Fortran reference — grep the
> actual subroutine text, don't recall it from memory.

---
