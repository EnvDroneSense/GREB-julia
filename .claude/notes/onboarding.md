# Onboarding & documentation accuracy

**Status:** Done — 2026-08-21 usability review; fixes landed. One item deferred.

A new-user usability pass: follow the documentation literally, from a clean
start, and record where it misleads. Everything below was reproduced by
running it, not inferred by reading.

---

## The bug this pass existed to find

`README.md`'s headline quick-start called `load_greb_jld2!(...)` and
**discarded its return value**, then called `greb_model!(RunSpec(), cfg)`
without `fields=`. Since `greb_model!` defaulted `fields = ClimateFields()`
(all zeros), the documented first-run experience silently simulated a zero
climatology:

| | global-mean Ts, control month 1 | model's own printout |
|---|---|---|
| README snippet | **233.15 K (−40.0 °C)** | `1970  -40.0  -40.0  -40.0` |
| Correct call | 276.64 K (**14.77 °C**) | `1970  14.77  29.82  3.69` |

It printed `✅ All GREB data loaded successfully`, returned 12 well-formed
monthly records, and `all(isfinite, Ts)` was `true`. Nothing signalled failure.

**Fix:** `ClimateFields` gained a `loaded::Bool` flag (set by
`load_greb_jld2!`), and `greb_model!` refuses unloaded fields unless the caller
passes `allow_uninitialized=true`.

**Why a flag rather than a hard error:** 14 call sites legitimately run
data-free — the `@compile_workload` precompile block (which must not require a
580 MB download) and 12 config/scenario-plumbing tests. A bare throw would
break precompilation and turn two `@test_throws` assertions green for the wrong
reason. Data-free runs stay possible; they just have to say so.

## Other findings, all fixed

- A **second instance of the same trap** in README §Loading Data, which taught
  the discard-the-return-value pattern in prose.
- `README.md` line 7 described the project as something that "runs in an
  interactive Pluto.jl notebook" — contradicting the docs, which correctly
  present a package with an optional notebook front-end.
- All user-facing data documentation was written around obtaining raw `.bin`
  files and running the converter. That is a maintainer path; users get the
  `.jld2` bundle. `DATA_README.md` and the converter are now labelled as such.
- `scripts/convert_greb_to_jld2.jl` defaulted its input directory to
  `Data/input`, which does not exist — the real layout is flat `Data/` plus
  `Data/solar_forcing_scenarios/`. The documented no-argument invocation
  therefore always failed. Default corrected, with an actionable error.
- `docs/src/tutorial.md` referenced `claude/BENCHMARKS.md`, which never existed.
- `docs/src/index.md` claimed Julia 1.9; `Project.toml` requires 1.10.
- README's project-structure tree omitted `benchmark/`, `archive/`, `Data/`,
  `greb_input_data/`, `DATA_README.md`, and `notebooks/PultoUI.jl`.
- `docs/Manifest.toml` was stale — it still recorded the package under its
  pre-rename name `GREB` with the same UUID (see `67680e0`), so the documented
  local docs build failed with "Refusing to add package". Removed; it is
  gitignored and regenerates.

## Checked and found accurate — left alone

Worth recording so the next pass does not redo it:

- `docs/src/switches.md` documents all 35 `PhysicsConfig` fields, every one of
  the ~45 symbols it names exists in `src/`, and its "Known Limitations"
  section is correct (`log_vapor_dmc` declared but never read; `log_conv` not
  exposed as a `create_experiment_config` keyword).
- README §Running the Model is correct, as are `examples/run_greb.jl`,
  `docs/src/tutorial.md`'s code, and `MonthlyRecord`'s documented field list.

## Deferred

The model prints unconditionally to stdout — emoji progress lines, per-year
diagnostics, an `@info` on hydrology setup — with no `verbose` switch. Callers
currently wrap runs in `redirect_stdout(devnull)`, which the test suite does
throughout. Fine for notebook use, wrong for a library called in a loop. Not
addressed in this pass; would be a small, self-contained change.

`examples/run_greb.jl`'s optional plot uses `@eval using Plots` inside a
`try`. `Plots` is not a dependency of the package environment, so on a machine
that has it available elsewhere Julia spends ~4.5 minutes precompiling it
before the `using` fails and the script correctly reports "Plots.jl not
available". The message is accurate and the run still succeeds; the wasted
precompile is the annoyance. Checking `Base.find_package("Plots") !== nothing`
first would avoid it. Not changed in this pass.

Deliberately **not** treated as a problem: the ~38 s first-call precompilation
cost is a chosen trade — heavy first run, fast subsequent ones.

## Prevention

`.claude/skills/docs-check/SKILL.md` exists because of this pass: it runs the
documentation's code against the real package and checks references, exports
and version claims. The −40 °C bug, the dead `BENCHMARKS.md` link and the
1.9/1.10 split were all mechanically detectable.
