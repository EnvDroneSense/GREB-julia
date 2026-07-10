# Workplan: Reorganize GREB into a Julia Package

**Status:** Draft for review — no repo changes will be made until you approve this plan.

**Goal:** Turn this repository into a standard Julia package with a `src/` module and a
`test/` subpackage, **without rewriting the GREB model logic** ("reorganize, don't
re-code").

---

## 1. What the repo is today

| File | Size | What it is |
|:--|:--|:--|
| `GREB_julia.jl` | 4303 lines / 143 cells | **Pluto notebook** holding the whole model: 32 functions, 4 structs (`PhysicsConfig`, `CirculationWorkspace`, `MonthlyAccumulator`, `TimeState`), many global `const`s — **plus** 68 markdown doc cells and 14 `@bind` interactive-UI cells. |
| `PultoUI.jl` | 45 KB | An **unrelated** copy of the PlutoUI.jl "Interactivity" featured demo notebook. Not GREB code. |
| `launch_pluto.jl` | 6 lines | Helper that activates the env and launches Pluto. |
| `Project.toml` | — | An **environment**, not a package: it has `[deps]` but **no** `name`, `uuid`, or `version`. |
| `Manifest.toml` | 158 KB | Pinned full dependency tree. |
| `README.md`, `DATA_README.md`, `LICENSE` | — | Docs + MIT license. |

There is **no** `src/`, no `test/`, and no package identity yet.

---

## 2. The core tension (needs your decision)

Your instruction is *"do not add or modify any existing code, just reorganize."* A Pluto
notebook, however, **cannot** be loaded as a package module as-is. Turning it into
`src/GREB.jl` unavoidably requires **mechanical** (not logical) changes:

- stripping Pluto cell markers (`# ╔═╡ …`) and the cell-order footer,
- dropping the notebook-only cells (`md"""…"""` docs, `@bind` widgets, the `@bind` macro),
- wrapping the definitions in `module GREB … end` with `using`/`export`,
- possibly re-ordering top-level `const`/`struct` cells so definitions precede use
  (Pluto stores cells in *display* order, not a valid script execution order).


Clarification: leave the notebook as is and copy into a `notebook` folder. Extract all functions from the notebook for src and all tests for test directory. Do not change the functions themselves (we will revise in the next section.)

**None of these touch the model's physics logic** — the functions, structs, and constants
are copied verbatim. But they are edits, so I want your explicit sign-off on the
interpretation. I read *"do not add or modify existing code"* as **"do not rewrite the GREB
model logic"**, while the module wrapper + test scaffold are the structural work the task
inherently requires.

I propose **Plan B** below and offer the ultra-conservative **Plan A** as the alternative.

### Plan A — Layout-only (truly zero code edits)
Create the standard *directory layout* and package metadata, but leave the model as a Pluto
notebook untouched. Result: the repo *looks* like a package but **`using GREB` would not
work** and there is nothing meaningful to test. Low value; included only for completeness.

### Plan B — Real package (recommended)
Extract the model definitions from the notebook into a proper `src/GREB.jl` module,
copying the physics **verbatim** (only the mechanical changes listed above). Keep the
Pluto notebook as an interactive front-end that does `using GREB`. This yields a package
you can `using`, precompile, and test — which is the point of the request.

The rest of this plan assumes **Plan B**.

---

## 3. Target structure (Plan B)

```
GREB.jl/
├── Project.toml          # becomes a PACKAGE manifest (name, uuid, version, [deps], [compat])
├── README.md             # updated install/usage snippet (docs only, no model code)
├── DATA_README.md        # unchanged
├── LICENSE               # unchanged
├── .gitignore            # NEW: ignore Manifest.toml, /docs/build, etc.
├── src/
│   └── GREB.jl           # NEW module: model code extracted verbatim from the notebook
├── test/
│   ├── Project.toml      # NEW: test-only deps (Test, …)
│   └── runtests.jl       # NEW: smoke/integration tests
├── notebooks/
│   ├── GREB_julia.jl     # the original Pluto notebook, MOVED here (UI + docs preserved)
│   └── launch_pluto.jl   # MOVED here
└── examples/  (optional) # a plain-Julia driver script, if wanted
```

Open questions on structure are collected in §6.

---

## 4. Step-by-step execution (after approval)

Each step is small and independently reviewable. I will keep the original notebook working.

1. **Package identity.** Rewrite `Project.toml` to a package manifest:
   - `name = "GREB"`, freshly generated `uuid`, `version = "0.1.0"`.
   - Keep the existing `[deps]` that the *module* needs. **Note:** `Plots` and `PlutoUI`
     are notebook/UI concerns, not model concerns — I'll propose moving `PlutoUI` out of the
     package deps (it belongs to the notebook env) and decide `Plots` per §6.
   - Add a `[compat]` section (at least `julia = "1.9"`, matching the README).

2. **Create `src/GREB.jl`.** Mechanically extract from `GREB_julia.jl`:
   - all `struct`/`mutable struct` definitions,
   - all top-level `const` definitions (grid dims, lookup tables, `MonthlyRecord`, etc.),
   - all 32 `function` definitions,
   wrapped in `module GREB` with the needed `using` (Statistics, NCDatasets, StaticArrays,
   LoopVectorization) and an `export` list. **No function body is altered.** Where Pluto's
   display order breaks top-level ordering, I reorder *cells only* so definitions precede
   use. I'll produce a short mapping (notebook cell → src location) for your review.

3. **Move the notebook** `GREB_julia.jl` → `notebooks/`, and adapt its dependency cell to
   `using GREB` so the interactive UI keeps working on top of the package. (This is the one
   notebook edit; the model cells become thin re-exports or are removed in favor of the
   package. Exact approach flagged in §6, item E.)

4. **Move** `launch_pluto.jl` → `notebooks/`.

5. **Decide on `PultoUI.jl`** (unrelated demo) — remove or archive; see §6 item B.

6. **Create the test package:**
   - `test/Project.toml` with `[deps]` = `Test` (+ anything a smoke test needs).
   - `test/runtests.jl`. Scope depends on §6 item D — either an empty scaffold or a minimal
     smoke test (load module, build a `PhysicsConfig`, run a very short integration, assert
     output shapes/types). Any test code is **new** and additive; it does not touch the model.

7. **`.gitignore`** + `Manifest.toml` handling — see §6 item C.

8. **Update `README.md`** install/usage to reflect the package layout (docs only).

9. **Verify:** run `julia --project -e 'using GREB'` (precompile check) and
   `julia --project=test test/runtests.jl`. Report results. No push/commit unless you ask.

---

## 5. Explicitly out of scope

- No changes to model physics, numerics, or algorithms.
- No performance tuning, no bug fixes (incl. the README's noted qflux issue).
- No new features, no CI setup, no docs site — unless you request them.

---

## 6. Decisions I need from you

| # | Decision | My recommendation |
|:--|:--|:--|
| **A** | Overall approach: **Plan A** (layout only) or **Plan B** (real package)? | **Plan B** |
| **B** | `PultoUI.jl` (unrelated PlutoUI demo): delete / move to `notebooks/examples/` / keep as-is? | **Delete** (it isn't GREB code) |
| **C** | `Manifest.toml`: gitignore & remove from VCS (standard for packages), or keep it committed? | **Gitignore & remove** |
| **D** | Tests: **empty scaffold** only, or a **minimal smoke test** I write (new, additive code)? | **Minimal smoke test** |
| **E** | Notebook after extraction: keep it as a full interactive UI running `using GREB`, or reduce it to a thin example? | **Keep interactive**, backed by `using GREB` |
| **F** | Package name `GREB` and version `0.1.0` OK? (matches the `GREB.jl` remote) | **Yes** |
| **G** | `Plots` / `PlutoUI` deps: move to the notebook environment rather than package deps? | **Move both out** of package `[deps]` |

---

## 7. Risks / notes

- **Cell ordering:** the biggest mechanical risk is top-level `const`/`struct` ordering when
  flattening the notebook into a module. I'll validate by actually precompiling (`using GREB`)
  before declaring done.
- **`Plots` dep is barely used** (one plotting cell, line ~3111) and belongs to the notebook,
  not the model — reinforces item G.
- **`BenchmarkTools`/`Profile`/`ProfileSVG`** are `using`-ed in one benchmark cell but are
  **not** in `Project.toml` deps today; I'll leave that benchmark cell in the notebook env,
  not the package.
- Everything is reversible via git; I'll work in small commits (only if/when you ask me to
  commit) and won't push without approval.

---

**Please review and tell me your choices for §6 (A–G), or edit this file directly. I will not
modify the repo until you approve.**
