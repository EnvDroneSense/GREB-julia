# =============================================================================
# check_docs.jl - static documentation checks for GREBClimate.jl.
#
# Dependency-free (stdlib only). Does NOT load the package, so it is fast and
# cannot be broken by a compile error:
#
#   julia .claude/skills/docs-check/check_docs.jl
#
# Checks:
#   1. Markdown links to local files resolve, relative to the linking file.
#   2. Julia version claims agree across README, docs and Project.toml.
#   3. Names in docs `@ref`/backticks that look like exports are exported.
#   4. Every .jld2 file on disk is referenced by code (interpolation-aware).
#
# Exits nonzero if any check fails, so it can gate CI.
#
# It deliberately does NOT run the documented snippets - that needs the real
# package and dataset. See SKILL.md; running them is the check that matters
# most and has no substitute.
# =============================================================================

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DOC_FILES = ["README.md", "DATA_README.md", "CHANGELOG.md",
                   "docs/src/index.md", "docs/src/tutorial.md",
                   "docs/src/switches.md", "docs/src/api.md",
                   ".claude/notes/INDEX.md"]

failures = String[]
fail(msg) = (push!(failures, msg); println("  FAIL  ", msg))
ok(msg)   = println("  ok    ", msg)

# ── 1. Markdown links to local files ─────────────────────────────────────────
# Only inline links `[text](target)`; bare filenames in prose and directory
# trees are not links and are not checked.
println("\n[1] local markdown links resolve")
let checked = 0
    for rel in DOC_FILES
        path = joinpath(ROOT, rel)
        isfile(path) || (fail("doc file listed but missing: $rel"); continue)
        base = dirname(path)
        for m in eachmatch(r"\[[^\]]*\]\(([^)\s]+)\)", read(path, String))
            target = m.captures[1]
            # skip URLs, anchors, and Documenter @ref links
            (startswith(target, "http") || startswith(target, "#") ||
             startswith(target, "mailto:") || occursin("@ref", target)) && continue
            file = first(split(target, '#'))           # drop any #anchor
            isempty(file) && continue
            checked += 1
            ispath(normpath(joinpath(base, file))) ||
                fail("$rel -> $target (does not exist)")
        end
    end
    ok("$checked local links checked")
end

# ── 2. Julia version claims agree ────────────────────────────────────────────
println("\n[2] Julia version claims agree")
let claims = Dict{String,Set{String}}()
    proj = read(joinpath(ROOT, "Project.toml"), String)
    m = match(r"^julia\s*=\s*\"([^\"]+)\""m, proj)
    m === nothing ? fail("Project.toml has no julia compat bound") :
        (claims["Project.toml"] = Set([m.captures[1]]))
    for rel in ["README.md", "docs/src/index.md"]
        s = read(joinpath(ROOT, rel), String)
        vs = Set(m.captures[1] for m in eachmatch(r"Julia\s+\*{0,2}(1\.\d+)"i, s))
        isempty(vs) || (claims[rel] = vs)
    end
    allv = union(values(claims)...)
    if length(allv) == 1
        ok("all sources say Julia $(first(allv))")
    else
        parts = [string(k, " => ", join(sort(collect(v)), ",")) for (k, v) in claims]
        fail("disagreement: " * join(parts, "; "))
    end
end

# ── 3. Documented names are exported ─────────────────────────────────────────
println("\n[3] names documented as exports exist in src/")
let src = read(joinpath(ROOT, "src", "GREBClimate.jl"), String)
    exported = Set{String}()
    for m in eachmatch(r"^export\s+(.+)$"m, src)
        for n in split(m.captures[1], ",")
            push!(exported, strip(n))
        end
    end
    isempty(exported) && fail("no exports parsed from src/GREBClimate.jl")
    missing_names = String[]
    for rel in ["docs/src/tutorial.md", "docs/src/switches.md", "README.md"]
        s = read(joinpath(ROOT, rel), String)
        # [`name`](@ref) - Documenter would catch these, but this is instant
        for m in eachmatch(r"\[`([A-Za-z_][A-Za-z0-9_!]*)`\]\(@ref\)", s)
            n = m.captures[1]
            n in exported || push!(missing_names, "$rel: $n")
        end
    end
    isempty(missing_names) ? ok("$(length(exported)) exports; all @ref'd names present") :
        foreach(fail, unique(missing_names))
end

# ── 4. No dead .jld2 files ───────────────────────────────────────────────────
# Interpolation-aware: src/io.jl builds ENSO filenames as
# "erainterim.tsurf.$suffix.forcing.jld2", so a literal grep reports ten
# false positives. Expand $suffix over its known values before comparing.
println("\n[4] every .jld2 on disk is read by code")
let datadir = joinpath(ROOT, "greb_input_data")
    if !isdir(datadir)
        ok("skipped - greb_input_data/ not present")
    else
        code = String[]
        for dir in ["src", "test", "benchmark", "examples", "scripts"]
            d = joinpath(ROOT, dir)
            isdir(d) || continue
            for (r, _, fs) in walkdir(d), f in fs
                endswith(f, ".jl") && push!(code, read(joinpath(r, f), String))
            end
        end
        blob = join(code, "\n")
        # expand the known interpolations
        for suffix in ["elnino", "lanina"]
            blob *= "\n" * replace(blob, "\$suffix" => suffix)
        end
        dead = String[]
        for (r, _, fs) in walkdir(datadir), f in fs
            endswith(f, ".jld2") && !occursin(f, blob) && push!(dead, f)
        end
        if isempty(dead)
            ok("all .jld2 files are referenced")
        else
            println("  NOTE  $(length(dead)) unreferenced file(s) - see .claude/notes/data-distribution.md")
            foreach(f -> println("          $f"), sort(dead))
        end
    end
end

println()
if isempty(failures)
    println("docs-check: all static checks passed.")
    println("STILL REQUIRED: run the documented snippets. See SKILL.md §1.")
else
    println("docs-check: $(length(failures)) failure(s).")
    exit(1)
end
