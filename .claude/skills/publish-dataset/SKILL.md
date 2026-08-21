---
name: publish-dataset
description: Build and publish the GREB JLD2 input dataset as a GitHub release asset, and reconcile DATA_SHA256 in src/data.jl. Use after regenerating greb_input_data/ from raw .bin files, when the DataDep download 404s or fails checksum verification, or when asked to cut/update a dataset release.
---

# Publish the dataset

The `.jld2` dataset (39 files, 439 MB) is distributed as a single `.tar.gz`
attached to a GitHub release and fetched by `greb_data_dir()` via DataDeps. Three
things must agree or users get a 404 or a checksum failure that only shows up on
a machine with no local copy:

| Thing | Where | Must match |
|:------|:------|:-----------|
| release tag | `DATA_RELEASE_TAG` in `src/data.jl` | the git release tag |
| asset filename | `DATA_ARCHIVE_NAME` in `src/data.jl` | the uploaded file's name, exactly |
| archive hash | `DATA_SHA256` in `src/data.jl` | the built archive's SHA256 |

`DATA_URL` is assembled from the first two, so it needs no separate edit.

## Steps

### 1. Build

```bash
julia --project=. tools/package_dataset.jl greb_input_data greb_input_data-v1.tar.gz
```

It validates the tree against the converter's `MODEL_FIELD_NAMES` allowlist
**before** building — so a stray field cannot silently add ~13 MB to every
user's download — then prints the file count, size and SHA256.

Expect `39 .jld2 files, 438.5 MB unpacked` and a 353 MB archive. A different
file count means the tree is wrong; fix that before publishing, not after.

### 2. Reconcile the hash

Compare the printed SHA256 with `DATA_SHA256` in `src/data.jl`.

- **Same** → the dataset is unchanged. Nothing to update.
- **Different** → the dataset contents changed. Update `DATA_SHA256` and the
  byte count in its docstring.

**A changed hash after a plain regeneration is not expected.** The archive is
built with `--sort=name --owner=0 --group=0 --numeric-owner --mtime=<fixed>`
and `gzip -n`, so it depends only on contents. If you regenerated from identical
`Data/` and the hash moved, something really did change — investigate rather
than pasting the new value.

> This is worth knowing because it bit us. The first version of the script
> omitted `--mtime`, so tar embedded each file's modification time and the hash
> depended on *when* the dataset was written. A regeneration produced a
> byte-identical tree and a 17-byte-different archive, which looked exactly like
> data corruption. Never remove the `--mtime` pin, and never change its value —
> published checksums stop reproducing if you do.

### 3. Publish

```bash
gh release create data-v1 greb_input_data-v1.tar.gz \
  --repo EnvDroneSense/GREBClimate.jl \
  --title "GREB input dataset v1" \
  --notes "JLD2 input dataset for GREBClimate.jl (39 files, 439 MB unpacked). Fetched automatically by greb_data_dir() via DataDeps."
```

Updating an existing release instead:

```bash
gh release upload data-v1 greb_input_data-v1.tar.gz --clobber --repo EnvDroneSense/GREBClimate.jl
```

**The filename must equal `DATA_ARCHIVE_NAME`.** `gh` uses the local basename as
the asset name, so building to a scratch path with a different name and
uploading that is the easy way to publish a working release with a dead URL.

If the dataset changes incompatibly rather than being regenerated, cut a *new*
tag (`data-v2`) and bump both `DATA_RELEASE_TAG` and `DATA_ARCHIVE_NAME` — don't
silently replace an asset older installs may have cached.

### 4. Verify

The static checks confirm the constants are self-consistent, not that the URL
resolves:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'   # includes the constants test
curl -sIL -o /dev/null -w '%{http_code}\n' \
  "$(julia --project=. -e 'using GREBClimate; print(GREBClimate.DATA_URL)')"
```

`200` means the asset is reachable. A real end-to-end check needs the DataDep
branch to actually run, which means no local dataset and no `GREB_DATA`:

```bash
# on a machine or container without greb_input_data/
DATADEPS_ALWAYS_ACCEPT=true julia --project=. -e 'using GREBClimate; println(greb_data_dir())'
```

**Do not simulate this by moving `greb_input_data/` aside.** A rename-and-restore
destroyed the directory once — 439 MB of files gone. Use a clean checkout, a
container, or accept the URL check plus the constants test.

## Failure modes and what they mean

| Symptom | Cause |
|:--------|:------|
| DataDeps reports a 404 | release or asset missing, or asset name ≠ `DATA_ARCHIVE_NAME` |
| DataDeps reports a hash mismatch | `DATA_SHA256` not updated, or the asset was rebuilt without `--mtime` |
| Process hangs with no output | consent prompt waiting on stdin — set `DATADEPS_ALWAYS_ACCEPT=true` |
| `package_dataset.jl` errors about extra files | the tree was built with `--all`; regenerate without it |
| `package_dataset.jl` errors about missing fields | conversion was incomplete, or `MODEL_FIELD_NAMES` gained an entry the converter cannot produce |

## Related

- `tools/convert_greb_to_jld2.jl` — builds `greb_input_data/` from raw `.bin`.
  Pass `--all` only to inspect extra fields, never to build a release tree.
- `.claude/notes/datadeps-context.md` — measured sizes, why `.tar.gz` and not
  `.tar.xz`, and the split-bundle option that was deliberately deferred.
