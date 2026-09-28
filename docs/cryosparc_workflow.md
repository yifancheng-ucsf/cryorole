# CryoSPARC CS Workflow

This guide shows the recommended cryoROLE 2.0 workflow for CryoSPARC `.cs`
metadata.

## Inputs

cryoROLE 2.0 treats `.cs` files as native inputs:

```text
ref_domain.cs
mov_domain.cs
```

The default CryoSPARC pose field is:

```text
alignments3D/pose
```

The pose is interpreted as a rotation vector / axis-angle representation and is
converted to cryoROLE internal active rotation matrices during normalization.

## Default Matching

CryoSPARC `.cs` inputs default to `uid` matching:

```bash
cryorole run --ref ref_domain.cs --mov mov_domain.cs
```

If matching reorders rows or drops unmatched particles, cryoROLE records that in
the run reports and manifest. Downstream selections and exports use the recorded
source-row provenance from the run bundle.

Public uid matching requires at least 50% overlap unless
`--allow-low-overlap` is explicit. Zero matches and duplicate uid values fail
before pose normalization.

## Row-Aligned Inputs

Use `--row-aligned` only when the two `.cs` files have already been prepared so
that row `N` in both files represents the same particle:

```bash
cryorole run \
  --ref aligned_ref.cs \
  --mov aligned_mov.cs \
  --row-aligned
```

In row-aligned mode, cryoROLE requires equal row counts and does not key-match,
reorder, or drop particles.

## Standard Workflow

Preflight and run:

```bash
cryorole preflight --ref ref_domain.cs --mov mov_domain.cs
cryorole run --ref ref_domain.cs --mov mov_domain.cs
```

Preflight validates `uid`, the `alignments3D/pose` vector shape and finite
values, matching coverage, and resources with no run-bundle side effects. The
equivalent check is `cryorole run ... --dry-run`.

Canonicalize:

```bash
cryorole canonicalize --run-dir cryorole_outputs
```

Visualize:

```bash
cryorole visualize --run-dir cryorole_outputs --space canonical
```

The default is two PNG triptychs for `sld_display >= 1`; use `--view 2d,1d`
or `--view 3d` only when those additional display products are needed.

Or explore linked Euler/RV projections locally and create a standard selection
only after explicit Confirm:

```bash
cryorole explore --run-dir cryorole_outputs --space canonical
```

Select:

```bash
cryorole select \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --space canonical \
  -c 13 0 14 \
  -r 6
```

The center values are usually chosen after inspecting the raw or canonical
visualizations.

Export:

```bash
cryorole export \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --domain both
```

## Select by CryoSPARC metadata

Use the run-time source field directly, without converting CS to STAR:

```bash
cryorole select --run-dir cryorole_outputs --selection-id ref_classes_01 --mode metadata --metadata-domain ref --metadata-column alignments3D/class --metadata-value 0,1
cryorole select --run-dir cryorole_outputs --selection-id ref_class --mode metadata --metadata-domain ref --metadata-column alignments3D/class --split-by-value
cryorole export --run-dir cryorole_outputs --selection-id ref_classes_01
```

The field must exist in the chosen source; class values are not renumbered.
Choose ref or mov explicitly. Selection uses the full parent landscape's
recorded source-row indices, including when run reordered or dropped unmatched
input rows. Source content must match its recorded SHA-256.

Supported CS fields are scalar integer (including uint64), boolean, and UTF-8
text. Text comparison is exact; empty strings are excluded and counted as
missing. Float/vector fields and invalid UTF-8 are rejected. Split writes one
standard selection per non-missing matched value, with at most 100 groups.
Name collisions fail before writing even with overwrite. See the
[CLI reference](cli_reference.md#cryosparc-metadata-fields) for value syntax.

## Outputs

The run bundle is written under:

```text
cryorole_outputs/
```

To choose a different run bundle, add `--output-dir RUN_DIR` to `cryorole run`
and use that same directory with downstream `--run-dir` commands.

Important outputs include:

```text
run_manifest.json
run_summary.json
run_report.md
data/raw_landscape.npz
data/raw_landscape.csv
data/match_table.csv
reports/
visualizations/quicklook/  # seven flat PNGs
canonical/
selections/
exports/
```

The run quick-look contains Euler/RV triptychs for all particles, `sld_raw >=
1`, and the top 40%, plus the full-data log-SLD distribution. Use `cryorole
visualize` for expanded display artifacts.

See `docs/output_files.md` for details.

`cryorole status --run-dir cryorole_outputs` derives progress and integrity from
the actual artifacts. `cryorole next --run-dir cryorole_outputs` prints exact
required, recommended, and optional commands; canonicalization remains optional.

## Export Notes

CryoSPARC export preserves the structured-array dtype and vector-valued fields
of the source metadata when writing selected subsets. Export does not reselect,
rematch, or write canonical/display coordinates back as physical source poses.

Export first verifies the source `.cs` SHA-256 recorded by `run`. A moved file
must be supplied explicitly with `--relocated-ref` / `--relocated-mov` and must
match byte-for-byte; a changed or replaced `.cs` file is rejected.

## Difference from 0.x Workflows

Older cryoROLE workflows often converted CryoSPARC `.cs` files to STAR before
analysis. cryoROLE 2.0 can read `.cs` inputs directly, so conversion is not the
default path.
