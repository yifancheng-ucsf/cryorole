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

Run:

```bash
cryorole run --ref ref_domain.cs --mov mov_domain.cs
```

Canonicalize:

```bash
cryorole canonicalize --run-dir cryorole_outputs
```

Visualize:

```bash
cryorole visualize --run-dir cryorole_outputs --space canonical
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
data/raw_landscape.npz
data/raw_landscape.csv
data/match_table.csv
reports/
visualizations/
canonical/
selections/
exports/
```

See `docs/output_files.md` for details.

## Export Notes

CryoSPARC export preserves the structured-array dtype and vector-valued fields
of the source metadata when writing selected subsets. Export does not reselect,
rematch, or write canonical/display coordinates back as physical source poses.

## Difference from 0.x Workflows

Older cryoROLE workflows often converted CryoSPARC `.cs` files to STAR before
analysis. cryoROLE 2.0 can read `.cs` inputs directly, so conversion is not the
default path.
