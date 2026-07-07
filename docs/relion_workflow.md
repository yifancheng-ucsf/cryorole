# RELION STAR Workflow

This guide shows the recommended cryoROLE 2.0 workflow for RELION STAR
metadata. It focuses on the public commands and the default matching behavior.

## Inputs

cryoROLE compares two independently refined domains from the same particle set.
For RELION inputs, each domain is provided as a STAR metadata file:

```text
ref_domain.star
mov_domain.star
```

RELION `Rot/Tilt/Psi` values are normalized by cryoROLE before relative
orientation computation. Source convention handling is part of the
normalization layer and should not be adjusted downstream.

## Default Matching

By default, `cryorole run` does not assume the two STAR files have the same row
order. It tries safe STAR identity matching in this order:

```text
_rlnTomoParticleName
_rlnImageName / rlnImageName
```

If matching reorders rows or drops unmatched particles, cryoROLE records that in
the run reports and manifest. If matching cannot be resolved safely, prepare
aligned inputs before running. Use `cryorole align --key ...` when the default
STAR identity keys are not sufficient and you need to provide an explicit
matching key.

## Standard Run

```bash
cryorole run --ref ref_domain.star --mov mov_domain.star
```

The default run bundle is written under:

```text
cryorole_outputs/
```

To choose a different run bundle, add `--output-dir RUN_DIR` to `cryorole run`
and use that same directory with downstream `--run-dir` commands.

Important outputs include:

```text
cryorole_outputs/run_manifest.json
cryorole_outputs/run_summary.json
cryorole_outputs/data/raw_landscape.npz
cryorole_outputs/data/raw_landscape.csv
cryorole_outputs/data/match_table.csv
cryorole_outputs/reports/
cryorole_outputs/visualizations/
```

## Row-Aligned Inputs

Use `--row-aligned` only when row `N` in both STAR files is known to represent
the same particle:

```bash
cryorole run \
  --ref aligned_ref.star \
  --mov aligned_mov.star \
  --row-aligned
```

In row-aligned mode, cryoROLE checks that the row counts match, pairs rows by
index, and records the row-aligned policy.

## Preparing Aligned STAR Files

When default matching is not enough, use `cryorole align` to prepare aligned
STAR files. The command can use safe automatic identity keys, or explicit keys
with `--key` when the automatic STAR keys are not enough:

```bash
cryorole align --ref ref_domain.star --mov mov_domain.star
```

Then run:

```bash
cryorole run \
  --ref alignments/default/aligned_ref.star \
  --mov alignments/default/aligned_mov.star \
  --row-aligned
```

`cryorole align` is a preparation step. It does not compute RO, SLD,
canonicalization, selection, or export.

## Canonicalize

```bash
cryorole canonicalize --run-dir cryorole_outputs
```

Canonicalization writes a separate canonical landscape under:

```text
cryorole_outputs/canonical/default/
```

It does not overwrite the raw landscape.

## Visualize

```bash
cryorole visualize --run-dir cryorole_outputs --space canonical
```

Visualization outputs are display-only. Display filters, ranges, and color
scales do not alter raw/canonical landscapes and do not create selections.

## Select

```bash
cryorole select \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --space canonical \
  -c 13 0 14 \
  -r 6
```

The default radius selection uses SO(3) geodesic distance when the center is
provided in Euler degrees.
The center values are usually chosen after inspecting the raw or canonical
visualizations.

## Export

```bash
cryorole export \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --domain both
```

Export uses the recorded `ref_source_row_id` and `mov_source_row_id`
provenance. It subsets the original source metadata and does not rewrite source
poses with canonical or display coordinates.
