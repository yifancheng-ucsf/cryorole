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

Default key matching requires at least 50% overlap. Lower overlap requires the
explicit audited `--allow-low-overlap` override; zero overlap and duplicate
particle identities are hard failures before RELION Euler normalization.

## Standard Run

```bash
cryorole preflight --ref ref_domain.star --mov mov_domain.star
cryorole run --ref ref_domain.star --mov mov_domain.star
```

Preflight validates STAR pose columns, resolves the RELION convention, checks
identity uniqueness and overlap, estimates resources, and prints the exact run
command. It has no run-bundle side effects. `cryorole run ... --dry-run` uses
the same service. Review `READY_WITH_WARNINGS`; fix `BLOCKED` inputs before RO
computation.

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
cryorole_outputs/run_report.md
cryorole_outputs/data/raw_landscape.npz
cryorole_outputs/data/raw_landscape.csv
cryorole_outputs/data/match_table.csv
cryorole_outputs/reports/
cryorole_outputs/visualizations/quicklook/  # seven flat PNGs
```

The run quick-look contains Euler/RV triptychs for all particles, `sld_raw >=
1`, and the top 40%, plus the full-data log-SLD distribution. Use `cryorole
visualize` for expanded display products.

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
The default is two PNG triptychs for `sld_display >= 1`; 1D and 3D views are
explicit through `--view`.

## Select

To choose a center interactively from linked Euler and rotation-vector views:

```bash
cryorole explore --run-dir cryorole_outputs --space canonical
```

Display filters and display sampling affect only the browser view. Evaluate and
Confirm use the full parent landscape and the Python SO(3) evaluator. Confirm
writes a normal selection; closing the browser before Confirm writes none.

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

The source STAR identity is also content-addressed. If a STAR file has moved,
use `--relocated-ref` or `--relocated-mov`; export continues only when the
streamed SHA-256 matches the run record. Changed or replaced STAR files are
rejected.

At any time, `cryorole status --run-dir cryorole_outputs` reports actual bundle
and artifact integrity, while `cryorole next --run-dir cryorole_outputs` prints
exact required/recommended/optional next commands.
