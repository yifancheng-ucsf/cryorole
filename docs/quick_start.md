# cryoROLE 2.0 Quick Start

This guide shows the smallest practical cryoROLE 2.0 workflow. For the stable
artifact and policy contract, see `docs/architecture.md`.

## Install

From a local checkout:

```bash
python -m pip install -e ".[test]"
cryorole --help
```

The installed command is `cryorole`.

## RELION STAR Workflow

Run the relative-orientation analysis:

```bash
cryorole run --ref ref_domain.star --mov mov_domain.star
```

By default, the run bundle is written under `cryorole_outputs/`. The run bundle
contains raw numeric arrays, CSV tables, reports, manifests, and default
quick-look visualizations.

Canonicalize the raw landscape:

```bash
cryorole canonicalize --run-dir cryorole_outputs
```

Visualize the canonical landscape:

```bash
cryorole visualize --run-dir cryorole_outputs --space canonical
```

Select particles around a canonical Euler center:

```bash
cryorole select \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --space canonical \
  -c 13 0 14 \
  -r 6
```

Export the selected source metadata:

```bash
cryorole export \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --domain both
```

## CryoSPARC CS Workflow

CryoSPARC `.cs` files are native inputs in cryoROLE 2.0:

```bash
cryorole run --ref ref_domain.cs --mov mov_domain.cs
cryorole canonicalize --run-dir cryorole_outputs
cryorole visualize --run-dir cryorole_outputs --space canonical
cryorole select --run-dir cryorole_outputs --selection-id state_1 --space canonical -c 13 0 14 -r 6
cryorole export --run-dir cryorole_outputs --selection-id state_1 --domain both
```

Default `.cs` matching uses `uid`.

## Choosing a Run Directory

The compact public `run` command writes to the default run bundle location
`cryorole_outputs/`. Downstream commands use that directory through
`--run-dir`.

If you have multiple runs, keep each run bundle in a separate directory and use
the matching `--run-dir` for canonicalization, visualization, selection, and
export.

## Matching and Row Alignment

Use default matching when possible:

- CryoSPARC `.cs`: `uid`
- RELION STAR: `_rlnTomoParticleName`, then `_rlnImageName` / `rlnImageName`
  when safe

Use `--row-aligned` only when you know the files are already particle-aligned:

```bash
cryorole run --ref aligned_ref.star --mov aligned_mov.star --row-aligned
```

In row-aligned mode, cryoROLE pairs row `N` with row `N` and requires equal row
counts.

If default matching is not safe, prepare aligned metadata first:

```bash
cryorole align --ref ref_domain.star --mov mov_domain.star
cryorole run \
  --ref alignments/default/aligned_ref.star \
  --mov alignments/default/aligned_mov.star \
  --row-aligned
```

Manual pre-alignment is also acceptable when the row-order assertion is true and
auditable.

## What the Main Artifacts Mean

- Raw landscape: direct relative-orientation result from matched input
  particles.
- Canonical landscape: a derived coordinate frame for easier comparison and
  inspection.
- Visualization: display-only figures and display tables.
- Selection: an explicit scientific decision artifact.
- Export: a source-metadata subset created from a selection.

Display filtering, such as plotting only high-density points, is not a
scientific selection. Use `cryorole select` when you intend to create a particle
set for export or reconstruction.

## Complete Workflow Template

Copy and edit the paths:

```bash
cryorole run --ref path/to/ref.star --mov path/to/mov.star

cryorole canonicalize --run-dir cryorole_outputs

cryorole visualize \
  --run-dir cryorole_outputs \
  --space canonical \
  --visual-id overview

cryorole select \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --space canonical \
  -c 13 0 14 \
  -r 6

cryorole export \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --domain both
```
