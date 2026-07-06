# cryoROLE

cryoROLE (cryo-EM Relative Orientation LandscapE) is a tool for describing
inter-domain rotational dynamics in single-particle cryo-EM. It is designed for 
complexes where two near-rigid domains, modules, or subcomplexes can be refined 
separately from the same particle set. For each matched particle, cryoROLE 
computes the relative orientation between the two domains and represents the 
particle population as a relative-orientation landscape.

The landscape can be used to inspect continuous domain motion, identify densely
sampled orientation states, select particles from a region of the landscape, and
export the selected source metadata for downstream RELION or CryoSPARC
reconstruction.

## Citation

If you use cryoROLE as a method, please cite the cryoROLE methodology preprint:

- Chengmin Li, Wooyoung Choi, Hao Wu, Yifan Cheng. **CryoROLE: describing large
  inter-domain rotation in single particle cryo-EM.** bioRxiv, 2026.
  [https://www.biorxiv.org/content/10.64898/2026.07.04.736454v1](https://www.biorxiv.org/content/10.64898/2026.07.04.736454v1)

The first application of cryoROLE was in the human fatty acid synthase study:

- Wooyoung Choi, Chengmin Li, Yifei Chen, YongQiang Wang, Yifan Cheng.
  **Structural dynamics of human fatty acid synthase in the condensing cycle.**
  Nature, 2025. [https://doi.org/10.1038/s41586-025-08782-w](https://doi.org/10.1038/s41586-025-08782-w)

## Status

This repository currently contains the cryoROLE 2.0.0a1 beta-preview workflow.
The command-line workflow is the supported public interface. MPI-enabled
large-scale execution and a GUI are in development.

## Installation

From a GitHub checkout:

```bash
git clone https://github.com/yifancheng-ucsf/cryorole.git
cd cryorole
python -m pip install .
cryorole --help
```


See [Installation](docs/installation.md) for conda and verification notes.

## Workflow

The public workflow is:

```text
align (optional) -> run RO compute -> canonicalize (optional) -> visualize -> select -> export
```

- `align` is only needed when STAR metadata cannot be safely matched by the
  default identity columns and you want to prepare row-aligned STAR files.
- `run` computes the raw relative-orientation (RO) landscape and writes a run
  bundle.
- `canonicalize` is optional. It derives a reusable coordinate frame that makes
  the landscape easier to inspect and compare; it does not overwrite the raw
  landscape.
- `visualize` makes display-only plots.
- `select` creates an explicit scientific particle-selection artifact.
- `export` subsets the original source metadata using the recorded selection.

## Quick Start

The examples below use the default run directory, `cryorole_outputs/`.

### 1. Decide whether your inputs are already aligned

If row `N` in both metadata files is known to be the same particle, you can use
row-aligned mode:

```bash
cryorole run --ref aligned_ref.star --mov aligned_mov.star --row-aligned
```

If you are not sure whether the files are already aligned, do not use
`--row-aligned`. Let cryoROLE resolve particle identity and report whether rows
were reordered or dropped:

```bash
cryorole run --ref REF.star --mov MOV.star
```

For CryoSPARC `.cs` inputs, default matching uses `uid`:

```bash
cryorole run --ref REF.cs --mov MOV.cs
```

For RELION STAR inputs, cryoROLE first tries `_rlnTomoParticleName`, then
`_rlnImageName` / `rlnImageName` when safe. If default matching is not possible,
prepare aligned STAR files:

```bash
cryorole align --ref REF.star --mov MOV.star
cryorole run \
  --ref alignments/default/aligned_ref.star \
  --mov alignments/default/aligned_mov.star \
  --row-aligned
```

### 2. Inspect the raw landscape

`cryorole run` writes default quick-look visualizations. You can also make
explicit raw visualizations:

```bash
cryorole visualize --run-dir cryorole_outputs --space raw
```

### 3. Optionally canonicalize

Canonicalization is optional. It rotates/re-expresses the landscape into a
consistent coordinate frame so that the main motion is easier to view and
compare. It is a derived coordinate frame and does not change the raw RO
landscape.

```bash
cryorole canonicalize --run-dir cryorole_outputs
```

Then visualize the canonical landscape:

```bash
cryorole visualize --run-dir cryorole_outputs --space canonical
```

### 4. Select particles

For a radius selection around a canonical Euler center:

```bash
cryorole select \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --space canonical \
  -c 13 0 14 \
  -r 6
```

By default, radius selection uses SO(3) geodesic distance, not simple Euclidean
distance in Euler-angle space.

### 5. Export selected metadata

```bash
cryorole export \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --domain both
```

Export subsets the original source metadata using recorded source-row
provenance. It does not rewrite source poses with canonical or display
coordinates.

## CLI Summary

Use `cryorole COMMAND --help` for the exact current options.

### `cryorole align`

Prepare RELION STAR files for `cryorole run --row-aligned`.

Common options:

| Option | Meaning |
|---|---|
| `--ref REF` | Reference STAR metadata. |
| `--mov MOV` | Moving-domain STAR metadata. |
| `--align-id ID` | Output id under `alignments/`; default `default`. |
| `--key COL [COL ...]` | Explicit STAR key columns. |
| `--float-tol COL=TOL` | Tolerance for numeric key columns; may repeat. |
| `--path-mode exact|basename|suffix:N` | Normalize path-like keys. |
| `--duplicate-policy exclude|first` | Duplicate key handling. |
| `--overwrite` | Replace artifacts for the same align id. |

### `cryorole run`

Compute the raw relative-orientation landscape.

Common options:

| Option | Meaning |
|---|---|
| `--ref REF` | Reference-domain pose metadata, STAR or CS. |
| `--mov MOV` | Moving-domain pose metadata, STAR or CS. |
| `--ref-domain NAME` | Provenance label for the reference domain. |
| `--mov-domain NAME` | Provenance label for the moving domain. |
| `--row-aligned` | Assert row `N` in ref and mov is the same particle. |
| `--no-visualize` | Skip default raw quick-look plots. |

### `cryorole canonicalize`

Optionally derive a canonical coordinate frame from a run bundle.

Common options:

| Option | Meaning |
|---|---|
| `--run-dir RUN` | Existing cryoROLE run bundle. |
| `--canonical-id ID` | Output id under `RUN/canonical/`; default `default`. |
| `--fit-top F` | Highest-`sld_raw` fraction used for frame fitting. |
| `--positive-side low|high` | Axis sign direction policy. |
| `--use-frame FRAME` | Apply an existing canonical frame. |
| `--no-visualize` | Skip canonical quick-look plots. |

### `cryorole visualize`

Render display-only views from raw, canonical, or selected landscapes.

Common options:

| Option | Meaning |
|---|---|
| `--run-dir RUN` | Existing run bundle. |
| `--space raw|canonical` | Landscape space to visualize. |
| `--selection-id ID` | Visualize rows from a selection. |
| `--visual-id ID` | Output id under `visualizations/`. |
| `--representation euler|rotvec|both` | Coordinate representation to plot. |
| `--top-fraction F` | Display-only high-density filter. |
| `--threshold VALUE` | Display-only SLD threshold. |
| `--range AXIS:LOWER:UPPER` | Display-only coordinate range. |
| `--formats png,svg,pdf` | Figure format list. |

### `cryorole select`

Create an explicit particle-selection artifact.

Common options:

| Option | Meaning |
|---|---|
| `--run-dir RUN` | Existing run bundle. |
| `--selection-id ID` | Selection id under `RUN/selections/`. |
| `--space raw|canonical` | Parent landscape space. |
| `--mode radius|threshold|range|random|metadata` | Selection mode. |
| `-c A B C`, `--center A B C` | Radius center; Euler degrees by default. |
| `-r DEG`, `--radius DEG` | Radius in degrees. |
| `--radius-rad RAD` | Radius in radians. |
| `--metric so3|rotvec` | Radius metric; default `so3`. |
| `--sld-min`, `--sld-max` | `sld_raw` threshold bounds. |
| `--metadata-domain ref|mov` | Metadata domain for metadata selection. |
| `--metadata-column COLUMN` | Source metadata column. |
| `--metadata-value VALUE[,VALUE...]` | Metadata values to include. |
| `--write-selected-landscape` | Also write a selected-derived landscape. |

### `cryorole export`

Export selected source metadata.

Common options:

| Option | Meaning |
|---|---|
| `--run-dir RUN` | Existing run bundle. |
| `--selection-id ID` | Selection id to export. |
| `--domain ref|mov|both` | Source metadata domain to export. |
| `--format auto|relion_star|cryosparc_cs|keys` | Export format. |
| `--output-dir DIR` | Advanced output override. |
| `--overwrite` | Replace export artifacts in the output directory. |

## Outputs

A cryoROLE run bundle records numeric data, tables, reports, visualizations,
selections, exports, and provenance:

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

`data/raw_landscape.npz` is the machine-readable source of truth for the raw
landscape. `data/raw_landscape.csv` is a flat table for inspection and external
tools. JSON files are reports, summaries, and manifests.

## Key Concepts

- Raw landscape: direct relative-orientation result from matched input
  particles.
- Canonical landscape: optional derived coordinate frame for inspection and
  comparison.
- Visualization filters: display-only choices; they do not create scientific
  selections.
- Selections: explicit decision artifacts that can be exported.
- `sld_raw`: the scientific SLD density field used by default policies.
- `sld_display`: a display-only density/color field.
- Export uses recorded source-row provenance and does not rewrite source poses.

## More Documentation

- [Installation](docs/installation.md)
- [Quick start](docs/quick_start.md)
- [RELION workflow](docs/relion_workflow.md)
- [CryoSPARC workflow](docs/cryosparc_workflow.md)
- [Output files](docs/output_files.md)
- [Migration from cryoROLE 0.x](docs/migration_from_0x.md)
- [Architecture contract](docs/architecture.md)
- [Roadmap](docs/roadmap.md)

## Notes for cryoROLE 0.x Users

The old multi-script workflow maps to the 2.0 CLI like this:

| cryoROLE 0.x command | cryoROLE 2.0 command |
|---|---|
| `orientation_analysis` | `cryorole run` |
| `landscape_projection` | `cryorole visualize` |
| `point_select` | `cryorole select` |
| `particle_backtrack` | `cryorole export` |

See [Migration from cryoROLE 0.x](docs/migration_from_0x.md) for details.

## License

cryoROLE is distributed under the BSD 3-Clause License.
