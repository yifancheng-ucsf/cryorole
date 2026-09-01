# cryoROLE

cryoROLE (cryo-EM Relative Orientation LandscapE) quantifies inter-domain rotational motion in single-particle cryo-EM by computing the per-particle relative orientation between two independently refined near-rigid domains.

It is designed for cryo-EM projects in which two domains, modules, or subcomplexes can be treated as near-rigid bodies and refined separately from the same particle set. For each matched particle, cryoROLE computes the relative orientation (RO) between a reference domain and a moving domain, represents the particle population as a relative-orientation landscape, and enables landscape-guided particle selection for downstream RELION or CryoSPARC reconstruction.

## What cryoROLE does

Given two per-domain pose metadata files, cryoROLE can:

1. match particles between the two domain refinements;
2. normalize RELION or CryoSPARC pose conventions into a common internal representation;
3. compute per-particle relative orientations,

   ```text
   RO = R_ref^-1 R_mov
   ```

4. estimate local sampling density in rotation-vector (RV) space;
5. visualize raw or motion-aligned RO landscapes;
6. select particles from defined regions of the landscape; and
7. export the selected source metadata for reconstruction or further refinement.

cryoROLE does not replace 3D classification, 3D variability analysis, cryoDRGN, 3DFlex, or other continuous-heterogeneity methods. It is a complementary, physically interpretable SO(3)-based analysis for cases where domain-specific poses are already available.

## Status

This repository contains the cryoROLE 2.0 command-line workflow. The package metadata version is currently `2.0.0a1`.

The command-line interface is the supported public interface. Core
Productionization and the first Workflow UX milestone are implemented:
side-effect-free preflight, artifact-derived status/next guidance, an offline
interactive explorer, and an optional conservative guide now sit on top of the
stable run-bundle, selection, and export contracts.

## Installation

We recommend installing cryoROLE in a dedicated conda environment, then installing the package from a GitHub checkout.

```bash
git clone https://github.com/yifancheng-ucsf/cryorole.git
cd cryorole

conda create -n cryorole python=3.10 -y
conda activate cryorole

python -m pip install .
cryorole --help
```

cryoROLE requires Python 3.9 or newer. Core Python dependencies are declared in
`pyproject.toml` and include `numpy`, `scipy`, `pandas`, `matplotlib`, and
`Pillow`.

For an editable development/test installation:

```bash
python -m pip install -e ".[test]"
python -m pytest
```

The repository also provides `environment.yml`, which creates an editable conda environment from the repository root:

```bash
conda env create -f environment.yml
conda activate cryorole
cryorole --help
```

See [Installation](docs/installation.md) for additional setup and verification notes.

## Input requirements

cryoROLE starts from two per-domain pose metadata files generated from the same particle set:

- one metadata file for the reference domain;
- one metadata file for the moving domain;
- RELION STAR (`.star`) or CryoSPARC (`.cs`) input;
- enough particle identity information to match the two pose tables.

The two domains should be approximately rigid over the range of motion being analyzed. If a domain undergoes major internal deformation, a single per-particle orientation may no longer be a meaningful descriptor for that domain.

### Particle matching

By default, cryoROLE matches particles by identity rather than assuming row order.

For CryoSPARC `.cs` inputs, the default identity key is `uid`.

For RELION STAR inputs, cryoROLE first tries `_rlnTomoParticleName` when available, then `_rlnImageName` / `rlnImageName` when the resolved column is unambiguous.

Use `--row-aligned` only when row `N` in the reference metadata is known to describe the same particle as row `N` in the moving-domain metadata. If safe automatic matching is not possible for STAR files, use `cryorole align --key ...` to prepare row-aligned inputs, or manually pre-align the metadata before running with `--row-aligned`.

## Workflow

The public workflow is:

```text
preflight -> run -> status/next -> [canonicalize] -> explore/visualize -> confirm/select -> export
```

| Step | Role |
|---|---|
| `preflight` / `run --dry-run` | Validates inputs, matching, pose schema, and resources without creating a run bundle. |
| `align` | Optional STAR pre-processing step that prepares row-aligned metadata when default matching is insufficient. |
| `run` | Computes the raw RO landscape and writes the run bundle. |
| `status` / `next` | Reads actual artifacts and proposes exact required, recommended, and optional commands. |
| `canonicalize` | Optionally derives a motion-aligned coordinate frame for interpretation and comparison. It does not overwrite the raw landscape. |
| `visualize` | Creates display-only plots and display tables from raw, canonical, or selected rows. |
| `explore` | Opens a localhost-only offline explorer; clicking/evaluation makes a draft, and Confirm writes a standard selection. |
| `select` | Creates an explicit, auditable particle-selection artifact. |
| `export` | Subsets the original source metadata using the recorded selection and source-row provenance. |
| `guide` | Conservatively combines preflight/status/next without silently making scientific decisions. |

A key design rule is that visualization filters are display-only. They do not create scientific selections and do not modify the raw or canonical landscape. To generate particles for downstream reconstruction, use `cryorole select` followed by `cryorole export`.

Each command writes provenance into the run bundle so selections and exports can be audited later.

## Quick start

The examples below use the default run directory, `cryorole_outputs/`, and assume no custom run directory was requested. To choose a different run bundle name, add `--output-dir RUN_DIR` to `cryorole run`, then pass that same directory to downstream commands with `--run-dir RUN_DIR`.

### RELION STAR inputs

If the reference and moving-domain STAR files contain safe particle identity columns, run:

```bash
cryorole preflight --ref ref_domain.star --mov mov_domain.star
cryorole run --ref ref_domain.star --mov mov_domain.star
```

If the STAR files cannot be matched safely by default, prepare aligned STAR files first:

```bash
cryorole align --ref ref_domain.star --mov mov_domain.star

cryorole run \
  --ref alignments/default/aligned_ref.star \
  --mov alignments/default/aligned_mov.star \
  --row-aligned
```

Use direct row-aligned mode only when you already know both files are in the same particle order:

```bash
cryorole run --ref aligned_ref.star --mov aligned_mov.star --row-aligned
```

### CryoSPARC `.cs` inputs

For CryoSPARC input, cryoROLE reads `.cs` files directly and matches particles by `uid`:

```bash
cryorole run --ref ref_domain.cs --mov mov_domain.cs --dry-run
cryorole run --ref ref_domain.cs --mov mov_domain.cs
```

`preflight` and `run --dry-run` use the same validation service as the real
run, create no run directory, and report `READY`, `READY_WITH_WARNINGS`, or
`BLOCKED`. Warnings return exit code 1 and blocked inputs return 2.

### Inspect the raw landscape

`cryorole run` writes seven flat quick-look PNGs: Euler/RV triptychs for all
particles, `sld_raw >= 1`, and the top 40%, plus a full-data log-SLD
distribution. `run_report.md` explains the outputs and warnings. Generate
an independent compact view with the explicit visualization command. Its
default is two PNG triptychs (Euler and RV) for `sld_display >= 1`:

```bash
cryorole visualize --run-dir cryorole_outputs --space raw
cryorole visualize --run-dir cryorole_outputs --space raw --view 2d,1d
cryorole visualize --run-dir cryorole_outputs --space raw --view 3d
```

The 1D option uses every row that passes the display filters. The 3D default is
a self-contained offline viewer; it is exploratory and cannot create a
Selection.

For local linked Euler/RV views and exact radius-selection previews:

```bash
cryorole explore --run-dir cryorole_outputs --space raw
```

The explorer binds only to `127.0.0.1`, uses packaged assets with no CDN, and
may downsample points only for display. Exact counts and Confirm always evaluate
the full parent landscape in Python using the same SO(3) evaluator as
`cryorole select`. No Selection is written until Confirm.

### Inspect a landscape interactively in ChimeraX

The standalone `cryorole_chimerax_viewer.py` script registers display-only
ChimeraX commands for cryoROLE landscape CSV files. Load it once in ChimeraX,
then open a raw or canonical landscape CSV:

```text
open /path/to/cryorole_chimerax_viewer.py
cryorole open /path/to/raw_landscape.csv
```

The viewer can switch between Euler and rotation-vector coordinates, color by
SLD, and apply threshold or top-fraction display filters. It does not perform
RO analysis, create scientific selections, or modify source metadata.

### Canonicalize the landscape

Canonicalization is optional. It re-expresses the landscape in a motion-aligned coordinate frame so that the dominant motion is easier to view and compare. It does not change the raw RO facts.

```bash
cryorole canonicalize --run-dir cryorole_outputs
cryorole visualize --run-dir cryorole_outputs --space canonical
```

### Test a moving-domain rotation

The standalone diagnostic script applies one extrinsic fixed-axis ZYX rotation
to every RO by right multiplication and inherits the parent SLD values:

```bash
python scripts/rotate_landscape.py \
  --input cryorole_outputs \
  --space raw \
  --rotation-euler 20 10 0 \
  --output-dir rotated_landscape

cryorole visualize --run-dir rotated_landscape --space raw
```

Use `--space canonical --canonical-id ID` to define the input rotation in an
existing canonical frame. The script writes a derived landscape bundle only;
it does not modify or export source STAR/CS poses. SLD is inherited rather than
recomputed so the same particles keep the same visualization colors.

### Export an offline animation bundle

`cryorole animate` implements SO(3) trajectory generation, synchronized
landscape PNG frames, and ChimeraX script export. The default `script-only`
mode does not require ChimeraX:

```bash
cryorole animate \
  --run-dir cryorole_outputs \
  --coordinate-set raw \
  --path-csv waypoints.csv \
  --path-space rv \
  --chimerax-session prepared_scene.cxs \
  --reference-model-id "#1" \
  --moving-model-id "#2" \
  --pivot 0 0 0 \
  --baseline-ro identity \
  --map-frame raw \
  --output-dir animation_output
```

Add `--render-mode execute --chimerax-bin PATH --no-encode` to request
validated structure and composite frames. Omit `--no-encode` and provide
`--ffmpeg-bin PATH --ffprobe-bin PATH` for validated H.264 MP4 output.
Animation display filters and legacy style are shared with
`cryorole visualize`; Linux execution uses ChimeraX offscreen mode and requires
an explicit renderer-completion status plus valid frames. Composition preserves
the full input frames and their aspect ratios; no frame set is automatically
deleted. Animation filters are display-only and never create selections.
Repeat either model-ID option to define disjoint rigid groups; all movers use
the same absolute trajectory delta from their own saved scene transforms.

### Select particles

For a radius selection around a canonical Euler center:

```bash
cryorole select \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --space canonical \
  -c 13 0 14 \
  -r 6
```

By default, radius selection uses SO(3) geodesic distance, not simple Euclidean distance in Euler-angle space.
The center values are usually chosen after inspecting the raw or canonical visualizations.

### Export selected metadata

```bash
cryorole export \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --domain both
```

Export subsets the original source metadata using recorded source-row provenance. It does not rewrite source poses with canonical or display coordinates.

## Common commands

Use `cryorole COMMAND --help` for the exact current options.

```bash
cryorole preflight --ref REF_METADATA --mov MOV_METADATA [--json REPORT.json]
cryorole run --ref REF_METADATA --mov MOV_METADATA --dry-run
cryorole align --ref REF.star --mov MOV.star
cryorole run --ref REF_METADATA --mov MOV_METADATA [--output-dir RUN_DIR]
cryorole status --run-dir RUN_DIR
cryorole next --run-dir RUN_DIR
cryorole guide --run-dir RUN_DIR --non-interactive
cryorole canonicalize --run-dir cryorole_outputs
cryorole visualize --run-dir cryorole_outputs --space canonical
cryorole explore --run-dir cryorole_outputs --space canonical --canonical-id default
cryorole animate --run-dir cryorole_outputs --path-csv PATH.csv --path-space rv --chimerax-session SCENE.cxs --reference-model-id "#1" --moving-model-id "#2" --pivot 0 0 0 --baseline-ro identity --map-frame raw --output-dir ANIMATION_DIR
cryorole select --run-dir cryorole_outputs --selection-id state_1 --space canonical -c A B C -r DEG
cryorole export --run-dir cryorole_outputs --selection-id state_1 --domain both
```

## Outputs

A standard cryoROLE run bundle records numeric arrays, flat tables, reports, visualizations, selections, exports, and provenance:

```text
cryorole_outputs/
  run_manifest.json
  run_summary.json
  run_report.md
  data/
    raw_landscape.npz
    raw_landscape.csv
    match_table.csv
  reports/
  visualizations/
  canonical/
  selections/
  exports/
```

If `cryorole run --output-dir my_run` is used, the same layout is written under `my_run/` instead.

Important user-facing artifacts include:

| Artifact | Meaning |
|---|---|
| `run_manifest.json` | Top-level provenance and artifact index for the run bundle. |
| `run_summary.json` | Human-readable summary of inputs, matching, policies, and outputs. |
| `run_report.md` | Concise guide to files, quick-look figures, warnings, and next commands. |
| `data/raw_landscape.npz` | Machine-readable source of truth for the raw RO landscape. |
| `data/raw_landscape.csv` | Flat table for inspection, plotting, and external tools. |
| `data/match_table.csv` | Matched-particle provenance linking reference and moving-domain source rows. |
| `canonical/<id>/canonical_landscape.csv` | Canonical landscape table with raw and canonical coordinates. |
| `visualizations/` | Display-only figures and display tables. |
| `selections/<selection_id>/` | Auditable particle-selection artifact. |
| `exports/<selection_id>/` | Selected source metadata for downstream reconstruction. |

JSON files are used for reports, summaries, manifests, and provenance. Full object-record landscape JSON is not the production-scale persistence format.

## Key concepts

| Concept | Meaning |
|---|---|
| Relative orientation (RO) | Per-particle inter-domain rotation, defined as `RO = R_ref^-1 R_mov`. |
| Reference domain | Domain whose frame is used as the reference for the relative orientation. |
| Moving domain | Domain expressed relative to the reference domain. |
| Raw landscape | Direct RO result from matched input particles. |
| Canonical landscape | Optional motion-aligned representation derived from the raw landscape. |
| RV space | Rotation-vector coordinate space used for statistics, density estimation, canonicalization, and geodesic reasoning. |
| Euler space | User-facing display coordinate. cryoROLE 2.0 uses extrinsic fixed-axis ZYX for public RO/RV-derived Euler output and records that convention in reports/manifests. |
| `sld_raw` | Scientific SLD (kNN-scaled local density) value used by default policies. |
| `sld_display` | Display-only SLD/color value. |
| Visualization filter | Display-only row/filter choice; does not create a selection. |
| Selection | Explicit particle subset that can be exported. |
| Export | Source metadata subset; source poses are not rewritten with display or canonical coordinates. |

## Documentation

- [Installation](docs/installation.md)
- [Quick start](docs/quick_start.md)
- [Workflow UX tutorial](docs/workflow_ux.md)
- [CLI reference](docs/cli_reference.md)
- [FAQ and troubleshooting](docs/faq.md)
- [RELION workflow](docs/relion_workflow.md)
- [CryoSPARC workflow](docs/cryosparc_workflow.md)
- [Output files](docs/output_files.md)
- [Migration from cryoROLE 0.x](docs/migration_from_0x.md)
- [Architecture contract](docs/architecture.md)
- [Animation export](docs/animation_export.md)
- [Roadmap](docs/roadmap.md)

## Notes for cryoROLE 0.x users

The old multi-script workflow maps to the 2.0 CLI as follows:

| cryoROLE 0.x command | cryoROLE 2.0 command |
|---|---|
| `orientation_analysis` | `cryorole run` |
| `landscape_projection` | `cryorole visualize` |
| `point_select` | `cryorole select` |
| `particle_backtrack` | `cryorole export` |

See [Migration from cryoROLE 0.x](docs/migration_from_0x.md) for details.

## Citation

If you use cryoROLE as a method, please cite the cryoROLE methodology preprint:

- Chengmin Li, Wooyoung Choi, Hao Wu, Yifan Cheng. **CryoROLE: describing large inter-domain rotation in single particle cryo-EM.** bioRxiv, 2026. [https://www.biorxiv.org/content/10.64898/2026.07.04.736454v1](https://www.biorxiv.org/content/10.64898/2026.07.04.736454v1)

The first application of cryoROLE was in the human fatty acid synthase study:

- Wooyoung Choi, Chengmin Li, Yifei Chen, YongQiang Wang, Yifan Cheng. **Structural dynamics of human fatty acid synthase in the condensing cycle.** Nature, 2025. [https://doi.org/10.1038/s41586-025-08782-w](https://doi.org/10.1038/s41586-025-08782-w)

## License

cryoROLE is distributed under the BSD 3-Clause License.
