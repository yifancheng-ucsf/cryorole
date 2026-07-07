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

The command-line interface is the supported public interface. The current development priority is hardening the 2.0 workflow for production-scale datasets while preserving stable run-bundle, selection, and export contracts.

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

cryoROLE requires Python 3.9 or newer. Core Python dependencies are declared in `pyproject.toml` and include `numpy`, `scipy`, `pandas`, and `matplotlib`.

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
prepare/match inputs -> run -> canonicalize -> visualize -> select -> export
```

| Step | Role |
|---|---|
| `align` | Optional STAR pre-processing step that prepares row-aligned metadata when default matching is insufficient. |
| `run` | Computes the raw RO landscape and writes the run bundle. |
| `canonicalize` | Optionally derives a motion-aligned coordinate frame for interpretation and comparison. It does not overwrite the raw landscape. |
| `visualize` | Creates display-only plots and display tables from raw, canonical, or selected rows. |
| `select` | Creates an explicit, auditable particle-selection artifact. |
| `export` | Subsets the original source metadata using the recorded selection and source-row provenance. |

A key design rule is that visualization filters are display-only. They do not create scientific selections and do not modify the raw or canonical landscape. To generate particles for downstream reconstruction, use `cryorole select` followed by `cryorole export`.

Each command writes provenance into the run bundle so selections and exports can be audited later.

## Quick start

The examples below use the default run directory, `cryorole_outputs/`, and assume no custom run directory was requested.

### RELION STAR inputs

If the reference and moving-domain STAR files contain safe particle identity columns, run:

```bash
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
cryorole run --ref ref_domain.cs --mov mov_domain.cs
```

### Inspect the raw landscape

`cryorole run` writes default quick-look visualizations. You can also generate explicit raw visualizations:

```bash
cryorole visualize --run-dir cryorole_outputs --space raw
```

### Canonicalize the landscape

Canonicalization is optional. It re-expresses the landscape in a motion-aligned coordinate frame so that the dominant motion is easier to view and compare. It does not change the raw RO facts.

```bash
cryorole canonicalize --run-dir cryorole_outputs
cryorole visualize --run-dir cryorole_outputs --space canonical
```

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
cryorole align --ref REF.star --mov MOV.star
cryorole run --ref REF_METADATA --mov MOV_METADATA
cryorole canonicalize --run-dir cryorole_outputs
cryorole visualize --run-dir cryorole_outputs --space canonical
cryorole select --run-dir cryorole_outputs --selection-id state_1 --space canonical -c A B C -r DEG
cryorole export --run-dir cryorole_outputs --selection-id state_1 --domain both
```

## Outputs

A standard cryoROLE run bundle records numeric arrays, flat tables, reports, visualizations, selections, exports, and provenance:

```text
cryorole_outputs/
  run_manifest.json
  run_summary.json
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

Important user-facing artifacts include:

| Artifact | Meaning |
|---|---|
| `run_manifest.json` | Top-level provenance and artifact index for the run bundle. |
| `run_summary.json` | Human-readable summary of inputs, matching, policies, and outputs. |
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
- [RELION workflow](docs/relion_workflow.md)
- [CryoSPARC workflow](docs/cryosparc_workflow.md)
- [Output files](docs/output_files.md)
- [Migration from cryoROLE 0.x](docs/migration_from_0x.md)

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
