# cryoROLE 2.0

**Understand inter-domain rotation, explore orientation landscapes, and export particle subsets for reconstruction.**

cryoROLE (cryo-EM Relative Orientation LandscapE) compares the per-particle orientations of two separately refined, approximately rigid domains. Start with CryoSPARC `.cs` or RELION `.star` metadata; obtain a relative-orientation landscape, plots, and selected metadata for downstream refinement or reconstruction.

The public interface is the `cryorole` command. The current package version is `2.0.0a1` (pre-release).

[Get started](#your-first-analysis) · [Visualize](#visualize-create-and-customize-plots) · [Select](#select-save-particle-subsets) · [Find outputs](#find-your-results) · [Help](#common-questions)

## Contents

- [What can I do with cryoROLE?](#what-can-i-do-with-cryorole)
- [Before you start](#before-you-start)
- [Install](#install)
- [Your first analysis](#your-first-analysis)
- [Explore the workflow](#explore-the-workflow)
- [Find your results](#find-your-results)
- [Common questions](#common-questions)
- [Advanced workflows and documentation](#advanced-workflows-and-documentation)
- [Citation](#citation) and [license](#license)

## What can I do with cryoROLE?

| Your task | Use |
| --- | --- |
| Compare two domains' orientations for each matched particle | `run` |
| View the population as 2D projections, 1D distributions, or interactive 3D | `visualize` |
| Express the landscape in a frame aligned with its dominant motion | `canonicalize` (optional) |
| Preview a region interactively and confirm a particle subset | `explore` |
| Save subsets by orientation, SLD, coordinate range, random sampling, or source metadata | `select` |
| Export the selected original metadata for CryoSPARC or RELION | `export` |

cryoROLE starts from existing domain-specific pose estimates. It does not perform the domain refinements or reconstruct maps. Its relative orientation is `RO = R_ref^-1 R_mov`. SLD summarizes local sampling density in the orientation landscape; higher SLD indicates more densely sampled regions under the recorded density policy.

## Before you start

Prepare two pose metadata files from separate refinements of approximately rigid domains of the **same particle population**. Choose one domain as **ref** (reference) and the other as **mov** (moving). Each matched particle must have a pose for both domains.

| Input | Required information | Default matching |
| --- | --- | --- |
| **CryoSPARC `.cs`** | `uid` and `alignments3D/pose` | Particle `uid`; no STAR conversion needed |
| **RELION `.star`** | Particle identities and `Rot/Tilt/Psi` pose columns | `_rlnTomoParticleName` first, otherwise an unambiguous `_rlnImageName` / `rlnImageName` column |

The files need not have the same row order when identity matching is available. Preflight checks the match and reports unmatched particles or other issues. Duplicate identities and zero matches block the normal run; overlap below 50% also blocks it unless explicitly overridden after review.

Use `--row-aligned` only when you know that corresponding rows describe the same particle and the row counts agree. It is an assertion about your inputs, not a general fix for matching errors. For STAR matching that needs explicit keys or preparation, see the [RELION workflow](docs/relion_workflow.md).

## Install

From a terminal with Git and Conda available:

```bash
git clone https://github.com/yifancheng-ucsf/cryorole.git
cd cryorole
conda create -n cryorole python=3.10 -y
conda activate cryorole
python -m pip install .
cryorole --help
```

The final command should list the available commands. Python 3.9 or newer is required; pip installs the core Python dependencies. Basic analysis and plotting do not require ChimeraX or FFmpeg.

See [Installation](docs/installation.md) for other installation routes and development setup. Commands below use single lines so they can be copied into Bash or PowerShell. Replace example file paths with your own and quote paths containing spaces.

## Your first analysis

This example uses **CryoSPARC `.cs` inputs** and saves the analysis in `my_run`. Run the steps in order. Keep your original input files available for later export.

### 1. Check your inputs

```bash
cryorole preflight --ref ref_domain.cs --mov mov_domain.cs
```

This checks the files and particle matching without creating an analysis directory. `READY` means the check passed; `READY_WITH_WARNINGS` means read and resolve or consciously accept the reported issues; `BLOCKED` means correct the inputs before running. Exit codes are 0, 1, and 2 respectively.

### 2. Compute the landscape

```bash
cryorole run --ref ref_domain.cs --mov mov_domain.cs --output-dir my_run
```

Open **`my_run/run_report.md`** first. It explains matching, diagnostics, output files, and next commands. Then open the PNGs in **`my_run/visualizations/quicklook/`**: Euler and rotation-vector projections for all particles, SLD ≥ 1, and the top 40% by SLD, plus a log-SLD distribution.

These are automatic previews. You already have a computed landscape; use [visualize](#visualize-create-and-customize-plots) when you want additional views or different display settings.

### 3. Save a particle subset

After inspecting your landscape, choose a region and give it a name:

```bash
cryorole select --run-dir my_run --selection-id region_01 --mode radius --center 0 0 0 --radius 15
```

**The center and radius are demonstration values, not recommended scientific cutoffs.** Replace them with values appropriate to your result. Here the center is in raw Euler degrees and the radius is 15 degrees of SO(3) rotation distance. The command reports the selected count and saves the subset under `my_run/selections/region_01/`.

For help choosing a center, use [explore](#explore-preview-and-confirm-interactively). For other selection rules, see [select modes](#select-save-particle-subsets).

### 4. Export the selected metadata

```bash
cryorole export --run-dir my_run --selection-id region_01
```

By default, both source domains are exported in their original metadata format. For these inputs, look for:

```text
my_run/exports/region_01/ref/selected_ref.cs
my_run/exports/region_01/mov/selected_mov.cs
```

These are metadata subsets for downstream CryoSPARC use. Export preserves source poses and does not modify the original files. See the [CryoSPARC workflow](docs/cryosparc_workflow.md) for more details.

### Using RELION STAR instead

Replace the first two commands with the following, then follow the same inspection, selection, and export steps. This is an alternative start: use a fresh `my_run` directory if you already ran the CryoSPARC example.

```bash
cryorole preflight --ref ref_domain.star --mov mov_domain.star
cryorole run --ref ref_domain.star --mov mov_domain.star --output-dir my_run
```

Automatic export will produce `.star` subsets. If matching is blocked (for example after signal subtraction or re-extraction changed the particle names), `preflight` lists candidate keys and prints an explicit `cryorole align` command; see the [RELION workflow](docs/relion_workflow.md#preparing-aligned-star-files).

## Explore the workflow

```text
preflight → run → inspect → select → export
                    ↳ optional canonicalize, visualize, or interactive explore
```

The examples below reuse `my_run`. Plot settings change what you see; a saved selection records which particles you choose.

### Run: compute once, inspect the report

`run` computes RO and SLD, writes the raw landscape and reports, and generates quick-look images. The default density metric is Euclidean kNN distance in rotation-vector space; SO(3) geodesic SLD is an explicit alternative through `--sld-metric so3_geodesic`.

Use `--output-dir` to name each analysis. Add `--no-visualize` to skip automatic preview images; the landscape and `run_report.md` are still written. Subsequent commands use `--run-dir` to read the saved analysis without repeating `run`.

### Canonicalize: optionally align the display frame

Canonicalization re-expresses the landscape in a motion-aligned frame. It is useful when the dominant motion is hard to interpret in raw coordinates, and is optional for selection and export.

```bash
cryorole canonicalize --run-dir my_run
cryorole visualize --run-dir my_run --space canonical --visual-id canonical_overview
```

Results are saved under `my_run/canonical/default/`, including a reusable `canonical_frame.json` and preview images. Raw results remain intact. Subsequent plotting or coordinate-based selection must use `--space canonical` to work in this frame.

| Control | Use |
| --- | --- |
| `--canonical-id NAME` | Name the frame/result; default is `default` |
| `--fit-top 0.4` | Fit using the highest-SLD fraction; default is 40% |
| `--positive-side low` or `high` | Choose the density-skew direction of the axes; default is `low` |
| `--use-frame PATH` | Apply an existing `canonical_frame.json` instead of fitting |

The fit fraction controls frame fitting, not which particles are retained. Canonical axis directions do not define an absolute biological clockwise/counterclockwise direction.

### Visualize: create and customize plots

Start with the default two PNG figures: three Euler projections and three rotation-vector projections, colored by **SLD**, for rows with `sld_display >= 1`.

```bash
cryorole visualize --run-dir my_run --visual-id overview
```

Open the files in `my_run/visualizations/raw/overview/`. Each visualization also writes `visualization_report.json` with its filters, counts, settings, and generated files.

**Choose your views:**

```bash
cryorole visualize --run-dir my_run --view 2d,1d --visual-id distributions
cryorole visualize --run-dir my_run --view 3d --visual-id interactive_3d
```

The second command writes `landscape_3d.html`: open it locally in a browser to inspect the point cloud, change representation, rotate, zoom, reset, and hover over points. It is self-contained and needs no network connection. See the [3D troubleshooting note](#why-does-an-older-3d-html-open-blank) for older blank viewers and the current validation limitation.

| Task | Options and behavior |
| --- | --- |
| Choose views | `--view 2d`, `1d`, `3d`, or a combination such as `2d,1d,3d` |
| Choose coordinates | `--representation euler`, `rotvec`, or `both` (default); Euler axes are degrees, RV axes radians |
| Choose raw/canonical | `--space raw` (default) or `canonical`; use `--canonical-id` for a named canonical result |
| Show all candidate rows | `--all` removes the default SLD threshold; 2D/3D point limits still apply |
| Filter by density | `--sld-threshold 1.5` or `--top-fraction 0.4`; these and `--all` are mutually exclusive |
| Filter by coordinates | Repeat `--range AXIS:LOWER:UPPER`, e.g. `--range alpha:-30:30` |
| Adjust only the viewport | Repeat `--axis-limit AXIS:LOWER:UPPER`; this does not filter rows |
| Style points | `--point-size 2 --opacity 0.5`; opacity ranges from 0 (transparent) to 1 (opaque) |
| Style colors | `--colormap`, `--vmin`, `--vmax`; default colormap is `rainbow_r` |
| Limit displayed points | `--max-points 50000`; deterministic sampling for 2D/3D only |
| Configure 1D | `--bins auto`, `--hist-mode percent` are defaults; `--kde` adds optional coordinate smoothing |
| Choose static output | `--format png` or, for example, `--format png,pdf`; `--3d-mode static` selects static 3D |
| Save another plot configuration | `--visual-id NAME`; default is `default`. Choose a new name or explicitly use `--overwrite` to replace that visualization |

For example, show all candidate rows in an Euler window with lighter points:

```bash
cryorole visualize --run-dir my_run --all --representation euler --range alpha:-30:30 --opacity 0.4 --point-size 2 --visual-id alpha_window
```

Or change the viewport without filtering the candidate rows:

```bash
cryorole visualize --run-dir my_run --all --representation euler --axis-limit alpha:-30:30 --visual-id alpha_viewport
```

1D distributions use **every row passing the display filters**, independently of the 2D/3D plotting sample. Optional KDE smooths a coordinate distribution; it is not an SO(3) density estimate. Color bars read `SLD`; reports retain the actual field (`sld_display` for public visualization). Color limits do not rewrite stored SLD.

### Explore: preview and confirm interactively

```bash
cryorole explore --run-dir my_run --space raw
```

Use the local browser page to inspect linked Euler/RV projections and evaluate a radius-selection draft. Review the exact count, then use **Confirm** and supply a new selection name to save it. Clicking or evaluating alone does not create a selection. Keep the local server running while using the page.

| Tool | What it saves |
| --- | --- |
| `visualize --view 3d` | A display-only HTML viewer; no selection |
| `explore` | A standard selection only after explicit Confirm |
| `select` | A standard selection from the supplied command parameters |

Explore serves packaged assets on `127.0.0.1` without a CDN. Display sampling affects only the preview. Exact evaluation and Confirm use the full parent landscape and the shared Python SO(3) selection evaluator. See the [interactive tutorial](docs/workflow_ux.md).

### Select: save particle subsets

Every selection needs a **user-chosen name** through `--selection-id`, for example `region_02` or `high_sld`. There is no default name. Results go under `my_run/selections/NAME/`; the command prints the count, location, and a suggested export command.

Use a different name for each subset. Existing names are protected; `--overwrite` explicitly replaces that selection's artifacts (the generated child names in metadata split mode). IDs must be names, not paths.

All modes evaluate the full parent landscape, independently of visualization filters or sampling. Coordinates default to `--space raw`. To select in an existing canonical frame, add `--space canonical` and, if needed, `--canonical-id NAME`. Canonical ID applies only to canonical space.

| Mode | Required controls | Useful choices and semantics |
| --- | --- | --- |
| `radius` (default) | `--center A B C` and exactly one of `--radius` / `--radius-rad` | Euler center in degrees by default; `--center-representation rotvec` takes radians. `--metric so3` is the default; `rotvec` explicitly requests Euclidean RV distance. Radius units are independent of center units. |
| `threshold` | At least one of `--sld-min` / `--sld-max` | Inclusive original `sld_raw` bounds, independent of display colors; both may be supplied |
| `range` | At least one `--range-bound AXIS:LOWER:UPPER` | Use Euler `alpha/beta/gamma` in degrees **or** RV `x/y/z` in radians; do not mix representations |
| `random` | `--fraction F`, where `0 < F <= 1` | `--seed` sets reproducibility; selects `ceil(F × parent row count)` rows without replacement |
| `metadata` | `--metadata-domain`, `--metadata-column`, and either `--metadata-value` or `--split-by-value` | Domain must be explicit (`ref` or `mov`); supports run-time CryoSPARC CS and RELION STAR metadata |

The numeric values below illustrate syntax. Choose bounds appropriate to your data.

**Around an orientation:** a raw Euler center in degrees, with a 6-degree SO(3) radius. This does not use Euclidean distance between Euler angles.

```bash
cryorole select --run-dir my_run --selection-id region_02 --mode radius --center 13 0 14 --radius 6
```

**By SLD:** keep rows with original SLD at least 2. Add `--sld-max` for an upper limit. Bounds include their endpoints.

```bash
cryorole select --run-dir my_run --selection-id high_sld --mode threshold --sld-min 2
```

**By coordinate range:** select an Euler alpha interval. Repeat `--range-bound` for multiple axes; their constraints are intersected. Blank endpoints are open (e.g. `alpha::20`). Euler lower > upper wraps across the periodic seam; RV bounds must be ordered. Repeating the same axis uses the last bound.

```bash
cryorole select --run-dir my_run --selection-id alpha_subset --mode range --range-bound alpha:-20:20
```

**A reproducible random subset:** choose 10% with a fixed seed. Omitting the seed uses fresh randomness and records a null seed.

```bash
cryorole select --run-dir my_run --selection-id sample_10pct --mode random --fraction 0.1 --seed 7
```

**By source metadata:** if the reference CryoSPARC file contains `alignments3D/class`, select classes 0 and 1, or create one selection per matched class value:

```bash
cryorole select --run-dir my_run --selection-id ref_classes_01 --mode metadata --metadata-domain ref --metadata-column alignments3D/class --metadata-value 0,1
cryorole select --run-dir my_run --selection-id ref_class --mode metadata --metadata-domain ref --metadata-column alignments3D/class --split-by-value
```

Use a field and values actually present in the run-time source; not every CS file contains a class field. Class numbers are used as stored, without renumbering. For a RELION run, use its column (for example `rlnClassNumber`) and class values instead. Split results use names such as `ref_class_0` with filename-safe value components.

CS supports scalar integers, booleans, and text. Integer comparison preserves full precision; booleans accept `true`, `false`, `1`, or `0`. UTF-8 text matches exactly, including leading zeros and spaces; empty strings are counted as missing and excluded. Float/vector fields and invalid UTF-8 are rejected. Commas separate requested values; escaping a comma inside one requested value is not supported.

CS metadata selection verifies the recorded source SHA-256 and requires valid source-row indices. It does not rematch particles or accept an external annotation file. CS splitting is limited to 100 distinct non-missing values; use explicit value selection for higher-cardinality fields. Child-name collisions or existing outputs are checked before writing; `--overwrite` does not resolve ambiguous names.

**Inspect and export a saved subset:** reuse the `region_01` selection from the first workflow. `--all` shows all selected candidates subject to the plotting point limit.

```bash
cryorole visualize --run-dir my_run --selection-id region_01 --all --visual-id review
```

Open `my_run/visualizations/selections/region_01/parent_raw/review/`. After inspection, export using the command in the next section. No extra selected-landscape file is required for this workflow.

Optionally, add `--write-selected-landscape` when creating a selection to save a derived landscape as well. It inherits parent SLD unless you explicitly also request `--recompute-sld`, which recomputes subset SLD and preserves parent SLD fields separately. View that artifact with `visualize --selection-id NAME --use-selected-landscape`. Recomputing SLD is not required for export.

For complete option details and mode-specific errors, run `cryorole select --help` or see the [CLI reference](docs/cli_reference.md#selection-modes-and-names).

### Export: use the subset downstream

For the first workflow's saved selection:

```bash
cryorole export --run-dir my_run --selection-id region_01
```

This is the same export as in the first workflow; run it once. If you deliberately regenerate an existing export, add `--overwrite`. To export another subset, substitute its selection ID.

The defaults are `--domain both` and `--format auto`. Choose `--domain ref` or `mov` for a single domain. Output is under `my_run/exports/SELECTION_ID/`, with a report and per-domain `.cs` or `.star` files.

Export subsets the original metadata using saved source-row provenance; it neither reselects particles nor writes canonical/display coordinates into source poses. Source content is verified against the run. If files moved, use the explicit relocation options described in the [FAQ](docs/faq.md).

## Find your results

Start with the report and images. Other directories appear as you run the corresponding commands:

```text
my_run/
  run_report.md                    Start here: results, diagnostics, next steps
  data/raw_landscape.csv           Per-particle table for inspection
  data/raw_landscape.npz           Numeric landscape used by cryoROLE
  visualizations/quicklook/        Automatic run previews
  visualizations/raw/overview/     Example named visualization
  canonical/default/              Optional canonical landscape and frame
  selections/region_01/            Saved subset and selection summary
  exports/region_01/               Export report and source metadata subsets
```

`run_manifest.json`, `run_summary.json`, and `reports/` record provenance and diagnostics. Keep the run bundle together for downstream commands. See [Output files](docs/output_files.md) for the complete layout.

## Common questions

### Why are fewer particles shown than reported by run?

Default `visualize` filters to `sld_display >= 1`; 2D/3D also use display point limits. `--all` removes the SLD threshold, while `--max-points` controls the plotting cap. Read `visualization_report.json` for candidate, filtered, and plotted counts. None of these settings removes particles from the parent landscape.

### Why did a display filter not create an exportable subset?

Display filters only control plots. Use `select`, or Confirm in `explore`, to save a named subset before export.

### Which coordinates and SLD should I use?

Raw coordinates are the direct RO result; canonical coordinates express that result in an optional motion-aligned frame. Choose the same space when reading a center from a plot and selecting around it. Public Euler coordinates use extrinsic fixed-axis ZYX in degrees; rotation-vector coordinates are in radians.

SLD color bars use one readable label, while artifacts retain `sld_raw` (scientific density) and `sld_display` (display field). Threshold selection uses `sld_raw`; color scaling does not change it. SLD is a sampling-density measure, not a direct particle-quality score.

### Why is selection-id required?

You may save many regions from one run. An explicit name makes them distinguishable and avoids silently replacing a default subset. Choose a new name for a new region; use `--overwrite` only to intentionally replace an existing one.

### Do RO coordinate coincidences mean duplicate particles?

No. The diagnostic counts rows sharing quantized RO rotation-vector coordinates (default grid step `1e-8 rad`). Dispersed pairs or triples are informational. Under the initial heuristic policy, groups of at least 10 contribute to the concentration count; a warning appears if their combined size reaches 100 rows **or** 1% of all landscape rows.

These thresholds are diagnostic heuristics, not scientific cutoffs. Coordinate coincidence alone establishes neither duplicate particle identity nor duplicate images. Existing identity/matching diagnostics handle their own evidence and failures. The RO diagnostic does not remove particles or change SLD.

### Why are the coordinates of my subtracted particles wrong?

RELION's Particle subtraction with recentring (`--center_x/y/z`) moves each box but does not update `_rlnCoordinateX/Y`, so the subtracted STAR places particles up to the projected recentring shift away from their true centre (a median of 92 px in our test data). Refinement is not affected. Re-extraction, polishing, distance-based duplicate removal and coordinate matching are affected. `cryorole preflight` warns about it, and `cryorole align --fix-subtract-coordinates Subtract/jobNNN/` writes a verified, corrected copy without changing your files. See the [RELION workflow](docs/relion_workflow.md#signal-subtraction-with-recentring-leaves-stale-coordinates).

### Why does an older 3D HTML open blank?

Older generated files can contain a JavaScript newline-escaping error. Regenerate the viewer with the updated code and a new `--visual-id`; existing HTML files are not repaired automatically. Generated-script and simulated interaction tests pass, but actual Chrome rendering and interaction acceptance for this fix remains pending. See [FAQ](docs/faq.md) for troubleshooting.

## Advanced workflows and documentation

| Need | Where to go |
| --- | --- |
| CryoSPARC input and metadata export | [CryoSPARC workflow](docs/cryosparc_workflow.md) |
| RELION input, matching, and optional `align` preparation | [RELION workflow](docs/relion_workflow.md) |
| Resume work, inspect status, or get next-step guidance | [Workflow tutorial: `status`, `next`, `guide`, `explore`](docs/workflow_ux.md) |
| Look up command arguments | [CLI reference](docs/cli_reference.md) and `cryorole COMMAND --help` |
| Inspect landscape CSVs in ChimeraX | [Standalone viewer](cryorole_chimerax_viewer.py): its opening documentation gives setup and commands |
| Render landscape/rigid-body movies or canonical camera views | [Animation and canonical-views guide](docs/animation_export.md); structure execution requires ChimeraX, MP4 encoding additionally requires FFmpeg/FFprobe |
| Apply a diagnostic rotation to a derived landscape | [Rotation script](scripts/rotate_landscape.py), with usage through `--help` |
| Resolve setup or workflow errors | [Installation](docs/installation.md) · [FAQ](docs/faq.md) |
| Inspect files or migrate an older workflow | [Output files](docs/output_files.md) · [0.x migration](docs/migration_from_0x.md) |

Animation represents an interpolated rigid-body rendering, not a new reconstruction at each trajectory position. Its frame, baseline, and pivot prerequisites are covered in the linked guide.

## Citation

If you use cryoROLE as a method, please cite the cryoROLE methodology preprint:

- Chengmin Li, Wooyoung Choi, Hao Wu, Yifan Cheng. **CryoROLE: describing large inter-domain rotation in single particle cryo-EM.** bioRxiv, 2026. [https://www.biorxiv.org/content/10.64898/2026.07.04.736454v1](https://www.biorxiv.org/content/10.64898/2026.07.04.736454v1)

The first application of cryoROLE was in the human fatty acid synthase study:

- Wooyoung Choi, Chengmin Li, Yifei Chen, YongQiang Wang, Yifan Cheng. **Structural dynamics of human fatty acid synthase in the condensing cycle.** Nature, 2025. [https://doi.org/10.1038/s41586-025-08782-w](https://doi.org/10.1038/s41586-025-08782-w)

## License

cryoROLE is distributed under the BSD 3-Clause License.
