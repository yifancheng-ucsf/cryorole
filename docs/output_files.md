# cryoROLE Output Files

cryoROLE 2.0 writes a run bundle: one directory that contains numeric arrays,
flat tables, reports, visualizations, selections, exports, and provenance. The
default run bundle is `cryorole_outputs/`; `cryorole run --output-dir RUN_DIR`
writes the same layout under `RUN_DIR/`.

## Run Bundle Overview

A typical run directory contains:

```text
run_manifest.json
run_summary.json
run_report.md
bundle_state.json
.cryorole_bundle_complete

data/
  raw_landscape.npz
  raw_landscape.csv
  match_table.csv

reports/
  *.json

visualizations/
canonical/
selections/
exports/
```

## Top-Level Files

`run_manifest.json`

Records the run inputs, policies, artifact paths, and provenance needed to
audit downstream analysis.

`run_summary.json`

Provides a concise human-readable summary of the run, including inputs,
matching, resolved policies, and important output paths.

`run_report.md`

Provides a short human guide to the run outputs, quick-look figures, warnings,
and next commands. JSON reports and the manifest remain authoritative.

### Input sanity

Every run reports the RO-angle distribution (median, 90th and 99th percentile)
in `run_report.md` under "Input sanity" and in `run_summary.json` under
`input_sanity`. Warnings appear only for patterns that usually mean an input
mistake; the thresholds are heuristics and are recorded in
`input_sanity.policy`. They never change the analysis.

| Code | Level | Trigger |
| --- | --- | --- |
| `SAME_INPUT_FILE` | strong warning (blocked in `preflight`) | `--ref` and `--mov` resolve to the same path or have the same SHA-256 |
| `IDENTICAL_POSES` | strong warning | at least 99% of particles have an RO angle below 1e-6 rad |
| `NEARLY_IDENTICAL_ORIENTATIONS` | warning | median RO angle below 1° and 99th percentile below 2° |

The summary also reports the fractions below 0.1°, 0.5°, 1° and 5° and the
histogram mode. The findings are prompts to inspect the inputs, not a
conclusion that the pairing is wrong: stable domains can have small relative
rotations.

### Alignment lineage (`--row-aligned` runs)

`run_summary.json` and `run_manifest.json` (`results.alignment_provenance`)
record the `cryorole align` lineage when an `align_report.json` next to the
inputs verifies. The record has these fields:
- `attached`;
- `lineage`: strategy, original paths and hashes, match table;
- `original_files`: `verified`, `unavailable` or `mismatch` for each original.

When the report cannot be verified, `attached` is `false` and `reason` says why.
The field is `null` when no report is found. This never affects the run.

If every particle has the same RO, the local density is undefined and `run`
stops before writing a landscape, with the input-sanity explanation as the
error message.

New summaries and manifests share a unique `run_id` and a ref/mov
`source_identities` mapping containing the original path, run-resolved absolute
path, source type, file size, mtime, streamed SHA-256, and source row count.
Selections, canonical reports/frames, and exports record their parent run ID.

`bundle_state.json` records transactional lifecycle history.
`.cryorole_bundle_complete` is written last; new transactional bundles without
both a completed manifest state and this marker are rejected downstream.

## Data Directory

`data/raw_landscape.npz`

The production machine-readable source of truth for the raw landscape. This is
the preferred input for downstream cryoROLE commands.

`data/raw_landscape.csv`

A user-facing flat table for inspection, plotting in external tools, and quick
checks. The CSV is derived from the stored landscape arrays.

`data/match_table.csv`

Records matched particle provenance, including source rows used for later
selection/export backtracking.

## Reports

`reports/*.json` files record import, identity, matching, density, and other
policy/report details. JSON reports and manifests are the audit layer.

Full object-record landscape JSON is debug-only and is not the production
persistence contract.

`cryorole preflight --json REPORT.json` may write a schema-1.0 report at an
explicit path outside a run bundle. The report contains source identities,
matching and pose-schema diagnostics, resource-estimate formulas/assumptions,
readiness, and exact commands. Preflight does not create the run directory.

If `cryorole guide --execute-run` launches a run, it may add
`reports/guide_history.json` to record that orchestration action. Status itself
is derived from actual artifacts and does not maintain a separate mutable
workflow-state file.

## Visualizations

`visualizations/` contains display-only figures, offline viewers, and reports. Display
filters, display ranges, downsampling, and color scaling do not alter raw or
canonical landscapes and do not create scientific selections.

The default `run` quick-look is intentionally minimal:

```text
visualizations/quicklook/
  all_euler_3view_projection.png
  all_rotvec_3view_projection.png
  sld_ge_1_euler_3view_projection.png
  sld_ge_1_rotvec_3view_projection.png
  top_40pct_euler_3view_projection.png
  top_40pct_rotvec_3view_projection.png
  sld_log_distribution.png
```

That directory contains only the seven PNG files. Range counts, top-cutoff tie
policy, full-data histogram diagnostics, and the shared color bound are
recorded in `run_summary.json` and the manifest artifact index. The bundle-root
`run_report.md` explains them. Use explicit `cryorole visualize` for expanded
display products.

Each explicit public visualization directory includes `visualization_report.json`.
The default also contains only `euler_3view_projection.png` and
`rotvec_3view_projection.png`. Optional 1D and 3D products are requested with
`--view` and remain flat in that directory. For NPZ inputs the report records
the full parent/input count, count after display filtering, independently
resolved per-view counts, random seed, and deterministic sampling policy. 1D
uses every filtered row; 2D/3D point limits never alter the scientific parent
landscape or create a selection.

## Canonical Landscapes

Canonical outputs live under:

```text
canonical/<canonical_id>/
```

Typical files include:

```text
canonical_landscape.npz
canonical_landscape.csv
canonicalization_report.json
canonicalize_summary.json
canonical_frame.json
canonical_frame.npz
```

Canonicalization derives a coordinate frame and writes new artifacts. It does
not overwrite raw landscape artifacts.

## Resolved options (`resolved_by`)

When `--run-dir`, `--canonical-id` or `--selection-id` is omitted and cryoROLE
fills it in (see `docs/cli_reference.md`, Common behaviour), the report written
by that command records how, for example:

```json
"resolved_by": {
  "run_dir": {"value": "cryorole_outputs", "resolved_by": "default_output_dir"},
  "canonical_id": {"value": "default", "resolved_by": "only_candidate"}
}
```

`resolved_by` is `explicit`, `current_directory`, `default_output_dir` or
`only_candidate`. It appears in `canonicalize_summary.json`,
`selection_summary.json`, `visualization_report.json`, the export
`export_report.json`, and the `status` / `next` JSON output.

## Selections

Selections live under:

```text
selections/<selection_id>/
```

Typical files include:

```text
selection.json
selection.csv
selected_particle_keys.csv
selected_landscape_rows.csv
selection_summary.json
```

A selection is a scientific decision artifact. It is not the same thing as a
visualization filter.

Random-mode selections always record the seed used: `random_seed` in
`selection.json` and `selection_summary.json`, with `random_seed_source`
(`user` for `--seed`, `generated` when it was omitted) in the summary.

CLI select and interactive Confirm use the same standard artifact writer for
these five files. `selected_landscape_rows.csv` retains ref/mov source-row IDs
for export backtracking; re-export consumes the existing Selection and does not
create a different selection schema.

An interactive explorer click or evaluation is only an in-memory draft and
writes nothing here. Explicit Confirm creates the normal files above and adds
interaction provenance to `selection.json`: timestamp, parent run ID,
landscape path/SHA-256/space, input and evaluated center, Euler convention,
SO(3) metric and radius, full/selected counts, and display filter/downsample
settings. Confirm refuses to overwrite an existing selection ID.

## Exports

Exports live under:

```text
exports/<selection_id>/
```

Export reads an explicit selection and writes source metadata subsets for the
requested domain (`ref`, `mov`, or `both`). Export does not reselect, rematch, or
rewrite source poses with display/canonical coordinates.

Before subset writing, export verifies the current source against the recorded
SHA-256. `export_report.json` records per-domain verification and relocation
status. Export also checks that the source particle table still has the row
count the run recorded (`source_row_count_verified` in each domain report) and
refuses to export on a mismatch. STAR subsets always come from the loop in the
`data_particles` block; other blocks (optics, tomograms) are copied verbatim.
A missing original can be replaced only with explicit
`--relocated-ref` / `--relocated-mov` whose hash matches. Legacy bundles without
hashes require `--allow-unverified-source`, and that decision is recorded.

## NPZ, CSV, and JSON Roles

- NPZ: machine-readable numeric arrays and the production landscape source of
  truth.
- CSV: user-facing flat tables.
- JSON: reports, summaries, manifests, policies, and provenance.

## Derived Coordinate Conventions

Rotation matrices are the internal source of truth. Rotation vectors,
quaternions, and Euler angles are derived representations.

Public RO/RV-derived Euler output uses extrinsic fixed-axis ZYX. CSV column names
may use compact labels such as `raw_ea_zyx_alpha_deg`; the reports and manifests
record the resolved Euler convention.

## SLD Fields

`sld_raw`

The scientific SLD density value used by default selection and canonicalization
policies.

`sld_unfloored`

SLD without the distance floor. A particle whose k nearest neighbours all have
exactly its RO (for example duplicated particles, symmetry expansion, or a rigid
subpopulation larger than k) has zero local distance and `+inf` here; the count
is `n_inf_sld_unfloored` in `reports/density_report.json` and is explained in
`run_report.md`. `sld_raw` stays finite because of the distance floor. NaN is
never a valid SLD value.

`sld_display`

A display-only density/color value. It must not be treated as the default
scientific density field.

`sld_display_is_outlier`

A display-policy diagnostic for tail-jump high-density outliers. Marked rows
remain present in the landscape unless a later explicit selection policy says
otherwise.

## Raw, Canonical, Visualization, and Selection

- Raw landscape: direct relative-orientation facts from the run.
- Canonical landscape: derived coordinate frame for inspection and comparison.
- Visualization: display-only renderings and display tables.
- Selection: explicit scientific particle subset with provenance for export.

Interactive display sampling is also display-only. The UI reports both the
displayed count and full candidate count; exact selection always evaluates the
full parent landscape.
