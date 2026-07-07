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

## Visualizations

`visualizations/` contains display-only figures and display tables. Display
filters, display ranges, downsampling, and color scaling do not alter raw or
canonical landscapes and do not create scientific selections.

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

## Exports

Exports live under:

```text
exports/<selection_id>/
```

Export reads an explicit selection and writes source metadata subsets for the
requested domain (`ref`, `mov`, or `both`). Export does not reselect, rematch, or
rewrite source poses with display/canonical coordinates.

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
