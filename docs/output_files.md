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
New identities also record `input_path_absolute`: the absolute input path before
resolving links. The post-read guard checks that this path still resolves to the
recorded source, including links in parent directories, and rechecks source
content. Retargeting a link during input loading aborts the run before publication.
Old identities without the additional field retain content checks; an original
absolute input path also permits checking its current target.
Selections, canonical reports/frames, and exports record their parent run ID.

Source convention policies record RELION Euler input in degrees and CryoSPARC
pose vectors in radians. A custom policy that contradicts those format units is
rejected before normalization.

STAR table readers remove surrounding single or double quotes from cell values
before particle matching and numeric pose conversion. Quoting alone does not
change particle identity. STAR subset export still copies the selected source
rows verbatim, including their original quotes and spacing.
Unquoted `#` comments after a complete row are ignored by readers; quoted
hashes and hashes inside a path remain literal. Subset export retains those
original row comments together with the selected source text.

`bundle_state.json` records transactional lifecycle history.
`.cryorole_bundle_complete` is written last; new transactional bundles without
both a completed manifest state and this marker are rejected downstream.

Publication rechecks whether the destination exists, so concurrent runs without
`--overwrite` cannot replace a run that finished while they were computing.
A short-lived sibling `.<run-name>.publish-lock` directory serializes cooperating
publishers. During overwrite, the previous bundle is retained under
`.<run-name>.previous-<run-id-prefix>` until the new bundle is published. A handled
publication error, including Ctrl-C, restores the previous bundle. A forced
process termination can leave a lock or backup; inspect those artifacts before
retrying. A backup cleanup error leaves the new completed bundle in place.
The recovery includes an interruption raised just after a filesystem rename
has completed. A new destination is removed if its first install is interrupted;
an overwrite restores the exact previous bundle. Failed recovery keeps the old
backup, prepared stage (if still present) and publication lock for manual recovery,
without moving or relabeling the prepared artifacts as a failed debug bundle.

## Data Directory

`data/raw_landscape.npz`

The production machine-readable source of truth for the raw landscape. This is
the preferred input for downstream cryoROLE commands.

`data/raw_landscape.csv`

A user-facing flat table for inspection, plotting in external tools, and quick
checks. The CSV is derived from the stored landscape arrays.

When a landscape CSV is read through the compatibility path, `particle_key`
remains exact text: leading zeros, large integer identifiers and strings such as
`NA` or `null` are preserved. Numeric and boolean columns retain their own types.

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
Before returning readiness, it rechecks input content and path resolution against
the recorded source identities; a change during inspection blocks the result.
Suggested run commands preserve the assessed domain labels, density metric and
neighbor count, backend, query batch size, CSV/preview choices and identity policy.
These checks apply at inspection boundaries; they are not immutable input snapshots.

Run, alignment and subtraction-correction output directories must be separate
from their input metadata and mapping files. An output directory that contains
an input file or its recorded path is refused, including with `--overwrite`.
The indexed STAR subset writer also refuses writing through a source alias or
hard link and validates integer source-row bounds before replacing an output.

Canonicalization, visualization, animation, canonical camera views and the
rotated-landscape diagnostic also refuse outputs that contain their consumed
landscape, frame or configuration inputs, including path aliases and a frame's
referenced NPZ. `--overwrite` does not authorize replacing these sources. Choose
a separate derived directory or canonical ID when reusing a frame. Animation
still keeps its failure manifest, logs and generated frames after a rendering
failure; these input checks run before output replacement.
Writers that update a derived directory in place also check existing file and
directory links, so a linked output cannot write through to an external source.

Alignment and subtraction-correction outputs are staged as one group and their
STAR inputs are rechecked before publication. A caught error or interruption
before commit preserves the previous output group; unrelated files survive an
explicit overwrite. Reports and suggested commands refer to the published paths.
Cooperating writers use a destination lock. Power loss or an unsuccessful rollback
can leave a lock and a hidden stage directory with backups for manual recovery;
retain those files until recovery. Cleanup failure after commit leaves the new
completed output group in place. These guarantees do not cover uncoordinated
external writers or transient input replacements restored between boundary checks.

Explicit alignment float tolerances must be finite and positive, and their
numeric key values must remain finite after scaling. For a cross-column key pair,
a tolerance named on either side applies to both columns in the same units;
conflicting tolerances are refused. The alignment key policy records the resolved
tolerances. Optics values and group IDs decode STAR quotes without changing the
stored source text; duplicate optics-group IDs are ambiguous and refused.
Subtraction correction with `--apply-to` requires unique image names among the
verified subtracted source rows as well as unique target names. Ambiguity stops
the operation; no first/last-occurrence rule or automatic deduplication is applied.
Repeatable `--float-tol` and `--identity-column` arguments retain every supplied
value in order; parsing another command with the same parser keeps empty defaults.

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
The standalone ChimeraX viewer reads UTF-8 CSV/TSV/semicolon tables with standard
quoted fields, including particle labels containing delimiters, quotes or line
breaks. Quoted labels do not shift the numeric coordinate/color columns. Legacy
whitespace tables and the existing nonfinite-row display filter remain supported.
The default also contains only `euler_3view_projection.png` and
`rotvec_3view_projection.png`. Optional 1D and 3D products are requested with
`--view` and remain flat in that directory. For NPZ inputs the report records
the full parent/input count, count after display filtering, independently
resolved per-view counts, random seed, and deterministic sampling policy. 1D
uses every filtered row; 2D/3D point limits never alter the scientific parent
landscape or create a selection.

An explicit GUI **Save figures** action writes the same PNG/report layout under
`visualizations/raw/<name>/` or `visualizations/canonical/<canonical-id>/<name>/`.
It stages the complete set, refuses existing names and publishes only after
rendering and parent-identity verification succeed. `gui_display_snapshot` in
the report records the run ID, parent landscape hash, equivalent command, and
that highlight overlays/mouse pan-zoom were excluded. The saved filter, colour
scale, coordinate representation and axis ranges use standard visualization
fields. No Selection or source-metadata artifact is written by saving figures
or loading a saved-selection overlay.

## Canonical Landscapes

Canonical outputs live under:

```text
canonical/<canonical_id>/
```

Typical files include:

```text
canonical_landscape.npz
canonicalization_report.json
canonicalize_summary.json
canonical_frame.json
canonical_frame.npz
canonical_landscape.csv        # only with --write-csv
```

`canonical_landscape.csv` is not written by default: it is large and
`canonical_landscape.npz` holds the same data. `--write-csv` writes it.

The canonical preview figures go to
`visualizations/canonical/<canonical_id>/<display set>/`, one folder per display
set (all particles, `sld_raw > 1.5`, and the top-SLD fit fraction). Each has the
Euler and rotation-vector 3-view projections and 3D views. `--projection-panels`
also writes every 3-view panel as its own PNG (`euler_alpha_beta.png`, …,
`rotvec_yz.png`; 18 extra figures in total).

Canonicalization derives a coordinate frame and writes new artifacts. It does
not overwrite raw landscape artifacts.

`canonicalization_report.json` flags near-degenerate PCA axes, including a fit
whose points lie on a line: its principal direction is determined, but the two
transverse axes are not unique. The canonical frame is still written with the
warning. Distinct planar principal axes do not have that two-direction ambiguity.

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

The selection files, export report and optional `selected_landscape/` are
prepared together before publication. A caught write failure or interruption
keeps the previous selection unchanged; a new selection is not exposed partly
written. Explicit overwrite retains unrelated files but replaces the managed
selection outputs. If the new request omits a selected-derived landscape, an
older `selected_landscape/` is removed so it cannot describe the previous subset.
Reports contain published paths. Directory publication coordinates cooperating
writers; an uncatchable process kill or failed recovery can leave a private
backup for manual recovery.

RELION metadata selection preserves exact integer categories, including adjacent
64-bit labels. Equivalent integral numeric spellings such as `01`, `1.0` and
`1e0` still match the same value. CryoSPARC metadata uses its recorded scalar
type and retains its existing type-specific comparison rules.

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

Export checks each selected particle's source-row mapping against the parent
run's raw landscape. If the Selection includes a companion CSV, its row IDs
must agree. Existing selections without the companion use the parent mapping
directly. A parent without the required mapping must be regenerated before
exporting source metadata. Both domains are prepared together before publication.

Before subset writing, export verifies the current source against the recorded
SHA-256 and rechecks it after preparing the export, immediately before publication.
A persistent source change during reading refuses publication and preserves an
existing export. The check binds the resolved or explicitly relocated export
source; the original run-time alias need not remain available.
The output directory must not contain a recorded absolute input entry, its
resolved target, or a currently used relocated input. This includes input
symlinks whose targets are outside the output and paths that enter the output
through parent-directory links. Explicit overwrite does not permit replacing
these input names. A separate output directory and ordinary overwrite remain
valid. Old relative-only provenance without a saved base is protected through
its recorded resolved source; export does not guess the historical cwd.
`export_report.json` records per-domain verification and relocation
status. Export also checks that the source particle table still has the row
count the run recorded (`source_row_count_verified` in each domain report) and
refuses to export on a mismatch. STAR subsets always come from the loop in the
`data_particles` block; other blocks (optics, tomograms) are copied verbatim.
A missing original can be replaced only with explicit
`--relocated-ref` / `--relocated-mov` whose hash matches. Legacy bundles without
hashes require `--allow-unverified-source`, and that decision is recorded.
For that explicit legacy override, `export_source_sha256` binds the current
read only; `sha256_verified` remains false because the old run recorded no hash.

Native CS and STAR subset writers require one-dimensional integer row IDs,
without duplicates or out-of-bounds values. CS particle tables must also be
one-dimensional. Output files, including an optional CS UID list, must be
distinct from each other and from source paths, symlinks and hard links. All
requested file destinations are checked before writing. Prepared files publish
as a recoverable group: a caught failure restores the previous files, and failed
recovery retains hidden stages, backups and locks. File readers can observe
individual renames; the guarantee is not a filesystem-wide atomic transaction.
Export directory publication likewise restores the previous directory after a
caught interruption. These checks do not guarantee immutable reads against
uncoordinated external writes or a source changed and restored between checks.

```text
exports/<selection_id>/
  export_report.json             # includes next_steps per domain
  NEXT_STEPS.txt                 # the same steps as plain text
  selected_particle_keys.txt
  selected_landscape_rows.csv
  ref/ (and/or mov/)
    selected_ref.star | selected_ref.cs   # verbatim subset of the run-time source
    selected_ref_uids.txt                 # CryoSPARC: selected uids
    selected_ref_<FILE>.cs                # CryoSPARC: subset of each other file of the job's group
    selected_ref.csg                      # CryoSPARC: result group for Import Result Group
    export_ref_report.json                # incl. cryosparc_group (source .csg, how it was found, slots)
```

The CryoSPARC result group mirrors the source job's own `.csg` slot by slot
(e.g. `alignments3D` and `ctf` from `selected_ref.cs`, `blob` from
`selected_ref_J75_passthrough_particles.cs`); see `docs/cli_reference.md`.

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
its RO (for example duplicated particles, symmetry expansion, or a rigid
subpopulation larger than k) has zero local distance and `+inf` here. "Zero"
means a mean kNN distance of at most 1e-12 rad, so that floating-point rounding,
which differs between platforms, does not change the result. The count
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
