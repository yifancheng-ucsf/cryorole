# Workflow UX Tutorial

For the same workflow in a desktop window, including saved-selection review
and standard figure export, see [Desktop GUI](gui.md). Offscreen checks cover
2D/1D and artifact behavior; native desktop/3D acceptance remains separate.

This tutorial follows the production workflow without hiding scientific
decisions. The compact path is:

```text
preflight -> run -> status/next -> optional canonicalize -> explore -> Confirm -> export
```

## 1. Understand the two inputs

Reference domain (`--ref`) means the domain whose frame defines RO. Moving domain
(`--mov`) means the domain expressed relative to it. Their metadata must
describe the same particle population. cryoROLE matches
particles by a safe identity field (`uid` for CryoSPARC; supported RELION
particle-name fields for STAR) unless `--row-aligned` is an explicit user
assertion.

The scientific definition remains:

```text
RO = R_ref^-1 R_mov
```

Changing which domain is ref versus mov changes the relative-orientation
question. The guide never chooses this for you.

## 2. Preflight without side effects

```bash
cryorole preflight --ref ref_domain.star --mov mov_domain.star
```

or:

```bash
cryorole run --ref ref_domain.star --mov mov_domain.star --dry-run
```

Both use the same service as production run and create no run bundle. The
result is:

- `READY` (exit 0): validation passed;
- `READY_WITH_WARNINGS` (exit 1): review matching/resource warnings;
- `BLOCKED` (exit 2): fix errors before computation.

Use `--json` or `--json preflight.json` for schema-1.0 structured output. It
contains input path/hash identity, source convention, pose-schema checks,
duplicate/collision examples, overlap counts, transparent memory/disk formulas
and assumptions, and exact commands. Resource numbers are estimates, not
guarantees.

In the desktop Inputs page, choosing another file immediately clears the old
file choice and readiness until the new inspection finishes. If several file
checks overlap, only the latest choice can update the run inputs. Changing an
input or parameter also invalidates an in-progress preflight result; Start run
waits for the current inputs and options to pass their own checks.

## 3. Run and inspect actual state

```bash
cryorole run --ref ref_domain.star --mov mov_domain.star --output-dir my_run
cryorole status --run-dir my_run
cryorole next --run-dir my_run
```

Run revalidates source identity and publishes the bundle transactionally.
Without `--overwrite`, publication refuses a destination created by another run
while computation was in progress. A simultaneous publisher for the same
destination must finish before retrying.
`status` derives state from the manifest, bundle marker, raw NPZ, source hashes,
canonical frames, visualizations, selections, and exports. `next` produces
exact commands labelled required, recommended, or optional. It does not treat
a visualization as a selection, and it does not require canonicalization.
Both exit with status 3 when the bundle is missing, failed or incomplete, so a
script can stop early. From inside `my_run` (or with the default
`cryorole_outputs`), `--run-dir` can be omitted; cryoROLE prints which bundle it
used.

## 4. Raw and canonical coordinates

The raw landscape is the direct relative-orientation result. A canonical
landscape is an optional, separately stored motion-aligned coordinate frame:

```bash
cryorole canonicalize --run-dir my_run
```

Canonicalization helps interpretation; it does not replace raw facts or write
display coordinates back to STAR/CS poses.

Derived output directories must be separate from their input landscapes, frames
and configuration files, even with `--overwrite`. When applying an existing
frame, choose another canonical ID or a separate output directory. Visualization
accepts `--opacity 0` through `1`; zero hides the points while retaining axes
and the color scale.

RV means rotation vector: its direction is the rotation axis and its magnitude
is the rotation angle in radians. EA means the public derived Euler display,
using recorded extrinsic fixed-axis ZYX. Euler coordinates are useful to read,
but scientific radius selection uses rotation-native SO(3) geodesic distance.

SLD is local sampling density. `sld_raw` is scientific data;
`sld_display` and display outlier flags are presentation aids. A high-SLD
display threshold is not a particle selection.

## 5. Explore locally

```bash
cryorole explore --run-dir my_run --space raw
```

For an existing canonical frame:

```bash
cryorole explore --run-dir my_run --space canonical --canonical-id default
```

The page opens from a server bound only to `127.0.0.1`, uses packaged assets,
and has no CDN dependency. It shows linked EA and RV projections. Click a point
to propose a center, choose a radius, then Evaluate.

The displayed cloud may use deterministic display downsampling for responsiveness.
The UI always reports displayed and full candidate counts separately. Display
threshold/top-fraction controls and downsampling affect presentation only.
Every exact count and Confirm evaluates the full parent landscape in Python via
the same SO(3) radius evaluator used by CLI selection.

An evaluated draft is not a Selection and writes no scientific artifact. You
can download draft/preview data for review, but export cannot consume it.
Python draft reports and display masks are detached copies: editing them does
not change the internally evaluated particles or the geometry used by Confirm.
Change the centre or radius by evaluating a new draft before confirming it.
Equivalent selection commands preserve the full run path, including spaces
or words that also appear in internal placeholder names.

## 6. Confirm a scientific Selection

Confirm in the explorer only after reviewing the center, evaluated RV/EA,
radius, metric, full candidate count, and selected count. Supply a new
selection ID. Confirm refuses overwrite and verifies the run ID and parent
landscape SHA-256 before writing:

```text
my_run/selections/<selection_id>/
```

The standard Selection records the parent run/landscape, input and evaluated
center, Euler convention, SO(3) metric, radius, counts, timestamp, and UI
display/filter/downsample provenance. The equivalent non-browser command is:

```bash
cryorole select --run-dir my_run --selection-id state_1 --space raw -c A B C -r DEG
```

For CLI selection, `--selection-id` is your chosen name, not an automatically
assigned identifier. Keep separate regions under separate names. An existing
name requires explicit `--overwrite`; ordinary CLI commands do not prompt for
missing names. Use `cryorole select --help` for mode-grouped parameters and
complete examples (radius, threshold, range, random, and metadata). See the
selection table in `docs/cli_reference.md` for units and dependencies. Selection
success guidance includes the saved name and an export command on stderr;
stdout remains the output directory.

Selection files and an optional selected-derived landscape are saved as one
group. A failed write keeps an existing selection unchanged. When overwriting,
request `--write-selected-landscape` again to create a landscape for the new
subset; otherwise the old derived landscape is removed. Unrelated files remain.
Metadata selection and split-by-value preserve distinct large integer labels
in RELION sources, including values above the exact-integer range of floats.

## 7. Export without reselecting

```bash
cryorole export --run-dir my_run --selection-id state_1 --domain both
```

Export checks the confirmed Selection against the parent run's source-row
mapping before publishing the requested subsets. It
verifies each source hash, does not rematch or reselect, does not modify source
files, and does not write canonical/display coordinates as physical poses.
Choose an output directory separate from the recorded input paths and any
explicitly relocated inputs. `--overwrite` can replace an earlier export but
cannot replace an input symlink inside that directory, even when the source
file itself is elsewhere. Input links outside the export directory are valid.

## 8. Optional guide

For a new analysis:

```bash
cryorole guide --ref ref_domain.star --mov mov_domain.star --non-interactive
```

To resume:

```bash
cryorole guide --run-dir my_run --non-interactive
```

Non-interactive guide prints a plan and never waits. `--execute-run` is an
explicit opt-in to execute a ready run. Guide never silently enables row
alignment/low overlap, chooses canonicalization or a neighborhood, creates or
overwrites a Selection, or exports metadata.

### Evaluated draft snapshots

Changing the browser's center or radius clears the evaluated draft and disables
Confirm until a new exact evaluation finishes. Display filters and representation
changes remain display-only. Picking considers only visible particles; top-fraction
filtering sorts once per redraw, and large display samples avoid argument-count
limits when calculating plot/color bounds.

Each evaluated result has a `draft_id`. HTTP Confirm must return that ID together
with the session/run/landscape identity; a missing or stale ID stops before writing.
Reload older Explore pages after upgrading. Concurrent evaluation publishes one
complete snapshot at a time, and confirmation retains the snapshot it approved even
if another evaluation finishes during source revalidation. The Python session's
optional `expected_draft_id` supplies the same check, and the desktop view always
uses it. Confirmed Selection interaction provenance records the evaluated draft ID.
Initial Explore top-fraction display filtering uses the shared cutoff policy and
retains tied values before deterministic sampling, matching other display paths.
