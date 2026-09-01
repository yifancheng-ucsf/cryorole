# Workflow UX Tutorial

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

## 3. Run and inspect actual state

```bash
cryorole run --ref ref_domain.star --mov mov_domain.star --output-dir my_run
cryorole status --run-dir my_run
cryorole next --run-dir my_run
```

Run revalidates source identity and publishes the bundle transactionally.
`status` derives state from the manifest, bundle marker, raw NPZ, source hashes,
canonical frames, visualizations, selections, and exports. `next` produces
exact commands labelled required, recommended, or optional. It does not treat
a visualization as a selection, and it does not require canonicalization.

## 4. Raw and canonical coordinates

The raw landscape is the direct relative-orientation result. A canonical
landscape is an optional, separately stored motion-aligned coordinate frame:

```bash
cryorole canonicalize --run-dir my_run
```

Canonicalization helps interpretation; it does not replace raw facts or write
display coordinates back to STAR/CS poses.

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

## 7. Export without reselecting

```bash
cryorole export --run-dir my_run --selection-id state_1 --domain both
```

Export consumes the confirmed Selection and recorded source-row provenance. It
verifies each source hash, does not rematch or reselect, does not modify source
files, and does not write canonical/display coordinates as physical poses.

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
