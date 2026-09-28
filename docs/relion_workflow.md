# RELION STAR Workflow

This guide shows the recommended cryoROLE 2.0 workflow for RELION STAR
metadata. It focuses on the public commands and the default matching behavior.

## Inputs

cryoROLE compares two independently refined domains from the same particle set.
For RELION inputs, each domain is provided as a STAR metadata file:

```text
ref_domain.star
mov_domain.star
```

RELION `Rot/Tilt/Psi` values are normalized by cryoROLE before relative
orientation computation. Source convention handling is part of the
normalization layer and should not be adjusted downstream.

## Default Matching

By default, `cryorole run` does not assume the two STAR files have the same row
order. It tries safe STAR identity matching in this order:

```text
_rlnTomoParticleName
_rlnImageName / rlnImageName
```

If matching reorders rows or drops unmatched particles, cryoROLE records that in
the run reports and manifest. If matching cannot be resolved safely, prepare
aligned inputs before running. Use `cryorole align --key ...` when the default
STAR identity keys are not sufficient and you need to provide an explicit
matching key.

Default key matching requires at least 50% overlap. Lower overlap requires the
explicit audited `--allow-low-overlap` override; zero overlap and duplicate
particle identities are hard failures before RELION Euler normalization.

## Standard Run

```bash
cryorole preflight --ref ref_domain.star --mov mov_domain.star
cryorole run --ref ref_domain.star --mov mov_domain.star
```

Preflight validates STAR pose columns, resolves the RELION convention, checks
identity uniqueness and overlap, estimates resources, and prints the exact run
command. It has no run-bundle side effects. `cryorole run ... --dry-run` uses
the same service. Review `READY_WITH_WARNINGS`; fix `BLOCKED` inputs before RO
computation.

The default run bundle is written under:

```text
cryorole_outputs/
```

To choose a different run bundle, add `--output-dir RUN_DIR` to `cryorole run`
and use that same directory with downstream `--run-dir` commands.

Important outputs include:

```text
cryorole_outputs/run_manifest.json
cryorole_outputs/run_summary.json
cryorole_outputs/run_report.md
cryorole_outputs/data/raw_landscape.npz
cryorole_outputs/data/raw_landscape.csv
cryorole_outputs/data/match_table.csv
cryorole_outputs/reports/
cryorole_outputs/visualizations/quicklook/  # seven flat PNGs
```

The run quick-look contains Euler/RV triptychs for all particles, `sld_raw >=
1`, and the top 40%, plus the full-data log-SLD distribution. Use `cryorole
visualize` for expanded display products.

## Row-Aligned Inputs

Use `--row-aligned` when you **know** that row `N` in both STAR files is the
same particle:

```bash
cryorole run \
  --ref aligned_ref.star \
  --mov aligned_mov.star \
  --row-aligned
```

`--row-aligned` is your assertion. cryoROLE checks that the row counts match and
that the inputs are valid, then pairs rows by index. It does not compare image
names, coordinates or any other column, because these legitimately change after
signal subtraction or re-extraction. Starting both refinements from the same
STAR does not by itself prove that the final rows still correspond: RELION
jobs, subsets and joins can reorder or drop particles.

If the files were written by `cryorole align`, the run also records where they
came from (see *Alignment lineage* below). That is provenance only and never a
requirement for running.

## Preparing Aligned STAR Files

| Situation | What to do |
| --- | --- |
| You know row *i* in both files is the same particle | `cryorole run --row-aligned` (no `align`) |
| Same particle names, different order or subsets | `cryorole run` (default matching) |
| Names changed (subtraction, re-extraction, moved project) | `cryorole preflight` → review the candidate table → run the printed `align` command → run the printed `run` command |

When default matching fails, `cryorole preflight` evaluates candidate keys on
your files and prints a table with the estimated overlap and ambiguity of each,
plus an explicit `align` command for the best unambiguous one. These are
suggestions: `align` repeats every check on all particles before writing.

```bash
cryorole preflight --ref consensus.star --mov body.star
# ...
#  * _rlnImageName = _rlnImageOriginalName: ~356280 matched (ref 100%, mov 100%)
# Suggested: cryorole align --ref consensus.star --mov body.star --key-pair _rlnImageName=_rlnImageOriginalName
```

### Strategies

All strategies are exact; none uses a distance tolerance unless you ask for one
with `--float-tol`.

- **Same-name keys** (`--key COL ...`): exact values, optionally with
  `--path-mode basename|suffix:N` for path-like columns. `preflight` warns when
  basename matching would merge distinct stacks.
- **Cross-column keys** (`--key-pair REF_COL=MOV_COL`, repeatable). After
  signal subtraction, `_rlnImageOriginalName` of the subtracted particles
  usually holds the `_rlnImageName` of the particles before subtraction.
  It records only the previous processing step.
- **Unchanged coordinates** (`--coordinate-match unchanged`): micrograph plus
  identical coordinates, e.g. subtraction *without* recentring.
- **Recentred re-extraction** (`--coordinate-match recentered-exact
  --recenter-shift X Y Z [--ref-angpix A]`): `--ref` is the STAR given to a
  RELION re-extraction with `--recenter`, and `--mov` is its output. `X Y Z` are
  the `--recenter_x/y/z` values (reference pixels, converted with `--ref_angpix`
  or else the particle pixel size). For every input row cryoROLE predicts the
  output integer coordinates and residual origins exactly as RELION does. It
  accepts a pair only if all of these hold:
  - the coordinates match;
  - the origins match (0.002 Å);
  - the angles are identical;
  - the assignment is one-to-one.

  Rows failing any check go to `unverified_*.star`. Duplicate picks that
  converged to identical metadata cannot be told apart and go to
  `ambiguous_*.star`. If fewer than 99 % of the remaining rows verify, nothing is
  written. `preflight` fills in the vector from the extraction job's `note.txt`
  when it finds one.
- **Chain** (`--via-extraction INPUT OUTPUT --recenter-shift X Y Z`): pairs
  `--ref` and `--mov` through an extraction job in three links:
  1. `--ref` ↔ `INPUT` by particle name;
  2. `INPUT` ↔ `OUTPUT` by recentered-exact;
  3. `OUTPUT` ↔ `--mov` by particle name.

  Only particles verified on every link are paired, and the report lists what
  each link kept.

### Outputs

`cryorole align` writes to `<directory of --ref>/cryorole_alignments/<align-id>/`
(or `--output-dir`):

| File | Contents |
| --- | --- |
| `aligned_ref.star`, `aligned_mov.star` | Matched rows, same order, copied verbatim |
| `ref_only.star`, `mov_only.star` | Rows without a partner |
| `duplicate_*.star`, `ambiguous_*.star`, `unverified_*.star` | Rows excluded because their key was repeated, another candidate was too close, or exact verification failed |
| `match_table.csv` | One row per aligned pair: source row IDs on both sides, plus verification columns for exact modes |
| `ambiguous_groups.csv` | Exact modes only: each group of rows that geometry cannot tell apart, with image names, micrographs, coordinates and half-set (`_rlnRandomSubset`) where present |
| `align_report.json` | Strategy, key, counts, coverage, suspected duplicate groups, warnings, SHA-256 of inputs and outputs, and the next command |

**Coverage.** Every report carries a `coverage` block. `full` means every row of
both files is paired. Otherwise the result is labelled a **matchable subset**, and
`align`, `run_report.md` and the run summary all say so. A landscape from a
matchable subset covers only the paired rows: it is neither the full data nor a
deduplicated dataset.

**Suspected duplicate groups.** In exact modes, rows whose predicted geometry
is identical cannot be paired unambiguously. cryoROLE never guesses a pairing
inside such a group and never removes particles. It lists the groups in
`ambiguous_groups.csv` and reports two different numbers:
- the rows excluded from geometric matching;
- the rows that would be in excess *if* each group is confirmed to be one
  physical particle.

When a unique particle-name key would pair every row, the report says so and
prints the `--key-pair` command for a full paired baseline.

It prints the exact next command:

```bash
cryorole run --ref …/cryorole_alignments/default/aligned_ref.star \
             --mov …/cryorole_alignments/default/aligned_mov.star --row-aligned
```

`cryorole align` is a preparation step. It does not compute RO, SLD,
canonicalization, selection, or export, and it never modifies your STAR files.

### Alignment lineage

When a `--row-aligned` run finds an `align_report.json` next to its inputs, it
attaches the recorded lineage to `run_summary.json` and `run_manifest.json`
(`alignment_provenance`). It does so only after verifying three things:
- the hashes of the two aligned files;
- the hash of `match_table.csv`;
- that the match table maps one aligned row to one source row on each side.

The current original files are reported separately as `verified`,
`unavailable` or `mismatch`. If the report cannot be verified, the run prints
`alignment provenance not attached: <reason>` and continues.

## Signal Subtraction With Recentring Leaves Stale Coordinates

When RELION's Particle subtraction recentres the box on a 3D point
(`--center_x/y/z`), it writes the new offset into `_rlnOriginX/YAngst` but does
not update `_rlnCoordinateX/Y` (confirmed in `src/particle_subtractor.cpp`). The
subtracted STAR therefore places each particle away from the true box centre by
the projected recentring shift. In our test data that is a median of 92 px
(77 Å) and up to 143 px.

Refinement and classification of the subtracted particles are not affected,
because they use the images and origins. Any step that uses micrograph
coordinates is affected:
- re-extraction from the subtracted STAR;
- Bayesian polishing;
- distance-based duplicate removal;
- coordinate-based matching (including `cryorole align`).

`cryorole preflight` warns (`STALE_SUBTRACTION_COORDINATES`) when an input sits
in a Subtract job whose `note.txt` shows recentring. To write a copy with
corrected coordinates:

```bash
cryorole align --fix-subtract-coordinates Subtract/job041/
# optionally transfer the correction to a later refinement of the subtracted particles:
cryorole align --fix-subtract-coordinates Subtract/job041/ --apply-to Refine3D/job042/run_data.star
```

cryoROLE reads the centre and the subtraction input from the job's `note.txt`
(or `--center`, `--subtract-input`, `--model-angpix`). It recomputes the shift
RELION applied, rounding in particle pixels as RELION does. It then checks that
this reproduces the subtracted origins and angles for at least 99 % of
particles. Only after that does it write `…_coords_corrected.star` under
`cryorole_alignments/`, with a `coordinate_correction_report.json`. The original
files are never changed.

Safeguards:

- **Only coordinate fields change.** Only `_rlnCoordinateX` and
  `_rlnCoordinateY` are rewritten, formatted to 6 decimals. Every other byte of
  each written row is unchanged, including the spacing between fields. One
  comment line naming the report is added at the top.
- **Rows are checked before correction.** A row is corrected only if its
  origins and angles verify and its coordinates still hold the stale value that
  RELION writes (equal to the subtraction input). Rows that are already
  corrected or hold unexpected coordinates are left out and listed in
  `excluded_rows.star`.
- **A second application is refused.** Running the correction on an
  already-corrected file stops with "already coordinate-corrected" and writes
  nothing.
- **`--apply-to` transfers by identity.** It joins the target rows to verified
  subtracted rows by `_rlnImageName`, which must be unique. It copies the
  corrected value computed from the subtraction input; the target's own angles
  and origins are never used and are left unchanged. A target that already holds
  the corrected coordinates is refused. Rows without a verified partner, or with
  coordinates that are neither stale nor corrected, are left out and listed.

## Canonicalize

```bash
cryorole canonicalize --run-dir cryorole_outputs
```

Canonicalization writes a separate canonical landscape under:

```text
cryorole_outputs/canonical/default/
```

It does not overwrite the raw landscape.

## Visualize

```bash
cryorole visualize --run-dir cryorole_outputs --space canonical
```

Visualization outputs are display-only. Display filters, ranges, and color
scales do not alter raw/canonical landscapes and do not create selections.
The default is two PNG triptychs for `sld_display >= 1`; 1D and 3D views are
explicit through `--view`.

## Select

To choose a center interactively from linked Euler and rotation-vector views:

```bash
cryorole explore --run-dir cryorole_outputs --space canonical
```

Display filters and display sampling affect only the browser view. Evaluate and
Confirm use the full parent landscape and the Python SO(3) evaluator. Confirm
writes a normal selection; closing the browser before Confirm writes none.

```bash
cryorole select \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --space canonical \
  -c 13 0 14 \
  -r 6
```

The default radius selection uses SO(3) geodesic distance when the center is
provided in Euler degrees.
The center values are usually chosen after inspecting the raw or canonical
visualizations.

## Export

```bash
cryorole export \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --domain both
```

Export uses the recorded `ref_source_row_id` and `mov_source_row_id`
provenance. It subsets the original source metadata and does not rewrite source
poses with canonical or display coordinates.

The source STAR identity is also content-addressed. If a STAR file has moved,
use `--relocated-ref` or `--relocated-mov`; export continues only when the
streamed SHA-256 matches the run record. Changed or replaced STAR files are
rejected.

At any time, `cryorole status --run-dir cryorole_outputs` reports actual bundle
and artifact integrity, while `cryorole next --run-dir cryorole_outputs` prints
exact required/recommended/optional next commands.
