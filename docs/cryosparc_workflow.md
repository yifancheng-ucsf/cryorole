# CryoSPARC CS Workflow

This guide shows the recommended cryoROLE 2.0 workflow for CryoSPARC `.cs`
metadata.

## Which output should I use?

Choose the **particle output containing 3D alignments** for each separately
refined domain. A local refinement, or another refinement aligned to the chosen
domain, can provide suitable poses; the actual requirement is `uid` plus
`alignments3D/pose`, not a particular job label. A picking job or 2D result without
3D poses is insufficient. Verify which domain was aligned by inspecting its map
and processing history.

In the GUI, **Job folder...** finds particle result groups with an `alignments3D`
slot and offers the `.cs` named by that slot. Without a usable group candidate,
it offers matching `*_particles.cs` files excluding passthrough files. If there
are several candidates, choose the intended output; cryoROLE does not guess the
latest iteration. **File...** lets you choose the pose `.cs` directly.

Keep the job's `.csg` and its associated passthrough files available for later
result-group export. A passthrough file alone is not the pose input unless it
actually contains the required pose field. Two domains with independently picked
particles need a validated association before this workflow: `uid` matching
identifies existing shared particles, not biological partners inferred by proximity.

For a GUI walkthrough, see the [README](../README.md#6-your-own-data).

## Inputs

cryoROLE 2.0 treats `.cs` files as native inputs:

```text
ref_domain.cs
mov_domain.cs
```

The default CryoSPARC pose field is:

```text
alignments3D/pose
```

The pose is an axis-angle rotation vector in radians. CryoSPARC defines it as
the rotation applied during image back-projection, which is the inverse
(transpose) of the reference-to-image rotation that RELION's
`_rlnAngleRot/Tilt/Psi` describe. cryoROLE therefore converts it with the same
mapping as pyem `csparc2star.py`:

```text
active_matrix = Rotation.from_rotvec(pose).as_matrix().T
```

A `.cs` input and the same particles converted to STAR with `csparc2star.py`
give identical relative orientations; a golden test against real
`csparc2star.py` output locks this. The rule is recorded as the CryoSPARC `ConventionPolicy` in
`run_manifest.json` and in the preflight report. Rotation vectors of any norm,
including the ones above pi that CryoSPARC sometimes writes, are valid.

Runs made with cryoROLE before this rule was introduced used
`Rotation.from_rotvec(pose).as_matrix()` without the transpose. They have the
same RO angles but different RO axes, so re-run any `.cs`-based analysis.

## Default Matching

CryoSPARC `.cs` inputs default to `uid` matching:

```bash
cryorole run --ref ref_domain.cs --mov mov_domain.cs
```

If matching reorders rows or drops unmatched particles, cryoROLE records that in
the run reports and manifest. Downstream selections and exports use the recorded
source-row provenance from the run bundle.

Public uid matching requires at least 50% overlap unless
`--allow-low-overlap` is explicit. Zero matches and duplicate uid values fail
before pose normalization.

## Row-Aligned Inputs

Use `--row-aligned` only when the two `.cs` files have already been prepared so
that row `N` in both files represents the same particle:

```bash
cryorole run \
  --ref aligned_ref.cs \
  --mov aligned_mov.cs \
  --row-aligned
```

In row-aligned mode, cryoROLE requires equal row counts and does not key-match,
reorder, or drop particles. No column (including `uid`) is compared: the row
order is the user's assertion.

## Standard Workflow

Preflight and run:

```bash
cryorole preflight --ref ref_domain.cs --mov mov_domain.cs
cryorole run --ref ref_domain.cs --mov mov_domain.cs
```

Preflight validates `uid`, the `alignments3D/pose` vector shape and finite
values, matching coverage, and resources with no run-bundle side effects. The
equivalent check is `cryorole run ... --dry-run`.

Canonicalize:

```bash
cryorole canonicalize --run-dir cryorole_outputs
```

Visualize:

```bash
cryorole visualize --run-dir cryorole_outputs --space canonical
```

The default is two PNG triptychs for `sld_display >= 1`; use `--view 2d,1d`
or `--view 3d` only when those additional display products are needed.

Or explore linked Euler/RV projections locally and create a standard selection
only after explicit Confirm:

```bash
cryorole explore --run-dir cryorole_outputs --space canonical
```

Select:

```bash
cryorole select \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --space canonical \
  -c 13 0 14 \
  -r 6
```

The center values are usually chosen after inspecting the raw or canonical
visualizations.

Export:

```bash
cryorole export \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --domain both
```

## Select by CryoSPARC metadata

Use the run-time source field directly, without converting CS to STAR:

```bash
cryorole select --run-dir cryorole_outputs --selection-id ref_classes_01 --mode metadata --metadata-domain ref --metadata-column alignments3D/class --metadata-value 0,1
cryorole select --run-dir cryorole_outputs --selection-id ref_class --mode metadata --metadata-domain ref --metadata-column alignments3D/class --split-by-value
cryorole export --run-dir cryorole_outputs --selection-id ref_classes_01
```

The field must exist in the chosen source; class values are not renumbered.
Choose ref or mov explicitly. Selection uses the full parent landscape's
recorded source-row indices, including when run reordered or dropped unmatched
input rows. Source content must match its recorded SHA-256.

Supported CS fields are scalar integer (including uint64), boolean, and UTF-8
text. Text comparison is exact; empty strings are excluded and counted as
missing. Float/vector fields and invalid UTF-8 are rejected. Split writes one
standard selection per non-missing matched value, with at most 100 groups.
Name collisions fail before writing even with overwrite. See the
[CLI reference](cli_reference.md#cryosparc-metadata-fields) for value syntax.

## Outputs

The run bundle is written under:

```text
cryorole_outputs/
```

To choose a different run bundle, add `--output-dir RUN_DIR` to `cryorole run`
and use that same directory with downstream `--run-dir` commands.

Important outputs include:

```text
run_manifest.json
run_summary.json
run_report.md
data/raw_landscape.npz
data/raw_landscape.csv
data/match_table.csv
reports/
visualizations/quicklook/  # seven flat PNGs
canonical/
selections/
exports/
```

The run quick-look contains Euler/RV triptychs for all particles, `sld_raw >=
1`, and the top 40%, plus the full-data log-SLD distribution. Use `cryorole
visualize` for expanded display artifacts.

See `docs/output_files.md` for details.

`cryorole status --run-dir cryorole_outputs` derives progress and integrity from
the actual artifacts. `cryorole next --run-dir cryorole_outputs` prints exact
required, recommended, and optional commands; canonicalization remains optional.

## Export Notes

CryoSPARC export preserves the structured-array dtype and vector-valued fields
of the source metadata when writing selected subsets. Export does not reselect,
rematch, or write canonical/display coordinates back as physical source poses.

Export first verifies the source `.cs` SHA-256 recorded by `run`. A moved file
must be supplied explicitly with `--relocated-ref` / `--relocated-mov` and must
match byte-for-byte; a changed or replaced `.cs` file is rejected.

### Bringing a selection back into CryoSPARC

```bash
cryorole export --run-dir cryorole_outputs --selection-id high_sld --domain ref \
  --cryosparc-group-ref /path/to/project/J75/J75_particles.csg
```

`--cryosparc-group-ref` is only needed when the run used a copy of the `.cs`
file kept under its original name (the group names its files by name); when
the run read the file from the job folder, the job's `.csg` is found
automatically. The group and every file it names are checked before anything
is written; if one is missing (for example a copy without its passthrough
file), the export stops and writes nothing, and `--no-cryosparc-group` writes
`selected_ref.cs` alone. The export writes `selected_ref.cs`, a subset of the job's
passthrough file, and `selected_ref.csg`, and prints how to finish:

1. Copy the three files into a folder inside the project directory, for example
   `<project>/exports/groups/cryorole_J75_high_sld_ref/`. Import Result Group
   rejects files outside the project. If the project belongs to the CryoSPARC
   user, copy as that user (`sudo -u cryosparc mkdir -p …` and `sudo -u cryosparc cp …`).
2. In CryoSPARC, run **Import Result Group** with the path of `selected_ref.csg`
   in that folder. The particles keep their original poses, CTF and half-set
   split, so a job such as Homogeneous Reconstruction Only can use them directly.

## Difference from 0.x Workflows

Older cryoROLE workflows often converted CryoSPARC `.cs` files to STAR before
analysis. cryoROLE 2.0 can read `.cs` inputs directly, so conversion is not the
default path.
