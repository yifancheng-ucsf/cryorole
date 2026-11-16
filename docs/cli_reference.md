# CLI Reference

Use `cryorole COMMAND --help` as the executable source of truth; most commands
end their help with copyable examples. This page summarizes the public Workflow
UX commands. The same operations are available from Python through
`cryorole.api` (see `docs/python_api.md`).

For complete argument lists, see [All command options](cli_options.md), generated from this release's help output. For worked plotting and selection examples, see [Command-line workflow](command_workflow.md).

## Common behaviour

`cryorole --version` prints the installed version.

**Errors.** Expected problems print one line, `cryorole: error: MESSAGE`,
followed where possible by `  what to do: REMEDY`, and exit with status 2.
Unexpected internal errors print a short report and exit with status 1; set
`CRYOROLE_DEBUG=1` to see the Python traceback instead. Ctrl-C exits with
status 130. `preflight` keeps its own readiness codes (below), and `status` /
`next` exit with 3 for an unusable bundle.

**Omitted `--run-dir`.** Commands that read a run bundle (`canonicalize`,
`visualize`, `select`, `export`, `explore`, `status`, `next`) use, in order:
the current directory if it is a run bundle, then `./cryorole_outputs` (the
default `run --output-dir`). Otherwise they stop and list the bundles found in
the current directory; they never pick one of several.

**Omitted `--canonical-id`** (`visualize`, `select`, `explore` with
`--space canonical`) uses the bundle's only canonical frame. With no frame the
command asks you to run `canonicalize` first; with several it lists them.
`canonicalize --canonical-id` still names the frame being *written* (default
`default`).

**Omitted `--selection-id`** (`export`, and `visualize --use-selected-landscape`)
uses the bundle's only selection, with the same rules. `select` always requires
an explicit `--selection-id`: a selection is a decision you name.

Every value filled in this way is printed to stderr (for example
`[cryorole] using --run-dir cryorole_outputs (./cryorole_outputs is the default
run output)`) and recorded under `resolved_by` in the report the command writes.
A raw-space `select` in a bundle that has canonical frames prints a note, unless
`--space raw` was given explicitly.

Older option spellings (`canonicalize --fit-top-fraction`, `animate
--sld-threshold`, `visualize --alpha`, …) remain accepted as hidden aliases.

## Preflight

```bash
cryorole preflight --ref REF --mov MOV [--output-dir RUN] [--row-aligned]
                   [--allow-low-overlap] [--sld-metric METRIC]
                   [--no-visualize] [--json [PATH]]
```

Performs source identity, convention, pose-schema, identity/matching, and
resource checks without creating `RUN`. `--json` with no path writes JSON to
stdout. Exit codes are 0 ready, 1 ready with warnings, and 2 blocked.

Preflight also checks the machine (`environment` in the JSON report): the
cryoROLE, Python, NumPy, SciPy, pandas and Matplotlib versions; whether the
output location can be written (blocking if not); whether `RUN` already exists
(a warning, because `run` then needs `--overwrite`); free disk against the
estimated bundle size (blocking below a 20 % margin); and available memory
against the estimated peak (a warning).

Selecting the same file for `--ref` and `--mov` (same path, or a copy with the
same SHA-256) is blocked before any matching. A formal `run` records the same
finding as a strong warning and always reports the RO-angle summary and any
input-sanity warnings (see `docs/output_files.md`, "Input sanity").

`--row-aligned` is the user's assertion that row N is the same particle in both
files. cryoROLE checks only that the row counts are equal and then pairs rows
by index. It deliberately does not compare image names, coordinates, or any
other column, because these can legitimately differ between the two
refinements (for example after RELION signal subtraction or re-extraction).

RELION STAR inputs: the particle table is the loop in the `data_particles`
block (RELION 3.1+). A file with more than one `data_particles` loop is
rejected rather than merged. `run`, metadata selection, and export all use this
same rule.

## Run dry-run

```bash
cryorole run --ref REF --mov MOV --dry-run [--json [PATH]]
```

Uses the same preflight service. `run --json` is valid only with `--dry-run`.

A formal run writes seven flat PNGs under `RUN/visualizations/quicklook/`:
Euler/RV triptychs for all particles, `sld_raw >= 1`, and the top 40%, plus the
full-data log-SLD distribution. `--no-visualize` skips them; every run still
writes concise `RUN/run_report.md`. Expanded display controls are owned by
`cryorole visualize`, not `run`.

## Align

```bash
cryorole align --ref REF --mov MOV [--key COL ... | --key-pair REF_COL=MOV_COL ...]
               [--coordinate-match unchanged|recentered-exact --recenter-shift X Y Z [--ref-angpix A]]
               [--via-extraction INPUT OUTPUT --recenter-shift X Y Z]
               [--output-dir DIR] [--align-id ID] [--overwrite]
cryorole align --fix-subtract-coordinates SUBTRACT_JOB [--apply-to STAR]
               [--subtract-input STAR] [--center X Y Z] [--model-angpix A]
```

Establishes particle correspondence when default matching fails, and writes
verbatim aligned STAR files with `match_table.csv` and `align_report.json`. The
default location is `<directory of --ref>/cryorole_alignments/<align-id>/`. It
prints the exact `cryorole run … --row-aligned` command to run next.
`--fix-subtract-coordinates` writes a verified copy of a recentred RELION
subtraction with corrected coordinates. See `docs/relion_workflow.md`.

When default STAR matching fails, `preflight` adds an `align_diagnosis` block
and prints a candidate table and a suggested `align` command. These are
suggestions only.

## Status and next

```bash
cryorole status [--run-dir RUN] [--json [PATH]]
cryorole next [--run-dir RUN] [--json [PATH]]
```

Status reads real artifacts and integrity evidence. Next returns exact
required/recommended/optional commands derived from that status. Both exit 0
for a completed (or legacy) bundle and 3 when the bundle is missing, failed or
incomplete, so scripts can test the result; the output is printed either way.

## Explore

```bash
cryorole explore --run-dir RUN [--space raw|canonical]
                 [--canonical-id ID] [--selection-id ID]
                 [--threshold SLD | --top-fraction FRACTION]
                 [--max-display-points N] [--colormap NAME]
                 [--port PORT] [--no-open]
```

Starts an offline loopback-only explorer. `--selection-id` overlays an existing
selection. Threshold/top fraction and point limits are display-only. The
default accessible colormap is `viridis`; legacy `rainbow_r` remains available.
The process serves until the page requests shutdown or Ctrl+C is used.

## GUI

```bash
cryorole gui [--run-dir RUN] [--no-gl]
```

Opens the desktop interface; `--run-dir` opens that run bundle in the View page
(otherwise use File > Open run, or make a new run on the Inputs, Parameters and
Run pages). It needs the optional extra:
`python -m pip install ".[gui]"` from the cryorole checkout (PySide6, pyqtgraph, PyOpenGL). On Linux, Qt 6.5
and later also need the system library `libxcb-cursor0`.

Before the window opens, a short check in a separate process tests that Qt can
open a window and that OpenGL works for the 3D view:

- if Qt cannot open a window, the command stops with the reason and what to
  install or change;
- if OpenGL is unusable, the GUI opens without the 3D view; the 2D and 1D views
  do not need OpenGL;
- with software OpenGL (for example a remote desktop), the 3D view works but is
  slow at large point counts.

On a workstation whose monitor is driven by an onboard chip while the NVIDIA
GPUs are only secondary devices, start the GUI with NVIDIA PRIME render offload
so the 3D view uses the GPU:

```bash
__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia cryorole gui
```

`--no-gl` skips the OpenGL check and the 3D view.

The View page shows the raw landscape or a canonical frame as three 2D
projections and 1D distributions, with Euler or rotation-vector coordinates;
"3D view..." opens the 3D points in their own window, which follows the same
settings.
Its defaults match `cryorole visualize`:

- display filter `sld_display >= 1`; alternatives are the top N % by SLD (like
  `--top-fraction`) or all particles (like `--all`);
- standard axis ranges (Euler ±180°, rotation vector ±3.14 rad) with equal axis
  units; alternatives are fitting the central 99 % of the displayed particles,
  fitting all of them, or custom ranges per axis (like `--axis-limit`). The
  chosen range applies to the 3D, 2D and 1D views;
- colour scale `vmax = min(max SLD, 100)` unless an SLD colour max is set (like
  `--vmax`), legacy colours.

The 3D view draws a deterministic sample of at most 100,000 points by default
(seed 0, adjustable as "3D max points"); the 2D projections and 1D distributions
always use every particle that passes the display filter. All of these controls
are display-only: they change no landscape, selection or export, and an axis
range never filters particles.

The Select page adds a panel next to the views. In radius mode, a click in the
2D or 3D view sets the centre to the nearest displayed particle (or type the
centre), and the exact count over the full landscape updates as the centre or
radius changes, using the same SO(3) evaluator as `cryorole select` in the space
shown. Selected particles are highlighted. The other modes are those of
`cryorole select`: SLD threshold (`sld_raw`), coordinate range, random fraction
and source metadata values. Confirm writes a standard selection under
`selections/NAME`, never overwrites one, and adds the equivalent `cryorole
select` command to the Commands panel. Display filters and the 3D sample never
change a selection.

The Inputs page takes the reference and moving body, each a `.cs` / `.star`
file or a job folder. In a CryoSPARC job folder the particle file named by the
job's particle result group (`.csg`) is used; in a RELION job folder,
`run_data.star`. When a folder holds several candidates, the GUI lists them and
you choose; it never guesses (for example between iterations). Both inputs must
come from the same software. The page then runs the same checks as
`cryorole preflight` and shows READY / READY_WITH_WARNINGS / BLOCKED with the
matching summary and the memory estimate. The Parameters page sets the output
folder and the `cryorole run` options (SLD metric, *k*, quick-look figures, the
raw CSV table), and shows the equivalent command. The Run page runs
`cryorole run` in a separate process with stage progress and Cancel (a cancelled
run publishes nothing); the new bundle then opens in the View page. The GUI never
overwrites a run. "Canonicalize..." on the View page runs
`cryorole canonicalize` and opens the new frame.

The Export page runs `cryorole export` for a selection and domain. For
CryoSPARC sources it shows the result group (`.csg`) found next to each source
`.cs`, or lets you choose the source job's particle group; after the export it
shows the files written and the steps to import them (copy commands and the
Import Result Group path). It never overwrites an existing export folder.

The desktop View page also offers **Save figures...** and a read-only saved
Selection overlay. Saving captures the current display settings and uses the
standard visualization renderer; it never overwrites an existing figure name.
Overlays do not create or confirm a Selection. See [Desktop GUI](gui.md) for
the complete workflow, installation and acceptance limits.

## Guide

```bash
cryorole guide --ref REF --mov MOV [--output-dir RUN]
               [--row-aligned] [--allow-low-overlap]
               [--non-interactive] [--execute-run] [--json [PATH]]

cryorole guide --run-dir RUN [--non-interactive] [--json [PATH]]
```

New mode wraps preflight; resume mode wraps status/next. `--execute-run` is an
explicit run-only action and never creates a selection or export.

## Existing scientific commands

```bash
cryorole run --ref REF --mov MOV [--output-dir RUN] [--overwrite]
cryorole canonicalize --run-dir RUN [--canonical-id ID] [--overwrite] [--write-csv] [--projection-panels]
cryorole visualize --run-dir RUN --space raw|canonical
# default: two PNG 3-view projections for sld_display >= 1
cryorole visualize --run-dir RUN --view 2d,1d
cryorole visualize --run-dir RUN --view 3d
cryorole select --run-dir RUN --selection-id ID --space SPACE -c A B C -r DEG
cryorole export --run-dir RUN --selection-id ID --domain ref|mov|both
```

`run --overwrite` replaces an existing bundle only after the new run has
completed and validated; `canonicalize --overwrite` replaces
`canonical/<canonical-id>/` only.

**Figure speed options** (`run` and `canonicalize`):

- `--jobs auto|N` renders independent figure sets in separate processes
  (`run`: the three quick-look subsets; `canonicalize`: the three display
  sets). `auto` (default) uses up to 3, limited by CPU cores and free memory,
  and stays at 1 below 100k particles. Output files are identical to `--jobs 1`;
  each extra process needs memory (about 200 MiB + 1 KB per plotted particle).
- `--point-style density` draws 2D projections as an image with one cell per
  screen pixel holding the SLD value that would be visible there, instead of one
  marker per point (`scatter`, default). About 15x faster at 500k points. The
  picture is near-identical; sparse regions look slightly thinner because a
  scatter marker (~1.4 px) also tints neighbouring pixels. 3D views keep markers.

`export` refuses to write a subset when the source particle table no longer has
the row count recorded by the run.

**Back into CryoSPARC.** For a CryoSPARC source, export also mirrors the source
job's own result group, so the subset can be brought back with **Import Result
Group**. It uses the job's particle `.csg` file (e.g. `J75/J75_particles.csg`):
`--cryosparc-group-ref` / `--cryosparc-group-mov PATH`, or the only `.csg` next
to the source `.cs` that names it. Every file the group names is subset: the
source `.cs` by the recorded source rows, the others (for example
`J75_passthrough_particles.cs`, which holds the particle images) by `uid` in the
same order; a missing or duplicated `uid` is an error. The new group,
`selected_<domain>.csg`, names those subsets. When no group is found (for
example a run on a copied `.cs` without its job folder), only the `.cs` subset
is written and the export says which `--cryosparc-group-*` file to pass. The
group is checked before the first file is written: when it cannot be mirrored
(several candidates, a missing passthrough file, a different row count), the
export stops without writing anything. `--no-cryosparc-group` writes the `.cs`
subsets alone.

Import Result Group only accepts files inside the CryoSPARC project directory.
The export report, `NEXT_STEPS.txt` and the command's stderr give the copy
commands (as the project's owner when you cannot write there, e.g.
`sudo -u cryosparc`) and the path to import. The import copies the files into
its own job folder.

See `docs/output_files.md` for the artifact layout and the scientific and
display distinctions.


## Visualization terminology

SLD color bars use `SLD`; reports retain `sld_raw` or `sld_display` as the
actual source field. Use `--opacity 0.4` for point opacity (0 transparent,
1 opaque); omission keeps the style default. The old `--alpha` remains a
hidden compatibility alias. Supplying both names is an error.

## Selection modes and names

`cryorole select --help` groups all controls by mode and provides complete
examples. Every command requires `--selection-id NAME` (and `--run-dir RUN`
unless it can be resolved as described under Common behaviour); choose a
name such as `region_01` or `high_sld`. Names cannot contain path separators or redirect output outside selections/.
No default is assigned. Reusing a name
requires explicit `--overwrite`; choose another name to keep both selections.

| Mode | Required mode inputs | Optional controls and interpretation |
| --- | --- | --- |
| `radius` | `--center A B C` and one of `--radius DEG`, `--radius-rad RAD` | `--center-representation euler|rotvec`, `--metric so3|rotvec`. Euler center: degrees; rotvec center: radians. Default metric: SO(3) geodesic. |
| `threshold` | `--sld-min`, `--sld-max`, or both | Inclusive `sld_raw` bounds; independent of display filters. |
| `range` | Repeatable `--range-bound AXIS:LOWER:UPPER` | Euler alpha/beta/gamma in degrees OR rotvec x/y/z in radians. Inclusive bounds, intersected across axes; blank sides are open. Euler reversed bounds wrap across the seam; rotvec reversed bounds are invalid. Last bound wins for a repeated axis. |
| `random` | `--fraction F`, with `0 < F <= 1` | Optional non-negative `--seed`. When omitted, a seed is generated; either way the seed used is recorded (`random_seed`, with `random_seed_source` `user` or `generated`), so rerunning with `--seed` reproduces the selection. Samples `ceil(F*N)` rows without replacement from all parent rows. |
| `metadata` | `--metadata-domain ref|mov`, `--metadata-column`, and one of `--metadata-value VALUES`, `--split-by-value` | Comma-separated values form a union. Uses run-time source rows; split writes standard child selections named `NAME_VALUE` with filename-safe components. |

`--space raw|canonical` chooses the parent coordinate space; `--canonical-id`
applies only to canonical space. Euler interpretation follows the recorded
parent convention. Options explicitly supplied for another mode are errors.

```bash
cryorole select --run-dir RUN --selection-id region_01 --mode radius --center 0 0 0 --radius 15
cryorole select --run-dir RUN --selection-id high_sld --mode threshold --sld-min 2
cryorole select --run-dir RUN --selection-id alpha_window --mode range --range-bound alpha:-20:20
cryorole select --run-dir RUN --selection-id sample_10pct --mode random --fraction 0.1 --seed 7
cryorole select --run-dir RUN --selection-id ref_classes --mode metadata --metadata-domain ref --metadata-column rlnClassNumber --metadata-value 1,3
cryorole select --run-dir RUN --selection-id ref_classes --mode metadata --metadata-domain ref --metadata-column rlnClassNumber --split-by-value
```

`--write-selected-landscape` writes an additional landscape for visualization.
`--recompute-sld` requires that option and preserves parent SLD separately.
Selection always evaluates the full parent, not the plotted sample. Success
guidance appears on stderr; stdout remains the output directory.

### CryoSPARC metadata fields

Metadata mode accepts native run-time CS sources as well as STAR. For a CS
source containing `alignments3D/class`:

```bash
cryorole select --run-dir RUN --selection-id cs_classes --mode metadata --metadata-domain ref --metadata-column alignments3D/class --metadata-value 0,1
cryorole select --run-dir RUN --selection-id cs_class --mode metadata --metadata-domain mov --metadata-column alignments3D/class --split-by-value
```

CS fields must be scalar integer, boolean, or text. Integers use exact decimal
matching without float conversion; booleans accept true/false/1/0. Byte strings
are decoded as strict UTF-8 and text matches exactly (including spaces and
leading zeros). Empty strings are excluded and counted as missing. Unsupported
float/vector/nested/object fields, invalid byte strings, invalid row indices,
and missing source hashes fail before writing. Commas separate values; no
comma-escaping syntax is provided.

CS split has a 100-group limit; use `--metadata-value` for explicit selections
if it is exceeded. All split child names are checked for case-insensitive
sanitization collisions and existing outputs before writing. STAR value
normalization and its existing uncapped group behavior remain unchanged.
Source hashes are verified when recorded; legacy STAR compatibility is
retained and recorded as unverified. CS requires a verified source identity.
