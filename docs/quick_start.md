# cryoROLE 2.0 Quick Start

This guide shows the smallest practical cryoROLE 2.0 workflow. For the
artifact layout, see `docs/output_files.md`; for every command option, see
`docs/cli_reference.md`.

## Install

From a local checkout:

```bash
python -m pip install -e ".[test]"
cryorole --help
```

The installed command is `cryorole`.

## RELION STAR Workflow

Preflight first, then run the relative-orientation analysis:

```bash
cryorole preflight --ref ref_domain.star --mov mov_domain.star
cryorole run --ref ref_domain.star --mov mov_domain.star
```

Preflight uses the same matching and pose-schema service as the real run, but
does not create a run bundle or modify the inputs. `READY_WITH_WARNINGS` is a
review gate (exit code 1); `BLOCKED` is a hard stop (exit code 2). The equivalent
run-shaped check is `cryorole run ... --dry-run`.

By default, the run bundle is written under `cryorole_outputs/`. The run bundle
contains raw numeric arrays, CSV tables, reports, manifests, and default
quick-look visualizations. The quick-look directory contains Euler/RV
triptychs for all particles, `sld_raw >= 1`, and the top 40%, plus one
full-data log-SLD distribution. `run_report.md` explains these files; use
explicit `cryorole visualize` for expanded display products.
Commands write provenance into the run bundle so selections and exports can be
audited later.

The public matcher fails below 50% overlap by default. If a scientifically
reviewed partial join must continue, add `--allow-low-overlap`; the override and
all coverage fractions are recorded. Zero matches and duplicate identities
always fail before pose normalization.

Canonicalize the raw landscape:

```bash
cryorole canonicalize --run-dir cryorole_outputs
```

Visualize the canonical landscape:

```bash
cryorole visualize --run-dir cryorole_outputs --space canonical
```

This defaults to two PNG three-view projections for `sld_display >= 1`. Add
`--view 2d,1d` for full-filtered-row coordinate distributions or `--view 3d`
for the self-contained offline display-only viewer.

Select particles around a canonical Euler center:

```bash
cryorole select \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --space canonical \
  -c 13 0 14 \
  -r 6
```

The center values can be chosen from static visualizations or from the offline
interactive explorer:

```bash
cryorole explore --run-dir cryorole_outputs --space canonical
```

Clicking and Evaluate create only a draft. Confirm asks for a selection ID and
writes the same standard selection artifact consumed by export. A display
filter or display downsample never limits the exact scientific candidate set.

Export the selected source metadata:

```bash
cryorole export \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --domain both
```

Export verifies ref/mov SHA-256 identities recorded by `run`. If an input moved,
provide `--relocated-ref NEW_REF` or `--relocated-mov NEW_MOV`; matching content
is required. For an old pre-hash bundle only, the advanced
`--allow-unverified-source` override is available and appears in the export
report.

## CryoSPARC CS Workflow

CryoSPARC `.cs` files are native inputs in cryoROLE 2.0:

```bash
cryorole preflight --ref ref_domain.cs --mov mov_domain.cs
cryorole run --ref ref_domain.cs --mov mov_domain.cs
cryorole status --run-dir cryorole_outputs
cryorole canonicalize --run-dir cryorole_outputs
cryorole explore --run-dir cryorole_outputs --space canonical
cryorole select --run-dir cryorole_outputs --selection-id state_1 --space canonical -c 13 0 14 -r 6
cryorole export --run-dir cryorole_outputs --selection-id state_1 --domain both
```

Default `.cs` matching uses `uid`.

## Choosing a Run Directory

By default, `cryorole run` writes to `cryorole_outputs/`. To choose a different
run bundle name, use `--output-dir`:

```bash
cryorole run --ref ref_domain.star --mov mov_domain.star --output-dir my_run
```

Downstream commands then use that directory through `--run-dir my_run`.

If you have multiple runs, keep each run bundle in a separate directory and use
the matching `--run-dir` for canonicalization, visualization, selection, and
export.

## Matching and Row Alignment

Use default matching when possible:

- CryoSPARC `.cs`: `uid`
- RELION STAR: `_rlnTomoParticleName`, then `_rlnImageName` / `rlnImageName`
  when safe

Use `--row-aligned` when you know that row `N` in both files is the same particle:

```bash
cryorole run --ref aligned_ref.star --mov aligned_mov.star --row-aligned
```

In row-aligned mode, cryoROLE pairs row `N` with row `N` and requires equal row
counts; it does not compare names or coordinates.

If default matching fails (for example after signal subtraction or
re-extraction changed the particle names), `cryorole preflight` lists candidate
keys and prints an explicit `cryorole align` command. `align` writes verified,
verbatim aligned files and prints the `run` command to use next:

```bash
cryorole preflight --ref ref_domain.star --mov mov_domain.star
cryorole align --ref ref_domain.star --mov mov_domain.star --key-pair _rlnImageName=_rlnImageOriginalName
cryorole run \
  --ref cryorole_alignments/default/aligned_ref.star \
  --mov cryorole_alignments/default/aligned_mov.star \
  --row-aligned
```

See `docs/relion_workflow.md` for all strategies, including recentred
re-extraction and the correction of stale coordinates after RELION subtraction
with recentring. Manual pre-alignment is also acceptable when the row-order
assertion is true and auditable.

## What the Main Artifacts Mean

- Raw landscape: direct relative-orientation result from matched input
  particles.
- Canonical landscape: a derived coordinate frame for easier comparison and
  inspection.
- Visualization: display-only figures and display tables.
- Selection: an explicit scientific decision artifact.
- Export: a source-metadata subset created from a selection.

Display filtering, such as plotting only high-density points, is not a
scientific selection. Use `cryorole select` when you intend to create a particle
set for export or reconstruction.

Use artifact-derived guidance at any point:

```bash
cryorole status --run-dir cryorole_outputs
cryorole next --run-dir cryorole_outputs
```

`next` labels commands required, recommended, or optional. Canonicalization is
optional and visualization never counts as a selection. `cryorole guide` wraps
the same preflight/status/next services; non-interactive mode never waits for
input or makes a scientific choice.

## Offline Animation Export

Prepare a ChimeraX session with stable reference and moving model IDs, then
write at least two RV waypoints:

```csv
label,rv_x_rad,rv_y_rad,rv_z_rad
start,0,0,0
end,0,0,0.5
```

Generate the Phase 1–4 script-only bundle without launching ChimeraX:

```bash
cryorole animate \
  --run-dir cryorole_outputs \
  --coordinate-set raw \
  --path-csv waypoints.csv \
  --path-space rv \
  --chimerax-session prepared_scene.cxs \
  --reference-model-id "#1" \
  --moving-model-id "#2" \
  --pivot 0 0 0 \
  --baseline-ro identity \
  --map-frame raw \
  --output-dir animation_output
```

The bundle contains `trajectory.csv`, landscape frames, a per-frame transform
artifact, generated ChimeraX scripts, logs, and a manifest. Add
`--render-mode execute --chimerax-bin PATH --no-encode` for validated
structure and composite frames. For MP4, also provide
`--ffmpeg-bin PATH --ffprobe-bin PATH`.
Repeat `--reference-model-id` or `--moving-model-id` for a disjoint rigid
group. All moving models follow the same absolute delta from their own saved
scene transforms; independent motion within the group is not supported.
Animation `--threshold`, `--top-fraction`, `--colormap`, `--vmin`, `--vmax`,
and `--range` use the same display-only policy as `cryorole visualize`. Linux
execute mode uses ChimeraX offscreen rendering and requires its completion
status artifact as well as valid PNG frames. MP4 output is validated as H.264,
`yuv420p`, fixed-size constant-fps video; source frames are retained.

## Complete Workflow Template

Copy and edit the paths:

```bash
cryorole preflight --ref path/to/ref.star --mov path/to/mov.star

cryorole run --ref path/to/ref.star --mov path/to/mov.star

cryorole status --run-dir cryorole_outputs

cryorole canonicalize --run-dir cryorole_outputs

cryorole explore \
  --run-dir cryorole_outputs \
  --space canonical

cryorole select \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --space canonical \
  -c 13 0 14 \
  -r 6

cryorole export \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --domain both
```
