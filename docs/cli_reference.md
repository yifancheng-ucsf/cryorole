# CLI Reference

Use `cryorole COMMAND --help` as the executable source of truth. This page
summarizes the public Workflow UX commands.

## Preflight

```bash
cryorole preflight --ref REF --mov MOV [--output-dir RUN] [--row-aligned]
                   [--allow-low-overlap] [--sld-metric METRIC]
                   [--no-visualize] [--json [PATH]]
```

Performs source identity, convention, pose-schema, identity/matching, and
resource checks without creating `RUN`. `--json` with no path writes JSON to
stdout. Exit codes are 0 ready, 1 ready with warnings, and 2 blocked.

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

## Status and next

```bash
cryorole status --run-dir RUN [--json [PATH]]
cryorole next --run-dir RUN [--json [PATH]]
```

Status reads real artifacts and integrity evidence. Next returns exact
required/recommended/optional commands derived from that status.

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
cryorole run --ref REF --mov MOV [--output-dir RUN]
cryorole canonicalize --run-dir RUN
cryorole visualize --run-dir RUN --space raw|canonical
# default: two PNG 3-view projections for sld_display >= 1
cryorole visualize --run-dir RUN --view 2d,1d
cryorole visualize --run-dir RUN --view 3d
cryorole select --run-dir RUN --selection-id ID --space SPACE -c A B C -r DEG
cryorole export --run-dir RUN --selection-id ID --domain ref|mov|both
```

See `docs/architecture.md` for the stable artifact and scientific policy
contract.


## Visualization terminology

SLD color bars use `SLD`; reports retain `sld_raw` or `sld_display` as the
actual source field. Use `--opacity 0.4` for point opacity (0 transparent,
1 opaque); omission keeps the style default. The old `--alpha` remains a
hidden compatibility alias. Supplying both names is an error.

## Selection modes and names

`cryorole select --help` groups all controls by mode and provides complete
examples. Every command requires `--run-dir RUN --selection-id NAME`; choose a
name such as `region_01` or `high_sld`. Names cannot contain path separators or redirect output outside selections/.
No default is assigned. Reusing a name
requires explicit `--overwrite`; choose another name to keep both selections.

| Mode | Required mode inputs | Optional controls and interpretation |
| --- | --- | --- |
| `radius` | `--center A B C` and one of `--radius DEG`, `--radius-rad RAD` | `--center-representation euler|rotvec`, `--metric so3|rotvec`. Euler center: degrees; rotvec center: radians. Default metric: SO(3) geodesic. |
| `threshold` | `--sld-min`, `--sld-max`, or both | Inclusive `sld_raw` bounds; independent of display filters. |
| `range` | Repeatable `--range-bound AXIS:LOWER:UPPER` | Euler alpha/beta/gamma in degrees OR rotvec x/y/z in radians. Inclusive bounds, intersected across axes; blank sides are open. Euler reversed bounds wrap across the seam; rotvec reversed bounds are invalid. Last bound wins for a repeated axis. |
| `random` | `--fraction F`, with `0 < F <= 1` | Optional non-negative `--seed`; omitted seed uses fresh randomness and is recorded as null. Samples `ceil(F*N)` rows without replacement from all parent rows. |
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
