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
