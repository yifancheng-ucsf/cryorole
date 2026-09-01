# FAQ and Troubleshooting

## Why did preflight return exit code 1?

`READY_WITH_WARNINGS` means validation completed but something needs review,
such as reordered/dropped matches or resource pressure. It intentionally has a
non-zero exit code so automation cannot mistake a warning gate for unconditional
success. Review the human or JSON report before running.

## Why is preflight BLOCKED?

Common causes are a missing/unsupported source, missing STAR pose columns,
missing or incorrectly shaped CryoSPARC `alignments3D/pose`, NaN/Inf poses,
duplicate particle keys, zero overlap, unsafe low overlap, or unequal row counts
under `--row-aligned`. Duplicate diagnostics include counts and bounded
collision examples.

## Does preflight create or overwrite a run?

No. `cryorole preflight` and `cryorole run --dry-run` create no run bundle and
do not modify the source files. A JSON report is written only when you provide
`--json PATH`.

## Is canonicalization required?

No. It is an optional derived coordinate frame. You may explore/select/export
from raw space. `cryorole next` labels canonicalization optional.

## Why do displayed and full candidate counts differ?

The explorer may deterministically downsample or filter for display. Those
controls do not change the scientific candidate universe. Evaluate and Confirm
always run against the exact full parent landscape.

## Did clicking a point create a Selection?

No. Clicking and Evaluate make an in-memory draft. Only explicit Confirm with a
new selection ID writes `selections/<id>/`. Closing the page before Confirm
writes no scientific selection artifact.

## What distance does radius selection use?

The public default is SO(3) geodesic distance. It is not naive Euclidean
distance between Euler triplets. Euler input is converted using the recorded
extrinsic fixed-axis ZYX policy, then evaluated rotation-natively.

## Can I expose the explorer on the network?

No public binding is supported. It binds only to `127.0.0.1`, serves fixed
packaged routes, validates a session token/run ID/landscape hash, caps request
size, and loads no CDN. Use the CLI for remote/headless workflows.

## Does explore require another Python package?

No. It uses the standard library plus cryoROLE's existing dependencies and a
normal local browser. Use `--no-open` if you want to open the printed loopback
URL yourself.

## What does status trust?

It inspects the actual manifest, transaction state and completion marker, raw
NPZ, source identities, canonical frames, visualizations, selections, and
exports. It does not trust a standalone editable workflow-state flag.

## Will guide make scientific decisions for me?

No. It will not infer ref/mov, silently assert row alignment or low overlap,
choose a canonical frame or radius, create/overwrite a Selection, or export.
Non-interactive mode never waits for input. Run execution requires explicit
`--execute-run`.

## Why did Confirm reject my session?

The run ID or parent landscape hash changed after the page loaded, or the
selection ID already exists. Reload from the current completed run and choose a
new ID. This prevents stale drafts and accidental overwrite.
