# Python API

`cryorole.api` runs the same services as the `cryorole` command line and
returns their typed results. Use it from notebooks and scripts; it is also the
layer a graphical interface builds on.

```python
from cryorole import api

check = api.preflight("ref.cs", "mov.cs")
print(check.report["readiness"])          # READY, READY_WITH_WARNINGS or BLOCKED

result = api.run("ref.cs", "mov.cs", output_dir="my_run")
print(result.output_dir, result.matched_count)

api.canonicalize("my_run")                # canonical/default/
api.select("my_run", selection_id="region_01", mode="radius",
           center=(0, 0, 0), radius=15)
api.visualize("my_run", space="canonical", view="2d,1d")
api.export("my_run", selection_id="region_01", domain="both")

api.status("my_run")["bundle_status"]    # completed / legacy / incomplete / failed / missing
api.next_actions("my_run")
```

## Functions

| Function | Returns | CLI equivalent |
| --- | --- | --- |
| `preflight(ref, mov, *, row_aligned=False, identity_mode=None, identity_columns=(), mapping_file=None, allow_low_overlap=False, **options)` | `PreflightResult` (`.report`, `.exit_code`) | `cryorole preflight` |
| `run(ref, mov, *, output_dir="cryorole_outputs", row_aligned=False, overwrite=False, progress=None, cancel=None, stream_progress=False, **options)` | `RunResult` (`.output_dir`, `.matched_count`, `.output_artifacts`) | `cryorole run` |
| `canonicalize(run_dir=None, *, canonical_id="default", **options)` | `CanonicalizeResult` (`.output_dir`, `.summary`) | `cryorole canonicalize` |
| `visualize(run_dir=None, *, space="raw", canonical_id=None, selection_id=None, **options)` | `VisualizationResult` (`.output_dir`, `.report`) | `cryorole visualize` |
| `select(run_dir=None, *, selection_id, mode="radius", space="raw", canonical_id=None, **options)` | `SelectResult` (`.output_dir`, `.selection_dirs`, `.selected_counts`) | `cryorole select` |
| `export(run_dir=None, *, selection_id=None, domain="both", format="auto", output_dir=None, overwrite=False, ...)` | export report (dict) | `cryorole export` |
| `align(ref, mov, **options)` | align report (dict) | `cryorole align` |
| `fix_subtract_coordinates(job, **options)` | correction report (dict) | `cryorole align --fix-subtract-coordinates` |
| `status(run_dir=None)`, `next_actions(run_dir=None)` | dict / list of dicts | `cryorole status`, `cryorole next` |

`**options` are the command's options with underscores instead of hyphens, for
example `k_neighbors=30`, `no_visualize=True`, `fit_top_fraction=0.5`,
`sld_min=2`, `metadata_domain="ref"`. An unknown name raises `CryoroleError`
listing the valid ones. For `align`, use the keyword names of
`cryorole.align.align_star_files` (`key_pairs=["_rlnImageName=_rlnImageOriginalName"]`,
`coordinate_match="recentered-exact"`, `recenter_shift=(-9, 22, -102)`, …).

Omitted `run_dir`, `canonical_id` and `selection_id` follow the CLI rules in
`docs/cli_reference.md` (Common behaviour): they are filled in only when exactly
one candidate exists, and the written report records how under `resolved_by`.
`select` always needs an explicit `selection_id`.

## Errors

Problems the user can act on raise `cryorole.api.CryoroleError`, a subclass of
`ValueError`, with:

* `code`: a stable identifier such as `input_not_found`, `input_invalid`,
  `output_exists`, `run_dir_unresolved`, `id_unresolved`, `blocked`,
  `cancelled`;
* `message`: what went wrong;
* `remedy`: what to do, when there is a clear next step;
* `details`: structured extras, for example the candidate list for
  `run_dir_unresolved`.

Some services still raise plain `ValueError` / `FileExistsError`;
`cryorole.errors.classify_exception(exc)` maps any exception to a
`CryoroleError` for display. Any other exception is a bug.

## Progress and cancellation

```python
from cryorole import api

token = api.CancelToken()

def on_progress(event: api.ProgressEvent) -> None:
    print(event.kind, event.stage, event.message, event.fraction)
    # a GUI would update a progress bar, or call token.cancel() from its Cancel button

result = api.run("ref.cs", "mov.cs", output_dir="my_run", progress=on_progress, cancel=token)
```

`ProgressEvent.kind` is `stage_start`, `stage_done`, `batch`, `info` or
`warning`. The cancel token is thread-safe and is checked at every stage
boundary and before the bundle is published. A cancelled run raises
`api.CancelledError` (code `cancelled`), removes its staging directory, and
publishes nothing; an existing bundle at `output_dir` is left untouched.

## Output and logging

The API never prints and never exits the interpreter. `run` writes no progress
text unless `stream_progress=True`. Warnings and notices (for example a
resolved `run_dir`) go to the standard `logging` logger named `cryorole`,
which has no output handler by default:

```python
import logging
logging.basicConfig(level=logging.INFO)   # show cryoROLE warnings and notices
```

The command line installs its own handler, which prints warnings as
`[cryorole] warning: …` on stderr.
