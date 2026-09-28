"""Typed Python API for cryoROLE (the layer a GUI or notebook builds on).

Every function here calls the same services as the ``cryorole`` CLI and
returns their typed results. The API never prints, never calls ``sys.exit``
and never parses ``argv``:

* problems are raised as :class:`cryorole.errors.CryoroleError`
  (``code``, ``message``, ``remedy``); other exceptions are bugs;
* user-relevant warnings and notices go to the ``cryorole`` logger;
* long operations accept ``progress`` (a callable receiving
  :class:`ProgressEvent`) and ``cancel`` (a :class:`CancelToken`), checked at
  stage boundaries. A cancelled ``run`` publishes nothing.

Omitted ``run_dir`` / ``canonical_id`` / ``selection_id`` follow the same
rules as the CLI (``cryorole.workflow.resolve``): a value is filled in only
when exactly one candidate exists, and how it was filled in is recorded in
the written report under ``resolved_by``.

Example::

    from cryorole import api

    check = api.preflight("ref.cs", "mov.cs")
    if check.report["readiness"] != "BLOCKED":
        result = api.run("ref.cs", "mov.cs", output_dir="my_run",
                         progress=lambda e: print(e.kind, e.message))
        api.canonicalize("my_run")
        api.select("my_run", selection_id="region_01", mode="radius",
                   center=(0, 0, 0), radius=15)
        api.export("my_run", selection_id="region_01")
"""

from cryorole.api._core import (
    align,
    canonicalize,
    export,
    fix_subtract_coordinates,
    next_actions,
    preflight,
    run,
    select,
    status,
    visualize,
)
from cryorole.errors import CancelledError, CryoroleError
from cryorole.workflows.progress import CancelToken, ProgressEvent

__all__ = [
    "CancelToken",
    "CancelledError",
    "CryoroleError",
    "ProgressEvent",
    "align",
    "canonicalize",
    "export",
    "fix_subtract_coordinates",
    "next_actions",
    "preflight",
    "run",
    "select",
    "status",
    "visualize",
]
