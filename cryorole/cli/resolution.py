"""Fill in omitted ``--run-dir`` / ``--canonical-id`` / ``--selection-id`` for CLI commands.

The rules live in ``cryorole.workflow.resolve`` (shared with ``cryorole.api``);
this adapter applies them to an argparse namespace, prints one notice per
implicit value to stderr, and stores ``args.resolved_by`` for the reports.
"""

from __future__ import annotations

import sys

from cryorole.workflow.resolve import (
    Resolved,
    canonical_frame_ids,
    resolve_canonical_id,
    resolve_run_dir,
    resolve_selection_id,
)


def _announce(resolved: Resolved) -> None:
    if resolved.implicit:
        print(f"[cryorole] {resolved.notice()}", file=sys.stderr)


def resolve_cli_ids(args, *, canonical: bool = False, selection: str | None = None) -> dict[str, dict[str, str]]:
    """Resolve ids on ``args`` in place.

    ``canonical``: resolve ``--canonical-id`` when ``--space canonical``.
    ``selection``: ``"required"`` always resolves ``--selection-id``;
    ``"if_selected_landscape"`` resolves it only with ``--use-selected-landscape``.
    """

    records: dict[str, dict[str, str]] = {}
    run = resolve_run_dir(getattr(args, "run_dir", None))
    _announce(run)
    args.run_dir = run.value
    records["run_dir"] = run.record()
    if canonical:
        if getattr(args, "space", "raw") == "canonical":
            frame = resolve_canonical_id(args.run_dir, getattr(args, "canonical_id", None))
            _announce(frame)
            args.canonical_id = frame.value
            records["canonical_id"] = frame.record()
        elif getattr(args, "canonical_id", None) is None:
            args.canonical_id = "default"
    wants_selection = selection == "required" or (
        selection == "if_selected_landscape" and getattr(args, "use_selected_landscape", False)
    )
    if wants_selection:
        chosen = resolve_selection_id(args.run_dir, getattr(args, "selection_id", None))
        _announce(chosen)
        args.selection_id = chosen.value
        records["selection_id"] = chosen.record()
    args.resolved_by = records
    return records


def note_raw_space_with_canonical_frames(args) -> None:
    """Tell the user a raw-space command ran although canonical frames exist."""

    if getattr(args, "space", "raw") != "raw" or "space" in getattr(args, "explicit_options", ()):
        return
    frames = canonical_frame_ids(args.run_dir)
    if frames:
        print(
            f"[cryorole] note: selecting in raw space; this bundle also has canonical frame(s) "
            f"{', '.join(frames)}. Add --space canonical to select in a canonical frame.",
            file=sys.stderr,
        )
