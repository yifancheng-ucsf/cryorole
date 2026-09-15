"""Thin CLI adapters for downstream domain services."""

from __future__ import annotations

import os
import shlex
import sys

from cryorole.canonicalize.service import CanonicalizeRequest, canonicalize
from cryorole.select import SelectRequest, create_selection
from cryorole.visualize import VisualizationRequest, visualize


def select_command(args) -> int:
    result = create_selection(SelectRequest.from_namespace(args))
    for path, count in zip(result.selection_dirs, result.selected_counts):
        print(f"[cryorole] selection {path.name!r}: {count} particles; saved to {path}", file=sys.stderr)
        print(f"Next: cryorole export --run-dir {_quote_argument(str(args.run_dir))} "
              f"--selection-id {_quote_argument(path.name)}", file=sys.stderr)
    print(str(result.output_dir))
    return 0


def _quote_argument(value: str) -> str:
    """Format copyable guidance for PowerShell on Windows and POSIX shells elsewhere."""
    return "'" + value.replace("'", "''") + "'" if os.name == "nt" else shlex.quote(value)


def visualize_command(args) -> int:
    result = visualize(VisualizationRequest.from_namespace(args))
    report = result.report
    counts = report.get("n_points_after_display_filter") or {}
    displayed = next(iter(counts.values()), 0) if isinstance(counts, dict) else counts
    print(
        "[cryorole] visualize: "
        f"views={','.join(report.get('requested_views') or [])}; "
        f"filter={report.get('display_filter_mode')}; displayed={displayed}; "
        f"files={len(report.get('generated_files') or {})}",
        file=sys.stderr,
    )
    print(str(result.output_dir))
    return 0


def canonicalize_command(args) -> int:
    return canonicalize(CanonicalizeRequest.from_namespace(args))


def animate_command(args) -> int:
    from cryorole.animation.workflow import run_animation

    result = run_animation(args)
    print(result.output_dir)
    return 0


def canonical_views_command(args) -> int:
    from cryorole.animation.canonical_views import run_canonical_views

    result = run_canonical_views(args)
    print(result.output_dir)
    return 0
