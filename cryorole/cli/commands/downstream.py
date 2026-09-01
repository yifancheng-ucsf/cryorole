"""Thin CLI adapters for downstream domain services."""

from __future__ import annotations

import sys

from cryorole.canonicalize.service import CanonicalizeRequest, canonicalize
from cryorole.select import SelectRequest, create_selection
from cryorole.visualize import VisualizationRequest, visualize


def select_command(args) -> int:
    result = create_selection(SelectRequest.from_namespace(args))
    print(str(result.output_dir))
    return 0


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
