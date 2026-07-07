"""Command-line frontend for cryoROLE workflow orchestration."""

from __future__ import annotations

import argparse
import ctypes
import csv
import json
from dataclasses import replace
from datetime import datetime
from datetime import timezone
import os
from pathlib import Path
import sys
from types import SimpleNamespace
from typing import Any
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from cryorole.cli.progress import ProgressReporter
from cryorole.core.density import compute_sld_display_values, compute_sld_values
from cryorole.core.euler_conventions import (
    CANONICAL_EULER_ANGLE_COLUMNS,
    DEFAULT_EULER_CONVENTION,
    EULER_CONVENTIONS,
    LEGACY_MISSING_EULER_CONVENTION_SOURCE,
    RAW_EULER_ANGLE_COLUMNS,
    resolve_euler_convention,
)
from cryorole.export import (
    export_selection_metadata_subset,
    export_selection,
    read_landscape,
    read_landscape_json,
    read_landscape_metadata,
    read_landscape_npz_arrays,
    read_selection_json,
    resolve_landscape_path,
    write_canonical_landscape_csv_from_npz,
    write_canonical_landscape_csv,
    write_json_artifact,
    write_landscape_json,
    write_landscape_npz,
    write_landscape_npz_arrays,
    write_landscape_visualizations,
    write_raw_landscape_csv,
    write_report_json,
)
from cryorole.export.serialization import to_json_safe
from cryorole.canonicalize import canonicalize_landscape, canonicalize_landscape_arrays
from cryorole.canonicalize.transforms import apply_canonical_transform, validate_canonical_transform
from cryorole.io.readers import read_relion_star
from cryorole.models.canonicalization_report import CanonicalizationReport
from cryorole.models.landscape_arrays import LandscapeArrays
from cryorole.models.policies import (
    CanonicalizationPolicy,
    ConventionPolicy,
    DensityPolicy,
    IdentityPolicy,
    MatchPolicy,
    RepresentationPolicy,
    RunManifestPolicy,
    SelectionExportPolicy,
    SelectionMetadataExportPolicy,
    SelectionPolicy,
)
from cryorole.models.landscape import Landscape
from cryorole.select import select_particles
from cryorole.workflows.pipeline_runner import PipelineRunner


_HIDDEN_HELP = argparse.SUPPRESS


def build_parser() -> argparse.ArgumentParser:
    """Build the cryoROLE CLI parser without embedding scientific logic."""

    parser = argparse.ArgumentParser(prog="cryorole")
    subparsers = parser.add_subparsers(dest="command", required=True)

    _add_align_parser(subparsers)
    _add_run_parser(subparsers)
    _add_visualize_parser(subparsers)
    _add_canonicalize_parser(subparsers)
    _add_select_parser(subparsers)
    _add_export_parser(subparsers)
    _add_manifest_parser(subparsers)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the cryoROLE CLI."""

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except NotImplementedError as exc:
        parser.exit(2, f"cryorole: {exc}\n")
    except (ValueError, FileExistsError, IsADirectoryError) as exc:
        parser.exit(2, f"cryorole: error: {exc}\n")
    return 0


def _add_align_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "align",
        help="Align STAR metadata before cryorole run --row-aligned.",
        description="Align STAR metadata before cryorole run --row-aligned.",
    )
    parser.add_argument("--ref", required=True, help="Reference STAR metadata path.")
    parser.add_argument("--mov", required=True, help="Moving STAR metadata path.")
    parser.add_argument("--align-id", default="default", help="Output id under alignments/. Default: default.")
    parser.add_argument("--key", nargs="+", help="Explicit STAR key columns for alignment.")
    parser.add_argument(
        "--float-tol",
        action="append",
        default=(),
        type=_parse_float_tolerance,
        metavar="COL=TOL",
        help="Numeric key tolerance, e.g. _rlnCoordinateX=0.1. May repeat.",
    )
    parser.add_argument(
        "--path-mode",
        default="exact",
        help="Path normalization for path-like keys: exact, basename, or suffix:N. Default: exact.",
    )
    parser.add_argument(
        "--duplicate-policy",
        choices=("exclude", "first"),
        default="exclude",
        help="Duplicate key handling. Default: exclude.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace artifacts for this align id only; source STAR files are unchanged.",
    )
    parser.set_defaults(handler=align_command)


def _add_run_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "run",
        help="Generate RO/SLD landscape artifacts without selection by default.",
    )
    parser.add_argument("--ref", required=True, help="Reference-domain pose metadata.")
    parser.add_argument("--mov", required=True, help="Moving-domain pose metadata.")
    parser.add_argument("--ref-domain", default="ref", help="Reference domain name.")
    parser.add_argument("--mov-domain", default="mov", help="Moving domain name.")
    parser.add_argument(
        "--row-aligned",
        action="store_true",
        help="Assert that ref and mov row N refer to the same particle.",
    )
    parser.add_argument(
        "--run-backend",
        choices=("auto", "array_native", "dataframe_compat"),
        default="auto",
        help=_HIDDEN_HELP,
    )
    parser.add_argument(
        "--identity-mode",
        choices=("cryosparc_uid", "relion_user_columns", "explicit_mapping_file", "row_aligned"),
        default=None,
        help=_HIDDEN_HELP,
    )
    parser.add_argument(
        "--identity-column",
        action="append",
        default=(),
        help=_HIDDEN_HELP,
    )
    parser.add_argument("--mapping-file", help=_HIDDEN_HELP)
    parser.add_argument("--k-neighbors", type=int, default=50, help=_HIDDEN_HELP)
    parser.add_argument(
        "--euler-convention",
        choices=EULER_CONVENTIONS,
        help=_HIDDEN_HELP,
    )
    parser.add_argument("--canonicalize", action="store_true", help=_HIDDEN_HELP)
    parser.add_argument(
        "--sign-rule",
        choices=("density_weighted_skewness", "largest_component_positive"),
        default="density_weighted_skewness",
        help=_HIDDEN_HELP,
    )
    parser.add_argument(
        "--skewness-positive-side",
        "--positive-side",
        dest="positive_side",
        choices=("high_density_skew", "low_density_skew"),
        default="high_density_skew",
        help=_HIDDEN_HELP,
    )
    parser.add_argument("--display-top-fraction", type=float, default=0.40, help=_HIDDEN_HELP)
    parser.add_argument("--display-sld-threshold", type=float, help=_HIDDEN_HELP)
    parser.add_argument(
        "--display-density-field",
        choices=("sld_display", "sld_raw"),
        default="sld_display",
        help=_HIDDEN_HELP,
    )
    _add_output_format_argument(parser, help_text=_HIDDEN_HELP)
    parser.add_argument("--max-points-2d", type=int, default=500000, help=_HIDDEN_HELP)
    parser.add_argument("--max-points-3d", type=int, default=50000, help=_HIDDEN_HELP)
    parser.add_argument("--random-seed", type=int, default=0, help=_HIDDEN_HELP)
    parser.add_argument("--write-projection-csvs", action="store_true", help=_HIDDEN_HELP)
    _add_visualization_style_arguments(parser, help_text=_HIDDEN_HELP)
    parser.add_argument(
        "--output-dir",
        default="cryorole_outputs",
        help="Run bundle output directory. Default: cryorole_outputs.",
    )
    parser.add_argument("--manifest-output", help=_HIDDEN_HELP)
    parser.add_argument("--overwrite", action="store_true", help=_HIDDEN_HELP)
    parser.add_argument("--no-visualize", action="store_true", help="Skip raw quick-look visualization.")
    parser.add_argument("--quiet", action="store_true", help=_HIDDEN_HELP)
    parser.add_argument("--verbose", action="store_true", help=_HIDDEN_HELP)
    parser.add_argument("--profile-time", action="store_true", help=_HIDDEN_HELP)
    parser.add_argument("--profile-memory", action="store_true", help=_HIDDEN_HELP)
    parser.add_argument(
        "--raw-csv-chunk-size",
        type=_positive_int_argument,
        default=100000,
        help=_HIDDEN_HELP,
    )
    parser.add_argument(
        "--density-query-batch-size",
        type=_positive_int_argument,
        default=100000,
        help=_HIDDEN_HELP,
    )
    parser.add_argument("--no-raw-csv", action="store_true", help=_HIDDEN_HELP)
    parser.add_argument(
        "--write-debug-json",
        action="store_true",
        help=_HIDDEN_HELP,
    )
    parser.set_defaults(handler=run_command)


def _add_visualize_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "visualize",
        help="Render display-only views from an existing run bundle.",
    )
    parser.add_argument("--run-dir", required=True, help="Existing cryoROLE run bundle directory.")
    parser.add_argument(
        "--space",
        choices=("raw", "canonical"),
        default="raw",
        help="Landscape space to visualize. Default: raw.",
    )
    parser.add_argument("--canonical-id", default="default")
    parser.add_argument(
        "--selection-id",
        help="Visualize only particles from run_dir/selections/SELECTION_ID.",
    )
    parser.add_argument(
        "--use-selected-landscape",
        action="store_true",
        help="Read selections/SELECTION_ID/selected_landscape/landscape.npz directly.",
    )
    parser.add_argument("--visual-id", default="default", help="Visualization output id. Default: default.")
    parser.add_argument(
        "--colormap",
        default="rainbow_r",
        help="Matplotlib colormap for display colors. Default: rainbow_r.",
    )
    parser.add_argument(
        "--representation",
        choices=("euler", "rotvec", "both"),
        default="both",
    )
    parser.add_argument("--top-fraction", type=float, help="Display the top fraction by sld_display.")
    parser.add_argument("--threshold", type=float, help="Display rows with sld_display >= threshold.")
    parser.add_argument(
        "--range",
        dest="range_bound",
        action="append",
        default=None,
        type=_parse_range_bound,
        metavar="AXIS:LOWER:UPPER",
        help="Display-only range filter, e.g. --range alpha:-30:30 or --range x:-0.4:0.4.",
    )
    parser.add_argument(
        "--formats",
        default="png",
        help="Comma-separated figure formats. Default: png. Supported: png,svg,pdf.",
    )
    parser.add_argument("--vmin", type=float, help="Display-only color minimum.")
    parser.add_argument("--vmax", type=float, help="Display-only color maximum.")
    parser.add_argument(
        "--bins",
        type=_positive_int_argument,
        default=72,
        help="1D histogram bin count. Default: 72.",
    )
    parser.add_argument(
        "--hist-mode",
        choices=("count", "percent"),
        default="count",
        help="1D histogram y-axis mode. Default: count.",
    )
    parser.add_argument(
        "--kde-bandwidth",
        default="scott",
        help="1D KDE bandwidth. Default: scott. Accepts scott, silverman, or a positive float.",
    )
    parser.add_argument("--xlim", type=_parse_float_pair_bound, help="1D x-axis limit as MIN:MAX.")
    parser.add_argument("--ylim", type=_parse_float_pair_bound, help="1D y-axis limit as MIN:MAX.")
    parser.add_argument("--overwrite", action="store_true")
    parser.set_defaults(handler=visualize_command)


def _add_canonicalize_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "canonicalize",
        help="Canonicalize an existing landscape explicitly.",
    )
    parser.add_argument("--landscape", help=_HIDDEN_HELP)
    parser.add_argument("--run-dir", help="Existing cryoROLE run bundle directory.")
    parser.add_argument("--canonical-id", default="default", help="Canonical output id. Default: default.")
    parser.add_argument(
        "--output-dir",
        dest="output_dir",
        help=_HIDDEN_HELP,
    )
    parser.add_argument(
        "--output",
        dest="output_dir",
        help=_HIDDEN_HELP,
    )
    parser.add_argument("--overwrite", action="store_true", help=_HIDDEN_HELP)
    parser.add_argument("--no-csv", action="store_true", help=_HIDDEN_HELP)
    parser.add_argument(
        "--euler-convention",
        choices=EULER_CONVENTIONS,
        help=_HIDDEN_HELP,
    )
    parser.add_argument(
        "--csv-chunk-size",
        type=int,
        default=100000,
        help=_HIDDEN_HELP,
    )
    parser.add_argument(
        "--fit-top",
        "--fit-top-fraction",
        dest="fit_top_fraction",
        type=float,
        default=None,
        help="Fraction of highest-sld_raw points used to fit canonical axes. Default: 0.40.",
    )
    parser.add_argument(
        "--profile-memory",
        action="store_true",
        help=_HIDDEN_HELP,
    )
    parser.add_argument(
        "--sign-rule",
        choices=("density_weighted_skewness", "largest_component_positive"),
        default="density_weighted_skewness",
        help=_HIDDEN_HELP,
    )
    parser.add_argument(
        "--positive-side",
        dest="positive_side",
        choices=("low", "high"),
        default=None,
        help="Canonical axis direction: low or high density skew. Default: low.",
    )
    parser.add_argument(
        "--skewness-positive-side",
        dest="legacy_positive_side",
        choices=("high_density_skew", "low_density_skew"),
        default=None,
        help=_HIDDEN_HELP,
    )
    parser.add_argument("--use-frame", help="Existing canonical_frame.json to apply instead of fitting.")
    parser.add_argument("--no-visualize", action="store_true", help="Skip canonical quick-look visualization.")
    parser.set_defaults(handler=canonicalize_command)


def _add_select_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "select",
        help="Create a Selection from an existing landscape after inspection.",
    )
    parser.add_argument("--run-dir", required=True, help="Existing cryoROLE run bundle directory.")
    parser.add_argument("--selection-id", required=True, help="Selection output id under RUN/selections/.")
    parser.add_argument(
        "--space",
        choices=("raw", "canonical"),
        default="raw",
        help="Parent landscape space. Default: raw.",
    )
    parser.add_argument(
        "--canonical-id",
        default="default",
        help="Canonical id when --space canonical. Default: default.",
    )
    parser.add_argument(
        "--mode",
        dest="selection_mode",
        choices=("radius", "threshold", "range", "random", "metadata"),
        default="radius",
        help="Selection mode. Default: radius.",
    )
    parser.add_argument("--sld-min", type=float, help="Lower sld_raw bound for threshold mode.")
    parser.add_argument("--sld-max", type=float, help="Upper sld_raw bound for threshold mode.")
    parser.add_argument("--fraction", type=float, help="Random selection fraction; 0 < F <= 1.")
    parser.add_argument("--seed", type=int, help="Optional random seed.")
    parser.add_argument("--metadata-domain", choices=("ref", "mov"), help="Source metadata domain.")
    parser.add_argument("--metadata-column", help="Source metadata column for metadata mode.")
    parser.add_argument("--metadata-value", help="Metadata value(s); comma-separated values are allowed.")
    parser.add_argument("--split-by-value", action="store_true", help="Write one child selection per metadata value.")
    parser.add_argument(
        "--center",
        "-c",
        nargs=3,
        type=float,
        metavar=("A", "B", "C"),
        help="Radius center coordinates; default representation is Euler degrees.",
    )
    parser.add_argument(
        "--center-representation",
        choices=("euler", "rotvec"),
        default="euler",
        help="Representation of --center. Default: euler.",
    )
    parser.add_argument("--radius", "-r", type=float, help="Radius in degrees for radius mode.")
    parser.add_argument("--radius-rad", type=float, help="Radius in radians for rotation selection.")
    parser.add_argument(
        "--metric",
        choices=("so3", "rotvec"),
        default="so3",
        help="Radius metric. Default: so3.",
    )
    parser.add_argument(
        "--range-bound",
        action="append",
        default=None,
        type=_parse_range_bound,
        metavar="AXIS:LOWER:UPPER",
        help="Coordinate range bound; use blank lower/upper for unconstrained sides.",
    )
    parser.add_argument(
        "--write-selected-landscape",
        action="store_true",
        help=(
            "Also write RUN/selections/ID/selected_landscape/ for "
            "visualize --selection-id ID --use-selected-landscape."
        ),
    )
    parser.add_argument(
        "--recompute-sld",
        action="store_true",
        help="Recompute SLD only in the selected-derived landscape; requires --write-selected-landscape.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing artifacts for this selection id only; raw/canonical landscapes are unchanged.",
    )
    parser.set_defaults(handler=select_command)


def _add_export_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "export",
        help="Export a selection from a run bundle; use --run-dir RUN --selection-id ID.",
    )
    parser.add_argument(
        "--run-dir",
        help="Existing cryoROLE run bundle directory; primary export input with --selection-id.",
    )
    parser.add_argument(
        "--selection-id",
        help="Selection ID under RUN/selections/; primary public input with --run-dir.",
    )
    parser.add_argument(
        "--selection",
        help="Advanced: direct path to a selection.json artifact.",
    )
    parser.add_argument(
        "--domain",
        choices=("ref", "mov", "both"),
        default="both",
        help="Source metadata domain to export: ref, mov, or both. Default: both.",
    )
    parser.add_argument(
        "--format",
        choices=("auto", "relion_star", "cryosparc_cs", "keys"),
        default="auto",
        help="Output format; auto resolves per source domain. Default: auto.",
    )
    parser.add_argument(
        "--output-dir",
        help="Advanced: override output directory. Default: RUN/exports/<selection_id>/.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace export artifacts in the output directory only; source, run, and selection files are unchanged.",
    )
    parser.set_defaults(handler=export_metadata_command)

    export_subparsers = parser.add_subparsers(dest="export_command")
    selection_parser = export_subparsers.add_parser(
        "selection",
        help="Compatibility alias for cryorole export.",
    )
    selection_parser.add_argument(
        "--selection",
        help="Advanced: direct path to a selection.json artifact.",
    )
    selection_parser.add_argument(
        "--run-dir",
        help="Existing cryoROLE run bundle directory; primary export input with --selection-id.",
    )
    selection_parser.add_argument(
        "--selection-id",
        help="Selection ID under RUN/selections/; primary public input with --run-dir.",
    )
    selection_parser.add_argument(
        "--domain",
        choices=("ref", "mov", "both"),
        default="both",
        help="Source metadata domain to export: ref, mov, or both. Default: both.",
    )
    selection_parser.add_argument(
        "--format",
        choices=("auto", "relion_star", "cryosparc_cs", "keys"),
        default="auto",
        help="Output format; auto resolves per source domain. Default: auto.",
    )
    selection_parser.add_argument(
        "--output-dir",
        help="Advanced: override output directory. Default: RUN/exports/<selection_id>/.",
    )
    selection_parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace export artifacts in the output directory only; source, run, and selection files are unchanged.",
    )
    selection_parser.set_defaults(handler=export_selection_command)


def _add_manifest_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "manifest",
        help="Write a provenance manifest from available workflow sections.",
    )
    parser.add_argument("--output", required=True, help="Manifest JSON output path.")
    parser.add_argument("--workflow-name")
    parser.add_argument("--command-string")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--compute-file-hashes", action="store_true")
    parser.add_argument("--hash-algorithm", default="sha256")
    parser.add_argument(
        "--include-selected-particle-keys",
        action="store_true",
        help="Include selected particle keys in the manifest.",
    )
    parser.set_defaults(handler=manifest_command)


def _add_visualization_style_arguments(parser, *, help_text: str | None = None) -> None:
    hidden = help_text == _HIDDEN_HELP
    parser.add_argument(
        "--visual-style",
        choices=("modern", "legacy", "paper"),
        default="modern",
        help=help_text or (
            "Display-only plotting style. legacy changes figure appearance only; "
            "use --display-filter-mode legacy_max_divisor for old max/3-like "
            "display filtering."
        ),
    )
    parser.add_argument("--color-map", help=help_text)
    parser.add_argument("--color-vmin", type=float, help=help_text)
    parser.add_argument("--color-vmax", type=float, help=help_text)
    parser.add_argument("--point-size", type=float, help=help_text)
    parser.add_argument("--point-alpha", type=float, help=help_text)
    parser.add_argument("--figure-width", type=float, help=help_text)
    parser.add_argument("--figure-height", type=float, help=help_text)
    parser.add_argument(
        "--colorbar-position",
        choices=("bottom", "right"),
        help=help_text,
    )
    parser.add_argument(
        "--sort-points-by-color",
        choices=("none", "ascending", "descending"),
        help=help_text,
    )
    parser.add_argument(
        "--axis-limit",
        action="append",
        default=None,
        type=_parse_range_bound,
        metavar=None if hidden else "AXIS:LOWER:UPPER",
        help=help_text or "Plot axis limit only; does not filter display rows.",
    )
    parser.add_argument(
        "--display-filter-mode",
        choices=("none", "top_fraction", "threshold", "legacy_max_divisor"),
        help=help_text or (
            "Display-only filter for visualization outputs; legacy_max_divisor "
            "approximates old max/3-like filtering using the selected SLD display field."
        ),
    )
    parser.add_argument("--display-max-divisor", type=float, default=3.0, help=help_text)
    parser.add_argument("--generate-histograms", action="store_true", help=help_text)
    parser.add_argument("--generate-axis-direction-map", action="store_true", help=help_text)


def run_command(args, *, runner: PipelineRunner | None = None) -> int:
    """Run read/normalize/match/RO/SLD orchestration and write minimal artifacts."""

    run_backend_resolved = _resolve_run_backend(args.run_backend)
    reporter = ProgressReporter(
        command="run",
        quiet=args.quiet,
        verbose=args.verbose,
        profile_time=args.profile_time,
        profile_memory=args.profile_memory,
    )
    reporter.sample_memory("start")
    reporter.info(
        (
            f"backend requested={args.run_backend}, resolved={run_backend_resolved}; "
            f"raw_csv_chunk_size={args.raw_csv_chunk_size}; "
            f"density_query_batch_size={args.density_query_batch_size}"
        ),
        verbose_only=True,
    )
    runner = runner or PipelineRunner()
    source_ref = _source_type_from_path(args.ref)
    source_mov = _source_type_from_path(args.mov)
    identity_policy_ref = _identity_policy_from_args(args, source_type=source_ref)
    identity_policy_mov = _identity_policy_from_args(args, source_type=source_mov)
    convention_policy_ref = _convention_policy_for_source(source_ref)
    convention_policy_mov = _convention_policy_for_source(source_mov)
    match_policy = MatchPolicy()
    density_policy = DensityPolicy(
        k_neighbors=args.k_neighbors,
        display_outlier_mode=getattr(args, "sld_display_outlier_mode", "tail_jump"),
        tail_search_fraction=getattr(args, "sld_tail_search_fraction", 0.01),
        tail_jump_factor=getattr(args, "sld_tail_jump_factor", 5.0),
        max_display_outlier_fraction=getattr(
            args,
            "sld_max_display_outlier_fraction",
            0.002,
        ),
    )
    euler_metadata = _resolve_run_euler_metadata(args)
    representation_policy = RepresentationPolicy(
        euler_convention=str(euler_metadata["euler_convention"]),
        scipy_euler_sequence=str(euler_metadata["scipy_euler_sequence"]),
        euler_degrees=True,
    )
    output_dir = _prepare_run_output_dir(args.output_dir, overwrite=args.overwrite)
    data_dir = output_dir / "data"
    reports_dir = output_dir / "reports"
    raw_visualization_dir = output_dir / "visualizations" / "raw_default"
    with reporter.timed_stage(
        "prepare_output_dirs",
        label="preparing run bundle directories",
        stage_number=1,
        stage_total=8,
        path=str(output_dir),
    ):
        for child in (
            data_dir,
            reports_dir,
            output_dir / "visualizations",
            output_dir / "canonical",
            output_dir / "selections",
            output_dir / "exports",
            output_dir / "debug",
        ):
            child.mkdir(parents=True, exist_ok=True)
    with reporter.timed_stage(
        "read_normalize_match",
        label="reading, normalizing, and matching metadata",
        stage_number=2,
        stage_total=8,
        detail=f"ref={args.ref}; mov={args.mov}",
    ) as stage:
        phase1 = runner.run_phase1(
            args.ref,
            args.mov,
            domain_a_name=args.ref_domain,
            domain_b_name=args.mov_domain,
            identity_policy_a=identity_policy_ref,
            identity_policy_b=identity_policy_mov,
            convention_policy_a=convention_policy_ref,
            convention_policy_b=convention_policy_mov,
            match_policy=match_policy,
        )
        stage["ref_row_count"] = len(phase1.pose_a.data)
        stage["mov_row_count"] = len(phase1.pose_b.data)
        stage["matched_count"] = _get_field(phase1.match_report, "matched_count")
    reporter.sample_memory("source_loaded")
    reporter.info(
        (
            f"matched particles: {_get_field(phase1.match_report, 'matched_count')} / "
            f"ref={len(phase1.pose_a.data)} / mov={len(phase1.pose_b.data)}"
        )
    )
    for warning in getattr(phase1.match_report, "warnings", ()):
        reporter.warning(str(warning))
    with reporter.timed_stage(
        "compute_ro",
        label="computing relative orientations",
        stage_number=3,
        stage_total=8,
    ):
        phase2 = runner.compute_ro_for_phase1_result(
            phase1,
            euler_sequence=representation_policy.scipy_euler_sequence,
        )
    reporter.sample_memory("ro_computed")
    with reporter.timed_stage(
        "compute_sld",
        label=f"computing SLD with k={density_policy.k_neighbors}",
        stage_number=4,
        stage_total=8,
        density_query_batch_size=args.density_query_batch_size,
    ):
        phase3 = runner.compute_density_for_ro_result(
            phase2,
            density_policy=density_policy,
        )
    if phase3.landscape.density_report is not None:
        for warning in phase3.landscape.density_report.warnings:
            reporter.warning(warning)
    reporter.sample_memory("density_computed")
    match_table_data = _match_table_data_for_phase1(phase1, phase3.landscape)
    raw_landscape = _landscape_with_match_rows(phase3.landscape, match_table_data)
    raw_active_policies = dict(raw_landscape.active_policies or {})
    raw_active_policies["representation_policy"] = representation_policy
    raw_landscape.active_policies = raw_active_policies
    with reporter.timed_stage(
        "write_npz",
        label="writing raw landscape NPZ",
        stage_number=5,
        stage_total=8,
        row_count=len(raw_landscape.data),
    ):
        raw_landscape_npz_path = write_landscape_npz(
            raw_landscape,
            data_dir / "raw_landscape.npz",
            overwrite=args.overwrite,
            artifact_type="raw_landscape",
        )
    reporter.sample_memory("raw_npz_written")
    raw_landscape_csv_path = None
    if args.no_raw_csv:
        reporter.stage_start(
            "raw CSV skipped by --no-raw-csv",
            stage_number=6,
            stage_total=8,
        )
        reporter._stages.append({"stage": "write_raw_csv", "elapsed_sec": 0.0, "skipped": True})
    else:
        with reporter.timed_stage(
            "write_raw_csv",
            label="writing raw landscape CSV",
            stage_number=6,
            stage_total=8,
            row_count=len(raw_landscape.data),
            raw_csv_chunk_size=args.raw_csv_chunk_size,
        ):
            raw_landscape_csv_path = write_raw_landscape_csv(
                raw_landscape,
                data_dir / "raw_landscape.csv",
                overwrite=args.overwrite,
                euler_sequence=representation_policy.scipy_euler_sequence,
            )
        reporter.sample_memory("raw_csv_written")
    with reporter.timed_stage(
        "write_core_reports",
        label="writing core reports",
        row_count=len(raw_landscape.data),
    ):
        match_table_path = _write_csv_artifact(
            match_table_data,
            data_dir / "match_table.csv",
            overwrite=args.overwrite,
        )
        density_report_path = write_report_json(
            raw_landscape.density_report,
            reports_dir / "density_report.json",
            overwrite=args.overwrite,
            artifact_type="density_report",
        )
        match_report_path = write_report_json(
            phase1.match_report,
            reports_dir / "match_report.json",
            overwrite=args.overwrite,
            artifact_type="match_report",
        )
        identity_ref_report_path = write_report_json(
            _identity_report_for_phase1(phase1, "identity_a", "ref"),
            reports_dir / "identity_ref_report.json",
            overwrite=args.overwrite,
            artifact_type="identity_report",
        )
        identity_mov_report_path = write_report_json(
            _identity_report_for_phase1(phase1, "identity_b", "mov"),
            reports_dir / "identity_mov_report.json",
            overwrite=args.overwrite,
            artifact_type="identity_report",
        )
        import_ref_report_path = write_json_artifact(
            _simple_import_report(args.ref, source_ref, len(phase1.pose_a.data)),
            reports_dir / "import_ref_report.json",
            overwrite=args.overwrite,
        )
        import_mov_report_path = write_json_artifact(
            _simple_import_report(args.mov, source_mov, len(phase1.pose_b.data)),
            reports_dir / "import_mov_report.json",
            overwrite=args.overwrite,
        )
    output_artifacts = {
        "raw_landscape_npz": str(raw_landscape_npz_path),
        "match_table_csv": str(match_table_path),
        "density_report_json": str(density_report_path),
        "match_report_json": str(match_report_path),
        "identity_ref_report_json": str(identity_ref_report_path),
        "identity_mov_report_json": str(identity_mov_report_path),
        "import_ref_report_json": str(import_ref_report_path),
        "import_mov_report_json": str(import_mov_report_path),
    }
    if raw_landscape_csv_path is not None:
        output_artifacts["raw_landscape_csv"] = str(raw_landscape_csv_path)
    if args.write_debug_json:
        landscape_debug_path = write_landscape_json(
            raw_landscape,
            output_dir / "debug" / "landscape_debug.json",
            overwrite=args.overwrite,
            artifact_type="landscape",
        )
        output_artifacts["landscape_debug_json"] = str(landscape_debug_path)
    canonicalization_policy = None
    canonical_landscape = None
    if args.canonicalize:
        with reporter.timed_stage("canonicalize", label="canonicalizing landscape"):
            canonical_output_dir = output_dir / "canonical" / "default"
            canonical_output_dir.mkdir(parents=True, exist_ok=True)
            canonicalization_policy = CanonicalizationPolicy(
                sign_rule=args.sign_rule,
                positive_side=args.positive_side,
            )
            canonical_input = SimpleNamespace(
                phase2=getattr(phase3, "phase2", None),
                landscape=raw_landscape,
            )
            canonical_result = runner.canonicalize_density_result(
                canonical_input,
                canonicalization_policy=canonicalization_policy,
            )
            canonical_landscape = canonical_result.landscape
            canonical_landscape_path = write_landscape_npz(
                canonical_landscape,
                canonical_output_dir / "canonical_landscape.npz",
                overwrite=args.overwrite,
                artifact_type="canonical_landscape",
            )
            canonical_landscape_csv_path = write_canonical_landscape_csv(
                canonical_landscape,
                canonical_output_dir / "canonical_landscape.csv",
                overwrite=args.overwrite,
                euler_sequence=representation_policy.scipy_euler_sequence,
            )
            canonicalization_report_path = write_report_json(
                canonical_landscape.canonicalization_report,
                canonical_output_dir / "canonicalization_report.json",
                overwrite=args.overwrite,
                artifact_type="canonicalization_report",
            )
        output_artifacts["canonical_landscape_npz"] = str(canonical_landscape_path)
        output_artifacts["canonical_landscape_csv"] = str(canonical_landscape_csv_path)
        output_artifacts["canonicalization_report_json"] = str(canonicalization_report_path)
        if args.write_debug_json:
            canonical_debug_path = write_landscape_json(
                canonical_landscape,
                output_dir / "debug" / "canonical_default_landscape_debug.json",
                overwrite=args.overwrite,
                artifact_type="canonical_landscape",
            )
            output_artifacts["canonical_landscape_debug_json"] = str(canonical_debug_path)

    raw_visualization_performed = False
    if args.no_visualize:
        reporter.stage_start(
            "raw visualization skipped by --no-visualize",
            stage_number=7,
            stage_total=8,
        )
        reporter._stages.append({"stage": "visualize_raw", "elapsed_sec": 0.0, "skipped": True})
    else:
        with reporter.timed_stage(
            "visualize_raw",
            label="writing raw visualization artifacts",
            stage_number=7,
            stage_total=8,
        ):
            visualization_report = _write_run_default_visualization_previews(
                args=args,
                raw_landscape=raw_landscape,
                raw_visualization_dir=raw_visualization_dir,
                representation_policy=representation_policy,
                euler_metadata=euler_metadata,
                raw_landscape_csv_path=raw_landscape_csv_path,
            )
        raw_visualization_performed = True
        output_artifacts["visualization_report_json"] = visualization_report["report_path"]
        for artifact_key, artifact_path in visualization_report["generated_files"].items():
            output_artifacts[f"raw_visualization_{artifact_key}"] = artifact_path

    manifest_path = Path(args.manifest_output) if args.manifest_output else output_dir / "run_manifest.json"
    run_summary_path = output_dir / "run_summary.json"
    timing_profile_path = reports_dir / "run_timing_profile.json"
    memory_profile_path = reports_dir / "run_memory_profile.json"
    if args.profile_time:
        output_artifacts["run_timing_profile_json"] = str(timing_profile_path)
    if args.profile_memory:
        output_artifacts["run_memory_profile_json"] = str(memory_profile_path)
    output_artifacts["run_summary_json"] = str(run_summary_path)
    output_artifacts["run_manifest_json"] = str(manifest_path)
    active_policies = {
        "identity_policy_ref": identity_policy_ref,
        "identity_policy_mov": identity_policy_mov,
        "convention_policy_ref": convention_policy_ref,
        "convention_policy_mov": convention_policy_mov,
        "match_policy": match_policy,
        "density_policy": density_policy,
        "representation_policy": representation_policy,
    }
    if canonicalization_policy is not None:
        active_policies["canonicalization_policy"] = canonicalization_policy
    with reporter.timed_stage(
        "manifest",
        label="writing run summary and manifest",
        stage_number=8,
        stage_total=8,
    ):
        write_json_artifact(
            _run_summary_payload(
                args=args,
                source_ref=source_ref,
                source_mov=source_mov,
                phase1=phase1,
                landscape=phase3.landscape,
                matched_count=_get_field(phase1.match_report, "matched_count"),
                k_neighbors=density_policy.k_neighbors,
                canonicalization_performed=args.canonicalize,
                output_artifacts=output_artifacts,
                run_backend_resolved=run_backend_resolved,
                raw_csv_performed=raw_landscape_csv_path is not None,
                raw_csv_backend=(
                    "dataframe_compat_full_table"
                    if raw_landscape_csv_path is not None
                    else "skipped"
                ),
                raw_visualization_performed=raw_visualization_performed,
                timing_profile_path=timing_profile_path if args.profile_time else None,
                memory_profile_path=memory_profile_path if args.profile_memory else None,
                euler_metadata=euler_metadata,
            ),
            run_summary_path,
            overwrite=args.overwrite,
        )
        runner.write_manifest(
            manifest_policy=RunManifestPolicy(
                output_path=manifest_path,
                overwrite=args.overwrite,
                workflow_name="run",
                schema_version="2.0",
            ),
            input_paths=(args.ref, args.mov),
            source_types={
                str(args.ref): source_ref,
                str(args.mov): source_mov,
            },
            row_counts={
                str(args.ref): len(phase1.pose_a.data),
                str(args.mov): len(phase1.pose_b.data),
            },
            active_policies=active_policies,
            match_report=phase1.match_report,
            density_report=phase3.landscape.density_report,
            canonicalization_report=(
                canonical_landscape.canonicalization_report
                if canonical_landscape is not None
                else None
            ),
            additional_results={"euler_metadata": euler_metadata},
            output_artifacts=output_artifacts,
        )
    if args.profile_time:
        write_json_artifact(
            reporter.timing_profile_payload(),
            timing_profile_path,
            overwrite=args.overwrite,
        )
    if args.profile_memory:
        reporter.sample_memory("completed")
        write_json_artifact(
            reporter.memory_profile_payload(),
            memory_profile_path,
            overwrite=args.overwrite,
        )
    reporter.info(f"completed successfully: {output_dir}")
    print(str(output_dir))
    return 0


def select_command(args) -> int:
    """Load a saved landscape, select particles, and export selection artifacts."""

    if getattr(args, "recompute_sld", False) and not getattr(args, "write_selected_landscape", False):
        raise ValueError("--recompute-sld requires --write-selected-landscape")
    landscape_source = args.run_dir
    landscape_metadata = read_landscape_metadata(
        landscape_source,
        space=args.space,
        canonical_id=args.canonical_id,
    )
    euler_metadata = _resolve_landscape_euler_metadata(
        args,
        columns=RAW_EULER_ANGLE_COLUMNS,
    )
    policy = _selection_policy_from_args(args, euler_metadata=euler_metadata)
    landscape = read_landscape(
        landscape_source,
        space=args.space,
        canonical_id=args.canonical_id,
    )

    source_metadata = None
    if policy.selection_mode in {"metadata_value", "metadata_group"}:
        source_metadata, metadata_source_file = _load_run_source_metadata_for_selection(
            args,
            policy=policy,
        )
        policy = replace(
            policy,
            metadata_source_file=str(metadata_source_file),
            metadata_source_row_id_field=f"{policy.metadata_domain}_source_row_id",
        )

    if policy.selection_mode == "metadata_group":
        selection_dirs = _write_metadata_group_selections(
            args=args,
            landscape=landscape,
            landscape_metadata=landscape_metadata,
            euler_metadata=euler_metadata,
            policy=policy,
            source_metadata=source_metadata,
        )
        print(str(_metadata_group_output_root(args)))
        return 0

    output_dir = _prepare_select_output_dir(
        _selection_output_dir(args, policy=policy),
        overwrite=args.overwrite,
    )
    selection = select_particles(landscape, policy=policy, source_metadata=source_metadata)
    _write_selection_outputs(
        args=args,
        landscape=landscape,
        landscape_metadata=landscape_metadata,
        euler_metadata=euler_metadata,
        policy=policy,
        selection=selection,
        output_dir=output_dir,
    )
    print(str(output_dir))
    return 0


def _write_selection_outputs(
    *,
    args,
    landscape: Landscape,
    landscape_metadata: dict[str, object],
    euler_metadata: dict[str, object],
    policy: SelectionPolicy,
    selection,
    output_dir: Path,
) -> None:
    export_report = export_selection(
        selection,
        policy=SelectionExportPolicy(
            output_dir=output_dir,
            overwrite=args.overwrite,
        ),
    )
    selection_csv_path = _write_selection_csv(
        selection,
        output_dir / "selection.csv",
        overwrite=args.overwrite,
    )
    selected_landscape_report = None
    if args.write_selected_landscape:
        selected_landscape_report = _write_selected_derived_landscape(
            parent_landscape=landscape,
            selection=selection,
            output_dir=output_dir / "selected_landscape",
            overwrite=args.overwrite,
            recompute_sld=bool(args.recompute_sld),
            density_policy=DensityPolicy(),
            coordinate_space=args.space,
            parent_landscape_metadata=landscape_metadata,
            euler_metadata=euler_metadata,
        )
        selected_rows_path = Path(
            selected_landscape_report["output_paths"]["selected_landscape_rows_csv"]
        )
    else:
        selected_rows_path = _write_selected_landscape_rows_csv(
            Landscape(
                data=_selected_rows_from_parent_landscape_allow_empty(
                    landscape,
                    selection.selected_particle_keys,
                ),
                canonical_transform=landscape.canonical_transform,
                active_policies=dict(landscape.active_policies or {}),
                density_report=None,
                canonicalization_report=landscape.canonicalization_report,
            ),
            output_dir / "selected_landscape_rows.csv",
            overwrite=args.overwrite,
        )
    write_json_artifact(
        {
            "artifact_type": "selection_summary",
            "schema_version": "1",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "selection_id": selection.selection_id,
            "overwrite": bool(args.overwrite),
            "mode": _public_selection_mode(selection.selection_mode),
            "run_dir": str(args.run_dir),
            "parent_space": args.space,
            "canonical_id": args.canonical_id if args.space == "canonical" else None,
            "parent_landscape_path": landscape_metadata["path"],
            "candidate_count": selection.total_count,
            "landscape_path": landscape_metadata["path"],
            "landscape_artifact_type": landscape_metadata["artifact_type"],
            "landscape_schema_version": landscape_metadata["schema_version"],
            "landscape_row_count": landscape_metadata["row_count"],
            "selection_mode": selection.selection_mode,
            "selection_basis": selection.selection_basis,
            "metric": selection.metric,
            "random_fraction": selection.random_fraction,
            "random_seed": selection.random_seed,
            "random_candidate_count": selection.random_candidate_count,
            "metadata_domain": selection.metadata_domain,
            "metadata_source_file": selection.metadata_source_file,
            "metadata_column": selection.metadata_column,
            "metadata_values": selection.metadata_values,
            "metadata_source_row_id_field": selection.metadata_source_row_id_field,
            "metadata_candidate_count": selection.metadata_candidate_count,
            "metadata_missing_count": selection.metadata_missing_count,
            "metadata_invalid_count": selection.metadata_invalid_count,
            "density_artifact_policy": selection.density_artifact_policy,
            "density_artifact_flag_field": selection.density_artifact_flag_field,
            "density_artifact_candidate_count_before": (
                selection.density_artifact_candidate_count_before
            ),
            "density_artifact_candidate_count_after": (
                selection.density_artifact_candidate_count_after
            ),
            "density_artifact_excluded_count": selection.density_artifact_excluded_count,
            "radius": selection.radius,
            "radius_unit": selection.radius_unit,
            "sld_min": selection.sld_min,
            "sld_max": selection.sld_max,
            "selected_count": selection.selected_count,
            "total_count": selection.total_count,
            "selection_policy": to_json_safe(policy),
            "policy": to_json_safe(policy),
            "euler_metadata": euler_metadata,
            **euler_metadata,
            "export_report_output_paths": to_json_safe(export_report.output_paths),
            "export_report_row_counts": to_json_safe(export_report.row_counts),
            "export_report": export_report,
            "selection_csv": str(selection_csv_path),
            "selected_particle_keys_path": to_json_safe(
                export_report.output_paths.get("selected_particle_keys_csv")
            ),
            "source_row_provenance_fields": [
                column
                for column in ("particle_key", "ref_source_row_id", "mov_source_row_id")
                if column in landscape.data.columns
            ],
            "warnings": [],
            "selected_landscape_rows_csv": str(selected_rows_path),
            "selected_landscape_written": selected_landscape_report is not None,
            "selected_landscape_report": selected_landscape_report,
        },
        output_dir / "selection_summary.json",
        overwrite=args.overwrite,
    )


def _write_metadata_group_selections(
    *,
    args,
    landscape: Landscape,
    landscape_metadata: dict[str, object],
    euler_metadata: dict[str, object],
    policy: SelectionPolicy,
    source_metadata: pd.DataFrame,
) -> list[Path]:
    if not policy.split_by_metadata:
        raise ValueError("metadata split selection requires --split-by-value")
    values = _metadata_group_values_for_landscape(
        landscape,
        source_metadata=source_metadata,
        policy=policy,
    )
    output_dirs: list[Path] = []
    for value in values:
        selection_id = _metadata_group_selection_id(policy, value)
        child_policy = replace(
            policy,
            selection_mode="metadata_value",
            metadata_values=(value,),
            selection_id=selection_id,
            split_by_metadata=False,
        )
        output_dir = _prepare_select_output_dir(
            _metadata_group_output_root(args) / selection_id,
            overwrite=args.overwrite,
        )
        selection = select_particles(
            landscape,
            policy=child_policy,
            source_metadata=source_metadata,
        )
        _write_selection_outputs(
            args=args,
            landscape=landscape,
            landscape_metadata=landscape_metadata,
            euler_metadata=euler_metadata,
            policy=child_policy,
            selection=selection,
            output_dir=output_dir,
        )
        output_dirs.append(output_dir)
    if not output_dirs:
        raise ValueError("metadata_group selection found no non-missing metadata values")
    return output_dirs


def _load_run_source_metadata_for_selection(
    args,
    *,
    policy: SelectionPolicy,
) -> tuple[pd.DataFrame, Path]:
    if not getattr(args, "run_dir", None):
        raise ValueError("metadata selection requires --run-dir")
    if policy.metadata_domain not in {"ref", "mov"}:
        raise ValueError("metadata selection requires --metadata-domain ref or mov")
    if not policy.metadata_column:
        raise ValueError("metadata selection requires --metadata-column")
    source_info = _resolve_run_source_info_for_metadata(Path(args.run_dir))
    domain_info = source_info.get(str(policy.metadata_domain), {})
    source_path_raw = domain_info.get("path")
    source_type = domain_info.get("source_type")
    if not source_path_raw:
        raise ValueError(f"Cannot locate original {policy.metadata_domain} source metadata path")
    source_path = Path(source_path_raw)
    if source_type not in {"relion", "relion_star", "star"} and source_path.suffix.lower() != ".star":
        raise ValueError("source metadata selection currently supports run-time RELION STAR files")
    star = read_relion_star(source_path)
    _resolve_cli_metadata_column(star.particles, str(policy.metadata_column))
    return star.particles, source_path


def _resolve_run_source_info_for_metadata(run_dir: Path) -> dict[str, dict[str, str | None]]:
    summary = _read_json_if_exists(run_dir / "run_summary.json") or {}
    manifest = _read_json_if_exists(run_dir / "run_manifest.json") or {}
    input_paths = _domain_string_mapping(summary.get("input_paths"))
    if not input_paths:
        input_paths = _manifest_input_paths_for_metadata(manifest)
    source_types = _domain_string_mapping(summary.get("source_types"))
    if not source_types:
        source_types = _manifest_source_types_for_metadata(manifest, input_paths)
    info: dict[str, dict[str, str | None]] = {}
    for domain in ("ref", "mov"):
        raw_path = input_paths.get(domain)
        if raw_path is None:
            info[domain] = {"path": None, "source_type": None}
            continue
        path = _resolve_run_relative_path(raw_path, run_dir)
        info[domain] = {
            "path": str(path),
            "source_type": source_types.get(domain) or _source_type_from_path(str(path)),
        }
    return info


def _domain_string_mapping(value: object) -> dict[str, str]:
    if not isinstance(value, dict):
        return {}
    return {
        domain: str(value[domain])
        for domain in ("ref", "mov")
        if value.get(domain) is not None
    }


def _manifest_input_paths_for_metadata(manifest: dict[str, object]) -> dict[str, str]:
    records = manifest.get("input_provenance", {}).get("inputs", [])
    if not isinstance(records, list) or len(records) < 2:
        return {}
    paths = []
    for record in records[:2]:
        if isinstance(record, dict) and record.get("path") is not None:
            paths.append(str(record["path"]))
    if len(paths) != 2:
        return {}
    return {"ref": paths[0], "mov": paths[1]}


def _manifest_source_types_for_metadata(
    manifest: dict[str, object],
    input_paths: dict[str, str],
) -> dict[str, str]:
    source_types = manifest.get("input_provenance", {}).get("source_types", {})
    if not isinstance(source_types, dict):
        return {}
    by_domain = {}
    for domain, path in input_paths.items():
        if path in source_types:
            by_domain[domain] = str(source_types[path])
    return by_domain


def _resolve_run_relative_path(path: str, run_dir: Path) -> Path:
    candidate = Path(path)
    if candidate.exists():
        return candidate
    relative = run_dir / candidate
    if relative.exists():
        return relative
    return candidate


def _resolve_cli_metadata_column(source_metadata: pd.DataFrame, requested: str) -> str:
    if requested in source_metadata.columns:
        return requested
    prefixed = requested if requested.startswith("_") else f"_{requested}"
    if prefixed in source_metadata.columns:
        return prefixed
    unprefixed = requested[1:] if requested.startswith("_") else requested
    if unprefixed in source_metadata.columns:
        return unprefixed
    raise ValueError(f"Source metadata missing column: {requested}")


def _metadata_group_values_for_landscape(
    landscape: Landscape,
    *,
    source_metadata: pd.DataFrame,
    policy: SelectionPolicy,
) -> tuple[str, ...]:
    row_field = policy.metadata_source_row_id_field or f"{policy.metadata_domain}_source_row_id"
    if row_field not in landscape.data.columns:
        raise ValueError(f"Landscape missing metadata source-row field: {row_field}")
    column = _resolve_cli_metadata_column(source_metadata, str(policy.metadata_column))
    values_by_row_id = {
        int(index): _cli_metadata_value_key(value)
        for index, value in source_metadata[column].items()
    }
    values: list[str] = []
    seen: set[str] = set()
    for source_row_id in landscape.data[row_field]:
        row_id = _coerce_cli_source_row_id(source_row_id)
        if row_id is None or row_id not in values_by_row_id:
            continue
        value = values_by_row_id[row_id]
        if value is None or value in seen:
            continue
        seen.add(value)
        values.append(value)
    return tuple(values)


def _coerce_cli_source_row_id(value) -> int | None:
    if isinstance(value, (bool, np.bool_)):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(numeric) or not numeric.is_integer():
        return None
    return int(numeric)


def _cli_metadata_value_key(value) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    try:
        numeric = float(text)
    except ValueError:
        return text
    if np.isfinite(numeric) and numeric.is_integer():
        return str(int(numeric))
    return text


def _metadata_group_selection_id(policy: SelectionPolicy, value: str) -> str:
    prefix = policy.selection_id
    if prefix:
        return f"{_sanitize_selection_id_component(prefix)}_{_sanitize_selection_id_component(value)}"
    column = _sanitize_selection_id_component(str(policy.metadata_column))
    domain = _sanitize_selection_id_component(str(policy.metadata_domain))
    return f"{domain}_{column}_{_sanitize_selection_id_component(value)}"


def _sanitize_selection_id_component(value: str) -> str:
    cleaned = "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in value)
    cleaned = cleaned.strip("_")
    return cleaned or "value"


def _metadata_group_output_root(args) -> Path:
    if args.run_dir:
        return Path(args.run_dir) / "selections"
    raise ValueError("select requires --run-dir")


def _selected_rows_from_parent_landscape_allow_empty(
    landscape: Landscape,
    selected_particle_keys: Sequence[object],
) -> pd.DataFrame:
    if selected_particle_keys:
        return _selected_rows_from_landscape(landscape, selected_particle_keys)
    columns = [
        column
        for column in ("particle_key", "ref_source_row_id", "mov_source_row_id")
        if column in landscape.data.columns
    ]
    selected = landscape.data.iloc[0:0].copy()
    if "sld_display_is_outlier" not in selected.columns:
        selected["sld_display_is_outlier"] = []
    return selected[columns + [column for column in selected.columns if column not in columns]]


def manifest_command(args) -> int:
    """Write a partial manifest through the workflow layer."""

    runner = PipelineRunner()
    result = runner.write_manifest(
        manifest_policy=RunManifestPolicy(
            output_path=args.output,
            overwrite=args.overwrite,
            workflow_name=args.workflow_name,
            command=args.command_string,
            compute_file_hashes=args.compute_file_hashes,
            hash_algorithm=args.hash_algorithm,
            include_selected_particle_keys=args.include_selected_particle_keys,
        )
    )
    print(result.manifest_report.output_path)
    return 0


def align_command(args) -> int:
    from cryorole.align import align_star_files

    report = align_star_files(
        ref=args.ref,
        mov=args.mov,
        align_id=args.align_id,
        key_columns=args.key,
        float_tolerances=dict(args.float_tol or ()),
        path_mode=args.path_mode,
        duplicate_policy=args.duplicate_policy,
        overwrite=args.overwrite,
    )
    print(report["output_dir"])
    return 0


def visualize_command(args) -> int:
    if not getattr(args, "run_dir", None):
        raise ValueError("visualize requires --run-dir")
    if getattr(args, "use_selected_landscape", False) and not getattr(args, "selection_id", None):
        raise ValueError("visualize --use-selected-landscape requires --selection-id")
    if args.top_fraction is not None and args.threshold is not None:
        raise ValueError("Use only one display density filter: --top-fraction or --threshold")
    if args.top_fraction is not None and not 0.0 < args.top_fraction <= 1.0:
        raise ValueError("--top-fraction must be > 0 and <= 1")
    if args.threshold is not None and not np.isfinite(args.threshold):
        raise ValueError("--threshold must be finite")
    kde_bandwidth = _parse_kde_bandwidth(args.kde_bandwidth)
    formats = _parse_formats(args.formats)
    _validate_visualize_color_field_args(args)
    landscape_source = args.run_dir
    output_dir = _visualization_output_dir(args)
    euler_metadata = _resolve_landscape_euler_metadata(
        args,
        columns=RAW_EULER_ANGLE_COLUMNS,
    )
    selection_metadata = None
    source_landscape_path = resolve_landscape_path(
        landscape_source,
        space=args.space,
        canonical_id=args.canonical_id,
    )
    if getattr(args, "use_selected_landscape", False):
        selected_landscape_path = _resolve_selected_landscape_path(
            Path(args.run_dir),
            args.selection_id,
        )
        landscape = read_landscape(selected_landscape_path)
        source_landscape_path = selected_landscape_path
        selection_metadata = _selected_landscape_visualization_metadata(
            Path(args.run_dir),
            args.selection_id,
            selected_landscape_path,
            selected_count=len(landscape.data),
        )
    else:
        landscape = read_landscape(
            landscape_source,
            space=args.space,
            canonical_id=args.canonical_id,
        )
    if args.selection_id and not getattr(args, "use_selected_landscape", False):
        landscape, selection_metadata = _filter_landscape_by_selection(args, landscape)
    _validate_visualize_color_field(landscape, "sld_display")
    range_bounds = dict(args.range_bound or ())
    axis_limits = _visualize_axis_limits_from_ranges(range_bounds)
    report = write_landscape_visualizations(
        landscape,
        output_dir,
        overwrite=args.overwrite,
        coordinate_source=_coordinate_source_for_landscape(args),
        representation=args.representation,
        color_field="sld_display",
        euler_convention=str(euler_metadata["euler_convention"]),
        euler_convention_source=str(euler_metadata["euler_convention_source"]),
        euler_degrees=True,
        display_top_fraction=args.top_fraction,
        display_sld_threshold=args.threshold,
        display_density_field="sld_display",
        range_bounds=range_bounds,
        formats=formats,
        max_points_2d=None,
        max_points_3d=50000,
        random_seed=0,
        write_projection_csvs=False,
        display_table_filename="display_table.csv",
        artifact_layout="run_bundle",
        selection_metadata=selection_metadata,
        visual_style="legacy",
        color_map=args.colormap,
        color_vmin=args.vmin,
        color_vmax=args.vmax,
        axis_limits=axis_limits,
        display_filter_mode=(
            "top_fraction"
            if args.top_fraction is not None
            else ("threshold" if args.threshold is not None else "none")
        ),
    )
    distribution_report = _write_visualization_1d_distributions(
        output_dir=Path(output_dir),
        display_table_path=Path(output_dir) / "display_table.csv",
        representation=args.representation,
        formats=formats,
        bins=args.bins,
        hist_mode=args.hist_mode,
        kde_bandwidth=kde_bandwidth,
        xlim=args.xlim,
        ylim=args.ylim,
        overwrite=args.overwrite,
    )
    report = _visualization_public_report(
        report,
        args=args,
        source_landscape_path=source_landscape_path,
        euler_metadata=euler_metadata,
        formats=formats,
        distribution_report=distribution_report,
    )
    write_json_artifact(
        report,
        Path(output_dir) / "visualization_report.json",
        overwrite=True,
    )
    print(str(Path(output_dir)))
    return 0


def canonicalize_command(args) -> int:
    memory_profiler = _MemoryProfiler(enabled=args.profile_memory)
    memory_profiler.sample("start")
    if not args.run_dir and not (args.landscape and args.output_dir):
        raise ValueError("canonicalize requires --run-dir")
    landscape_source = _landscape_source_from_args(args)
    _validate_canonicalize_frame_conflicts(args)
    if args.csv_chunk_size <= 0:
        raise ValueError("--csv-chunk-size must be positive")
    fit_top_fraction = _canonicalize_fit_top_fraction(args)
    if not 0.0 < fit_top_fraction <= 1.0:
        raise ValueError("--fit-top-fraction must be > 0 and <= 1")
    positive_side = _canonicalize_positive_side(args)
    output_dir = _prepare_canonicalize_output_dir(
        _canonicalize_output_dir(args),
        overwrite=args.overwrite,
    )
    memory_profiler.sample("output_dir_prepared")
    source_path = resolve_landscape_path(landscape_source)
    source_metadata = read_landscape_metadata(source_path)
    euler_metadata = _canonicalize_euler_metadata()
    memory_profiler.sample("source_metadata_read")
    frame_source_metadata = _read_canonical_frame(args.use_frame) if args.use_frame else None
    policy = _canonicalization_policy_for_command(
        fit_top_fraction=fit_top_fraction,
        positive_side=positive_side,
        frame_metadata=frame_source_metadata,
    )
    frame_applied = frame_source_metadata is not None

    if source_path.suffix.lower() == ".npz":
        canonicalization_backend = "array_native"
        arrays = read_landscape_npz_arrays(source_path)
        memory_profiler.sample("array_landscape_loaded")
        if frame_source_metadata is None:
            canonical_arrays, canonicalization_report = canonicalize_landscape_arrays(
                arrays,
                policy=policy,
            )
        else:
            canonical_arrays = _apply_canonical_frame_arrays(
                arrays,
                frame_source_metadata["canonical_transform"],
            )
            canonicalization_report = _frame_applied_canonicalization_report(
                policy=policy,
                row_count=canonical_arrays.n_points,
                frame_metadata=frame_source_metadata,
            )
        memory_profiler.sample("array_canonicalization_complete")
        output_row_count = canonical_arrays.n_points
        canonical_landscape_path = write_landscape_npz_arrays(
            canonical_arrays,
            output_dir / "canonical_landscape.npz",
            overwrite=args.overwrite,
            artifact_type="canonical_landscape",
        )
        memory_profiler.sample("canonical_npz_written")
        canonical_landscape_csv_path = None
        csv_backend = "skipped"
        if not args.no_csv:
            canonical_landscape_csv_path = write_canonical_landscape_csv_from_npz(
                canonical_landscape_path,
                output_dir / "canonical_landscape.csv",
                overwrite=args.overwrite,
                chunk_size=args.csv_chunk_size,
                euler_sequence=str(euler_metadata["scipy_euler_sequence"]),
            )
            csv_backend = "array_native_chunked"
            memory_profiler.sample("canonical_csv_written")
        visualization_landscape = _landscape_from_arrays_for_visualization(
            canonical_arrays,
            canonicalization_report=canonicalization_report,
        )
        canonical_transform = canonical_arrays.canonical_transform
    else:
        canonicalization_backend = "dataframe_compat"
        landscape = read_landscape(source_path)
        memory_profiler.sample("dataframe_landscape_loaded")
        if frame_source_metadata is None:
            canonical_landscape = canonicalize_landscape(landscape, policy=policy)
        else:
            canonical_landscape = _apply_canonical_frame_landscape(
                landscape,
                frame_source_metadata["canonical_transform"],
                policy=policy,
                frame_metadata=frame_source_metadata,
            )
        memory_profiler.sample("dataframe_canonicalization_complete")
        canonicalization_report = canonical_landscape.canonicalization_report
        output_row_count = len(canonical_landscape.data)
        canonical_landscape_path = write_landscape_npz(
            canonical_landscape,
            output_dir / "canonical_landscape.npz",
            overwrite=args.overwrite,
            artifact_type="canonical_landscape",
        )
        memory_profiler.sample("canonical_npz_written")
        canonical_landscape_csv_path = None
        csv_backend = "skipped"
        if not args.no_csv:
            canonical_landscape_csv_path = write_canonical_landscape_csv(
                canonical_landscape,
                output_dir / "canonical_landscape.csv",
                overwrite=args.overwrite,
                euler_sequence=str(euler_metadata["scipy_euler_sequence"]),
            )
            csv_backend = "dataframe_compat"
            memory_profiler.sample("canonical_csv_written")
        visualization_landscape = canonical_landscape
        canonical_transform = canonical_landscape.canonical_transform

    if canonical_transform is None:
        raise ValueError("canonicalization did not produce canonical_transform")
    frame_paths = _write_canonical_frame_artifacts(
        output_dir=output_dir,
        canonical_transform=canonical_transform,
        policy=policy,
        args=args,
        euler_metadata=euler_metadata,
        source_frame_metadata=frame_source_metadata,
        overwrite=args.overwrite,
    )
    canonicalization_report_payload = _canonicalization_report_payload(
        canonicalization_report,
        frame_paths=frame_paths,
        frame_applied=frame_applied,
        frame_source_metadata=frame_source_metadata,
    )
    canonicalization_report_path = write_report_json(
        canonicalization_report_payload,
        output_dir / "canonicalization_report.json",
        overwrite=args.overwrite,
        artifact_type="canonicalization_report",
    )
    memory_profiler.sample("canonicalization_report_written")
    output_artifacts = {
        "canonical_landscape_npz": str(canonical_landscape_path),
        "canonical_landscape_csv": (
            str(canonical_landscape_csv_path)
            if canonical_landscape_csv_path is not None
            else None
        ),
        "canonicalization_report_json": str(canonicalization_report_path),
        "canonical_frame_json": frame_paths["canonical_frame_json"],
        "canonical_frame_npz": frame_paths["canonical_frame_npz"],
    }
    memory_profile_path = output_dir / "canonical_memory_profile.json"
    if args.profile_memory:
        output_artifacts["memory_profile_json"] = str(memory_profile_path)
    visualization_performed = False
    visualization_report_path = None
    if args.no_visualize:
        pass
    else:
        visualization_report = _write_canonical_default_visualization(
            args=args,
            canonical_landscape=visualization_landscape,
            canonical_landscape_csv_path=canonical_landscape_csv_path,
            euler_metadata=euler_metadata,
            fit_top_fraction=policy.fit_top_fraction,
            frame_source_metadata=frame_source_metadata,
            overwrite=args.overwrite,
        )
        visualization_performed = True
        visualization_report_path = visualization_report["report_path"]
        output_artifacts["canonical_visualization_report_json"] = str(visualization_report_path)

    canonicalize_summary = {
            "artifact_type": "canonicalize_summary",
            "schema_version": "1",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "source_landscape_path": source_metadata["path"],
            "source_run_dir": args.run_dir,
            "source_landscape_artifact_type": source_metadata["artifact_type"],
            "source_landscape_schema_version": source_metadata["schema_version"],
            "source_landscape_row_count": source_metadata["row_count"],
            "canonicalization_policy": to_json_safe(policy),
            "fit_top_fraction": policy.fit_top_fraction,
            "canonicalization_performed": not frame_applied,
            "frame_applied": frame_applied,
            "frame_source_path": (
                frame_source_metadata["frame_source_path"] if frame_source_metadata else None
            ),
            "frame_metadata": (
                frame_source_metadata["metadata"] if frame_source_metadata else None
            ),
            "canonicalization_backend": canonicalization_backend,
            "canonical_landscape_npz": str(canonical_landscape_path),
            "csv_performed": canonical_landscape_csv_path is not None,
            "csv_backend": csv_backend,
            "csv_chunk_size": args.csv_chunk_size,
            "canonical_landscape_csv": (
                str(canonical_landscape_csv_path)
                if canonical_landscape_csv_path is not None
                else None
            ),
            "canonical_frame_json": frame_paths["canonical_frame_json"],
            "canonical_frame_npz": frame_paths["canonical_frame_npz"],
            "visualization_performed": visualization_performed,
            "visualization_report": visualization_report_path,
            "memory_profile_performed": args.profile_memory,
            "memory_profile": str(memory_profile_path) if args.profile_memory else None,
            "selection_performed": False,
            "updated_star_cs_export_performed": False,
            "composite_map_export_performed": False,
            "future_optional_exports": (
                "updated STAR/CS metadata and composite map export are planned "
                "future derived outputs and were not performed."
            ),
            "output_row_count": output_row_count,
            "output_artifacts": output_artifacts,
        }
    canonicalize_summary.update(euler_metadata)
    write_json_artifact(
        canonicalize_summary,
        output_dir / "canonicalize_summary.json",
        overwrite=args.overwrite,
    )
    memory_profiler.sample("canonicalize_summary_written")
    if args.profile_memory:
        memory_profiler.write(memory_profile_path, overwrite=args.overwrite)
    print(str(output_dir))
    return 0


CANONICAL_FRAME_TRANSFORM_DIRECTION = "canonical_rv = raw_rv @ canonical_transform"
CANONICAL_FRAME_COORDINATE_SPACE = "rotvec_ro_radians"


def _canonicalize_fit_top_fraction(args) -> float:
    return 0.40 if args.fit_top_fraction is None else float(args.fit_top_fraction)


def _canonicalize_positive_side(args) -> str:
    if getattr(args, "legacy_positive_side", None):
        return str(args.legacy_positive_side)
    if args.positive_side == "high":
        return "high_density_skew"
    return "low_density_skew"


def _validate_canonicalize_frame_conflicts(args) -> None:
    if not getattr(args, "use_frame", None):
        return
    conflicts = []
    if args.fit_top_fraction is not None:
        conflicts.append("--fit-top/--fit-top-fraction")
    if args.positive_side is not None:
        conflicts.append("--positive-side")
    if getattr(args, "legacy_positive_side", None) is not None:
        conflicts.append("--skewness-positive-side")
    if conflicts:
        raise ValueError("--use-frame conflicts with explicit fitting controls: " + ", ".join(conflicts))


def _canonicalize_euler_metadata() -> dict[str, object]:
    return resolve_euler_convention(
        DEFAULT_EULER_CONVENTION,
        source="public_default",
    ).metadata(euler_angle_columns=CANONICAL_EULER_ANGLE_COLUMNS)


def _canonicalization_policy_for_command(
    *,
    fit_top_fraction: float,
    positive_side: str,
    frame_metadata: dict[str, Any] | None,
) -> CanonicalizationPolicy:
    if frame_metadata is not None:
        metadata = frame_metadata["metadata"]
        fit_top_fraction = float(metadata.get("fit_top_fraction", fit_top_fraction))
        positive_side = str(metadata.get("positive_side", positive_side))
    return CanonicalizationPolicy(
        fit_top_fraction=fit_top_fraction,
        sign_rule="density_weighted_skewness",
        positive_side=positive_side,
    )


def _read_canonical_frame(path: str | None) -> dict[str, Any] | None:
    if not path:
        return None
    frame_path = Path(path)
    if not frame_path.exists():
        raise ValueError(f"Canonical frame does not exist: {frame_path}")
    if frame_path.suffix.lower() != ".json":
        raise ValueError("canonicalize --use-frame expects canonical_frame.json")
    with frame_path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    if not isinstance(metadata, dict):
        raise ValueError("canonical_frame.json must contain a JSON object")
    transform = _canonical_transform_from_frame_metadata(metadata, frame_path)
    transform_direction = metadata.get("transform_direction")
    if transform_direction != CANONICAL_FRAME_TRANSFORM_DIRECTION:
        raise ValueError(
            "canonical_frame transform_direction must be "
            f"{CANONICAL_FRAME_TRANSFORM_DIRECTION!r}"
        )
    coordinate_space = metadata.get("coordinate_space")
    if coordinate_space != CANONICAL_FRAME_COORDINATE_SPACE:
        raise ValueError(
            "canonical_frame coordinate_space must be "
            f"{CANONICAL_FRAME_COORDINATE_SPACE!r}"
        )
    if metadata.get("euler_convention") not in {None, DEFAULT_EULER_CONVENTION}:
        raise ValueError(f"canonical_frame euler_convention must be {DEFAULT_EULER_CONVENTION!r}")
    return {
        "frame_source_path": str(frame_path),
        "canonical_transform": transform,
        "metadata": metadata,
    }


def _canonical_transform_from_frame_metadata(metadata: dict[str, Any], frame_path: Path) -> np.ndarray:
    if "canonical_transform" in metadata:
        transform = np.asarray(metadata["canonical_transform"], dtype=float)
    else:
        npz_value = metadata.get("canonical_frame_npz")
        if not npz_value:
            raise ValueError("canonical_frame.json requires canonical_transform or canonical_frame_npz")
        npz_path = Path(str(npz_value))
        if not npz_path.is_absolute():
            npz_path = frame_path.parent / npz_path
        with np.load(npz_path, allow_pickle=False) as payload:
            if "canonical_transform" not in payload:
                raise ValueError("canonical_frame.npz missing canonical_transform")
            transform = np.asarray(payload["canonical_transform"], dtype=float)
    validate_canonical_transform(transform)
    return transform


def _apply_canonical_frame_arrays(
    arrays: LandscapeArrays,
    canonical_transform: np.ndarray,
) -> LandscapeArrays:
    canonical_coordinates = apply_canonical_transform(
        arrays.coordinates_analysis,
        canonical_transform,
    )
    return LandscapeArrays(
        particle_key=arrays.particle_key,
        coordinates_analysis=arrays.coordinates_analysis,
        coordinates_display=arrays.coordinates_display,
        sld_unfloored=arrays.sld_unfloored,
        sld_raw=arrays.sld_raw,
        sld_display=arrays.sld_display,
        sld_display_is_outlier=arrays.sld_display_is_outlier,
        sld_was_floored=arrays.sld_was_floored,
        sld_local_k_mean=arrays.sld_local_k_mean,
        sld_effective_local_k_mean=arrays.sld_effective_local_k_mean,
        sld_distance_floor=arrays.sld_distance_floor,
        ref_source_row_id=arrays.ref_source_row_id,
        mov_source_row_id=arrays.mov_source_row_id,
        coordinates_canonical=canonical_coordinates,
        canonical_transform=canonical_transform,
    )


def _apply_canonical_frame_landscape(
    landscape: Landscape,
    canonical_transform: np.ndarray,
    *,
    policy: CanonicalizationPolicy,
    frame_metadata: dict[str, Any],
) -> Landscape:
    data = landscape.data.copy(deep=True)
    coordinates = np.vstack([np.asarray(value, dtype=float) for value in data["coordinates_analysis"]])
    canonical_coordinates = apply_canonical_transform(coordinates, canonical_transform)
    data["coordinates_canonical"] = [row.copy() for row in canonical_coordinates]
    active_policies = dict(landscape.active_policies or {})
    active_policies["canonicalization_policy"] = policy
    return Landscape(
        data=data,
        canonical_transform=canonical_transform,
        active_policies=active_policies,
        density_report=landscape.density_report,
        canonicalization_report=_frame_applied_canonicalization_report(
            policy=policy,
            row_count=len(data),
            frame_metadata=frame_metadata,
        ),
    )


def _frame_applied_canonicalization_report(
    *,
    policy: CanonicalizationPolicy,
    row_count: int,
    frame_metadata: dict[str, Any],
) -> dict[str, Any]:
    source_metadata = frame_metadata["metadata"]
    return {
        "fit_subset": "imported_frame",
        "n_fit_points": 0,
        "density_support_field": source_metadata.get("density_support_field", "sld_raw"),
        "singular_values": (),
        "explained_variance_ratios": (),
        "axis_degeneracy_detected": False,
        "fit_top_fraction": policy.fit_top_fraction,
        "axis_assignment": source_metadata.get("axis_assignment", policy.axis_assignment),
        "pca_axis_order": (),
        "assigned_coordinate_names": ("alpha", "beta", "gamma"),
        "assigned_rotvec_columns": ("z", "y", "x"),
        "sign_rule": source_metadata.get("sign_rule", policy.sign_rule),
        "positive_side": source_metadata.get("positive_side", policy.positive_side),
        "sign_weight_field": policy.sign_weight_field,
        "axis_weighted_skewness": (),
        "flipped_axes": (),
        "ambiguous_sign_axes": (),
        "handedness_rule": policy.handedness_rule,
        "handedness_adjusted_axes": (),
        "transform_direction": CANONICAL_FRAME_TRANSFORM_DIRECTION,
        "warnings": (),
        "frame_applied": True,
        "frame_source_path": frame_metadata["frame_source_path"],
        "input_row_count": int(row_count),
    }


def _write_canonical_frame_artifacts(
    *,
    output_dir: Path,
    canonical_transform: np.ndarray,
    policy: CanonicalizationPolicy,
    args,
    euler_metadata: dict[str, object],
    source_frame_metadata: dict[str, Any] | None,
    overwrite: bool,
) -> dict[str, str]:
    frame_npz_path = output_dir / "canonical_frame.npz"
    if frame_npz_path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {frame_npz_path}")
    np.savez(
        frame_npz_path,
        schema_version=np.asarray("1"),
        canonical_transform=np.asarray(canonical_transform, dtype=float),
    )
    frame_json_path = output_dir / "canonical_frame.json"
    metadata = {
        "artifact_type": "canonical_frame",
        "schema_version": "1",
        "canonical_transform": np.asarray(canonical_transform, dtype=float),
        "canonical_frame_npz": str(frame_npz_path),
        "transform_direction": CANONICAL_FRAME_TRANSFORM_DIRECTION,
        "coordinate_space": CANONICAL_FRAME_COORDINATE_SPACE,
        "axis_assignment": policy.axis_assignment,
        "sign_rule": policy.sign_rule,
        "positive_side": policy.positive_side,
        "fit_top_fraction": policy.fit_top_fraction,
        "density_support_field": policy.density_support_field,
        "euler_convention": euler_metadata["euler_convention"],
        "source_run_dir": args.run_dir,
        "source_canonical_id": args.canonical_id,
        "frame_applied": source_frame_metadata is not None,
        "frame_source_path": (
            source_frame_metadata["frame_source_path"] if source_frame_metadata else None
        ),
    }
    if source_frame_metadata is not None:
        metadata["source_frame_metadata"] = source_frame_metadata["metadata"]
    write_json_artifact(metadata, frame_json_path, overwrite=overwrite)
    return {
        "canonical_frame_json": str(frame_json_path),
        "canonical_frame_npz": str(frame_npz_path),
    }


def _canonicalization_report_payload(
    report: Any,
    *,
    frame_paths: dict[str, str],
    frame_applied: bool,
    frame_source_metadata: dict[str, Any] | None,
) -> dict[str, Any]:
    payload = dict(to_json_safe(report))
    payload.update(
        {
            "frame_applied": frame_applied,
            "frame_source_path": (
                frame_source_metadata["frame_source_path"] if frame_source_metadata else None
            ),
            "frame_metadata": (
                frame_source_metadata["metadata"] if frame_source_metadata else None
            ),
            "canonical_frame_json": frame_paths["canonical_frame_json"],
            "canonical_frame_npz": frame_paths["canonical_frame_npz"],
        }
    )
    return payload


def _landscape_from_arrays_for_visualization(
    arrays: LandscapeArrays,
    *,
    canonicalization_report: Any,
) -> Landscape:
    data = pd.DataFrame(
        {
            "particle_key": arrays.particle_key.astype(str),
            "coordinates_analysis": list(arrays.coordinates_analysis),
            "coordinates_display": list(
                arrays.coordinates_display
                if arrays.coordinates_display is not None
                else arrays.coordinates_analysis
            ),
            "coordinates_canonical": list(arrays.coordinates_canonical),
            "sld_unfloored": arrays.sld_unfloored,
            "sld_raw": arrays.sld_raw,
            "sld_display": arrays.sld_display,
            "sld_display_is_outlier": arrays.sld_display_is_outlier,
            "sld_was_floored": arrays.sld_was_floored,
            "sld_local_k_mean": arrays.sld_local_k_mean,
            "sld_effective_local_k_mean": arrays.sld_effective_local_k_mean,
            "sld_distance_floor": arrays.sld_distance_floor,
        }
    )
    if arrays.ref_source_row_id is not None:
        data["ref_source_row_id"] = arrays.ref_source_row_id
    if arrays.mov_source_row_id is not None:
        data["mov_source_row_id"] = arrays.mov_source_row_id
    return Landscape(
        data=data,
        canonical_transform=arrays.canonical_transform,
        canonicalization_report=canonicalization_report,
    )


def _write_canonical_default_visualization(
    *,
    args,
    canonical_landscape: Landscape,
    canonical_landscape_csv_path: Path | None,
    euler_metadata: dict[str, object],
    fit_top_fraction: float,
    frame_source_metadata: dict[str, Any] | None,
    overwrite: bool,
) -> dict[str, Any]:
    output_dir = _canonical_visualization_output_dir(args)
    top_group = f"filter_particles_by_top_sld_{_fit_fraction_pct_label(fit_top_fraction)}pct"
    warnings: list[str] = []
    if frame_source_metadata is not None and "fit_top_fraction" not in frame_source_metadata["metadata"]:
        warnings.append(
            "canonical frame metadata missing fit_top_fraction; "
            "using current canonicalize default/policy for top-SLD preview"
        )

    preview_specs = (
        {
            "label": "all_particles",
            "description": "all canonical landscape rows",
            "display_filter_mode": "none",
            "display_sld_threshold": None,
            "display_top_fraction": None,
        },
        {
            "label": "filter_particles_by_sld_gt_1p5",
            "description": "display-only preview of particles with sld_raw > 1.5",
            "display_filter_mode": "threshold",
            "display_sld_threshold": 1.5,
            "display_top_fraction": None,
        },
        {
            "label": top_group,
            "description": "display-only preview of the top-sld_raw canonical fit support",
            "display_filter_mode": "top_fraction",
            "display_sld_threshold": None,
            "display_top_fraction": fit_top_fraction,
        },
    )

    preview_reports: dict[str, dict[str, Any]] = {}
    generated_files: dict[str, str] = {}
    full_landscape_table = (
        str(canonical_landscape_csv_path) if canonical_landscape_csv_path is not None else None
    )
    for spec in preview_specs:
        label = str(spec["label"])
        max_displayed_sld = _preview_max_displayed_sld(
            canonical_landscape,
            display_density_field="sld_raw",
            display_sld_threshold=spec["display_sld_threshold"],
        )
        resolved_vmax = (
            min(max_displayed_sld, 100.0) if max_displayed_sld is not None else None
        )
        color_vmax = resolved_vmax
        if color_vmax is None and spec["display_sld_threshold"] is not None:
            color_vmax = float(spec["display_sld_threshold"])
        selection_metadata: dict[str, Any] = {
            "preview_group": label,
            "preview_description": spec["description"],
            "preview_only": True,
            "not_a_selection": True,
            "display_filter_mode": spec["display_filter_mode"],
            "display_vmax_cap": 100.0,
            "display_vmax_policy": "min(max_displayed_sld, 100)",
            "max_displayed_sld": max_displayed_sld,
            "resolved_vmax": resolved_vmax,
        }
        if spec["display_sld_threshold"] is not None:
            selection_metadata["display_sld_threshold"] = spec["display_sld_threshold"]
            selection_metadata["display_filter_label"] = "sld_raw > 1.5"
        if spec["display_top_fraction"] is not None:
            selection_metadata["top_sld_fraction"] = spec["display_top_fraction"]
            selection_metadata["fit_support_preview"] = True
        report = write_landscape_visualizations(
            canonical_landscape,
            output_dir / label,
            overwrite=overwrite,
            coordinate_source="canonical",
            representation="both",
            color_field="sld_raw",
            euler_convention=DEFAULT_EULER_CONVENTION,
            euler_convention_source=str(euler_metadata["euler_convention_source"]),
            display_top_fraction=spec["display_top_fraction"],
            display_sld_threshold=spec["display_sld_threshold"],
            display_density_field="sld_raw",
            formats=("png",),
            max_points_2d=500000,
            max_points_3d=50000,
            random_seed=0,
            write_projection_csvs=False,
            full_landscape_table=full_landscape_table,
            display_table_filename="display_table.csv",
            artifact_layout="run_bundle",
            visual_style="legacy",
            color_map="rainbow_r",
            color_vmax=color_vmax,
            display_filter_mode=str(spec["display_filter_mode"]),
            selection_metadata=selection_metadata,
        )
        preview_reports[label] = report
        generated_files[f"{label}_visualization_report_json"] = report["report_path"]
        for artifact_key, artifact_path in report["generated_files"].items():
            generated_files[f"{label}_{artifact_key}"] = artifact_path

    aggregate_report = {
        "artifact_type": "canonical_default_visualization_report",
        "schema_version": "1",
        "status": "ok",
        "canonical_id": args.canonical_id,
        "preview_only": True,
        "not_a_selection": True,
        "visual_style": "legacy_rainbow",
        "color_map": "rainbow_r",
        "color_field": "sld_raw",
        "display_vmax_cap": 100.0,
        "display_vmax_policy": "min(max_displayed_sld, 100)",
        "fit_top_fraction": fit_top_fraction,
        "top_sld_preview_group": top_group,
        "preview_groups": [str(spec["label"]) for spec in preview_specs],
        "preview_reports": {
            label: report["report_path"] for label, report in preview_reports.items()
        },
        "generated_files": generated_files,
        "warnings": warnings,
    }
    report_path = write_json_artifact(
        aggregate_report,
        output_dir / "visualization_report.json",
        overwrite=overwrite,
    )
    aggregate_report["report_path"] = str(report_path)
    return aggregate_report


def _canonical_visualization_output_dir(args) -> Path:
    if args.run_dir:
        return Path(args.run_dir) / "visualizations" / "canonical" / args.canonical_id
    return Path(args.output_dir) / "visualizations"


def _fit_fraction_pct_label(fraction: float) -> str:
    percent = float(fraction) * 100.0
    rounded = round(percent)
    if np.isclose(percent, rounded):
        return str(int(rounded))
    return f"{percent:.3g}".replace(".", "p")


def _write_visualization_1d_distributions(
    *,
    output_dir: Path,
    display_table_path: Path,
    representation: str,
    formats: tuple[str, ...],
    bins: int,
    hist_mode: str,
    kde_bandwidth: str | float,
    xlim: tuple[float, float] | None,
    ylim: tuple[float, float] | None,
    overwrite: bool,
) -> dict[str, Any]:
    if hist_mode not in {"count", "percent"}:
        raise ValueError("--hist-mode must be count or percent")
    if bins <= 0:
        raise ValueError("--bins must be positive")
    if xlim is not None and xlim[0] >= xlim[1]:
        raise ValueError("--xlim requires MIN < MAX")
    if ylim is not None and ylim[0] >= ylim[1]:
        raise ValueError("--ylim requires MIN < MAX")

    table = pd.read_csv(display_table_path)
    distribution_dir = output_dir / "distributions_1d"
    distribution_dir.mkdir(parents=True, exist_ok=True)
    axes = _distribution_1d_axes(representation, table)
    expected_paths = [
        distribution_dir / "distribution_1d_report.json",
        distribution_dir / "distribution_1d_stats.csv",
    ]
    for axis in axes:
        for fmt in formats:
            expected_paths.append(distribution_dir / f"{axis['name']}_distribution.{fmt}")
    if not overwrite:
        for path in expected_paths:
            if path.exists():
                raise FileExistsError(f"Output path already exists: {path}")

    generated_files: dict[str, str] = {}
    stats_rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    for axis in axes:
        values = pd.to_numeric(table[axis["column"]], errors="coerce").to_numpy(dtype=float)
        paths = [distribution_dir / f"{axis['name']}_distribution.{fmt}" for fmt in formats]
        kde_status = _write_distribution_1d_figure(
            values,
            paths=paths,
            axis_label=str(axis["label"]),
            bins=bins,
            hist_mode=hist_mode,
            kde_bandwidth=kde_bandwidth,
            xlim=xlim,
            ylim=ylim,
        )
        if kde_status != "ok":
            warnings.append(f"{axis['name']}: {kde_status}")
        for path in paths:
            generated_files[f"{axis['name']}_{path.suffix.lstrip('.')}"] = str(path)
        stats_rows.append(_distribution_1d_stats(axis["name"], values, kde_status))

    stats_csv = distribution_dir / "distribution_1d_stats.csv"
    pd.DataFrame(stats_rows).to_csv(stats_csv, index=False)
    report = {
        "artifact_type": "distribution_1d_report",
        "schema_version": "1",
        "status": "ok",
        "representation": representation,
        "generated_axes": [axis["name"] for axis in axes],
        "bins": bins,
        "hist_mode": hist_mode,
        "kde_bandwidth": kde_bandwidth,
        "xlim": list(xlim) if xlim is not None else None,
        "ylim": list(ylim) if ylim is not None else None,
        "input_row_count": int(len(table)),
        "displayed_row_count": int(len(table)),
        "formats": list(formats),
        "generated_files": generated_files,
        "stats_csv": str(stats_csv),
        "warnings": warnings,
    }
    report_path = write_json_artifact(
        report,
        distribution_dir / "distribution_1d_report.json",
        overwrite=True,
    )
    report["report_path"] = str(report_path)
    return report


def _distribution_1d_axes(representation: str, table: pd.DataFrame) -> list[dict[str, str]]:
    specs: list[dict[str, str]] = []
    if representation in {"both", "euler"}:
        specs.extend(
            [
                {"name": "euler_alpha", "column": "euler_alpha", "label": "alpha (deg)"},
                {"name": "euler_beta", "column": "euler_beta", "label": "beta (deg)"},
                {"name": "euler_gamma", "column": "euler_gamma", "label": "gamma (deg)"},
            ]
        )
    if representation in {"both", "rotvec"}:
        specs.extend(
            [
                {"name": "rotvec_x", "column": "rotvec_x", "label": "rotvec x (rad)"},
                {"name": "rotvec_y", "column": "rotvec_y", "label": "rotvec y (rad)"},
                {"name": "rotvec_z", "column": "rotvec_z", "label": "rotvec z (rad)"},
            ]
        )
    missing = [spec["column"] for spec in specs if spec["column"] not in table.columns]
    if missing:
        raise ValueError(f"display_table.csv missing 1D distribution columns: {missing}")
    return specs


def _write_distribution_1d_figure(
    values: np.ndarray,
    *,
    paths: list[Path],
    axis_label: str,
    bins: int,
    hist_mode: str,
    kde_bandwidth: str | float,
    xlim: tuple[float, float] | None,
    ylim: tuple[float, float] | None,
) -> str:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde

    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    kde_status = "ok"
    if finite.size:
        weights = None
        density = False
        if hist_mode == "percent":
            weights = np.full(finite.size, 100.0 / finite.size)
        counts, edges, _ = ax.hist(
            finite,
            bins=bins,
            range=xlim,
            weights=weights,
            density=density,
            color="silver",
            edgecolor="black",
            alpha=0.55,
        )
        try:
            if finite.size < 2 or np.isclose(float(np.nanstd(finite)), 0.0):
                raise ValueError("not enough spread for KDE")
            xmin, xmax = xlim if xlim is not None else (float(np.min(finite)), float(np.max(finite)))
            if np.isclose(xmin, xmax):
                xmin -= 1.0
                xmax += 1.0
            grid = np.linspace(float(xmin), float(xmax), 512)
            kde = gaussian_kde(finite, bw_method=kde_bandwidth)
            kde_values = kde(grid)
            bin_width = float(edges[1] - edges[0]) if len(edges) > 1 else (float(xmax) - float(xmin)) / bins
            if hist_mode == "percent":
                kde_values = kde_values * 100.0 * bin_width
            else:
                kde_values = kde_values * finite.size * bin_width
            ax.plot(grid, kde_values, color="black", linewidth=2)
        except Exception as exc:
            kde_status = f"skipped: {exc}"
    else:
        kde_status = "skipped: no finite values"
        ax.hist([], bins=bins, color="silver", edgecolor="black", alpha=0.55)
    ax.set_xlabel(axis_label)
    ax.set_ylabel("Count" if hist_mode == "count" else "Percent")
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    for path in paths:
        fig.savefig(path)
    plt.close(fig)
    return kde_status


def _distribution_1d_stats(axis_name: str, values: np.ndarray, kde_status: str) -> dict[str, Any]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    row: dict[str, Any] = {
        "axis": axis_name,
        "n": int(finite.size),
        "kde_status": kde_status,
    }
    if not finite.size:
        row.update({key: None for key in ("min", "max", "median", "q05", "q25", "q75", "q95")})
        return row
    row.update(
        {
            "min": float(np.min(finite)),
            "max": float(np.max(finite)),
            "median": float(np.median(finite)),
            "q05": float(np.quantile(finite, 0.05)),
            "q25": float(np.quantile(finite, 0.25)),
            "q75": float(np.quantile(finite, 0.75)),
            "q95": float(np.quantile(finite, 0.95)),
        }
    )
    return row


def export_selection_command(args) -> int:
    if getattr(args, "run_dir", None):
        return export_metadata_command(args)
    return _export_selection_artifact_command(args)


def _export_selection_artifact_command(args) -> int:
    _validate_export_selection_inputs(args)
    selection = read_selection_json(args.selection)
    output_dir = _prepare_selection_export_output_dir(
        args.output_dir,
        overwrite=args.overwrite,
    )
    export_selection(
        selection,
        policy=SelectionExportPolicy(
            output_dir=output_dir,
            overwrite=args.overwrite,
        ),
    )
    print(str(output_dir))
    return 0


def export_metadata_command(args) -> int:
    """Export source-format metadata subsets for an existing Selection."""

    if getattr(args, "selection", None) and not getattr(args, "run_dir", None):
        return _export_selection_artifact_command(args)
    _validate_export_selection_inputs(args)
    selection_path = _resolve_export_selection_path(args)
    selection = read_selection_json(selection_path)
    output_dir = Path(args.output_dir) if getattr(args, "output_dir", None) else None
    report = export_selection_metadata_subset(
        selection,
        policy=SelectionMetadataExportPolicy(
            output_dir=output_dir,
            overwrite=args.overwrite,
            domain=args.domain,
            format=args.format,
            run_dir=args.run_dir,
        ),
        selection_path=selection_path,
    )
    outputs = report.get("outputs", {})
    selected_keys_path = outputs.get("selected_particle_keys_txt")
    if selected_keys_path:
        print(str(Path(selected_keys_path).parent))
    else:
        print(str(output_dir or Path(args.run_dir) / "exports" / selection.selection_id))
    return 0


def _resolve_export_selection_path(args) -> Path:
    if getattr(args, "selection", None):
        path = Path(args.selection)
        if not path.exists():
            raise ValueError(
                f"Selection JSON file does not exist: {path}. "
                "Provide --run-dir RUN --selection-id ID, or advanced "
                "--selection PATH/to/selection.json."
            )
        return path
    path = Path(args.run_dir) / "selections" / args.selection_id / "selection.json"
    if not path.exists():
        raise ValueError(
            f"Selection JSON file does not exist: {path}. "
            "Provide --run-dir RUN --selection-id ID, or advanced "
            "--selection PATH/to/selection.json."
        )
    return path


def _validate_export_selection_inputs(args) -> None:
    has_selection = bool(getattr(args, "selection", None))
    has_selection_id = bool(getattr(args, "selection_id", None))
    has_run_dir = bool(getattr(args, "run_dir", None))
    has_output_dir = bool(getattr(args, "output_dir", None))
    compact_hint = (
        "Provide --run-dir RUN --selection-id ID, or advanced "
        "--selection PATH/to/selection.json."
    )
    if has_selection and has_selection_id:
        raise ValueError("Use either --selection or --selection-id, not both.")
    if has_selection:
        if not has_run_dir and not has_output_dir:
            raise ValueError(
                "Direct --selection export requires --output-dir, or use "
                "--run-dir RUN --selection-id ID."
            )
        return
    if not has_run_dir:
        raise ValueError(compact_hint)
    if not has_selection_id:
        raise ValueError(compact_hint)


def _landscape_source_from_args(args):
    if getattr(args, "run_dir", None):
        return args.run_dir
    if getattr(args, "landscape", None):
        return args.landscape
    raise ValueError("Provide --run-dir or --landscape")


def _visualization_output_dir(args) -> Path:
    run_dir = Path(args.run_dir)
    visual_id = getattr(args, "visual_id", None) or "default"
    if getattr(args, "selection_id", None):
        if getattr(args, "use_selected_landscape", False):
            return run_dir / "visualizations" / "selections" / args.selection_id / "selected_landscape" / visual_id
        if args.space == "canonical":
            return (
                run_dir
                / "visualizations"
                / "selections"
                / args.selection_id
                / "parent_canonical"
                / args.canonical_id
                / visual_id
            )
        return run_dir / "visualizations" / "selections" / args.selection_id / "parent_raw" / visual_id
    if args.space == "canonical":
        return run_dir / "visualizations" / "canonical" / args.canonical_id / visual_id
    return run_dir / "visualizations" / "raw" / visual_id


def _validate_visualize_color_field_args(args) -> None:
    if args.vmin is not None and not np.isfinite(args.vmin):
        raise ValueError("--vmin must be finite")
    if args.vmax is not None and not np.isfinite(args.vmax):
        raise ValueError("--vmax must be finite")
    if args.vmin is not None and args.vmax is not None and args.vmin >= args.vmax:
        raise ValueError("--vmin must be less than --vmax")


def _validate_visualize_color_field(landscape: Landscape, color_field: str) -> None:
    if color_field not in landscape.data.columns:
        raise ValueError(f"visualize requires landscape color field: {color_field}")
    values = pd.to_numeric(landscape.data[color_field], errors="coerce")
    if values.notna().sum() == 0:
        raise ValueError(f"visualize color field must be numeric or boolean-like: {color_field}")


def _visualize_axis_limits_from_ranges(
    range_bounds: dict[str, tuple[float | None, float | None]],
) -> dict[str, tuple[float, float]] | None:
    limits: dict[str, tuple[float, float]] = {}
    for axis, bounds in range_bounds.items():
        lower, upper = bounds
        if lower is None or upper is None or lower >= upper:
            continue
        limits[axis] = (float(lower), float(upper))
    return limits or None


def _parse_kde_bandwidth(value: str | float | None) -> str | float:
    if value is None:
        return "scott"
    if isinstance(value, (int, float)):
        parsed = float(value)
        if parsed <= 0 or not np.isfinite(parsed):
            raise ValueError("--kde-bandwidth float must be positive")
        return parsed
    normalized = str(value).strip().lower()
    if normalized in {"scott", "silverman"}:
        return normalized
    try:
        parsed = float(normalized)
    except ValueError as exc:
        raise ValueError("--kde-bandwidth must be scott, silverman, or a positive float") from exc
    if parsed <= 0 or not np.isfinite(parsed):
        raise ValueError("--kde-bandwidth float must be positive")
    return parsed


def _visualization_public_report(
    report: dict[str, Any],
    *,
    args,
    source_landscape_path: Path,
    euler_metadata: dict[str, object],
    formats: tuple[str, ...],
    distribution_report: dict[str, Any],
) -> dict[str, Any]:
    report = dict(report)
    report.update(
        {
            "artifact_type": "visualization_report",
            "source_landscape_path": str(source_landscape_path),
            "space": args.space,
            "canonical_id": args.canonical_id if args.space == "canonical" else None,
            "visual_id": args.visual_id,
            "inherited_euler_convention": euler_metadata["euler_convention"],
            "inherited_euler_convention_source": euler_metadata["euler_convention_source"],
            "display_density_field": "sld_display",
            "colormap": args.colormap,
            "top_fraction": args.top_fraction,
            "threshold": args.threshold,
            "formats": list(formats),
            "vmin": args.vmin,
            "vmax": args.vmax,
            "distribution_1d_policy": {
                "enabled": True,
                "bins": args.bins,
                "hist_mode": args.hist_mode,
                "kde_bandwidth": distribution_report["kde_bandwidth"],
                "xlim": list(args.xlim) if args.xlim is not None else None,
                "ylim": list(args.ylim) if args.ylim is not None else None,
            },
            "distribution_1d_report": distribution_report["report_path"],
        }
    )
    generated_files = dict(report.get("generated_files") or {})
    generated_files["distribution_1d_report_json"] = distribution_report["report_path"]
    generated_files["distribution_1d_stats_csv"] = distribution_report["stats_csv"]
    for key, value in distribution_report["generated_files"].items():
        generated_files[f"distribution_1d_{key}"] = value
    report["generated_files"] = generated_files
    return report


def _write_selected_derived_landscape(
    *,
    parent_landscape: Landscape,
    selection,
    output_dir: Path,
    overwrite: bool,
    recompute_sld: bool,
    density_policy: DensityPolicy,
    coordinate_space: str,
    parent_landscape_metadata: dict[str, object],
    euler_metadata: dict[str, object],
) -> dict[str, object]:
    selected_data = _selected_rows_from_landscape(parent_landscape, selection.selected_particle_keys)
    density_source = "recomputed_on_selection" if recompute_sld else "parent_landscape"
    effective_k = None
    display_diagnostics: dict[str, object] = {}
    if recompute_sld:
        selected_data = _with_parent_sld_fields(selected_data)
        coordinates = _density_coordinates_for_selected_space(selected_data, coordinate_space)
        if len(coordinates) < 2:
            raise ValueError(
                "SLD recomputation on a selected-derived landscape requires at least "
                "two selected particles"
            )
        sld = compute_sld_values(
            coordinates,
            k_neighbors=density_policy.k_neighbors,
            distance_floor_fraction=density_policy.distance_floor_fraction,
        )
        display = compute_sld_display_values(
            sld["sld_raw"],
            mode=density_policy.display_normalization_mode,
            outlier_mode=density_policy.display_outlier_mode,
            tail_search_fraction=density_policy.tail_search_fraction,
            tail_jump_factor=density_policy.tail_jump_factor,
            max_display_outlier_fraction=density_policy.max_display_outlier_fraction,
        )
        selected_data["sld_unfloored"] = sld["sld_unfloored"]
        selected_data["sld_raw"] = sld["sld_raw"]
        selected_data["sld_display"] = display["sld_display"]
        selected_data["sld_display_is_outlier"] = display["sld_display_is_outlier"]
        selected_data["sld_was_floored"] = sld["sld_was_floored"]
        selected_data["sld_local_k_mean"] = sld["sld_local_k_mean"]
        selected_data["sld_effective_local_k_mean"] = sld["sld_effective_local_k_mean"]
        selected_data["sld_distance_floor"] = sld["sld_distance_floor"]
        effective_k = min(density_policy.k_neighbors, len(selected_data) - 1)
        display_diagnostics = {
            key: value
            for key, value in display.items()
            if key != "sld_display" and key != "sld_display_is_outlier"
        }

    selected_landscape = Landscape(
        data=selected_data,
        canonical_transform=parent_landscape.canonical_transform,
        active_policies=dict(parent_landscape.active_policies or {}),
        density_report=None,
        canonicalization_report=parent_landscape.canonicalization_report,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    npz_path = write_landscape_npz(
        selected_landscape,
        output_dir / "landscape.npz",
        overwrite=overwrite,
        artifact_type="selected_landscape",
    )
    csv_path = _write_selected_landscape_csv(
        selected_landscape,
        output_dir / "landscape.csv",
        overwrite=overwrite,
        euler_sequence=str(euler_metadata["scipy_euler_sequence"]),
    )
    rows_path = _write_selected_landscape_rows_csv(
        selected_landscape,
        output_dir.parent / "selected_landscape_rows.csv",
        overwrite=overwrite,
    )
    report = {
        "artifact_type": "selected_landscape_report",
        "schema_version": "1",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "selection_id": selection.selection_id,
        "parent_landscape_path": parent_landscape_metadata.get("path"),
        "parent_landscape_artifact_type": parent_landscape_metadata.get("artifact_type"),
        "parent_landscape_schema_version": parent_landscape_metadata.get("schema_version"),
        "parent_landscape_row_count": parent_landscape_metadata.get("row_count"),
        "selected_count": selection.selected_count,
        "coordinate_space": coordinate_space,
        "density_source": density_source,
        "sld_recomputed": recompute_sld,
        "requested_k_neighbors": density_policy.k_neighbors if recompute_sld else None,
        "effective_k_neighbors": effective_k,
        "distance_floor_fraction": (
            density_policy.distance_floor_fraction if recompute_sld else None
        ),
        "display_outlier_policy": (
            {
                "display_outlier_mode": density_policy.display_outlier_mode,
                "tail_search_fraction": density_policy.tail_search_fraction,
                "tail_jump_factor": density_policy.tail_jump_factor,
                "max_display_outlier_fraction": density_policy.max_display_outlier_fraction,
            }
            if recompute_sld
            else None
        ),
        "display_diagnostics": display_diagnostics,
        "parent_sld_fields_preserved": recompute_sld,
        "euler_metadata": euler_metadata,
        "output_paths": {
            "landscape_npz": str(npz_path),
            "landscape_csv": str(csv_path),
            "landscape_report_json": str(output_dir / "landscape_report.json"),
            "selected_landscape_rows_csv": str(rows_path),
        },
    }
    report_path = write_json_artifact(
        report,
        output_dir / "landscape_report.json",
        overwrite=overwrite,
    )
    report["output_paths"]["landscape_report_json"] = str(report_path)
    return report


def _selected_rows_from_landscape(
    landscape: Landscape,
    selected_particle_keys: Sequence[object],
) -> pd.DataFrame:
    if "particle_key" not in landscape.data.columns:
        raise ValueError("Landscape must contain particle_key to write selected landscape")
    data = landscape.data.copy(deep=True)
    key_to_indices: dict[str, list[int]] = {}
    for row_index, particle_key in enumerate(data["particle_key"].map(str)):
        key_to_indices.setdefault(particle_key, []).append(row_index)
    selected_indices = []
    for particle_key in selected_particle_keys:
        key = str(particle_key)
        if key not in key_to_indices:
            raise ValueError(
                "Selection contains particle_key values absent from the parent landscape: "
                f"{key}"
            )
        selected_indices.extend(key_to_indices[key])
    if not selected_indices:
        raise ValueError("Selected-derived landscape cannot be empty")
    selected = data.iloc[selected_indices].reset_index(drop=True)
    if "sld_display_is_outlier" not in selected.columns:
        selected["sld_display_is_outlier"] = False
    return selected


def _with_parent_sld_fields(data: pd.DataFrame) -> pd.DataFrame:
    output = data.copy(deep=True)
    for field in (
        "sld_unfloored",
        "sld_raw",
        "sld_display",
        "sld_display_is_outlier",
        "sld_was_floored",
        "sld_local_k_mean",
        "sld_effective_local_k_mean",
        "sld_distance_floor",
    ):
        if field in output.columns:
            output[f"parent_{field}"] = output[field]
    return output


def _density_coordinates_for_selected_space(data: pd.DataFrame, coordinate_space: str) -> np.ndarray:
    if coordinate_space == "canonical":
        if "coordinates_canonical" not in data.columns:
            raise ValueError("canonical selected SLD recomputation requires coordinates_canonical")
        column = "coordinates_canonical"
    else:
        column = "coordinates_analysis"
    return np.vstack([np.asarray(value, dtype=float) for value in data[column]])


def _write_selected_landscape_csv(
    landscape: Landscape,
    path: Path,
    *,
    overwrite: bool,
    euler_sequence: str,
) -> Path:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    data = landscape.data
    raw = np.vstack([np.asarray(value, dtype=float) for value in data["coordinates_analysis"]])
    raw_euler = Rotation.from_rotvec(raw).as_euler(euler_sequence, degrees=True)
    output = pd.DataFrame(
        {
            "particle_key": data["particle_key"].map(str),
            "ref_source_row_id": data.get("ref_source_row_id", -1),
            "mov_source_row_id": data.get("mov_source_row_id", -1),
            "raw_rv_x_rad": raw[:, 0],
            "raw_rv_y_rad": raw[:, 1],
            "raw_rv_z_rad": raw[:, 2],
            "raw_angle_deg": np.degrees(np.linalg.norm(raw, axis=1)),
            "raw_ea_zyx_alpha_deg": raw_euler[:, 0],
            "raw_ea_zyx_beta_deg": raw_euler[:, 1],
            "raw_ea_zyx_gamma_deg": raw_euler[:, 2],
            "sld_unfloored": data["sld_unfloored"],
            "sld_raw": data["sld_raw"],
            "sld_display": data["sld_display"],
            "sld_display_is_outlier": data["sld_display_is_outlier"],
            "sld_was_floored": data["sld_was_floored"],
            "sld_local_k_mean": data["sld_local_k_mean"],
            "sld_effective_local_k_mean": data["sld_effective_local_k_mean"],
            "sld_distance_floor": data["sld_distance_floor"],
        }
    )
    if "coordinates_canonical" in data.columns:
        canonical = np.vstack(
            [np.asarray(value, dtype=float) for value in data["coordinates_canonical"]]
        )
        canonical_euler = Rotation.from_rotvec(canonical).as_euler(
            euler_sequence,
            degrees=True,
        )
        output["canonical_rv_x_rad"] = canonical[:, 0]
        output["canonical_rv_y_rad"] = canonical[:, 1]
        output["canonical_rv_z_rad"] = canonical[:, 2]
        output["canonical_ea_zyx_alpha_deg"] = canonical_euler[:, 0]
        output["canonical_ea_zyx_beta_deg"] = canonical_euler[:, 1]
        output["canonical_ea_zyx_gamma_deg"] = canonical_euler[:, 2]
    for column in (
        "parent_sld_unfloored",
        "parent_sld_raw",
        "parent_sld_display",
        "parent_sld_display_is_outlier",
        "parent_sld_was_floored",
        "parent_sld_local_k_mean",
        "parent_sld_effective_local_k_mean",
        "parent_sld_distance_floor",
    ):
        if column in data.columns:
            output[column] = data[column]
    output.to_csv(path, index=False)
    return path


def _write_selected_landscape_rows_csv(
    landscape: Landscape,
    path: Path,
    *,
    overwrite: bool,
) -> Path:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        column
        for column in ("particle_key", "ref_source_row_id", "mov_source_row_id")
        if column in landscape.data.columns
    ]
    landscape.data[columns].to_csv(path, index=False)
    return path


def _write_selection_csv(selection, path: Path, *, overwrite: bool) -> Path:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = pd.DataFrame(
        {
            "selection_rank": range(1, len(selection.selected_particle_keys) + 1),
            "particle_key": list(selection.selected_particle_keys),
        }
    )
    rows.to_csv(path, index=False)
    return path


def _resolve_selected_landscape_path(run_dir: Path, selection_id: str) -> Path:
    selection_dir = run_dir / "selections" / selection_id / "selected_landscape"
    npz_path = selection_dir / "landscape.npz"
    csv_path = selection_dir / "landscape.csv"
    if npz_path.exists():
        return npz_path
    if csv_path.exists():
        return csv_path
    raise ValueError(
        "Selected-derived landscape not found; expected landscape.npz or "
        f"landscape.csv under {selection_dir}"
    )


def _selected_landscape_visualization_metadata(
    run_dir: Path,
    selection_id: str,
    selected_landscape_path: Path,
    *,
    selected_count: int,
) -> dict[str, object]:
    report_path = run_dir / "selections" / selection_id / "selected_landscape" / "landscape_report.json"
    report = _read_json_if_exists(report_path) or {}
    output_paths = report.get("output_paths") if isinstance(report, dict) else None
    return {
        "selection_id": selection_id,
        "selected_count": selected_count,
        "total_count": report.get("parent_landscape_row_count"),
        "selection_filter_applied": False,
        "selected_landscape_used": True,
        "selected_landscape_path": str(selected_landscape_path),
        "selected_landscape_report_path": str(report_path) if report_path.exists() else None,
        "selected_landscape_parent_path": report.get("parent_landscape_path"),
        "parent_landscape_path": report.get("parent_landscape_path"),
        "selected_landscape_density_source": report.get("density_source"),
        "density_source": report.get("density_source"),
        "selected_landscape_sld_recomputed": report.get("sld_recomputed"),
        "sld_recomputed": report.get("sld_recomputed"),
        "selected_landscape_output_paths": output_paths if isinstance(output_paths, dict) else None,
    }


def _filter_landscape_by_selection(args, landscape: Landscape) -> tuple[Landscape, dict[str, object]]:
    if not getattr(args, "run_dir", None):
        raise ValueError("visualize --selection-id requires --run-dir")
    selection_keys, selection_metadata = _read_visualization_selection_keys(
        Path(args.run_dir),
        args.selection_id,
        landscape_count=len(landscape.data),
    )
    if "particle_key" not in landscape.data.columns:
        raise ValueError("Landscape must contain particle_key for --selection-id filtering")

    selected_key_strings = [str(key) for key in selection_keys]
    selected_key_set = set(selected_key_strings)
    landscape_key_strings = landscape.data["particle_key"].map(str)
    missing_keys = sorted(selected_key_set.difference(set(landscape_key_strings)))
    if missing_keys:
        preview = ", ".join(missing_keys[:5])
        suffix = "" if len(missing_keys) <= 5 else f", ... ({len(missing_keys)} total)"
        raise ValueError(
            "Selection contains particle_key values absent from the requested landscape: "
            f"{preview}{suffix}"
        )

    filtered_data = landscape.data.loc[landscape_key_strings.isin(selected_key_set)].copy(deep=True)
    if filtered_data.empty:
        raise ValueError("Selection filtering produced an empty landscape")
    filtered_landscape = Landscape(
        data=filtered_data.reset_index(drop=True),
        canonical_transform=landscape.canonical_transform,
        active_policies=landscape.active_policies,
        density_report=None,
        canonicalization_report=landscape.canonicalization_report,
    )
    selection_metadata.update(
        {
            "selection_filter_applied": True,
            "selected_landscape_used": False,
            "landscape_count_before_selection_filter": int(len(landscape.data)),
            "landscape_count_after_selection_filter": int(len(filtered_data)),
            "missing_selected_key_count": 0,
        }
    )
    return filtered_landscape, selection_metadata


def _read_visualization_selection_keys(
    run_dir: Path,
    selection_id: str,
    *,
    landscape_count: int,
) -> tuple[list[object], dict[str, object]]:
    selection_dir = run_dir / "selections" / selection_id
    if not selection_dir.exists():
        raise ValueError(f"Selection directory does not exist: {selection_dir}")
    if not selection_dir.is_dir():
        raise ValueError(f"Selection path is not a directory: {selection_dir}")

    selected_keys_csv = selection_dir / "selected_particle_keys.csv"
    selection_json = selection_dir / "selection.json"
    if selected_keys_csv.exists():
        keys = _read_selected_particle_keys_csv(selected_keys_csv)
        metadata: dict[str, object] = {
            "selection_id": selection_id,
            "selected_count": int(len(keys)),
            "total_count": int(landscape_count),
            "selection_source_path": str(selected_keys_csv),
            "selection_count_source": "selected_particle_keys_csv",
            "selection_total_count_source": "landscape_row_count",
        }
        if selection_json.exists():
            counts = _read_selection_json_counts(selection_json)
            if counts.get("selected_count") is not None:
                metadata["selected_count"] = int(counts["selected_count"])
                if len(keys) != counts["selected_count"]:
                    raise ValueError(
                        "selected_particle_keys.csv row count does not match "
                        "selection.json selected_count"
                    )
            if counts.get("total_count") is not None:
                metadata["total_count"] = int(counts["total_count"])
                metadata["selection_total_count_source"] = "selection_json"
            if counts.get("selection_id") not in (None, selection_id):
                raise ValueError(
                    "selection.json selection_id does not match requested --selection-id"
                )
        return keys, metadata

    if selection_json.exists():
        selection = read_selection_json(selection_json)
        return list(selection.selected_particle_keys), {
            "selection_id": selection_id,
            "selected_count": int(selection.selected_count),
            "total_count": int(selection.total_count),
            "selection_source_path": str(selection_json),
            "selection_filter_applied": True,
            "selection_count_source": "selection_json",
            "selection_total_count_source": "selection_json",
        }

    raise ValueError(
        "Selection artifacts not found; expected selected_particle_keys.csv "
        f"or selection.json under {selection_dir}"
    )


def _read_selected_particle_keys_csv(path: Path) -> list[str]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Selected particle keys CSV is empty: {path}")
        if "particle_key" not in reader.fieldnames:
            raise ValueError("selected_particle_keys.csv must contain a particle_key column")
        keys = [row["particle_key"] for row in reader if row.get("particle_key") not in (None, "")]
    if not keys:
        raise ValueError(f"Selected particle keys CSV contains no particle keys: {path}")
    return keys


def _read_selection_json_counts(path: Path) -> dict[str, object]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Malformed selection JSON: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError("Selection JSON must be an object")

    counts: dict[str, object] = {"selection_id": payload.get("selection_id")}
    for field_name in ("selected_count", "total_count"):
        value = payload.get(field_name)
        if value is None:
            counts[field_name] = None
            continue
        if not isinstance(value, int) or value < 0:
            raise ValueError(f"Selection JSON {field_name} must be a non-negative integer")
        counts[field_name] = value
    selected_count = counts.get("selected_count")
    total_count = counts.get("total_count")
    if (
        isinstance(selected_count, int)
        and isinstance(total_count, int)
        and selected_count > total_count
    ):
        raise ValueError("Selection JSON selected_count must be <= total_count")
    return counts


def _canonicalize_output_dir(args) -> Path:
    if args.output_dir:
        return Path(args.output_dir)
    return Path(args.run_dir) / "canonical" / args.canonical_id


def _selection_output_dir(args, *, policy: SelectionPolicy | None = None) -> Path:
    if args.run_dir:
        selection_id = policy.selection_id if policy is not None else args.selection_id
        return Path(args.run_dir) / "selections" / str(selection_id)
    raise ValueError("select requires --run-dir")


def _coordinate_source_for_landscape(args) -> str:
    if args.run_dir and args.space == "raw":
        return "analysis"
    if args.run_dir and args.space == "canonical":
        return "canonical"
    return args.coordinate_source


def _resolve_run_euler_metadata(args) -> dict[str, object]:
    source = "cli_override" if args.euler_convention else "cli_default"
    resolved = resolve_euler_convention(args.euler_convention, source=source)
    return resolved.metadata(euler_angle_columns=RAW_EULER_ANGLE_COLUMNS)


def _resolve_landscape_euler_metadata(
    args,
    *,
    columns=RAW_EULER_ANGLE_COLUMNS,
    explicit_sequence: str | None = None,
    warn_on_legacy_missing: bool = True,
) -> dict[str, object]:
    if getattr(args, "euler_convention", None):
        resolved = resolve_euler_convention(args.euler_convention, source="cli_override")
        return resolved.metadata(euler_angle_columns=columns)
    if explicit_sequence:
        resolved = resolve_euler_convention(
            scipy_euler_sequence=explicit_sequence,
            source="cli_override",
        )
        return resolved.metadata(euler_angle_columns=columns)
    inherited = _read_parent_euler_metadata(args, columns=columns)
    if inherited is not None:
        return inherited
    resolved = resolve_euler_convention(
        DEFAULT_EULER_CONVENTION,
        source=LEGACY_MISSING_EULER_CONVENTION_SOURCE,
    )
    metadata = resolved.metadata(euler_angle_columns=columns)
    if warn_on_legacy_missing:
        print(
            "[cryorole] warning: parent landscape lacks Euler convention metadata; "
            f"defaulting to {DEFAULT_EULER_CONVENTION} for derived Euler coordinates.",
            file=sys.stderr,
        )
    return metadata


def _read_parent_euler_metadata(args, *, columns) -> dict[str, object] | None:
    run_dir_raw = getattr(args, "run_dir", None)
    if not run_dir_raw:
        return None
    run_dir = Path(run_dir_raw)
    payloads: list[dict[str, object]] = []
    if getattr(args, "space", "raw") == "canonical":
        canonical_id = getattr(args, "canonical_id", "default")
        payload = _read_json_if_exists(
            run_dir / "canonical" / canonical_id / "canonicalize_summary.json"
        )
        if payload is not None:
            payloads.append(payload)
    run_summary = _read_json_if_exists(run_dir / "run_summary.json")
    if run_summary is not None:
        payloads.append(run_summary)
    for payload in payloads:
        convention = payload.get("euler_convention")
        sequence = payload.get("scipy_euler_sequence")
        if convention or sequence:
            resolved = resolve_euler_convention(
                convention if isinstance(convention, str) else None,
                scipy_euler_sequence=sequence if isinstance(sequence, str) else None,
                source="inherited",
            )
            return resolved.metadata(euler_angle_columns=columns)
    return None


def _read_json_if_exists(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        return None
    return payload


def _identity_policy_from_args(args, *, source_type: str | None = None) -> IdentityPolicy | None:
    if getattr(args, "row_aligned", False):
        return IdentityPolicy(identity_mode="row_aligned")
    if args.identity_mode is None:
        if source_type == "cryosparc":
            return IdentityPolicy.cryosparc_uid()
        if source_type == "relion":
            return IdentityPolicy.relion_image_name()
        return None
    if args.identity_mode == "cryosparc_uid":
        return IdentityPolicy.cryosparc_uid()
    if args.identity_mode == "relion_user_columns":
        if not args.identity_column:
            raise ValueError("relion_user_columns identity requires --identity-column")
        return IdentityPolicy.relion_user_columns(tuple(args.identity_column))
    if args.identity_mode == "explicit_mapping_file":
        if not args.mapping_file:
            raise ValueError("explicit_mapping_file identity requires --mapping-file")
        return IdentityPolicy(
            identity_mode="explicit_mapping_file",
            mapping_file=args.mapping_file,
        )
    if args.identity_mode == "row_aligned":
        return IdentityPolicy(identity_mode="row_aligned")
    raise ValueError(f"Unsupported identity mode: {args.identity_mode}")


def _convention_policy_for_source(source_type: str) -> ConventionPolicy | None:
    if source_type == "relion":
        return ConventionPolicy.relion_default()
    return None


def _run_summary_payload(
    *,
    args,
    source_ref: str,
    source_mov: str,
    phase1,
    landscape,
    matched_count: int,
    k_neighbors: int,
    canonicalization_performed: bool,
    output_artifacts: dict[str, str],
    run_backend_resolved: str,
    raw_csv_performed: bool,
    raw_csv_backend: str,
    raw_visualization_performed: bool,
    timing_profile_path: Path | None,
    memory_profile_path: Path | None,
    euler_metadata: dict[str, object],
) -> dict[str, Any]:
    payload = {
        "artifact_type": "run_summary",
        "schema_version": "1",
        "input_paths": {"ref": args.ref, "mov": args.mov},
        "source_types": {"ref": source_ref, "mov": source_mov},
        "ref_domain": args.ref_domain,
        "mov_domain": args.mov_domain,
        "row_counts": {
            "ref": len(phase1.pose_a.data),
            "mov": len(phase1.pose_b.data),
        },
        "matched_count": int(matched_count),
        "match_key": _get_field(phase1.match_report, "match_key", None),
        "ref_row_count": _get_field(phase1.match_report, "ref_row_count", len(phase1.pose_a.data)),
        "mov_row_count": _get_field(phase1.match_report, "mov_row_count", len(phase1.pose_b.data)),
        "matched_row_count": _get_field(phase1.match_report, "matched_row_count", matched_count),
        "dropped_ref_only_count": _get_field(
            phase1.match_report,
            "dropped_ref_only_count",
            0,
        ),
        "dropped_mov_only_count": _get_field(
            phase1.match_report,
            "dropped_mov_only_count",
            0,
        ),
        "matched_rows_reordered": _get_field(
            phase1.match_report,
            "matched_rows_reordered",
            False,
        ),
        "match_warnings": list(_get_field(phase1.match_report, "warnings", ())),
        "landscape_row_count": len(landscape.data),
        "k_neighbors": k_neighbors,
        "canonicalization_performed": canonicalization_performed,
        "selection_performed": False,
        "run_backend_requested": args.run_backend,
        "run_backend_resolved": run_backend_resolved,
        "raw_landscape_npz": output_artifacts.get("raw_landscape_npz"),
        "raw_csv_performed": raw_csv_performed,
        "raw_csv_backend": raw_csv_backend,
        "raw_csv_chunk_size": args.raw_csv_chunk_size,
        "raw_visualization_performed": raw_visualization_performed,
        "raw_visualization_style": "legacy_rainbow",
        "raw_visualization_preview_groups": [
            "all_particles",
            "filter_particles_by_sld_gt_1p5",
        ],
        "raw_visualization_sld_preview_cutoff": 1.5,
        "raw_visualization_display_vmax_cap": 100.0,
        "density_backend": "current_dataframe_compat",
        "density_query_batch_size": args.density_query_batch_size,
        "timing_profile_performed": args.profile_time,
        "timing_profile": str(timing_profile_path) if timing_profile_path is not None else None,
        "memory_profile_performed": args.profile_memory,
        "memory_profile": str(memory_profile_path) if memory_profile_path is not None else None,
        "display_filtering": {
            "display_filter_mode": args.display_filter_mode,
            "display_top_fraction": (
                None if args.display_sld_threshold is not None else args.display_top_fraction
            ),
            "display_sld_threshold": args.display_sld_threshold,
            "display_density_field": args.display_density_field,
            "display_max_divisor": args.display_max_divisor,
            "visual_style": args.visual_style,
            "sld_display_mode": "identity",
            "sld_display_outlier_mode": getattr(
                args,
                "sld_display_outlier_mode",
                "tail_jump",
            ),
            "sld_tail_search_fraction": getattr(args, "sld_tail_search_fraction", 0.01),
            "sld_tail_jump_factor": getattr(args, "sld_tail_jump_factor", 5.0),
            "sld_max_display_outlier_fraction": getattr(
                args,
                "sld_max_display_outlier_fraction",
                0.002,
            ),
        },
        "output_artifacts": dict(output_artifacts),
    }
    payload.update(euler_metadata)
    return payload


def _get_field(value, field_name: str, default=None):
    if isinstance(value, dict):
        return value.get(field_name, default)
    return getattr(value, field_name, default)


def _landscape_with_match_rows(landscape, match_table: Any):
    data = landscape.data.copy(deep=True)
    if isinstance(match_table, pd.DataFrame):
        match_frame = match_table
    else:
        match_frame = pd.DataFrame(match_table)
    if {"particle_key", "domain_a_row", "domain_b_row"}.issubset(match_frame.columns):
        indexed = match_frame.set_index("particle_key")
        data["ref_source_row_id"] = (
            data["particle_key"].map(indexed["domain_a_row"]).fillna(-1).astype(int)
        )
        data["mov_source_row_id"] = (
            data["particle_key"].map(indexed["domain_b_row"]).fillna(-1).astype(int)
        )
    return Landscape(
        data=data,
        canonical_transform=landscape.canonical_transform,
        active_policies=landscape.active_policies,
        density_report=landscape.density_report,
        canonicalization_report=landscape.canonicalization_report,
    )


def _match_table_data_for_phase1(phase1, landscape) -> pd.DataFrame:
    match_table = getattr(phase1, "match_table", None)
    if match_table is not None and hasattr(match_table, "data"):
        return match_table.data
    particle_keys = list(landscape.data["particle_key"])
    return pd.DataFrame(
        {
            "particle_key": particle_keys,
            "domain_a_row": range(len(particle_keys)),
            "domain_b_row": range(len(particle_keys)),
            "match_status": "matched",
        }
    )


def _identity_report_for_phase1(phase1, attribute: str, label: str):
    identity = getattr(phase1, attribute, None)
    if identity is not None and hasattr(identity, "report"):
        return identity.report
    return {
        "identity_mode": "unavailable_in_test_runner",
        "identity_columns": (),
        "column_normalization_rules": {},
        "unique_rate": 1.0,
        "duplicate_count": 0,
        "collision_examples": (),
        "status": f"not_reported_for_{label}",
    }


def _write_csv_artifact(data: pd.DataFrame, path: Path, *, overwrite: bool) -> Path:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(path, index=False)
    return path


def _simple_import_report(path: str, source_type: str, row_count: int) -> dict[str, Any]:
    return {
        "artifact_type": "import_report",
        "schema_version": "1",
        "input_path": path,
        "source_type": source_type,
        "row_count": int(row_count),
        "status": "ok",
    }


def _write_run_default_visualization_previews(
    *,
    args,
    raw_landscape: Landscape,
    raw_visualization_dir: Path,
    representation_policy: RepresentationPolicy,
    euler_metadata: dict[str, object],
    raw_landscape_csv_path: Path | None,
) -> dict[str, Any]:
    formats = _parse_formats(args.formats)
    full_landscape_table = "data/raw_landscape.csv" if raw_landscape_csv_path is not None else None
    preview_specs = (
        {
            "label": "all_particles",
            "display_sld_threshold": None,
            "display_filter_mode": "none",
            "description": "all matched particles",
        },
        {
            "label": "filter_particles_by_sld_gt_1p5",
            "display_sld_threshold": 1.5,
            "display_filter_mode": "threshold",
            "description": "display-only preview of particles with sld_raw > 1.5",
        },
    )
    preview_reports: dict[str, dict[str, Any]] = {}
    generated_files: dict[str, str] = {}
    for spec in preview_specs:
        preview_dir = raw_visualization_dir / str(spec["label"])
        max_displayed_sld = _preview_max_displayed_sld(
            raw_landscape,
            display_density_field="sld_raw",
            display_sld_threshold=spec["display_sld_threshold"],
        )
        resolved_vmax = (
            min(max_displayed_sld, 100.0) if max_displayed_sld is not None else None
        )
        report = write_landscape_visualizations(
            raw_landscape,
            preview_dir,
            overwrite=args.overwrite,
            coordinate_source="analysis",
            representation="both",
            color_field="sld_raw",
            euler_convention=representation_policy.euler_convention,
            euler_convention_source=str(euler_metadata["euler_convention_source"]),
            display_top_fraction=None,
            display_sld_threshold=spec["display_sld_threshold"],
            display_density_field="sld_raw",
            formats=formats,
            max_points_2d=args.max_points_2d,
            max_points_3d=args.max_points_3d,
            random_seed=args.random_seed,
            write_projection_csvs=args.write_projection_csvs,
            full_landscape_table=full_landscape_table,
            display_table_filename="display_table.csv",
            artifact_layout="run_bundle",
            selection_metadata={
                "preview_group": spec["label"],
                "preview_description": spec["description"],
                "preview_only": True,
                "not_a_selection": True,
                "sld_preview_cutoff": spec["display_sld_threshold"],
                "display_vmax_cap": 100.0,
                "max_displayed_sld": max_displayed_sld,
                "resolved_vmax": resolved_vmax,
            },
            **_run_preview_visualization_style_kwargs(
                args,
                display_filter_mode=str(spec["display_filter_mode"]),
                resolved_vmax=resolved_vmax,
                display_sld_threshold=spec["display_sld_threshold"],
            ),
        )
        preview_reports[str(spec["label"])] = report
        generated_files[f'{spec["label"]}_visualization_report_json'] = report["report_path"]
        for artifact_key, artifact_path in report["generated_files"].items():
            generated_files[f'{spec["label"]}_{artifact_key}'] = artifact_path

    aggregate_report = {
        "artifact_type": "run_default_visualization_report",
        "schema_version": "1",
        "status": "ok",
        "preview_only": True,
        "not_a_selection": True,
        "visual_style": "legacy_rainbow",
        "color_map": "rainbow_r",
        "color_field": "sld_raw",
        "sld_preview_cutoff": 1.5,
        "display_vmax_cap": 100.0,
        "display_vmax_policy": "min(max_displayed_sld, 100)",
        "preview_groups": [str(spec["label"]) for spec in preview_specs],
        "preview_reports": {
            label: report["report_path"] for label, report in preview_reports.items()
        },
        "generated_files": generated_files,
    }
    report_path = write_json_artifact(
        aggregate_report,
        raw_visualization_dir / "visualization_report.json",
        overwrite=args.overwrite,
    )
    aggregate_report["report_path"] = str(report_path)
    return aggregate_report


def _preview_max_displayed_sld(
    landscape: Landscape,
    *,
    display_density_field: str,
    display_sld_threshold: float | None,
) -> float | None:
    values = np.asarray(landscape.data[display_density_field], dtype=float)
    mask = np.isfinite(values)
    if display_sld_threshold is not None:
        mask &= values > display_sld_threshold
    displayed = values[mask]
    if not displayed.size:
        return None
    return float(np.max(displayed))


def _visualization_style_kwargs(args) -> dict[str, Any]:
    return {
        "visual_style": args.visual_style,
        "color_map": args.color_map,
        "color_vmin": args.color_vmin,
        "color_vmax": args.color_vmax,
        "point_size": args.point_size,
        "point_alpha": args.point_alpha,
        "figure_width": args.figure_width,
        "figure_height": args.figure_height,
        "colorbar_position": args.colorbar_position,
        "sort_points_by_color": args.sort_points_by_color,
        "axis_limits": dict(args.axis_limit or ()),
        "display_filter_mode": args.display_filter_mode,
        "display_max_divisor": args.display_max_divisor,
        "generate_histograms": args.generate_histograms,
        "generate_axis_direction_map": args.generate_axis_direction_map,
    }


def _run_preview_visualization_style_kwargs(
    args,
    *,
    display_filter_mode: str,
    resolved_vmax: float | None,
    display_sld_threshold: float | None = None,
) -> dict[str, Any]:
    color_vmax = args.color_vmax if args.color_vmax is not None else resolved_vmax
    if color_vmax is None and display_sld_threshold is not None:
        color_vmax = display_sld_threshold
    kwargs = _visualization_style_kwargs(args)
    kwargs.update(
        {
            "visual_style": "legacy",
            "color_map": args.color_map or "rainbow_r",
            "color_vmax": color_vmax,
            "display_filter_mode": display_filter_mode,
        }
    )
    return kwargs


def _prepare_run_output_dir(output_dir: str, *, overwrite: bool) -> Path:
    path = Path(output_dir)
    if path.exists() and not path.is_dir():
        raise ValueError(f"Run output path exists and is not a directory: {path}")
    if path.exists() and not overwrite:
        raise FileExistsError(f"Run output directory already exists: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def _prepare_select_output_dir(output_dir: str, *, overwrite: bool) -> Path:
    path = Path(output_dir)
    if path.exists() and not path.is_dir():
        raise ValueError(f"Selection output path exists and is not a directory: {path}")
    if path.exists() and not overwrite:
        raise FileExistsError(f"Selection output directory already exists: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def _prepare_canonicalize_output_dir(output_dir: str, *, overwrite: bool) -> Path:
    path = Path(output_dir)
    if path.exists() and not path.is_dir():
        raise ValueError(f"Canonicalize output path exists and is not a directory: {path}")
    if path.exists() and not overwrite:
        raise FileExistsError(f"Canonicalize output directory already exists: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def _prepare_selection_export_output_dir(output_dir: str, *, overwrite: bool) -> Path:
    path = Path(output_dir)
    if path.exists() and not path.is_dir():
        raise ValueError(f"Selection export output path exists and is not a directory: {path}")
    if path.exists() and not overwrite:
        raise FileExistsError(f"Selection export output directory already exists: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def _resolve_run_backend(requested: str) -> str:
    if requested == "array_native":
        raise ValueError(
            "--run-backend array_native is not implemented for cryorole run yet; "
            "use --run-backend auto or dataframe_compat"
        )
    if requested in {"auto", "dataframe_compat"}:
        return "dataframe_compat"
    raise ValueError(f"Unsupported run backend: {requested}")


def _source_type_from_path(path: str) -> str:
    suffix = Path(path).suffix.lower()
    if suffix == ".star":
        return "relion"
    if suffix == ".cs":
        return "cryosparc"
    return "unknown"


def _selection_policy_from_args(
    args,
    *,
    euler_metadata: dict[str, object] | None = None,
) -> SelectionPolicy:
    radius, radius_unit = _resolve_radius_args(args)
    euler_metadata = euler_metadata or resolve_euler_convention().metadata()
    scipy_sequence = str(euler_metadata["scipy_euler_sequence"])
    euler_convention = str(euler_metadata["euler_convention"])
    internal_mode = _internal_selection_mode(args)
    coordinate_source = _selection_coordinate_source_from_space(args.space)
    range_representation = _range_representation_from_bounds(args.range_bound or ())
    parent_metadata = {
        "euler_metadata": dict(euler_metadata),
        "space": args.space,
        "canonical_id": args.canonical_id if args.space == "canonical" else None,
        "run_dir": str(args.run_dir),
    }
    return SelectionPolicy(
        selection_mode=internal_mode,
        density_support_field="sld_raw",
        density_artifact_policy="include_all",
        top_fraction=0.40,
        random_fraction=getattr(args, "fraction", None),
        random_seed=getattr(args, "seed", None),
        metadata_domain=getattr(args, "metadata_domain", None),
        metadata_column=getattr(args, "metadata_column", None),
        metadata_values=_metadata_values_from_args(args),
        split_by_metadata=bool(getattr(args, "split_by_value", False)),
        threshold=getattr(args, "sld_min", None),
        sld_min=getattr(args, "sld_min", None),
        sld_max=getattr(args, "sld_max", None),
        center_input=tuple(args.center) if args.center is not None else None,
        center_input_representation=args.center_representation,
        center_input_space="evaluation",
        center_euler_sequence=scipy_sequence,
        center_euler_convention=euler_convention,
        center_scipy_euler_sequence=scipy_sequence,
        center_degrees=True,
        evaluation_space=coordinate_source,
        metric=_internal_selection_metric(args.metric),
        radius=radius,
        radius_unit=radius_unit,
        range_coordinate_source=coordinate_source,
        range_representation=range_representation,
        range_euler_sequence=scipy_sequence,
        range_euler_convention=euler_convention,
        range_scipy_euler_sequence=scipy_sequence,
        range_degrees=True,
        range_bounds=dict(args.range_bound or ()),
        parent_landscape_metadata=parent_metadata,
        selection_id=args.selection_id,
    )


def _metadata_values_from_args(args) -> tuple[str, ...]:
    values: list[str] = []
    single_value = getattr(args, "metadata_value", None)
    if single_value is not None:
        values.extend(
            normalized
            for part in str(single_value).split(",")
            if (normalized := _cli_metadata_value_key(part)) is not None
        )
    return tuple(values)


def _internal_selection_mode(args) -> str:
    mode = args.selection_mode
    if mode == "radius":
        return "radius_around_center"
    if mode == "threshold":
        return "threshold_by_density"
    if mode == "range":
        return "range_by_coordinates"
    if mode == "random":
        return "random"
    if mode == "metadata":
        if getattr(args, "split_by_value", False):
            if getattr(args, "metadata_value", None) is not None:
                raise ValueError("--split-by-value cannot be used with --metadata-value")
            return "metadata_group"
        return "metadata_value"
    raise ValueError(f"Unsupported selection mode: {mode}")


def _public_selection_mode(internal_mode: str) -> str:
    return {
        "radius_around_center": "radius",
        "threshold_by_density": "threshold",
        "range_by_coordinates": "range",
        "random": "random",
        "metadata_value": "metadata",
        "metadata_group": "metadata",
    }.get(internal_mode, internal_mode)


def _selection_coordinate_source_from_space(space: str) -> str:
    if space == "raw":
        return "analysis"
    if space == "canonical":
        return "canonical"
    raise ValueError(f"Unsupported selection space: {space}")


def _internal_selection_metric(metric: str) -> str:
    if metric == "so3":
        return "so3_geodesic"
    if metric == "rotvec":
        return "rotvec_euclidean"
    raise ValueError(f"Unsupported selection metric: {metric}")


def _range_representation_from_bounds(
    range_bounds: Sequence[tuple[str, tuple[float | None, float | None]]],
) -> str:
    axes = {axis for axis, _bounds in range_bounds}
    euler_axes = {"alpha", "beta", "gamma"}
    rotvec_axes = {"x", "y", "z"}
    if not axes:
        return "euler"
    if axes <= euler_axes:
        return "euler"
    if axes <= rotvec_axes:
        return "rotvec"
    unsupported = sorted(axes - euler_axes - rotvec_axes)
    if unsupported:
        raise ValueError(f"Unsupported range axis names: {unsupported}")
    raise ValueError("Do not mix Euler and rotvec range axes in one selection")


def _default_metadata_value_selection_id(policy: SelectionPolicy) -> str:
    value_part = "_".join(_sanitize_selection_id_component(value) for value in policy.metadata_values)
    if not value_part:
        value_part = "value"
    return (
        f"{_sanitize_selection_id_component(str(policy.metadata_domain))}_"
        f"{_sanitize_selection_id_component(str(policy.metadata_column))}_"
        f"{value_part}"
    )


def _resolve_radius_args(args) -> tuple[float | None, str | None]:
    provided = [
        name
        for name, value in (
            ("--radius", args.radius),
            ("--radius-rad", args.radius_rad),
        )
        if value is not None
    ]
    if len(provided) > 1:
        raise ValueError("Use only one radius argument: --radius or --radius-rad")
    if args.radius_rad is not None:
        return args.radius_rad, "radians"
    if args.radius is not None:
        return args.radius, "degrees"
    return None, None


def _add_output_format_argument(parser, *, help_text: str | None = None) -> None:
    parser.add_argument(
        "--formats",
        "--output-formats",
        dest="formats",
        default="png",
        help=help_text or "Comma-separated figure formats. Default: png. Supported: png,pdf,svg.",
    )


def _parse_float_tuple(value: str) -> tuple[float, ...]:
    try:
        parsed = tuple(float(part.strip()) for part in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected comma-separated numeric values") from exc
    if len(parsed) != 3:
        raise argparse.ArgumentTypeError("expected exactly three comma-separated values")
    return parsed


def _parse_float_tolerance(value: str) -> tuple[str, float]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected COL=TOL")
    column, tolerance_raw = value.split("=", 1)
    if not column:
        raise argparse.ArgumentTypeError("column name must be non-empty")
    try:
        tolerance = float(tolerance_raw)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("tolerance must be numeric") from exc
    if tolerance <= 0:
        raise argparse.ArgumentTypeError("tolerance must be > 0")
    return column, tolerance


def _positive_int_argument(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def _parse_formats(value: str) -> tuple[str, ...]:
    formats = tuple(part.strip().lower().lstrip(".") for part in value.split(",") if part.strip())
    if not formats:
        raise ValueError("--formats/--output-formats must include at least one format")
    return formats


def _parse_float_pair_bound(value: str) -> tuple[float, float]:
    parts = value.split(":")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("expected MIN:MAX")
    try:
        lower = float(parts[0])
        upper = float(parts[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError("bounds must be numeric") from exc
    if lower >= upper:
        raise argparse.ArgumentTypeError("MIN must be less than MAX")
    return lower, upper


def _parse_range_bound(
    value: str,
) -> tuple[str, tuple[float | None, float | None]]:
    parts = value.split(":")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("expected AXIS:LOWER:UPPER")
    axis, lower_raw, upper_raw = parts
    if not axis:
        raise argparse.ArgumentTypeError("range axis name must be non-empty")
    try:
        lower = None if lower_raw == "" else float(lower_raw)
        upper = None if upper_raw == "" else float(upper_raw)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("range bounds must be numeric or blank") from exc
    return axis, (lower, upper)


class _MemoryProfiler:
    """Small RSS sampler for optional CLI diagnostics."""

    def __init__(self, *, enabled: bool) -> None:
        self.enabled = enabled
        self.samples: list[dict[str, Any]] = []
        self.backend: str | None = None

    def sample(self, stage: str) -> None:
        if not self.enabled:
            return
        rss_bytes, backend = _current_rss_bytes()
        if self.backend is None:
            self.backend = backend
        start_rss = self.samples[0]["rss_bytes"] if self.samples else rss_bytes
        self.samples.append(
            {
                "stage": stage,
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "rss_bytes": rss_bytes,
                "rss_mib": rss_bytes / (1024 * 1024),
                "delta_from_start_bytes": rss_bytes - start_rss,
                "delta_from_start_mib": (rss_bytes - start_rss) / (1024 * 1024),
            }
        )

    def write(self, path: Path, *, overwrite: bool) -> Path:
        peak_rss = max((sample["rss_bytes"] for sample in self.samples), default=None)
        payload = {
            "artifact_type": "canonical_memory_profile",
            "schema_version": "1",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "metric": "rss",
            "units": "bytes",
            "backend": self.backend,
            "sample_count": len(self.samples),
            "peak_rss_bytes": peak_rss,
            "peak_rss_mib": (
                peak_rss / (1024 * 1024) if peak_rss is not None else None
            ),
            "samples": self.samples,
        }
        return write_json_artifact(payload, path, overwrite=overwrite)


def _current_rss_bytes() -> tuple[int, str]:
    try:
        import psutil  # type: ignore

        return int(psutil.Process().memory_info().rss), "psutil"
    except Exception:
        pass
    if os.name == "nt":
        return _current_rss_bytes_windows(), "windows_psapi"
    if sys.platform.startswith("linux"):
        return _current_rss_bytes_linux_proc(), "linux_proc_statm"
    raise RuntimeError(
        "RSS memory profiling is not available on this platform without psutil"
    )


def _current_rss_bytes_linux_proc() -> int:
    with Path("/proc/self/statm").open("r", encoding="utf-8") as handle:
        fields = handle.read().split()
    if len(fields) < 2:
        raise RuntimeError("Could not read RSS from /proc/self/statm")
    page_size = os.sysconf("SC_PAGE_SIZE")
    return int(fields[1]) * int(page_size)


def _current_rss_bytes_windows() -> int:
    class PROCESS_MEMORY_COUNTERS(ctypes.Structure):
        _fields_ = [
            ("cb", ctypes.c_ulong),
            ("PageFaultCount", ctypes.c_ulong),
            ("PeakWorkingSetSize", ctypes.c_size_t),
            ("WorkingSetSize", ctypes.c_size_t),
            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
            ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
            ("PagefileUsage", ctypes.c_size_t),
            ("PeakPagefileUsage", ctypes.c_size_t),
        ]

    counters = PROCESS_MEMORY_COUNTERS()
    counters.cb = ctypes.sizeof(PROCESS_MEMORY_COUNTERS)
    process = ctypes.windll.kernel32.GetCurrentProcess()
    ok = ctypes.windll.psapi.GetProcessMemoryInfo(
        process,
        ctypes.byref(counters),
        counters.cb,
    )
    if not ok:
        raise RuntimeError("Could not read RSS from Windows Process Status API")
    return int(counters.WorkingSetSize)


if __name__ == "__main__":
    raise SystemExit(main())
