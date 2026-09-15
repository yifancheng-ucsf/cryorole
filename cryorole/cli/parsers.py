"""Argparse-only command registration for the public CLI."""

from __future__ import annotations

import argparse
from typing import Callable, Mapping

from cryorole.core.density import SLD_METRICS
from cryorole.core.euler_conventions import EULER_CONVENTIONS

_HIDDEN_HELP = argparse.SUPPRESS

class _CommandParser(argparse.ArgumentParser):
    def error(self, message):
        if self.prog.endswith(" select") and "required" in message and "--selection-id" in message:
            message += (
                "\nChoose a name for this selection: add --selection-id region_01."
                "\nEach name identifies a saved subset under RUN/selections/NAME/."
            )
        super().error(message)


class _StoreExplicit(argparse.Action):
    """Retain user intent for options whose defaults apply to only one mode."""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)
        namespace.explicit_options = (*getattr(namespace, "explicit_options", ()), self.dest)


def build_command_parser(handlers: Mapping[str, Callable[..., int]]) -> argparse.ArgumentParser:
    """Build the cryoROLE CLI parser without embedding scientific logic."""

    parser = _CommandParser(prog="cryorole")
    subparsers = parser.add_subparsers(dest="command", required=True)

    _add_align_parser(subparsers, handlers)
    _add_preflight_parser(subparsers, handlers)
    _add_run_parser(subparsers, handlers)
    _add_status_parser(subparsers, handlers)
    _add_next_parser(subparsers, handlers)
    _add_guide_parser(subparsers, handlers)
    _add_explore_parser(subparsers, handlers)
    _add_visualize_parser(subparsers, handlers)
    _add_animate_parser(subparsers, handlers)
    _add_canonical_views_parser(subparsers, handlers)
    _add_canonicalize_parser(subparsers, handlers)
    _add_select_parser(subparsers, handlers)
    _add_export_parser(subparsers, handlers)
    _add_manifest_parser(subparsers, handlers)
    return parser


def _add_align_parser(subparsers, handlers) -> None:
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
    parser.set_defaults(handler=handlers["align"])


def _add_run_parser(subparsers, handlers) -> None:
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
        "--allow-low-overlap",
        action="store_true",
        help="Explicitly allow key-match coverage below the public 50%% safety threshold.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run the shared preflight only; do not create a run bundle.",
    )
    parser.add_argument(
        "--json",
        nargs="?",
        const="-",
        dest="preflight_json",
        metavar="PATH",
        help="With --dry-run, write the structured preflight report (default: stdout).",
    )
    parser.add_argument(
        "--sld-metric",
        choices=SLD_METRICS,
        default="rotvec_euclidean",
        help=(
            "SLD kNN distance metric. Default: rotvec_euclidean; "
            "so3_geodesic is explicit opt-in."
        ),
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
        default="low_density_skew",
        help=_HIDDEN_HELP,
    )
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
    parser.set_defaults(handler=handlers["run"])


def _add_preflight_parser(subparsers, handlers) -> None:
    parser = subparsers.add_parser(
        "preflight",
        help="Check inputs, matching, policies, and resources without creating a run.",
    )
    parser.add_argument("--ref", required=True, help="Reference-domain pose metadata.")
    parser.add_argument("--mov", required=True, help="Moving-domain pose metadata.")
    parser.add_argument("--ref-domain", default="ref", help="Reference domain label.")
    parser.add_argument("--mov-domain", default="mov", help="Moving domain label.")
    parser.add_argument("--output-dir", default="cryorole_outputs", help="Prospective run directory.")
    parser.add_argument("--row-aligned", action="store_true", help="Assert row N matches row N.")
    parser.add_argument(
        "--allow-low-overlap",
        action="store_true",
        help="Explicitly allow overlap below the public 50%% threshold.",
    )
    parser.add_argument(
        "--sld-metric", choices=SLD_METRICS, default="rotvec_euclidean",
        help="Prospective SLD metric. Default: rotvec_euclidean.",
    )
    parser.add_argument("--no-visualize", action="store_true", help="Estimate a run without previews.")
    parser.add_argument(
        "--json", nargs="?", const="-", dest="preflight_json", metavar="PATH",
        help="Write JSON to PATH, or stdout when PATH is omitted.",
    )
    parser.add_argument("--k-neighbors", type=int, default=50, help=_HIDDEN_HELP)
    parser.add_argument(
        "--density-query-batch-size", type=int, default=25000, help=_HIDDEN_HELP
    )
    parser.add_argument(
        "--run-backend", choices=("auto", "array_native", "dataframe_compat"),
        default="auto", help=_HIDDEN_HELP,
    )
    parser.add_argument("--identity-mode", default=None, help=_HIDDEN_HELP)
    parser.add_argument("--identity-column", action="append", default=(), help=_HIDDEN_HELP)
    parser.add_argument("--mapping-file", help=_HIDDEN_HELP)
    parser.set_defaults(handler=handlers["preflight"])


def _add_status_parser(subparsers, handlers) -> None:
    parser = subparsers.add_parser("status", help="Inspect actual run-bundle artifacts and integrity.")
    parser.add_argument("--run-dir", required=True, help="Run bundle to inspect.")
    parser.add_argument("--json", nargs="?", const="-", dest="json_output", metavar="PATH")
    parser.set_defaults(handler=handlers["status"])


def _add_next_parser(subparsers, handlers) -> None:
    parser = subparsers.add_parser("next", help="Recommend conservative next commands from run artifacts.")
    parser.add_argument("--run-dir", required=True, help="Run bundle to inspect.")
    parser.add_argument("--json", nargs="?", const="-", dest="json_output", metavar="PATH")
    parser.set_defaults(handler=handlers["next"])


def _add_guide_parser(subparsers, handlers) -> None:
    parser = subparsers.add_parser("guide", help="Plan or resume the cryoROLE workflow safely.")
    parser.add_argument("--ref", help="Reference metadata for a new workflow.")
    parser.add_argument("--mov", help="Moving metadata for a new workflow.")
    parser.add_argument("--run-dir", help="Existing run bundle to resume.")
    parser.add_argument("--output-dir", default="cryorole_outputs", help="Prospective output for a new run.")
    parser.add_argument("--row-aligned", action="store_true", help="Explicit row-aligned assertion.")
    parser.add_argument("--allow-low-overlap", action="store_true", help="Explicit low-overlap override.")
    parser.add_argument("--non-interactive", action="store_true", help="Print a plan and never wait for input.")
    parser.add_argument("--execute-run", action="store_true", help="Explicitly execute a ready run; never creates selections.")
    parser.add_argument("--json", nargs="?", const="-", dest="json_output", metavar="PATH")
    parser.set_defaults(handler=handlers["guide"])


def _add_explore_parser(subparsers, handlers) -> None:
    parser = subparsers.add_parser("explore", help="Explore draft radius selections in an offline local browser.")
    parser.add_argument("--run-dir", required=True, help="Completed run bundle.")
    parser.add_argument("--space", choices=("raw", "canonical"), default="raw")
    parser.add_argument("--canonical-id", default="default")
    parser.add_argument("--selection-id", help="Optional existing selection overlay.")
    filters = parser.add_mutually_exclusive_group()
    filters.add_argument("--threshold", type=float, help="Initial display-only SLD threshold.")
    filters.add_argument("--top-fraction", type=float, help="Initial display-only top SLD fraction.")
    parser.add_argument("--max-display-points", type=_positive_int_argument, default=50000)
    parser.add_argument("--colormap", default="viridis", help="Accessible default display colormap; legacy rainbow_r remains available.")
    parser.add_argument("--port", type=int, default=0, help="Loopback port; 0 selects an available port.")
    parser.add_argument("--no-open", action="store_true", help="Do not open the default browser automatically.")
    parser.set_defaults(handler=handlers["explore"])


def _add_visualize_parser(subparsers, handlers) -> None:
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
        "--view",
        help="Comma-separated views to generate: 2d,1d,3d. Default: 2d.",
    )
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
    filters = parser.add_mutually_exclusive_group()
    filters.add_argument("--top-fraction", type=float, help="Display the top fraction by SLD.")
    filters.add_argument(
        "--sld-threshold",
        dest="threshold",
        type=float,
        help="Display rows with sld_display >= T. Default: 1.",
    )
    filters.add_argument(
        "--all",
        dest="all_particles",
        action="store_true",
        help="Display all candidate particles instead of the default SLD threshold.",
    )
    parser.add_argument("--threshold", dest="threshold", type=float, help=_HIDDEN_HELP)
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
        "--format",
        dest="formats",
        action="append",
        help="Static figure format; repeat or use commas. Default: png.",
    )
    parser.add_argument("--formats", dest="formats", action="append", help=_HIDDEN_HELP)
    parser.add_argument("--vmin", type=float, help="Display-only color minimum.")
    parser.add_argument("--vmax", type=float, help="Display-only color maximum.")
    parser.add_argument("--point-size", type=float, help="Display-only scatter point size.")
    opacity = parser.add_mutually_exclusive_group()
    opacity.add_argument("--opacity", dest="alpha", type=float,
                         help="Display-only point opacity: 0 transparent, 1 opaque. Uses the style default when omitted.")
    opacity.add_argument("--alpha", dest="alpha", type=float, help=_HIDDEN_HELP)
    parser.add_argument(
        "--axis-limit",
        dest="axis_limit",
        action="append",
        default=None,
        type=_parse_range_bound,
        metavar="AXIS:LOWER:UPPER",
        help="Viewport limit only; may repeat and does not filter rows.",
    )
    parser.add_argument(
        "--max-points",
        type=_positive_int_argument,
        help="Override the deterministic point cap for 2D/3D views.",
    )
    parser.add_argument(
        "--bins",
        default="auto",
        help="1D histogram bins: auto or a positive integer. Default: auto.",
    )
    parser.add_argument(
        "--hist-mode",
        choices=("count", "percent"),
        default="percent",
        help="1D histogram y-axis mode. Default: percent.",
    )
    parser.add_argument("--kde", action="store_true", help="Add coordinate KDE curves to 1D views.")
    parser.add_argument(
        "--kde-bandwidth",
        help="KDE bandwidth when --kde is used: scott, silverman, or a positive float.",
    )
    parser.add_argument(
        "--3d-mode",
        dest="three_d_mode",
        choices=("interactive", "static"),
        default="interactive",
        help="3D output mode. Default: interactive offline HTML.",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.set_defaults(handler=handlers["visualize"])


def _add_animate_parser(subparsers, handlers) -> None:
    parser = subparsers.add_parser(
        "animate",
        help="Export Phase 1-4 trajectory, rendered frames, composites, and optional MP4.",
    )
    parser.add_argument("--run-dir", required=True, help="Existing cryoROLE run bundle.")
    parser.add_argument(
        "--coordinate-set",
        choices=("raw", "canonical"),
        default="raw",
        help="Waypoint and landscape coordinate set. Default: raw.",
    )
    parser.add_argument("--canonical-id", default="default")
    parser.add_argument("--path-csv", required=True, help="EA or RV waypoint CSV.")
    parser.add_argument("--path-space", choices=("ea", "rv"), required=True)
    parser.add_argument("--chimerax-session", required=True, help="Preconfigured .cxs session.")
    parser.add_argument(
        "--secondary-chimerax-session",
        help="Optional second .cxs with the same model transforms and a different view.",
    )
    parser.add_argument(
        "--tertiary-chimerax-session",
        help=(
            "Optional third .cxs with the same model transforms and a different "
            "view; requires --secondary-chimerax-session."
        ),
    )
    parser.add_argument(
        "--reference-model-id",
        action="append",
        required=True,
        help="Exact stationary ChimeraX model ID; repeat for a rigid reference group.",
    )
    parser.add_argument(
        "--moving-model-id",
        action="append",
        required=True,
        help="Exact ChimeraX model ID; repeat for a rigid moving group.",
    )
    parser.add_argument(
        "--pivot",
        nargs=3,
        type=float,
        required=True,
        metavar=("X", "Y", "Z"),
        help="Rotation pivot in ChimeraX scene coordinates.",
    )
    parser.add_argument(
        "--baseline-ro",
        required=True,
        help="Saved moving-model RO assertion: identity, first-waypoint, ea:A,B,G, or rv:X,Y,Z.",
    )
    parser.add_argument(
        "--map-frame",
        choices=("raw", "canonical", "explicit"),
        required=True,
    )
    parser.add_argument("--map-frame-transform")
    parser.add_argument(
        "--euler-convention",
        choices=("auto",) + EULER_CONVENTIONS,
        default="auto",
    )
    parser.add_argument("--frames-per-segment", type=_positive_int_argument, default=30)
    parser.add_argument("--fps", type=float, default=30.0)
    parser.add_argument("--hold-frames", type=int, default=0)
    parser.add_argument("--reverse", action="store_true")
    parser.add_argument("--ping-pong", action="store_true")
    display_filter = parser.add_mutually_exclusive_group()
    display_filter.add_argument(
        "--threshold",
        "--sld-threshold",
        dest="threshold",
        type=float,
        help="Display rows with sld_display >= threshold; --sld-threshold is an alias.",
    )
    display_filter.add_argument("--top-fraction", type=float)
    parser.add_argument("--colormap", default="rainbow_r")
    parser.add_argument("--vmin", type=float)
    parser.add_argument("--vmax", type=float)
    parser.add_argument(
        "--range",
        dest="range_bound",
        action="append",
        default=None,
        type=_parse_range_bound,
        metavar="AXIS:LOWER:UPPER",
    )
    parser.add_argument("--point-size", type=float, default=1.0)
    parser.add_argument("--landscape-alpha", type=float, default=1.0)
    parser.add_argument("--trail-frames", type=int, default=0)
    parser.add_argument("--show-current-values", action="store_true")
    parser.add_argument("--axis-limits", help="JSON file defining all three panel limits.")
    parser.add_argument(
        "--render-mode",
        choices=("script-only", "execute"),
        default="script-only",
    )
    parser.add_argument("--chimerax-bin")
    parser.add_argument("--canvas-width", type=_positive_int_argument, default=1800)
    parser.add_argument("--canvas-height", type=_positive_int_argument, default=600)
    parser.add_argument("--structure-width", type=_positive_int_argument, default=900)
    parser.add_argument("--structure-height", type=_positive_int_argument, default=900)
    structure_horizontal_crop = parser.add_mutually_exclusive_group()
    structure_horizontal_crop.add_argument(
        "--structure-horizontal-crop",
        dest="structure_horizontal_crop",
        type=float,
        metavar="FRACTION",
        help="Crop this fraction from both horizontal sides of each multi-view frame.",
    )
    structure_horizontal_crop.add_argument(
        "--dual-structure-horizontal-crop",
        dest="structure_horizontal_crop",
        type=float,
        metavar="FRACTION",
        help=(
            "Compatibility alias for --structure-horizontal-crop."
        ),
    )
    parser.set_defaults(structure_horizontal_crop=0.0)
    parser.add_argument(
        "--structure-vertical-crop",
        type=float,
        default=0.0,
        metavar="FRACTION",
        help="Crop this fraction from both vertical sides of each multi-view frame.",
    )
    parser.add_argument("--composite-width", type=_positive_int_argument, default=1920)
    parser.add_argument("--composite-height", type=_positive_int_argument, default=1080)
    parser.add_argument(
        "--layout",
        choices=("stacked", "side_by_side"),
        default="stacked",
        help="Composite layout. Default: stacked.",
    )
    parser.add_argument("--background-color", default="#ffffff")
    parser.add_argument(
        "--no-encode",
        action="store_true",
        help="Write validated composite PNGs without encoding MP4.",
    )
    parser.add_argument("--ffmpeg-bin")
    parser.add_argument("--ffprobe-bin")
    parser.add_argument("--crf", type=int, default=18)
    parser.add_argument("--movie-name", default="animation.mp4")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.set_defaults(handler=handlers["animate"])


def _add_canonical_views_parser(subparsers, handlers) -> None:
    parser = subparsers.add_parser(
        "canonical-views",
        help="Export camera-only ChimeraX views along the canonical axes.",
    )
    parser.add_argument("--run-dir", required=True, help="Existing cryoROLE run bundle.")
    parser.add_argument("--canonical-id", default="default")
    parser.add_argument("--chimerax-session", required=True, help="Preconfigured .cxs session.")
    parser.add_argument(
        "--map-frame",
        choices=("raw", "explicit"),
        required=True,
        help="Scene basis of the input session.",
    )
    parser.add_argument("--map-frame-transform")
    parser.add_argument(
        "--render-mode",
        choices=("script-only", "execute"),
        default="script-only",
    )
    parser.add_argument("--chimerax-bin")
    parser.add_argument("--width", type=_positive_int_argument, default=900)
    parser.add_argument("--height", type=_positive_int_argument, default=900)
    parser.add_argument(
        "--save-sessions",
        action="store_true",
        help="Also save one camera-only .cxs session for each canonical axis.",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.set_defaults(handler=handlers["canonical_views"])


def _add_canonicalize_parser(subparsers, handlers) -> None:
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
    parser.set_defaults(handler=handlers["canonicalize"])


def _add_select_parser(subparsers, handlers) -> None:
    parser = subparsers.add_parser(
        "select", help="Create a named scientific Selection from a saved landscape.",
        description=("Save a scientific subset for inspection and export. Choose a selection name and a mode. "
                     "Selection evaluates the full parent landscape, independently of display sampling or filters."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples (replace RUN with your run directory):
  cryorole select --run-dir RUN --selection-id region_01 --mode radius --center 0 0 0 --radius 15
  cryorole select --run-dir RUN --selection-id high_sld --mode threshold --sld-min 2
  cryorole select --run-dir RUN --selection-id alpha_window --mode range --range-bound alpha:-20:20
  cryorole select --run-dir RUN --selection-id sample_10pct --mode random --fraction 0.1 --seed 7
  cryorole select --run-dir RUN --selection-id cs_classes --mode metadata --metadata-domain ref --metadata-column alignments3D/class --metadata-value 0,1
  cryorole select --run-dir RUN --selection-id ref_classes --mode metadata --metadata-domain ref --metadata-column rlnClassNumber --metadata-value 1,3
  cryorole select --run-dir RUN --selection-id ref_classes --mode metadata --metadata-domain ref --metadata-column rlnClassNumber --split-by-value""",
    )
    common = parser.add_argument_group("Common inputs and naming")
    common.add_argument("--run-dir", required=True, help="Existing cryoROLE run bundle directory.")
    common.add_argument("--selection-id", required=True, metavar="NAME",
                        help="Required user-chosen name, e.g. region_01. Saved under RUN/selections/NAME/. No default name.")
    common.add_argument("--space", choices=("raw", "canonical"), default="raw",
                        help="Parent coordinate space. Default: raw.")
    common.add_argument("--canonical-id", default="default", action=_StoreExplicit,
                        help="Canonical landscape ID; only with --space canonical. Default: default.")
    common.add_argument("--mode", dest="selection_mode",
                        choices=("radius", "threshold", "range", "random", "metadata"), default="radius",
                        help="radius: nearby orientations; threshold: SLD bounds; range: coordinate bounds; "
                             "random: random subset; metadata: source values. Default: radius.")

    radius = parser.add_argument_group("radius mode", "Select near a center using SO(3) geodesic distance by default. "
                                       "Requires --center and exactly one radius option. Euler conventions come from the parent landscape.")
    radius.add_argument("--center", "-c", nargs=3, type=float, metavar=("A", "B", "C"),
                        help="Center in the chosen --space: Euler degrees or rotvec radians.")
    radius.add_argument("--center-representation", choices=("euler", "rotvec"), default="euler", action=_StoreExplicit,
                        help="Representation of --center. Default: euler.")
    radius.add_argument("--radius", "-r", type=float, help="Radius in degrees; exclusive with --radius-rad.")
    radius.add_argument("--radius-rad", type=float, help="Radius in radians; independent of center units.")
    radius.add_argument("--metric", choices=("so3", "rotvec"), default="so3", action=_StoreExplicit,
                        help="so3: rotation geodesic; rotvec: Euclidean rotation-vector distance. Default: so3. Neither uses Euler Euclidean distance.")

    threshold = parser.add_argument_group("threshold mode", "Select by sld_raw, including both boundaries. "
                                          "Requires at least one bound; display colors and filters do not control selection.")
    threshold.add_argument("--sld-min", type=float, help="Inclusive lower SLD bound (>= 0).")
    threshold.add_argument("--sld-max", type=float, help="Inclusive upper SLD bound (>= lower bound).")

    ranges = parser.add_argument_group("range mode", "Select the intersection of coordinate-axis bounds. "
                                       "Use alpha/beta/gamma in degrees OR x/y/z in radians; do not mix representations. "
                                       "Euler lower > upper wraps across the periodic seam; rotvec bounds must be ordered.")
    ranges.add_argument("--range-bound", action="append", default=None, type=_parse_range_bound,
                        metavar="AXIS:LOWER:UPPER", help="Required; repeat for multiple axes. Bounds are inclusive. "
                        "Blank sides are open, e.g. alpha::20. At least one side must be constrained. Repeat of the same axis uses the last bound.")

    random = parser.add_argument_group("random mode", "Sample without replacement from all parent rows; "
                                       "ceil(fraction * row count) rows are selected. This creates a Selection, unlike visualization sampling.")
    random.add_argument("--fraction", type=float, help="Required fraction, 0 < F <= 1.")
    random.add_argument("--seed", type=int, help="Non-negative random seed for reproducibility. Omitted: fresh randomness; recorded seed is null.")

    metadata = parser.add_argument_group("metadata mode", "Use run-time CryoSPARC CS or RELION STAR metadata and recorded source-row provenance. "
                                         "CS supports scalar integers, booleans, and UTF-8 text; empty strings are excluded. "
                                         "Requires domain, column, and either values or split. No rematching or external annotation file.")
    metadata.add_argument("--metadata-domain", choices=("ref", "mov"), help="Required source domain; never inferred.")
    metadata.add_argument("--metadata-column", help="Required source particle column, e.g. alignments3D/class (CS) or rlnClassNumber (STAR).")
    metadata.add_argument("--metadata-value", help="Values to include (union), e.g. 0,1. CS uses typed exact matching; booleans accept true/false/1/0. Exclusive with --split-by-value.")
    metadata.add_argument("--split-by-value", action="store_true",
                          help="One standard Selection per matched value: NAME_VALUE, with filename-safe components. CS limit: 100 groups; collisions fail before writing. Exclusive with --metadata-value.")

    output = parser.add_argument_group("Output controls")
    output.add_argument("--write-selected-landscape", action="store_true",
                        help="Also write selected_landscape/ for visualize --selection-id NAME --use-selected-landscape.")
    output.add_argument("--recompute-sld", action="store_true",
                        help="Requires --write-selected-landscape. Recompute subset SLD and preserve parent SLD fields separately.")
    output.add_argument("--overwrite", action="store_true",
                        help="Explicitly replace this selection ID only (generated child IDs in split mode); parent landscapes and unrelated selections remain unchanged.")
    parser.set_defaults(handler=handlers["select"])


def _add_export_parser(subparsers, handlers) -> None:
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
    _add_export_source_verification_arguments(parser)
    parser.set_defaults(handler=handlers["export_metadata"])

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
    _add_export_source_verification_arguments(selection_parser)
    selection_parser.set_defaults(handler=handlers["export"])


def _add_export_source_verification_arguments(parser) -> None:
    parser.add_argument(
        "--relocated-ref",
        help="Explicit relocated ref source; content must match the recorded SHA-256.",
    )
    parser.add_argument(
        "--relocated-mov",
        help="Explicit relocated mov source; content must match the recorded SHA-256.",
    )
    parser.add_argument(
        "--allow-unverified-source",
        action="store_true",
        help="Advanced legacy override for bundles that predate source hashes; recorded in export report.",
    )


def _add_manifest_parser(subparsers, handlers) -> None:
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
    parser.set_defaults(handler=handlers["manifest"])


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


def _add_output_format_argument(parser, *, help_text: str | None = None) -> None:
    parser.add_argument(
        "--formats",
        "--output-formats",
        dest="formats",
        default="png",
        help=help_text or "Comma-separated figure formats. Default: png. Supported: png,pdf,svg.",
    )


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
