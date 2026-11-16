"""Display model for viewing one run landscape in the desktop GUI.

Everything here is display-only. It applies the same policy as
``cryorole visualize`` (display filter, deterministic sample with seed 0,
colour scale ``vmax = min(max sld_display, 100)``, legacy style, 1D bins),
so the screen matches the standard figures. It wraps an
:class:`~cryorole.interactive.session.ExploreSession`, which owns the
landscape arrays, the run identity and the exact SO(3) selection evaluator;
the full landscape is loaded once.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
import shlex
from typing import Any, Mapping

import numpy as np
from scipy.spatial.transform import Rotation

from cryorole.core.display_policy import (
    LEGACY_DISPLAY_STYLE,
    ColorScale,
    resolve_display_indices,
    sort_display_indices,
    visible_value_raster,
)
from cryorole.errors import CryoroleError
from cryorole.export.selection_export import read_selection_json
from cryorole.export import write_json_artifact
from cryorole.interactive.session import ExploreSession
from cryorole.io.writers.publication import ArtifactDirectory
from cryorole.select.evaluator import center_to_rotvec
from cryorole.visualize.service import (
    VisualizationRequest,
    VisualizationResult,
    _visualization_output_dir,
    _deterministic_subset_sample,
    _resolve_visualize_color_scale,
    _resolved_bin_count,
    visualize,
)
from cryorole.workflow.resolve import canonical_frame_ids, selection_ids

REPRESENTATIONS = ("euler", "rotvec")
DEFAULT_MAX_POINTS_3D = 100_000  # v2 plan D5; counts and selections always use every particle
DISPLAY_SEED = 0
_AXES = {
    "euler": (("alpha", "alpha (deg)"), ("beta", "beta (deg)"), ("gamma", "gamma (deg)")),
    "rotvec": (("x", "rotvec x (rad)"), ("y", "rotvec y (rad)"), ("z", "rotvec z (rad)")),
}
PROJECTIONS = ((0, 1), (1, 2), (0, 2))  # same panels as the 2D figures
GUI_INTERACTION_SOURCE = "cryorole_gui"
_SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


@dataclass(frozen=True)
class RadiusDraft:
    """An evaluated radius selection that is not yet a scientific Selection (nothing is written)."""

    token: str
    center_input: tuple[float, float, float]
    representation: str  # of center_input: "euler" (degrees) or "rotvec" (radians)
    center_rv: tuple[float, float, float]
    radius_deg: float
    selected_mask: np.ndarray  # exact, over the full landscape; display filters never apply
    selected_count: int
    total_count: int
    command_template: str  # `cryorole select ...` with SELECTION_ID as the name

    @property
    def selected_fraction(self) -> float:
        return self.selected_count / self.total_count if self.total_count else 0.0

    def command(self, selection_id: str) -> str:
        parts = shlex.split(self.command_template)
        parts[parts.index('--selection-id') + 1] = selection_id
        return shlex.join(parts)


@dataclass(frozen=True)
class DisplaySettings:
    """Display-only controls; none of them changes counts, selections or exports."""

    threshold: float | None = None
    top_fraction: float | None = None
    max_points_3d: int = DEFAULT_MAX_POINTS_3D
    color_vmax: float | None = None  # like `cryorole visualize --vmax`; None keeps min(max SLD, 100)


@dataclass(frozen=True)
class DisplayFrame:
    """Rows to show for one set of :class:`DisplaySettings`."""

    settings: DisplaySettings
    filtered_indices: np.ndarray  # every row passing the display filter (2D images, 1D histograms)
    sample_indices: np.ndarray  # deterministic 3D sample, in draw order (ascending SLD)
    color_scale: ColorScale
    vmin: float  # resolved colour range actually used
    vmax: float
    full_count: int

    @property
    def filtered_count(self) -> int:
        return int(self.filtered_indices.size)

    @property
    def sample_count(self) -> int:
        return int(self.sample_indices.size)

    @property
    def downsampled(self) -> bool:
        return self.sample_count < self.filtered_count


@dataclass(frozen=True)
class SavedSelectionOverlay:
    """A read-only, key-matched overlay. It is never a draft or a new Selection."""

    selection_id: str
    selected_mask: np.ndarray
    selected_count: int
    total_count: int


class LandscapeView:
    """One landscape (raw or one canonical frame) of a completed run bundle."""

    def __init__(
        self,
        run_dir: str | Path,
        *,
        space: str = "raw",
        canonical_id: str | None = None,
        interaction_source: str = GUI_INTERACTION_SOURCE,
    ) -> None:
        if space not in ("raw", "canonical"):
            raise CryoroleError("input_invalid", f"space must be 'raw' or 'canonical', not {space!r}.")
        if space == "canonical" and not canonical_id:
            raise CryoroleError("input_invalid", "A canonical view needs a canonical_id.")
        folder = Path(run_dir).expanduser()
        if not folder.is_dir():
            raise CryoroleError(
                "input_not_found",
                f"Run folder not found: {folder}",
                "Check the path, or choose the folder that `cryorole run` created.",
            )
        if not (folder / "run_manifest.json").is_file():
            raise CryoroleError(
                "input_invalid",
                f"{folder} is not a cryoROLE run bundle (it has no run_manifest.json).",
                "Choose the folder that `cryorole run` created; it contains run_manifest.json and data/.",
            )
        self.session = ExploreSession(
            run_dir,
            space=space,
            canonical_id=canonical_id or "default",
            max_display_points=1,  # the browser sample is not used here
            interaction_source=interaction_source,
        )
        self._draft_token: str | None = None
        self.space = space
        self.canonical_id = canonical_id if space == "canonical" else None
        self.canonical_ids = canonical_frame_ids(self.session.run_dir)
        self.style = dict(LEGACY_DISPLAY_STYLE)
        self._euler: np.ndarray | None = None

    # --- identity --------------------------------------------------------------

    @property
    def run_dir(self) -> Path:
        return self.session.run_dir

    @property
    def run_id(self) -> str:
        return self.session.run_id

    @property
    def n_points(self) -> int:
        return int(self.session.arrays.n_points)

    @property
    def euler_convention(self) -> str:
        return self.session.euler_convention

    @property
    def sld_display(self) -> np.ndarray:
        return self.session.arrays.sld_display

    # --- display policy --------------------------------------------------------

    default_max_points_3d = DEFAULT_MAX_POINTS_3D
    display_seed = DISPLAY_SEED
    projections = PROJECTIONS

    def apply(
        self,
        *,
        threshold: float | None = None,
        top_fraction: float | None = None,
        max_points_3d: int = DEFAULT_MAX_POINTS_3D,
        color_vmax: float | None = None,
    ) -> DisplayFrame:
        """Resolve filtered rows, the 3D sample and the colour range (display-only)."""

        settings = DisplaySettings(
            threshold=threshold, top_fraction=top_fraction, max_points_3d=max_points_3d, color_vmax=color_vmax
        )
        if settings.max_points_3d < 1:
            raise CryoroleError("input_invalid", "The 3D point cap must be at least 1.")
        try:
            filtered = resolve_display_indices(
                self.sld_display, threshold=settings.threshold, top_fraction=settings.top_fraction
            )
        except ValueError as exc:
            raise CryoroleError("input_invalid", str(exc)) from exc
        if not filtered.size:
            raise CryoroleError(
                "display_empty",
                "No particles pass the display filter.",
                "Lower the SLD threshold or raise the top fraction.",
            )
        try:
            scale = _resolve_visualize_color_scale(
                self.sld_display,
                display_outlier_mask=self.session.arrays.sld_display_is_outlier,
                color_vmin=None,
                color_vmax=settings.color_vmax,
                display_threshold=settings.threshold,
            )
        except ValueError as exc:
            raise CryoroleError("input_invalid", str(exc), "Raise the SLD colour max above the threshold.") from exc
        sample = _deterministic_subset_sample(
            filtered, max_points=settings.max_points_3d, random_seed=DISPLAY_SEED
        )
        sample = sort_display_indices(sample, self.sld_display, order=self.style["sort_points_by_color"])
        vmin = scale.vmin if scale.vmin is not None else float(np.min(self.sld_display[filtered]))
        vmax = float(scale.vmax) if scale.vmax is not None else float(np.max(self.sld_display[filtered]))
        if vmax <= vmin:
            vmax = vmin + max(abs(vmin) * 0.01, 0.1)
        return DisplayFrame(
            settings=settings,
            filtered_indices=filtered,
            sample_indices=sample,
            color_scale=scale,
            vmin=float(vmin),
            vmax=float(vmax),
            full_count=self.n_points,
        )

    # --- coordinates -----------------------------------------------------------

    def coordinates(self, representation: str, indices: np.ndarray | None = None) -> np.ndarray:
        """``(n, 3)`` coordinates: Euler angles in degrees (recorded convention) or rotation vectors."""

        if representation == "rotvec":
            values = self.session.coordinates
        elif representation == "euler":
            if self._euler is None:
                self._euler = Rotation.from_rotvec(self.session.coordinates).as_euler(
                    self.session.euler_sequence, degrees=True
                )
            values = self._euler
        else:
            raise CryoroleError("input_invalid", f"representation must be one of {REPRESENTATIONS}.")
        return values if indices is None else values[indices]

    @staticmethod
    def axis_names(representation: str) -> tuple[str, str, str]:
        return tuple(name for name, _ in _AXES[representation])

    @staticmethod
    def axis_labels(representation: str) -> tuple[str, str, str]:
        return tuple(label for _, label in _AXES[representation])

    def axis_limits(self, representation: str) -> tuple[tuple[float, float], ...]:
        limits = self.style["axis_limits"]
        return tuple(tuple(limits[name]) for name in self.axis_names(representation))

    # --- 2D and 1D -------------------------------------------------------------

    def projection_image(
        self,
        frame: DisplayFrame,
        representation: str,
        projection: tuple[int, int],
        *,
        x_range: tuple[float, float],
        y_range: tuple[float, float],
        width: int,
        height: int,
        rows: np.ndarray | None = None,
    ) -> np.ndarray:
        """Visible ``sld_display`` per screen pixel over every filtered row (NaN where empty).

        ``rows`` restricts the image to a subset of rows (e.g. a selection highlight).
        """

        rows = frame.filtered_indices if rows is None else rows
        coords = self.coordinates(representation, rows)
        a, b = projection
        return self.raster(
            coords[:, a], coords[:, b], self.sld_display[rows],
            x_range=x_range, y_range=y_range, width=width, height=height,
        )

    def raster(self, x, y, values, *, x_range, y_range, width: int, height: int) -> np.ndarray:
        """``visible_value_raster`` with this view's draw order (for callers that cache coordinates)."""

        return visible_value_raster(
            x, y, values, x_range=x_range, y_range=y_range, width=width, height=height,
            order=self.style["sort_points_by_color"],
        )

    def histogram(
        self,
        frame: DisplayFrame,
        representation: str,
        axis: int,
        *,
        rows: np.ndarray | None = None,
        edges: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Histogram of one coordinate, in percent of every filtered row (``visualize --view 1d`` bins).

        ``rows`` and ``edges`` draw a subset (e.g. a selection) on the same bins and scale.
        """

        values = self.coordinates(representation, frame.filtered_indices if rows is None else rows)[:, axis]
        return self.histogram_of(values, frame, edges=edges)

    def histogram_of(self, values, frame: DisplayFrame, *, edges: np.ndarray | None = None):
        """Percent-of-filtered-rows histogram of given coordinate values (``visualize`` bins unless ``edges``)."""

        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        bins = edges if edges is not None else _resolved_bin_count(values, "auto")
        counts, edges = np.histogram(values, bins=bins)
        return counts * (100.0 / max(1, frame.filtered_count)), edges

    def sld_histogram(self, frame: DisplayFrame) -> tuple[np.ndarray, np.ndarray]:
        """Percent histogram of log10(sld_display) over every filtered row."""

        values = self.sld_display[frame.filtered_indices]
        values = np.log10(values[np.isfinite(values) & (values > 0)])
        counts, edges = np.histogram(values, bins=_resolved_bin_count(values, "auto"))
        return counts * (100.0 / max(1, values.size)), edges

    # --- radius selection (exact, full landscape) --------------------------------

    def represent(self, center_rv, representation: str) -> tuple[float, float, float]:
        """A rotation vector as Euler degrees (recorded convention) or as itself."""

        rv = np.asarray(center_rv, dtype=float).reshape(1, 3)
        if representation == "rotvec":
            return tuple(float(v) for v in rv[0])
        if representation == "euler":
            euler = Rotation.from_rotvec(rv).as_euler(self.session.euler_sequence, degrees=True)
            return tuple(float(v) for v in euler[0])
        raise CryoroleError("input_invalid", f"representation must be one of {REPRESENTATIONS}.")

    def rotvec_of(self, center, representation: str) -> tuple[float, float, float]:
        """Inverse of :meth:`represent`: a centre given in ``representation`` as a rotation vector
        (the selection evaluator's own conversion, with the run's recorded Euler convention)."""

        if representation not in REPRESENTATIONS:
            raise CryoroleError("input_invalid", f"representation must be one of {REPRESENTATIONS}.")
        rv = center_to_rotvec(list(center), representation=representation,
                              scipy_euler_sequence=self.session.euler_sequence, euler_degrees=True)
        return tuple(float(v) for v in np.asarray(rv, dtype=float).reshape(3))

    def evaluate_radius(self, center, *, representation: str, radius_deg: float) -> RadiusDraft:
        """Exact SO(3)-geodesic radius selection over the full landscape (the CLI's evaluator)."""

        try:
            draft = self.session.evaluate(center=list(center), representation=representation, radius_deg=radius_deg)
            selected_mask = self.session.draft_selected_mask(expected_draft_id=draft["draft_id"])
        except ValueError as exc:
            raise CryoroleError("input_invalid", str(exc), "Use a finite centre and a positive radius.") from exc
        token = str(draft["draft_id"])
        self._draft_token = token
        return RadiusDraft(
            token=token,
            center_input=tuple(float(v) for v in center),
            representation=representation,
            center_rv=tuple(float(v) for v in draft["center_rv"]),
            radius_deg=float(radius_deg),
            selected_mask=selected_mask,
            selected_count=int(draft["selected_count"]),
            total_count=int(draft["full_candidate_count"]),
            command_template=str(draft["command_template"]),
        )

    def confirm(self, selection_id: str, draft: RadiusDraft, *, frame: DisplayFrame | None = None) -> dict:
        """Write ``draft`` as a standard Selection named ``selection_id``; never overwrites.

        ``frame``: the display shown when it was confirmed; its filter and point counts are recorded as
        the interaction provenance (display only; the selection is exact over the full landscape).
        """

        if not _SAFE_ID.fullmatch(selection_id or ""):
            raise CryoroleError(
                "input_invalid",
                f"{selection_id!r} is not a valid selection name.",
                "Use letters, numbers, '.', '_' or '-', starting with a letter or number (e.g. region_01).",
            )
        if draft.token != self._draft_token:
            raise CryoroleError(
                "input_invalid",
                "The selection changed after this count was shown.",
                "Wait for the count to update, then confirm again.",
            )
        if (self.run_dir / "selections" / selection_id).exists():
            raise CryoroleError(
                "output_exists",
                f"A selection named {selection_id} already exists in this run.",
                "Choose another name; the GUI never overwrites a selection.",
            )
        try:
            return self.session.confirm(
                selection_id, expected_draft_id=draft.token,
                display_provenance=None if frame is None else self.display_provenance(frame),
            )
        except ValueError as exc:
            raise CryoroleError(
                "input_invalid", str(exc), "Evaluate again; reopen the run if its landscape changed on disk."
            ) from exc

    def display_provenance(self, frame: DisplayFrame) -> dict:
        """How ``frame`` was shown: its display filter, the 2D/1D rows and the 3D sample."""

        settings = frame.settings
        if settings.threshold is not None:
            display_filter = {"mode": "threshold", "field": "sld_display", "threshold": float(settings.threshold)}
        elif settings.top_fraction is not None:
            display_filter = {"mode": "top_fraction", "field": "sld_display",
                              "top_fraction": float(settings.top_fraction)}
        else:
            display_filter = {"mode": "all", "field": "sld_display"}
        return {
            "displayed_point_count": frame.filtered_count,
            "displayed_3d_sample_count": int(frame.sample_indices.size),
            "sampling_policy": "all_rows_passing_the_display_filter_2d_1d;deterministic_3d_sample",
            "display_filter": display_filter,
        }

    def selection_ids(self) -> list[str]:
        return selection_ids(self.run_dir)

    def saved_selection_overlay(self, selection_id: str) -> SavedSelectionOverlay:
        """Read a standard Selection from this run and map keys into the current space."""

        if not _SAFE_ID.fullmatch(selection_id or ""):
            raise CryoroleError("input_invalid", "Choose a selection from this run.")
        self._verify_current_landscape()
        selection = read_selection_json(self.run_dir / "selections" / selection_id / "selection.json")
        if selection.selection_id != selection_id or selection.parent_run_id != self.run_id:
            raise CryoroleError(
                "input_invalid", "The saved selection does not identify this run.",
                "Open its original run or create a new verified selection; legacy selections without a run ID are not overlaid.",
            )
        if selection.total_count != self.n_points:
            raise CryoroleError("input_invalid", "The saved selection's parent count differs from this landscape.")
        keys = np.asarray([str(key) for key in selection.selected_particle_keys])
        parent_keys = np.asarray(self.session.arrays.particle_key).astype(str)
        mask = np.isin(parent_keys, keys)
        if int(mask.sum()) != selection.selected_count or len(set(keys.tolist())) != len(keys):
            raise CryoroleError("input_invalid", "Saved selection keys are missing or duplicated in this landscape.")
        mask.setflags(write=False)
        return SavedSelectionOverlay(selection_id, mask, int(mask.sum()), self.n_points)

    def _verify_current_landscape(self) -> None:
        self.session._validate_identity(self.run_id, self.session.landscape_sha256, rehash=True)

    def save_figures(
        self, visual_id: str, *, frame: DisplayFrame, representation: str,
        axis_ranges: tuple[tuple[float, float], ...], views: str = "2d,1d",
    ) -> VisualizationResult:
        """Save standard figures of a captured display state; publish only on success.

        Uses the public renderer. Draft/saved-selection highlights are deliberately
        excluded. Existing figure sets are never overwritten.
        """

        return self.save_figures_from_request(self.figure_request(
            visual_id, frame=frame, representation=representation, axis_ranges=axis_ranges, views=views,
        ))

    def figure_request(
        self, visual_id: str, *, frame: DisplayFrame, representation: str,
        axis_ranges: tuple[tuple[float, float], ...], views: str = "2d,1d",
    ) -> dict[str, Any]:
        """A checked, JSON-safe description of a figure set for a captured display state.

        The GUI builds it on its own thread and saves the figures in a child process
        (:meth:`save_figures_from_request`), so no plotting runs beside the Qt event loop.
        """

        if not _SAFE_ID.fullmatch(visual_id or ""):
            raise CryoroleError("input_invalid", "Use a new figure name containing letters, numbers, '.', '_' or '-'.")
        if representation not in REPRESENTATIONS or views not in ("2d", "2d,1d"):
            raise CryoroleError("input_invalid", "Choose 2D projections, optionally with 1D distributions.")
        if frame.full_count != self.n_points or len(axis_ranges) != 3:
            raise CryoroleError("input_invalid", "The display snapshot does not belong to this landscape.")
        settings = frame.settings
        return {
            "run_dir": str(self.run_dir), "space": self.space, "canonical_id": self.canonical_id,
            "run_id": self.run_id, "landscape_sha256": self.session.landscape_sha256, "n_points": self.n_points,
            "visual_id": visual_id, "representation": representation, "views": views,
            "threshold": None if settings.threshold is None else float(settings.threshold),
            "top_fraction": None if settings.top_fraction is None else float(settings.top_fraction),
            "vmin": float(frame.vmin), "vmax": float(frame.vmax),
            "axis_ranges": [[float(lo), float(hi)] for lo, hi in axis_ranges],
        }

    def save_figures_from_request(self, figure: Mapping[str, Any]) -> VisualizationResult:
        """Save the figure set that :meth:`figure_request` described, for this same landscape only."""

        same = (
            figure.get("run_id") == self.run_id and figure.get("landscape_sha256") == self.session.landscape_sha256
            and figure.get("space") == self.space and figure.get("canonical_id") == self.canonical_id
            and figure.get("n_points") == self.n_points
        )
        if not same:
            raise CryoroleError("input_invalid", "The display snapshot does not belong to this landscape.")
        visual_id, representation, views = figure["visual_id"], figure["representation"], figure["views"]
        if not _SAFE_ID.fullmatch(visual_id or ""):
            raise CryoroleError("input_invalid", "Use a new figure name containing letters, numbers, '.', '_' or '-'.")
        if representation not in REPRESENTATIONS or views not in ("2d", "2d,1d"):
            raise CryoroleError("input_invalid", "Choose 2D projections, optionally with 1D distributions.")
        axis_ranges = tuple(tuple(bounds) for bounds in figure["axis_ranges"])
        if len(axis_ranges) != 3:
            raise CryoroleError("input_invalid", "The display snapshot does not belong to this landscape.")
        threshold, top_fraction = figure["threshold"], figure["top_fraction"]
        vmin, vmax = float(figure["vmin"]), float(figure["vmax"])
        self._verify_current_landscape()
        axes = self.axis_names(representation)
        limits = tuple((axis, tuple(bounds)) for axis, bounds in zip(axes, axis_ranges))
        request = VisualizationRequest(
            run_dir=str(self.run_dir), space=self.space, canonical_id=self.canonical_id or "default",
            visual_id=visual_id, view=views, representation=representation, formats=("png",),
            threshold=threshold, top_fraction=top_fraction,
            all_particles=threshold is None and top_fraction is None,
            vmax=vmax, vmin=vmin, axis_limit=limits, max_points=self.n_points,
            resolved_by={"run_dir": {"value": str(self.run_dir), "resolved_by": "explicit"}},
        )
        output = _visualization_output_dir(request)
        command = ["cryorole", "visualize", "--run-dir", str(self.run_dir), "--space", self.space,
                   "--visual-id", visual_id, "--view", views, "--representation", representation,
                   "--formats", "png", "--vmin", str(vmin), "--vmax", str(vmax),
                   "--max-points", str(self.n_points)]
        if self.canonical_id:
            command += ["--canonical-id", self.canonical_id]
        if threshold is not None:
            command += ["--threshold", str(float(threshold))]
        elif top_fraction is not None:
            command += ["--top-fraction", str(float(top_fraction))]
        else:
            command += ["--all"]
        for axis, (lo, hi) in limits:
            command += ["--axis-limit", f"{axis}:{float(lo)}:{float(hi)}"]
        with ArtifactDirectory(output, inputs=(self.session.landscape_path, self.run_dir / "run_manifest.json"),
                               label="Figure") as publication:
            rendered = visualize(request, _output_dir=publication.path)
            self._verify_current_landscape()
            report = publication.publicize(rendered.report)
            report["gui_display_snapshot"] = {
                "source": "cryorole_gui", "run_id": self.run_id,
                "landscape_sha256": self.session.landscape_sha256,
                "highlight_overlays_included": False, "mouse_pan_zoom_included": False,
                "equivalent_cli_command": shlex.join(command),
            }
            write_json_artifact(report, publication.path / "visualization_report.json", overwrite=True)
        return VisualizationResult(Path(output).resolve(), report)

    def saved_selections(self) -> list[tuple[str, int | None]]:
        """(name, selected count) of the run's selections, from their small summaries only.

        Cheap enough to call after every confirm: it reads no selection.json and hashes nothing.
        """

        entries = []
        for name in self.selection_ids():
            summary = self.run_dir / "selections" / name / "selection_summary.json"
            try:
                count = json.loads(summary.read_text(encoding="utf-8")).get("selected_count")
            except (OSError, ValueError):
                count = None
            entries.append((name, count if isinstance(count, int) else None))
        return entries

