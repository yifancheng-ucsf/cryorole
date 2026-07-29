"""Three-panel landscape projection animation renderer."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cryorole.animation.schemas import TrajectoryFrame
from cryorole.core.display_policy import (
    ColorScale,
    EULER_AXIS_NAMES,
    EULER_PROJECTIONS,
    LEGACY_DISPLAY_STYLE,
    coordinate_range_mask,
    resolve_color_scale,
    sort_display_indices,
)


PANEL_SPECS = tuple(
    (
        name,
        x_index,
        y_index,
        f"{EULER_AXIS_NAMES[x_index]} (deg)",
        f"{EULER_AXIS_NAMES[y_index]} (deg)",
    )
    for name, x_index, y_index in EULER_PROJECTIONS
)

PROJECTION_LAYOUT_VERSION = "3"
PROJECTION_WSPACE = 0.035
MARKER_DIAMETER_PIXELS = 13.0
CURRENT_VALUE_FONT_PIXELS = 17.0
ANNOTATION_POLICY = "coordinate_only_once_when_enabled"
SHARED_COLORBAR_POLICY = "single_horizontal_centered"


def resolve_axis_limits(
    landscape_euler: np.ndarray,
    trajectory_euler: np.ndarray,
    *,
    explicit: Mapping[str, Sequence[float]] | None,
    range_bounds: Mapping[str, tuple[float | None, float | None]] | None = None,
) -> tuple[dict[str, tuple[float, float, float, float]], tuple[str, ...]]:
    """Resolve the legacy Euler viewport and reject clipped trajectory points."""

    _validate_euler_array(landscape_euler, "landscape_euler")
    trajectory = _validate_euler_array(trajectory_euler, "trajectory_euler")
    if explicit is not None and range_bounds:
        raise ValueError("Use only one explicit viewport policy: --axis-limits or --range")
    limits: dict[str, tuple[float, float, float, float]] = {}
    axis_viewport = {
        axis: tuple(LEGACY_DISPLAY_STYLE["axis_limits"][axis])
        for axis in EULER_AXIS_NAMES
    }
    for axis, bounds in dict(range_bounds or {}).items():
        if axis not in EULER_AXIS_NAMES:
            raise ValueError(f"Unsupported animation Euler range axis: {axis}")
        lower, upper = bounds
        if lower is not None and upper is not None and lower < upper:
            axis_viewport[axis] = (float(lower), float(upper))
    if range_bounds and not coordinate_range_mask(
        trajectory,
        axis_names=EULER_AXIS_NAMES,
        range_bounds=range_bounds,
        wraparound=True,
    ).all():
        raise ValueError("Trajectory lies outside explicit display range")
    if explicit is not None and set(explicit) != {spec[0] for spec in PANEL_SPECS}:
        raise ValueError("Explicit axis limits must define alpha_beta, beta_gamma, alpha_gamma")
    for panel, x_index, y_index, _, _ in PANEL_SPECS:
        trajectory_x = trajectory[:, x_index]
        trajectory_y = trajectory[:, y_index]
        if explicit is not None:
            values = np.asarray(explicit[panel], dtype=float)
            if values.shape != (4,) or not np.isfinite(values).all():
                raise ValueError(f"Explicit axis limits for {panel} must be four finite values")
            xmin, xmax, ymin, ymax = map(float, values)
            if xmin >= xmax or ymin >= ymax:
                raise ValueError(f"Explicit axis limits for {panel} require lower < upper")
            if (
                np.any(trajectory_x < xmin)
                or np.any(trajectory_x > xmax)
                or np.any(trajectory_y < ymin)
                or np.any(trajectory_y > ymax)
            ):
                raise ValueError(f"Trajectory lies outside explicit axis limits for {panel}")
            limits[panel] = (xmin, xmax, ymin, ymax)
            continue
        xmin, xmax = axis_viewport[EULER_AXIS_NAMES[x_index]]
        ymin, ymax = axis_viewport[EULER_AXIS_NAMES[y_index]]
        if (
            np.any(trajectory_x < xmin)
            or np.any(trajectory_x > xmax)
            or np.any(trajectory_y < ymin)
            or np.any(trajectory_y > ymax)
        ):
            raise ValueError(f"Trajectory lies outside explicit display range for {panel}")
        limits[panel] = (xmin, xmax, ymin, ymax)
    return limits, ()


class ProjectionRenderer:
    """Persist static landscape artists and update only dynamic overlays."""

    def __init__(
        self,
        *,
        landscape_euler: np.ndarray,
        sld_display: np.ndarray,
        sld_display_is_outlier: np.ndarray | None = None,
        trajectory_frames: Sequence[TrajectoryFrame],
        scipy_euler_sequence: str,
        axis_limits: Mapping[str, Sequence[float]],
        canvas_width: int,
        canvas_height: int,
        point_size: float,
        landscape_alpha: float,
        colormap: str = "rainbow_r",
        color_vmin: float | None = None,
        color_vmax: float | None = None,
        display_threshold: float | None = None,
        resolved_color_scale: ColorScale | None = None,
        trail_frames: int,
        show_current_values: bool,
        composite_scale: float = 1.0,
    ) -> None:
        self.landscape = _validate_euler_array(landscape_euler, "landscape_euler")
        self.sld = np.asarray(sld_display, dtype=float)
        self.sld_field = "sld_display"
        if self.sld.shape != (len(self.landscape),) or not np.isfinite(self.sld).all():
            raise ValueError("sld_display must be finite and match landscape rows")
        self.sld_display_is_outlier = (
            np.zeros(self.sld.shape, dtype=bool)
            if sld_display_is_outlier is None
            else np.asarray(sld_display_is_outlier, dtype=bool)
        )
        if self.sld_display_is_outlier.shape != self.sld.shape:
            raise ValueError("sld_display_is_outlier must match landscape rows")
        self.frames = tuple(trajectory_frames)
        if not self.frames:
            raise ValueError("trajectory must contain at least one frame")
        self.trajectory_euler = np.vstack(
            [
                frame.rotation.as_euler(scipy_euler_sequence, degrees=True)
                for frame in self.frames
            ]
        )
        if canvas_width <= 0 or canvas_height <= 0:
            raise ValueError("canvas dimensions must be positive")
        if not np.isfinite(point_size) or point_size <= 0:
            raise ValueError("point_size must be finite and positive")
        if not np.isfinite(landscape_alpha) or not 0 < landscape_alpha <= 1:
            raise ValueError("landscape_alpha must be in (0, 1]")
        if trail_frames < 0:
            raise ValueError("trail_frames must be non-negative")
        if not np.isfinite(composite_scale) or composite_scale <= 0:
            raise ValueError("composite_scale must be finite and positive")
        self.trail_frames = trail_frames
        self.show_current_values = show_current_values
        self.panel_specs = PANEL_SPECS
        self.static_source_row_count = len(self.landscape)
        self.annotation_policy = ANNOTATION_POLICY
        self.shared_colorbar_policy = SHARED_COLORBAR_POLICY
        self.colorbar_label = "SLD"
        self.projection_layout_version = PROJECTION_LAYOUT_VERSION
        dpi = 100
        self.resolved_marker_diameter_pixels = MARKER_DIAMETER_PIXELS
        self.resolved_current_value_font_size_pixels = CURRENT_VALUE_FONT_PIXELS
        marker_diameter_points = (
            MARKER_DIAMETER_PIXELS / composite_scale * 72.0 / dpi
        )
        self.marker_size_points2 = marker_diameter_points**2
        self.current_value_font_size_points = (
            CURRENT_VALUE_FONT_PIXELS / composite_scale * 72.0 / dpi
        )
        self.figure, axes = plt.subplots(
            1,
            3,
            figsize=(canvas_width / dpi, canvas_height / dpi),
            dpi=dpi,
        )
        left_margin = min(0.09, max(0.045, 42.0 / canvas_width))
        self.figure.subplots_adjust(
            left=left_margin,
            right=0.99,
            bottom=0.25,
            top=0.90,
            wspace=PROJECTION_WSPACE,
        )
        self.axes = tuple(np.atleast_1d(axes))
        self.markers = []
        self.trails = []
        self.static_scatter_count = 0
        color_scale = resolved_color_scale or resolve_color_scale(
            self.sld,
            display_outlier_mask=self.sld_display_is_outlier,
            visual_style="legacy",
            color_vmin=color_vmin,
            color_vmax=color_vmax,
            display_threshold=display_threshold,
        )
        self.resolved_vmin = color_scale.vmin
        self.resolved_vmax = color_scale.vmax
        self.color_vmin_source = color_scale.vmin_source
        self.color_vmax_source = color_scale.vmax_source
        draw_indices = sort_display_indices(
            np.arange(len(self.sld), dtype=int),
            self.sld,
            order="ascending",
        )
        static_scatter = None
        for axis, (panel, x_index, y_index, xlabel, ylabel) in zip(
            self.axes, PANEL_SPECS
        ):
            static_scatter = axis.scatter(
                self.landscape[draw_indices, x_index],
                self.landscape[draw_indices, y_index],
                c=self.sld[draw_indices],
                cmap=colormap,
                vmin=self.resolved_vmin,
                vmax=self.resolved_vmax,
                s=point_size,
                alpha=landscape_alpha,
                linewidths=0,
                rasterized=True,
            )
            self.static_scatter_count += 1
            xmin, xmax, ymin, ymax = map(float, axis_limits[panel])
            axis.set_xlim(xmin, xmax)
            axis.set_ylim(ymin, ymax)
            axis.set_xlabel(xlabel)
            axis.set_ylabel(ylabel)
            axis.set_aspect("equal", adjustable="box")
            trail, = axis.plot([], [], color="white", linewidth=1.5, alpha=0.65)
            marker = axis.scatter(
                [],
                [],
                s=self.marker_size_points2,
                c="#111111",
                edgecolors="white",
                linewidths=1.5,
                zorder=5,
            )
            self.trails.append(trail)
            self.markers.append(marker)
        first_panel = self.axes[0].get_position(original=True)
        self.annotation = self.figure.text(
            first_panel.x0,
            0.965,
            "",
            va="top",
            ha="left",
            fontsize=self.current_value_font_size_points,
            color="black",
            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
        )
        self.annotations = (self.annotation,)
        colorbar_axis = self.figure.add_axes((0.28, 0.075, 0.44, 0.035))
        self.colorbar = self.figure.colorbar(
            static_scatter,
            cax=colorbar_axis,
            label=self.colorbar_label,
            orientation="horizontal",
        )
        self.figure.canvas.draw()
        self.panel_rectangles = tuple(
            _axes_layout_rectangle(
                axis,
                canvas_width=canvas_width,
                canvas_height=canvas_height,
            )
            for axis in self.axes
        )
        panel_gaps = tuple(
            right[0] - (left[0] + left[2])
            for left, right in zip(
                self.panel_rectangles,
                self.panel_rectangles[1:],
            )
        )
        self.projection_gap_pixels = max(panel_gaps, default=0)
        self.projection_gap_fraction = self.projection_gap_pixels / canvas_width
        if self.projection_gap_fraction > 0.0125:
            raise RuntimeError(
                "Resolved projection panel gap exceeds 1.25% of canvas width"
            )
        self._static_background = self.figure.canvas.copy_from_bbox(self.figure.bbox)
        self.static_background_draw_count = 1

    def render(self, output_dir: str | Path) -> tuple[Path, ...]:
        """Write one fixed-size PNG for each trajectory row."""

        directory = Path(output_dir)
        directory.mkdir(parents=True, exist_ok=True)
        outputs: list[Path] = []
        for index in range(len(self.frames)):
            self.figure.canvas.restore_region(self._static_background)
            start = max(0, index - self.trail_frames)
            trail_values = self.trajectory_euler[start : index + 1]
            current = self.trajectory_euler[index]
            for panel_index, (_, x_index, y_index, _, _) in enumerate(PANEL_SPECS):
                self.trails[panel_index].set_data(
                    trail_values[:, x_index],
                    trail_values[:, y_index],
                )
                self.markers[panel_index].set_offsets(
                    np.asarray([[current[x_index], current[y_index]]])
                )
                self.axes[panel_index].draw_artist(self.trails[panel_index])
                self.axes[panel_index].draw_artist(self.markers[panel_index])
            self.annotation.set_text(
                (
                    f"\u03b1={current[0]:.2f} "
                    f"\u03b2={current[1]:.2f} "
                    f"\u03b3={current[2]:.2f}"
                )
                if self.show_current_values
                else ""
            )
            self.figure.draw_artist(self.annotation)
            output = directory / f"frame_{index:06d}.png"
            self.figure.canvas.blit(self.figure.bbox)
            plt.imsave(output, np.asarray(self.figure.canvas.buffer_rgba()))
            outputs.append(output)
        plt.close(self.figure)
        return tuple(outputs)


def _validate_euler_array(values: np.ndarray, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n, 3)")
    if len(array) == 0 or not np.isfinite(array).all():
        raise ValueError(f"{name} must be non-empty and finite")
    return array


def _axes_layout_rectangle(
    axis,
    *,
    canvas_width: int,
    canvas_height: int,
) -> tuple[int, int, int, int]:
    """Return a subplot cell as a top-left-origin pixel rectangle."""

    bounds = axis.get_position(original=True)
    x = int(round(bounds.x0 * canvas_width))
    y = int(round((1.0 - bounds.y1) * canvas_height))
    width = int(round(bounds.width * canvas_width))
    height = int(round(bounds.height * canvas_height))
    return x, y, width, height
