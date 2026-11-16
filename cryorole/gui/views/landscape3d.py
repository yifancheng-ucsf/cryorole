"""3D landscape view (pyqtgraph.opengl). Only created after the OpenGL probe succeeded."""

from __future__ import annotations

import numpy as np
import pyqtgraph.opengl as gl
from OpenGL import GL
from PySide6.QtCore import QPointF, Qt, Signal
from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget

POINT_SIZE_PX = 2.0
_AXIS_COLORS = ((0.80, 0.15, 0.15, 1.0), (0.15, 0.55, 0.15, 1.0), (0.15, 0.25, 0.80, 1.0))
_AXIS_HTML = ("#cc2626", "#268c26", "#2640cc")
DIM_RGBA_F = (0.55, 0.55, 0.55, 0.10)
CLICK_TOLERANCE_PX = 4  # a press and release closer than this is a click, not a rotation
PICK_RADIUS_PX = 12


class _ClickableGLView(gl.GLViewWidget):
    clicked = Signal(float, float)  # widget pixel coordinates

    def mousePressEvent(self, event) -> None:  # noqa: N802 - Qt override
        self._press = QPointF(event.position())
        super().mousePressEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802 - Qt override
        super().mouseReleaseEvent(event)
        press = getattr(self, "_press", None)
        if press is not None and event.button() == Qt.MouseButton.LeftButton:
            moved = event.position() - press
            if abs(moved.x()) + abs(moved.y()) <= CLICK_TOLERANCE_PX:
                self.clicked.emit(float(event.position().x()), float(event.position().y()))
        self._press = None


class Landscape3DView(QWidget):
    """Display sample as opaque points, scaled into a unit box around the data."""

    point_clicked = Signal(int)  # index into the points last passed to set_points

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.view = _ClickableGLView()
        self.view.setBackgroundColor("w")
        self.view.setCameraPosition(distance=3.2, elevation=22, azimuth=45)
        self.scatter = gl.GLScatterPlotItem(
            pos=np.zeros((1, 3), dtype=np.float32), size=POINT_SIZE_PX, pxMode=True, glOptions="opaque"
        )
        self.view.addItem(self.scatter)
        # While a selection is shown, unselected points are drawn faint and translucent after the
        # opaque selected points, so a selection inside the cloud stays visible.
        self.dim_scatter = gl.GLScatterPlotItem(
            pos=np.zeros((0, 3), dtype=np.float32), size=POINT_SIZE_PX, pxMode=True,
            color=DIM_RGBA_F, glOptions="translucent",
        )
        self.view.addItem(self.dim_scatter)
        self.marker = gl.GLScatterPlotItem(
            pos=np.zeros((0, 3), dtype=np.float32), size=14.0, pxMode=True, color=(0.0, 0.0, 0.0, 1.0)
        )
        self.marker.setGLOptions({GL.GL_DEPTH_TEST: False, GL.GL_BLEND: False})  # always on top
        self.view.addItem(self.marker)
        self.view.clicked.connect(self._clicked)
        self._positions = np.zeros((0, 3), dtype=np.float32)
        self.axis_lines = []
        for axis in range(3):
            start, end = np.zeros(3), np.zeros(3)
            start[axis], end[axis] = -1.0, 1.0
            line = gl.GLLinePlotItem(pos=np.array([start, end]), color=_AXIS_COLORS[axis], width=1.5, antialias=True)
            self.view.addItem(line)
            self.axis_lines.append(line)
        # Axis key as plain text: GL text items do not render the same on every platform.
        self.axis_key = QLabel()
        self.axis_key.setObjectName("axis_key")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.view, 1)
        layout.addWidget(self.axis_key)
        self.center = np.zeros(3)
        self.half_extent = 1.0

    def set_points(
        self,
        coordinates: np.ndarray,
        rgba: np.ndarray,
        labels: tuple[str, str, str],
        box: tuple[tuple[float, float], tuple[float, float], tuple[float, float]] | None = None,
        dimmed: np.ndarray | None = None,
    ) -> None:
        """``coordinates`` (n, 3) in data units; ``rgba`` (n, 4) uint8 in draw order.

        ``box`` (per-axis ranges) re-centres and re-scales the view; None keeps the current one.
        Points outside the box are still drawn. ``dimmed`` (bool per point) draws those points faint.
        """

        if box is not None:
            low = np.array([lo for lo, _ in box], dtype=float)
            high = np.array([hi for _, hi in box], dtype=float)
            self.center = (low + high) / 2.0
            self.half_extent = float(max(np.max(high - low) / 2.0, 1e-9))
            self.axis_key.setText(
                "   ".join(
                    f"<span style='color:{color}'>&#9632;</span> {label} [{lo:.3g}, {hi:.3g}]"
                    for color, label, (lo, hi) in zip(_AXIS_HTML, labels, box)
                )
            )
            self.view.setCameraPosition(distance=3.2, elevation=22, azimuth=45)
        positions = ((coordinates - self.center) / self.half_extent).astype(np.float32)
        self._positions = positions
        colors = rgba.astype(np.float32) / 255.0
        if dimmed is None or not np.any(dimmed):
            self.scatter.setData(pos=positions, color=colors)
            self.dim_scatter.setData(pos=np.zeros((0, 3), dtype=np.float32))
        else:
            keep = ~dimmed
            self.scatter.setData(pos=positions[keep], color=colors[keep])
            self.dim_scatter.setData(pos=positions[dimmed], color=DIM_RGBA_F)

    def set_colors(self, rgba: np.ndarray) -> None:
        self.scatter.setData(color=(rgba.astype(np.float32) / 255.0))

    def set_marker(self, coordinate: np.ndarray | None) -> None:
        """Centre marker at one point in data units, or none."""

        if coordinate is None:
            self.marker.setData(pos=np.zeros((0, 3), dtype=np.float32))
            return
        position = ((np.asarray(coordinate, dtype=float) - self.center) / self.half_extent).astype(np.float32)
        self.marker.setData(pos=position.reshape(1, 3))

    def pick(self, x: float, y: float, radius_px: float = PICK_RADIUS_PX) -> int | None:
        """Index of the drawn point nearest to widget pixel (x, y), if within ``radius_px``."""

        if not len(self._positions):
            return None
        viewport = self.view.getViewport()
        matrix = self.view.projectionMatrix(viewport, viewport) * self.view.viewMatrix()
        mvp = np.array(matrix.data(), dtype=float).reshape(4, 4).T  # QMatrix4x4.data() is column-major
        clip = np.c_[self._positions, np.ones(len(self._positions))] @ mvp.T
        in_front = clip[:, 3] > 1e-9
        if not np.any(in_front):
            return None
        ndc = clip[:, :2] / np.where(in_front, clip[:, 3], 1.0)[:, None]
        sx = (ndc[:, 0] + 1.0) / 2.0 * self.view.width()
        sy = (1.0 - ndc[:, 1]) / 2.0 * self.view.height()
        distance = np.hypot(sx - x, sy - y)
        distance[~in_front] = np.inf
        best = int(np.argmin(distance))
        return best if distance[best] <= radius_px else None

    def _clicked(self, x: float, y: float) -> None:
        index = self.pick(x, y)
        if index is not None:
            self.point_clicked.emit(index)

