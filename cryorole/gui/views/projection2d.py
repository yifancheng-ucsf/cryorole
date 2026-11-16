"""Three 2D projections drawn as density images, re-rasterised on every pan or zoom.

Each image holds one cell per screen pixel with the SLD of the point a scatter
would show there (``core.display_policy.visible_value_raster``), so it costs the
same for 10k or 1M points and needs no OpenGL.

With equal axis units (the standard figure style) each panel's plot area is
sized to the aspect ratio of its axis ranges, so the requested ranges are shown
exactly instead of one axis being stretched past them.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import QHBoxLayout, QWidget

# render(panel_index, x_range, y_range, width_px, height_px) -> (height, width, 4) uint8 RGBA
RenderFunction = Callable[[int, tuple[float, float], tuple[float, float], int, int], np.ndarray]
Range = tuple[float, float]
REDRAW_DELAY_MS = 15


class _AspectBox(QWidget):
    """Holds one PlotWidget; with ``aspect`` set, sizes it so the plot area has that width/height ratio."""

    def __init__(self, plot_widget: pg.PlotWidget, parent=None) -> None:
        super().__init__(parent)
        self.plot_widget = plot_widget
        plot_widget.setParent(self)
        self.aspect: float | None = None

    def set_aspect(self, aspect: float | None) -> None:
        self.aspect = aspect
        self._place()

    def resizeEvent(self, event) -> None:  # noqa: N802 - Qt override
        super().resizeEvent(event)
        self._place()

    def _place(self) -> None:
        width, height = self.width(), self.height()
        if self.aspect is None or width < 10 or height < 10:
            self.plot_widget.setGeometry(0, 0, width, height)
            return
        # Axis labels and ticks take a roughly fixed number of pixels around the plot area.
        box = self.plot_widget.getPlotItem().getViewBox()
        deco_w = max(40, self.plot_widget.width() - int(box.width())) if box.width() > 0 else 70
        deco_h = max(30, self.plot_widget.height() - int(box.height())) if box.height() > 0 else 55
        area_w, area_h = max(1, width - deco_w), max(1, height - deco_h)
        if area_w / area_h > self.aspect:
            area_w = int(area_h * self.aspect)
        else:
            area_h = int(area_w / self.aspect)
        total_w, total_h = area_w + deco_w, area_h + deco_h
        self.plot_widget.setGeometry((width - total_w) // 2, (height - total_h) // 2, total_w, total_h)


class ProjectionPanels(QWidget):
    point_clicked = Signal(int, float, float)  # panel index, x, y in data coordinates

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.plots: list[pg.PlotItem] = []
        self.images: list[pg.ImageItem] = []
        self.markers: list[pg.ScatterPlotItem] = []
        self.boxes: list[_AspectBox] = []
        self._timers: list[QTimer] = []
        self._render: RenderFunction | None = None
        self._ranges: list[tuple[Range, Range]] = []
        self.equal_units = True
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        for index in range(3):
            widget = pg.PlotWidget(background="w")
            plot = widget.getPlotItem()
            plot.showGrid(x=True, y=True, alpha=0.2)
            image = pg.ImageItem(axisOrder="row-major")
            plot.addItem(image)
            marker = pg.ScatterPlotItem(symbol="+", size=22, pen=pg.mkPen("k", width=2.5), brush=None)
            marker.setZValue(10)
            plot.addItem(marker)
            widget.scene().sigMouseClicked.connect(lambda event, i=index: self._clicked(i, event))
            timer = QTimer(self)
            timer.setSingleShot(True)
            timer.setInterval(REDRAW_DELAY_MS)
            timer.timeout.connect(lambda i=index: self._redraw(i))
            plot.getViewBox().sigRangeChanged.connect(lambda *_, t=timer: t.start())
            plot.getViewBox().sigResized.connect(lambda *_, t=timer: t.start())
            box = _AspectBox(widget)
            layout.addWidget(box, 1)
            self.plots.append(plot)
            self.images.append(image)
            self.markers.append(marker)
            self.boxes.append(box)
            self._timers.append(timer)

    def set_source(
        self,
        render: RenderFunction,
        labels: list[tuple[str, str]],
        ranges: list[tuple[Range, Range]],
        *,
        equal_units: bool,
        reset_view: bool,
    ) -> None:
        """Show new data; ``ranges`` (x, y) per panel are what Reset view returns to."""

        self._render = render
        self._ranges = ranges
        self.equal_units = equal_units
        for index, plot in enumerate(self.plots):
            plot.setLabel("bottom", labels[index][0])
            plot.setLabel("left", labels[index][1])
        if reset_view:
            self.reset_view()
        for index in range(3):
            self._redraw(index)

    def reset_view(self) -> None:
        for plot, box, ((x0, x1), (y0, y1)) in zip(self.plots, self.boxes, self._ranges):
            view = plot.getViewBox()
            view.setAspectLocked(False)
            box.set_aspect((x1 - x0) / (y1 - y0) if self.equal_units else None)
            view.setRange(xRange=(x0, x1), yRange=(y0, y1), padding=0)
            if self.equal_units:
                view.setAspectLocked(True)  # zooming keeps equal units

    def redraw_all(self) -> None:
        for index in range(len(self.plots)):
            self._redraw(index)

    def set_markers(self, points: list[tuple[float, float]] | None) -> None:
        """One centre marker per panel (data coordinates), or none."""

        for index, marker in enumerate(self.markers):
            marker.setData(pos=[] if points is None else [points[index]])

    def pixel_size(self, index: int) -> tuple[float, float]:
        """Data units per screen pixel along x and y of one panel."""

        size = self.plots[index].getViewBox().viewPixelSize()
        return float(size[0]), float(size[1])

    def _clicked(self, index: int, event) -> None:
        if event.button() != Qt.MouseButton.LeftButton or event.double():
            return
        box = self.plots[index].getViewBox()
        if not box.sceneBoundingRect().contains(event.scenePos()):
            return
        point = box.mapSceneToView(event.scenePos())
        self.point_clicked.emit(index, float(point.x()), float(point.y()))

    def _redraw(self, index: int) -> None:
        if self._render is None:
            return
        box = self.plots[index].getViewBox()
        (x0, x1), (y0, y1) = box.viewRange()
        ratio = self.devicePixelRatioF()
        width = max(1, int(box.width() * ratio))
        height = max(1, int(box.height() * ratio))
        rgba = self._render(index, (x0, x1), (y0, y1), width, height)
        image = self.images[index]
        image.setImage(rgba, autoLevels=False)
        image.setRect(x0, y0, x1 - x0, y1 - y0)
