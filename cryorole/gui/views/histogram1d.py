"""1D marginal distributions and the SLD distribution, over every filtered row."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg

_BRUSH = (192, 192, 192, 170)
_OVERLAY_PEN = (20, 60, 200)


class HistogramPanels(pg.GraphicsLayoutWidget):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setBackground("w")
        self.coordinate_plots = [self.addPlot(row=0, col=index) for index in range(3)]
        self.sld_plot = self.addPlot(row=1, col=0, colspan=3)
        for plot in (*self.coordinate_plots, self.sld_plot):
            plot.showGrid(x=True, y=True, alpha=0.25)
            plot.setLabel("left", "Particles (%)")
        self.sld_plot.setLabel("bottom", "log10(SLD)")
        self._overlays: list = []

    def set_data(
        self,
        coordinate_histograms: list[tuple[np.ndarray, np.ndarray]],
        labels: tuple[str, str, str],
        sld_histogram: tuple[np.ndarray, np.ndarray],
        x_ranges: list[tuple[float, float]],
    ) -> None:
        """``x_ranges`` are the coordinate axis ranges shared with the 2D and 3D views."""

        self._overlays = []  # _draw clears each plot, overlays included
        for plot, (counts, edges), label, x_range in zip(self.coordinate_plots, coordinate_histograms, labels, x_ranges):
            _draw(plot, counts, edges)
            plot.setLabel("bottom", label)
            plot.setXRange(*x_range, padding=0)
        _draw(self.sld_plot, *sld_histogram)
        self.sld_plot.enableAutoRange()


    def set_overlay(self, histograms: list[tuple[np.ndarray, np.ndarray]] | None) -> None:
        """Outline of a selection's distribution on the same bins and scale, or none."""

        for plot, item in self._overlays:
            plot.removeItem(item)
        self._overlays = []
        if histograms is None:
            return
        for plot, (counts, edges) in zip(self.coordinate_plots, histograms):
            item = plot.plot(edges, counts, stepMode="center", pen=pg.mkPen(_OVERLAY_PEN, width=2))
            self._overlays.append((plot, item))


def _draw(plot: pg.PlotItem, counts: np.ndarray, edges: np.ndarray) -> None:
    plot.clear()
    plot.plot(edges, counts, stepMode="center", fillLevel=0, brush=_BRUSH, pen=pg.mkPen("k", width=1))
    plot.enableAutoRange(axis="y")
