"""Horizontal SLD colour legend."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QColor, QLinearGradient, QPainter
from PySide6.QtWidgets import QSizePolicy, QWidget


class ColorLegend(QWidget):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._lut: np.ndarray | None = None
        self._vmin = 0.0
        self._vmax = 1.0
        self.setMinimumSize(220, 34)
        self.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)

    def set_scale(self, lut: np.ndarray, vmin: float, vmax: float) -> None:
        self._lut, self._vmin, self._vmax = lut, vmin, vmax
        self.setToolTip(f"SLD colour range {vmin:.3g} to {vmax:.3g} (display only; vmax capped at 100)")
        self.update()

    def label_text(self) -> str:
        return f"SLD {self._vmin:.3g} - {self._vmax:.3g}"

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt override
        if self._lut is None:
            return
        painter = QPainter(self)
        bar = QRectF(4, 2, self.width() - 8, 14)
        gradient = QLinearGradient(bar.topLeft(), bar.topRight())
        steps = len(self._lut) - 1
        for index in range(0, len(self._lut), 8):
            r, g, b, a = (int(value) for value in self._lut[index])
            gradient.setColorAt(index / steps, QColor(r, g, b, a))
        r, g, b, a = (int(value) for value in self._lut[-1])
        gradient.setColorAt(1.0, QColor(r, g, b, a))
        painter.fillRect(bar, gradient)
        painter.drawText(QRectF(4, 17, self.width() - 8, 16), Qt.AlignmentFlag.AlignLeft, f"{self._vmin:.3g}")
        painter.drawText(QRectF(4, 17, self.width() - 8, 16), Qt.AlignmentFlag.AlignHCenter, "SLD")
        painter.drawText(QRectF(4, 17, self.width() - 8, 16), Qt.AlignmentFlag.AlignRight, f"{self._vmax:.3g}")
        painter.end()
