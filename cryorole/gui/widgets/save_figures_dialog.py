"""Choose a new name for standard display figures of a captured view."""

from __future__ import annotations

from PySide6.QtCore import QRegularExpression
from PySide6.QtGui import QRegularExpressionValidator
from PySide6.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QFormLayout, QLabel, QLineEdit, QVBoxLayout,
)


class SaveFiguresDialog(QDialog):
    def __init__(self, description: str, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Save figures")
        self.setObjectName("save_figures_dialog")
        self.name_edit = QLineEdit("view_01")
        self.name_edit.setObjectName("figure_name")
        self.name_edit.setValidator(QRegularExpressionValidator(QRegularExpression(r"[A-Za-z0-9][A-Za-z0-9._-]*")))
        self.views_combo = QComboBox()
        self.views_combo.addItem("2D projections and 1D distributions", "2d,1d")
        self.views_combo.addItem("2D projections only", "2d")
        form = QFormLayout()
        form.addRow("New figure name", self.name_edit)
        form.addRow("Figures", self.views_combo)
        intro = QLabel(description)
        intro.setWordWrap(True)
        note = QLabel(
            "Saves PNGs and a report using the display filter, colour scale and Axis range settings. "
            "Use Custom ranges to save a zoom. Highlight overlays and mouse pan/zoom are not included. "
            "No selection is created. Existing figure sets are kept."
        )
        note.setWordWrap(True)
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.name_edit.textChanged.connect(self._validate)
        layout = QVBoxLayout(self)
        layout.addWidget(intro)
        layout.addLayout(form)
        layout.addWidget(note)
        layout.addWidget(self.buttons)
        self.resize(510, 230)
        self._validate()

    def _validate(self) -> None:
        self.buttons.button(QDialogButtonBox.StandardButton.Save).setEnabled(self.name_edit.hasAcceptableInput())

    @property
    def visual_id(self) -> str:
        return self.name_edit.text()

    @property
    def views(self) -> str:
        return self.views_combo.currentData()
