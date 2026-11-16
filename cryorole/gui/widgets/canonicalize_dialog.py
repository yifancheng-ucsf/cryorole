"""Options for `cryorole canonicalize`, with its defaults."""

from __future__ import annotations

import re

from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QLineEdit,
)

DEFAULT_FIT_TOP_PERCENT = 40.0  # cryorole canonicalize --fit-top 0.4
_SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class CanonicalizeDialog(QDialog):
    def __init__(self, parent=None, *, existing: list[str] = ()) -> None:
        super().__init__(parent)
        self.setWindowTitle("Canonicalize")
        self.existing = set(existing)
        name = "default" if "default" not in self.existing else self._free_name()
        self.id_edit = QLineEdit(name)
        self.fit_spin = QDoubleSpinBox()
        self.fit_spin.setRange(1.0, 100.0)
        self.fit_spin.setDecimals(1)
        self.fit_spin.setSuffix(" %")
        self.fit_spin.setValue(DEFAULT_FIT_TOP_PERCENT)
        self.fit_spin.setToolTip("Fit the frame to this top fraction of particles by sld_raw (--fit-top).")
        self.visualize_check = QCheckBox("Write preview figures (PNG)")
        self.visualize_check.setChecked(True)
        self.message = QLabel()
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        form = QFormLayout(self)
        form.addRow(QLabel("Fits a motion-aligned frame; raw artifacts are not changed."))
        form.addRow("Frame name", self.id_edit)
        form.addRow("Fit to top", self.fit_spin)
        form.addRow(self.visualize_check)
        form.addRow(self.message)
        form.addRow(self.buttons)
        self.id_edit.textChanged.connect(lambda _: self._check())
        self._check()

    def _free_name(self) -> str:
        index = 2
        while f"frame{index}" in self.existing:
            index += 1
        return f"frame{index}"

    def _check(self) -> None:
        name = self.id_edit.text()
        problem = None
        if not _SAFE_ID.fullmatch(name):
            problem = "Use letters, numbers, '.', '_' or '-'."
        elif name in self.existing:
            problem = f"A frame named {name} exists; choose another name (the GUI never overwrites)."
        self.message.setText(problem or "")
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(problem is None)

    def values(self) -> dict:
        return {"canonical_id": self.id_edit.text(), "fit_top_fraction": self.fit_spin.value() / 100.0,
                "visualize": self.visualize_check.isChecked()}
