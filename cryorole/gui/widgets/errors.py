"""Show any failure as message + remedy, the same text the CLI prints."""

from __future__ import annotations

from PySide6.QtWidgets import QMessageBox, QWidget

from cryorole.errors import classify_exception
from cryorole.logs import warn_user


def show_error(parent: QWidget, exc: BaseException, *, title: str = "cryoROLE") -> QMessageBox:
    """Open a non-blocking warning box; the message also goes to the log panel."""

    error = classify_exception(exc)
    warn_user(error.message + (f" What to do: {error.remedy}" if error.remedy else ""))
    box = QMessageBox(QMessageBox.Icon.Warning, title, error.message, parent=parent)
    if error.remedy:
        box.setInformativeText(error.remedy)
    box.setObjectName(f"error_{error.code}")
    box.open()
    return box
