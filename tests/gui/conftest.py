"""GUI tests run headless; Qt-dependent modules are skipped without the [gui] extra."""

from __future__ import annotations

import importlib.util
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# Pure tests (probe interpretation, launch checks) run everywhere. Tests that need
# Qt are named test_gui_qt_*.py and are collected only when the extra is installed.
if any(importlib.util.find_spec(name) is None for name in ("PySide6", "pyqtgraph", "OpenGL", "pytestqt")):
    collect_ignore_glob = ["test_gui_qt_*.py"]
