"""GUI tests run headless; Qt-dependent modules are skipped without the [gui] extra."""

from __future__ import annotations

import importlib.util
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# Pure tests (probe interpretation, launch checks) run everywhere. Tests that need
# Qt are named test_gui_qt_*.py and are collected only when the extra is installed.
if any(importlib.util.find_spec(name) is None for name in ("PySide6", "pyqtgraph", "OpenGL", "pytestqt")):
    collect_ignore_glob = ["test_gui_qt_*.py"]


if "collect_ignore_glob" not in globals():
    import gc

    import pytest

    @pytest.fixture(autouse=True)
    def _collect_qt_garbage_between_tests(qapp):
        """Run the garbage collector between Qt tests, not inside a later test's Qt call.

        With PySide6 6.11 on Python 3.11, a collection that starts inside a signal
        disconnect (pyqtgraph does this while building a new window) and finalizes the
        previous test's widgets can deadlock; CI then hits the job timeout.
        """
        yield
        qapp.processEvents()
        gc.collect()
        qapp.processEvents()
