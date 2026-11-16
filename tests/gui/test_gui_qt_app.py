"""Main window and the real probe child, headless (needs the [gui] extra and pytest-qt)."""

from __future__ import annotations

import logging
import time

from cryorole.gui.app import PAGES, MainWindow
from cryorole.gui.glprobe import ProbeResult, probe

NO_GL = ProbeResult(platform_ok=True, gl_status="failed", platform_name="xcb",
                    reason="Could not initialize GLX", remedy="Start the GUI with PRIME offload")


def test_window_lists_the_workflow_pages_and_switches(qtbot):
    window = MainWindow(NO_GL)
    qtbot.addWidget(window)
    assert [window.nav.item(i).text() for i in range(window.nav.count())] == list(PAGES)
    window.nav.setCurrentRow(3)
    assert window.pages.currentWidget().objectName() == "page_view"


def test_status_bar_shows_why_3d_is_off(qtbot):
    window = MainWindow(NO_GL)
    qtbot.addWidget(window)
    assert not window.gl_available
    assert window.gl_label.text() == "3D view off: Could not initialize GLX"
    assert window.gl_label.toolTip() == "Start the GUI with PRIME offload"


def test_cryorole_log_messages_reach_the_log_panel_until_close(qtbot):
    window = MainWindow(NO_GL)
    qtbot.addWidget(window)
    logger = logging.getLogger("cryorole")
    logger.setLevel(logging.INFO)
    logger.warning("matched 116,492 particles")
    qtbot.waitUntil(lambda: "matched 116,492 particles" in window.log_view.toPlainText(), timeout=2000)
    window.close()
    assert window._log_handler not in logger.handlers


def test_real_probe_child_is_fast_and_safe_offscreen():
    start = time.monotonic()
    result = probe()
    assert time.monotonic() - start < 15
    assert result.platform_ok
    assert result.platform_name == "offscreen"
    assert result.gl_status == "failed"
    assert "no OpenGL widgets" in result.reason
