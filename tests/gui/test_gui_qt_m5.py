"""User workflows for M5, with real standard artifacts and controlled late jobs."""

from __future__ import annotations

import hashlib
import threading

import numpy as np
from PySide6.QtWidgets import QDialogButtonBox, QMessageBox

from cryorole import api
from cryorole.gui import process as process_module
from cryorole.gui.app import MainWindow
from cryorole.gui.pages.view import ViewPage
from cryorole.gui.widgets.save_figures_dialog import SaveFiguresDialog

from test_gui_qt_view import NO_GL, _two_runs, page as _page_fixture, run_dir as _run_fixture

page = _page_fixture
run_dir = _run_fixture


def test_saved_overlay_is_read_only_and_does_not_enable_confirm(page, qtbot, run_dir):
    api.select(run_dir, selection_id="m5_review", mode="random", fraction=0.2, seed=14)
    source = run_dir / "selections/m5_review/selection.json"
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    page.refresh_overlays()
    with qtbot.waitSignal(page.overlay_loaded, timeout=5000):
        page.overlay_combo.setCurrentIndex(page.overlay_combo.findData("m5_review"))
    assert page._saved_overlay.selected_count == 60
    assert "read only" in page.overlay_status.text()
    assert page.highlight_mask is None
    assert np.array_equal(page._display_highlight(), page._saved_overlay.selected_mask)
    assert page.selection_panel.draft is None
    assert not page.selection_panel.confirm_button.isEnabled()
    page.filter_value.setValue(1e8)
    page.refresh(reset_view=False)
    assert "0 visible" in page.overlay_status.text()
    assert page._saved_overlay.selected_count == 60
    assert hashlib.sha256(source.read_bytes()).hexdigest() == digest


def test_space_switch_keeps_saved_keys_but_a_new_run_clears_them(qtbot, tmp_path):
    a, b = _two_runs(tmp_path)
    api.canonicalize(a, no_visualize=True)
    api.select(a, selection_id="review", mode="random", fraction=0.25, seed=8)
    widget = ViewPage(gl_available=False)
    qtbot.addWidget(widget)
    widget.set_landscape(api.open_landscape(a))
    with qtbot.waitSignal(widget.overlay_loaded, timeout=5000):
        widget.overlay_combo.setCurrentIndex(widget.overlay_combo.findData("review"))
    before = widget._saved_overlay.selected_mask.copy()
    with qtbot.waitSignal(widget.overlay_loaded, timeout=5000):
        widget.open_run(a, space="canonical", canonical_id="default")
    assert np.array_equal(widget._saved_overlay.selected_mask, before)
    with qtbot.waitSignal(widget.landscape_loaded, timeout=5000):
        widget.open_run(b)
    assert widget._saved_overlay is None and widget.overlay_combo.currentData() is None
    assert widget._display_highlight() is None


def test_clearing_overlay_while_loading_drops_late_result(page, qtbot, monkeypatch):
    started, release = threading.Event(), threading.Event()
    old = page.landscape.saved_selection_overlay

    def delayed(name):
        started.set()
        release.wait(5)
        return old(name)

    api.select(page.landscape.run_dir, selection_id="m5_late", mode="random", fraction=0.1, seed=9)
    monkeypatch.setattr(page.landscape, "saved_selection_overlay", delayed)
    page.refresh_overlays()
    page.overlay_combo.setCurrentIndex(page.overlay_combo.findData("m5_late"))
    qtbot.waitUntil(started.is_set)
    page.overlay_combo.setCurrentIndex(0)
    release.set()
    qtbot.waitUntil(lambda: not page._overlay_runner.busy, timeout=5000)
    assert page._saved_overlay is None and page._display_highlight() is None


def test_save_figures_flushes_pending_filter_and_records_a_standard_command(page, qtbot):
    page.filter_value.setValue(0.0)  # still inside the display debounce
    page.set_custom_ranges([(-60, 60), (-35, 35), (-70, 70)])
    with qtbot.waitSignal(page.figures_saved, timeout=30000) as completed:
        page.save_figures("m5_view", views="2d")
        assert page.saving_figures and not page.save_figures_button.isEnabled()
    result = completed.args[0]
    assert (result.output_dir / "visualization_report.json").is_file()
    assert result.report["threshold"] == 0.0
    assert "--visual-id m5_view" in result.report["gui_display_snapshot"]["equivalent_cli_command"]
    page.save_figures("m5_view", views="2d")
    qtbot.waitUntil(lambda: page.findChild(object, "error_output_exists") is not None, timeout=5000)
    assert not page.saving_figures


def test_saving_captures_the_old_run_while_a_new_run_is_opened(qtbot, tmp_path, monkeypatch):
    a, b = _two_runs(tmp_path)
    window = MainWindow(NO_GL)
    qtbot.addWidget(window)
    window.show()
    widget = window.view_page
    widget.set_landscape(api.open_landscape(a))
    started, release = tmp_path / "worker_started", tmp_path / "worker_release"
    actual = process_module.module_command

    def held(module, *args):  # the figure worker child waits for the test before it saves
        python, flag, code, *rest = actual(module, *args)
        hold = (f"import os, time\nopen({str(started)!r}, 'w').close()\n"
                f"while not os.path.exists({str(release)!r}):\n    time.sleep(0.05)\n")
        return [python, flag, hold + code, *rest]

    monkeypatch.setattr(process_module, "module_command", held)
    widget.save_figures("captured", views="2d")
    qtbot.waitUntil(started.exists, timeout=10000)
    with qtbot.waitSignal(widget.landscape_loaded, timeout=5000):
        widget.open_run(b)
    messages = []
    monkeypatch.setattr(QMessageBox, "information", lambda *_args: messages.append(_args[-1]))
    window.close()
    assert window.isVisible() and messages  # saving is allowed to finish, not orphaned
    with qtbot.waitSignal(widget.figures_saved, timeout=30000) as completed:
        release.touch()
    assert completed.args[0].output_dir == a.resolve() / "visualizations/raw/captured"
    assert not (b / "visualizations/raw/captured").exists()
    assert widget.landscape.run_dir == b.resolve()
    assert str(a.resolve()) in window.commands_view.toPlainText()


def test_save_dialog_rejects_unsafe_names_and_explains_its_scope(qtbot):
    dialog = SaveFiguresDialog("raw · Euler")
    qtbot.addWidget(dialog)
    button = dialog.buttons.button(QDialogButtonBox.StandardButton.Save)
    dialog.name_edit.setText("../replace")
    assert not button.isEnabled()
    dialog.name_edit.setText("figure_02")
    assert button.isEnabled()
    dialog.reject()
    assert dialog.result() == 0
