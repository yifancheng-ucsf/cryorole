"""Inputs -> preflight -> Run -> View, and Canonicalize, headless."""

from __future__ import annotations

import sys
import time
from pathlib import Path

import pytest

from cryorole.gui.app import PAGES, MainWindow
from cryorole.gui.glprobe import ProbeResult

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_gui_worker import _inputs  # noqa: E402

NO_GL = ProbeResult(platform_ok=True, gl_status="failed", platform_name="offscreen", reason="no OpenGL")


@pytest.fixture
def window(qtbot):
    win = MainWindow(NO_GL)
    qtbot.addWidget(win)
    win.resize(1400, 900)
    win.show()
    return win


def _set_inputs(qtbot, window, ref, mov):
    setup = window.run_setup
    window.inputs_page.boxes["ref"].set_path(ref)
    window.inputs_page.boxes["mov"].set_path(mov)
    qtbot.waitUntil(lambda: setup.inputs_ready, timeout=10000)
    qtbot.waitUntil(lambda: setup.preflight is not None or setup.preflight_error is not None, timeout=30000)


def test_inputs_preflight_run_and_open_in_view(window, qtbot, tmp_path):
    ref, mov = _inputs(tmp_path)
    window.parameters_page.output_edit.setText(str(tmp_path / "gui_run"))
    window.parameters_page.k_spin.setValue(10)
    window.parameters_page.visualize_check.setChecked(False)
    window.parameters_page._changed()
    _set_inputs(qtbot, window, ref, mov)
    setup = window.run_setup
    assert setup.readiness in ("READY", "READY_WITH_WARNINGS")
    text = window.inputs_page.summary.text()
    assert "Matched <b>300</b> particles by <code>uid</code>" in text
    command = window.parameters_page.command.text()
    assert command.startswith("cryorole run --ref") and "--k-neighbors 10" in command and "--no-visualize" in command
    run_page = window.run_page
    assert run_page.start_button.isEnabled()
    with qtbot.waitSignal(window.view_page.landscape_loaded, timeout=60000):
        run_page.start_button.click()
    assert (tmp_path / "gui_run" / "run_manifest.json").is_file()
    assert window.pages.currentWidget() is window.view_page
    assert window.view_page.landscape.n_points == 300
    assert "cryorole run --ref" in window.commands_view.toPlainText()
    assert run_page.progress.value() == 1000
    assert "Done: 300 particles" in run_page.stage_label.text()
    qtbot.waitUntil(lambda: not run_page.worker.running, timeout=10000)
    window.nav.setCurrentRow(PAGES.index("Run"))
    assert not run_page.start_button.isEnabled()  # the output folder exists now
    assert "already exists" in run_page.summary.text()


def test_mixed_software_is_blocked_before_preflight(window, qtbot, tmp_path):
    ref, _mov = _inputs(tmp_path)
    star = tmp_path / "mov.star"
    star.write_text("data_particles\n\nloop_\n_rlnImageName #1\n_rlnAngleRot #2\n", encoding="utf-8")
    window.inputs_page.boxes["ref"].set_path(ref)
    window.inputs_page.boxes["mov"].set_path(star)
    qtbot.waitUntil(lambda: window.run_setup.chosen["mov"] is not None and window.run_setup.chosen["ref"] is not None,
                    timeout=10000)
    assert window.run_setup.mixed_software and not window.run_setup.inputs_ready
    assert "same software" in window.inputs_page.readiness.text()
    assert not window.run_page.start_button.isEnabled()


def test_canonicalize_from_the_view_page_opens_the_new_frame(window, qtbot, tmp_path):
    from cryorole import api

    ref, mov = _inputs(tmp_path)
    api.run(ref, mov, output_dir=tmp_path / "run", k_neighbors=10, no_visualize=True)
    with qtbot.waitSignal(window.view_page.landscape_loaded, timeout=10000):
        window.open_run(tmp_path / "run")
    page = window.view_page
    with qtbot.waitSignal(page.landscape_loaded, timeout=60000):
        page.canonicalize(canonical_id="fit30", fit_top_fraction=0.3, visualize=False)
    assert page.landscape.space == "canonical" and page.landscape.canonical_id == "fit30"
    assert "Canonical: fit30" in [page.space_combo.itemText(i) for i in range(page.space_combo.count())]
    assert "cryorole canonicalize --run-dir" in window.commands_view.toPlainText()


def test_canonicalize_dialog_refuses_existing_names(qtbot):
    from cryorole.gui.widgets.canonicalize_dialog import CanonicalizeDialog

    dialog = CanonicalizeDialog(existing=["default"])
    qtbot.addWidget(dialog)
    assert dialog.id_edit.text() == "frame2"
    ok = dialog.buttons.button(dialog.buttons.StandardButton.Ok)
    dialog.id_edit.setText("default")
    assert not ok.isEnabled() and "never overwrites" in dialog.message.text()
    dialog.id_edit.setText("fit50")
    assert ok.isEnabled()
    assert dialog.values() == {"canonical_id": "fit50", "fit_top_fraction": 0.4, "visualize": True}


def _ready_run(window, qtbot, tmp_path):
    ref, mov = _inputs(tmp_path)
    window.parameters_page.output_edit.setText(str(tmp_path / "gui_run"))
    window.parameters_page.k_spin.setValue(10)
    window.parameters_page.visualize_check.setChecked(False)
    window.parameters_page._changed()
    _set_inputs(qtbot, window, ref, mov)
    assert window.run_page.start_button.isEnabled()
    return window.run_page


def test_cancel_button_stops_real_worker_without_publishing(window, qtbot, tmp_path):
    page = _ready_run(window, qtbot, tmp_path)
    with qtbot.waitSignal(page.worker.cancelled, timeout=15000):
        page.start_button.click()
        page.cancel_button.click()
    qtbot.waitUntil(lambda: not page.worker.running, timeout=5000)
    assert not (tmp_path / "gui_run").exists()
    assert not list(tmp_path.glob(".gui_run.preparing-*"))
    assert "Cancelled" in page.stage_label.text()


def test_stop_now_and_worker_crash_report_actionable_failure(window, qtbot, tmp_path, monkeypatch):
    from cryorole.gui import process
    from cryorole.gui.pages.run_setup import CANCEL_GRACE_S

    page = _ready_run(window, qtbot, tmp_path)
    monkeypatch.setattr(process, "module_command", lambda _module: [
        sys.executable, "-c", "import sys,time; print('test worker waiting',file=sys.stderr,flush=True); time.sleep(60)",
    ])
    page.start_button.click()
    qtbot.waitUntil(lambda: page.worker.running)
    page.cancel_button.click()
    page.cancel_requested_at = time.monotonic() - CANCEL_GRACE_S - 1
    page._tick()  # exercise the grace-period UI without a 20-second test sleep
    assert not page.stop_button.isHidden()
    with qtbot.waitSignal(page.worker.failed, timeout=5000) as failed:
        page.stop_button.click()
    assert failed.args[0]["code"] == "worker_stopped"
    assert "before finishing" in failed.args[0]["message"] and failed.args[0]["remedy"]
    assert not (tmp_path / "gui_run").exists()
    assert page.findChild(object, "error_worker_stopped") is not None


def test_canonicalize_failure_restores_controls_and_keeps_the_raw_run(window, qtbot, tmp_path):
    from cryorole import api

    ref, mov = _inputs(tmp_path)
    run = tmp_path / "run"
    api.run(ref, mov, output_dir=run, k_neighbors=10, no_visualize=True)
    with qtbot.waitSignal(window.view_page.landscape_loaded, timeout=5000):
        window.open_run(run)
    page = window.view_page
    with qtbot.waitSignal(page._canonicalizer.failed, timeout=10000):
        page.canonicalize(canonical_id="invalid", fit_top_fraction=2.0, visualize=False)
    assert page.landscape.space == "raw" and page.canonicalize_button.isEnabled()
    assert not page.canonicalizing
    assert not (run / "canonical/invalid").exists()
