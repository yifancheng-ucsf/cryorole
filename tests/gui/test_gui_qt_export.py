"""Export page, headless, on a synthetic CryoSPARC project (see tests/export/test_cryosparc_group_export.py)."""

from __future__ import annotations

import sys
import threading
from pathlib import Path

import pytest

from cryorole.gui.app import PAGES, MainWindow
from cryorole.gui.glprobe import ProbeResult

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "export"))
from test_cryosparc_group_export import build_project  # noqa: E402

NO_GL = ProbeResult(platform_ok=True, gl_status="failed", platform_name="offscreen", reason="no OpenGL")


@pytest.fixture
def project(tmp_path):
    return build_project(tmp_path)


@pytest.fixture
def window(qtbot, project):
    win = MainWindow(NO_GL)
    qtbot.addWidget(win)
    win.resize(1400, 900)
    win.show()
    with qtbot.waitSignal(win.view_page.landscape_loaded, timeout=10000):
        win.open_run(project["run"])
    win.nav.setCurrentRow(PAGES.index("Export"))
    page = win.export_page
    qtbot.waitUntil(lambda: not page._preview_runner.busy and bool(page.preview), timeout=10000)
    return win


def test_export_page_lists_selections_sources_and_found_groups(window, project):
    page = window.export_page
    assert window.pages.currentWidget() is page
    assert page.selection_id == "dense"
    assert page.output_edit.text() == str(project["run"].resolve() / "exports" / "dense")
    assert "J5_particles.csg (found next to the source)" in page.group_labels["ref"].text()
    assert "Choose the source job's particle .csg" in page.group_labels["mov"].text()
    assert page.command_label.text().startswith("cryorole export --run-dir")


def test_export_writes_the_group_and_shows_import_steps(window, qtbot, project):
    page = window.export_page
    page.domain_combo.setCurrentIndex(page.domain_combo.findData("ref"))
    with qtbot.waitSignal(page.exported, timeout=20000) as signal:
        page.export_button.click()
    report = signal.args[0]
    out = Path(report["output_dir"])
    assert (out / "ref" / "selected_ref.csg").is_file()
    assert (out / "ref" / "selected_ref_J5_passthrough_particles.cs").is_file()
    text = page.result_view.toPlainText()
    assert "Import Result Group" in text and str(project["root"] / "exports" / "groups") in text
    assert "cryorole export --run-dir" in window.commands_view.toPlainText()
    page.export_button.click()  # same folder again: refused, nothing overwritten
    qtbot.waitUntil(lambda: page.findChild(object, "error_output_exists") is not None, timeout=5000)


def test_a_group_of_another_job_is_refused_and_choices_survive_navigation(window, qtbot, project, tmp_path):
    page = window.export_page
    page.set_group("mov", str(project["ref"].parent / "J5_particles.csg"))  # J5's group for the J6 source
    window.nav.setCurrentRow(PAGES.index("View"))
    window.nav.setCurrentRow(PAGES.index("Export"))
    assert page.groups["mov"].endswith("J5_particles.csg")
    assert "--cryosparc-group-mov" in page.command_label.text()
    page.output_edit.setText(str(tmp_path / "wrong_group"))
    page.export_button.click()
    qtbot.waitUntil(lambda: page.findChild(object, "error_input_invalid") is not None, timeout=20000)
    box = page.findChild(object, "error_input_invalid")
    assert "does not name the run's source file" in box.text()


def test_output_folder_and_preview_belong_to_the_open_run(window, qtbot, project, tmp_path):
    page = window.export_page
    page.output_edit.setText(str(tmp_path / "typed_for_run_a"))
    page.output_edit.textEdited.emit(page.output_edit.text())
    (tmp_path / "b").mkdir()
    other = build_project(tmp_path / "b")
    with qtbot.waitSignal(window.view_page.landscape_loaded, timeout=10000):
        window.open_run(other["run"])
    window.nav.setCurrentRow(PAGES.index("Export"))
    assert page.output_edit.text() == str(other["run"].resolve() / "exports" / "dense")
    qtbot.waitUntil(lambda: not page._preview_runner.busy, timeout=10000)
    shown = page.source_labels["ref"].text()
    assert str(other["ref"]) in shown
    qtbot.wait(100)
    assert page.source_labels["ref"].text() == shown


def test_a_group_the_export_cannot_use_is_not_offered(window, qtbot, project):
    page = window.export_page
    (project["ref"].parent / "J5_passthrough_particles.cs").unlink()  # the .csg names a file that is gone
    window.nav.setCurrentRow(PAGES.index("View"))
    window.nav.setCurrentRow(PAGES.index("Export"))
    qtbot.waitUntil(lambda: not page._preview_runner.busy, timeout=10000)
    label = page.group_labels["ref"].text()
    assert "found next to the source" not in label
    assert "would stop before writing" in label


def test_late_preview_cannot_replace_the_new_run(window, qtbot, project, tmp_path, monkeypatch):
    from cryorole import api

    page = window.export_page
    started, release = threading.Event(), threading.Event()
    actual = api.export_preview
    old = project["run"].resolve()

    def slow_preview(run):
        if run == old:
            started.set()
            release.wait(5)
        return actual(run)

    monkeypatch.setattr(api, "export_preview", slow_preview)
    page.on_shown()
    qtbot.waitUntil(started.is_set)
    (tmp_path / "new").mkdir()
    other = build_project(tmp_path / "new")
    with qtbot.waitSignal(window.view_page.landscape_loaded, timeout=10000):
        window.open_run(other["run"])
    release.set()
    qtbot.waitUntil(lambda: not page._preview_runner.busy, timeout=10000)
    assert str(other["ref"]) in page.source_labels["ref"].text()
    assert str(project["ref"]) not in page.source_labels["ref"].text()


def test_preview_failure_is_shown_with_the_actionable_error(window, qtbot, monkeypatch):
    from cryorole import api
    from cryorole.errors import CryoroleError

    def broken(_run):
        raise CryoroleError("input_invalid", "Source information is unavailable", "Reopen a complete run")

    monkeypatch.setattr(api, "export_preview", broken)
    window.export_page.on_shown()
    qtbot.waitUntil(lambda: window.export_page.findChild(object, "error_input_invalid") is not None)
    box = window.export_page.findChild(object, "error_input_invalid")
    assert "Source information is unavailable" in box.text()
    assert "Reopen a complete run" in box.informativeText()
