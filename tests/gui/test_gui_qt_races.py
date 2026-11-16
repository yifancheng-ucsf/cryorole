"""Overlapping background answers in the GUI: only the answer to the latest request is shown."""

from __future__ import annotations

import time
from pathlib import Path
from types import SimpleNamespace

from cryorole.gui.pages import run_setup as rs


def _inspection(path: Path):
    candidate = SimpleNamespace(path=path, software="cryosparc", reason="r", job=None, group=None, particle_count=1)
    return SimpleNamespace(chosen=candidate, candidates=(candidate,), note=None)


def test_the_input_box_keeps_the_answer_for_its_latest_path(qtbot, monkeypatch, tmp_path):
    def inspect(path):
        time.sleep(0.8 if "slow" in str(path) else 0.05)  # e.g. a large job folder on NFS
        return _inspection(Path(path) / "particles.cs")

    monkeypatch.setattr(rs.api, "inspect_input", inspect)
    monkeypatch.setattr(rs.api, "preflight", lambda ref, mov, **o: SimpleNamespace(report={"readiness": "READY"}))
    setup = rs.RunSetup()
    setup.output_dir = tmp_path / "out"
    page = rs.InputsPage(setup)
    qtbot.addWidget(page)
    box = page.boxes["ref"]
    box.set_path("/data/J_slow")
    qtbot.wait(100)
    box.set_path("/data/J_fast")  # the user corrects the choice while the first is still being inspected
    assert setup.chosen["ref"] is None  # no input until the new path is inspected
    qtbot.waitUntil(lambda: setup.chosen["ref"] is not None, timeout=5000)
    qtbot.wait(1000)  # the slow answer arrives now and must be dropped
    assert box.path_edit.text() == "/data/J_fast"
    assert setup.chosen["ref"] == Path("/data/J_fast/particles.cs")


def test_a_preflight_answer_for_older_inputs_never_enables_start(qtbot, monkeypatch, tmp_path):
    def preflight(ref, mov, **options):
        time.sleep(0.6)
        return SimpleNamespace(report={"readiness": "BLOCKED" if "bad" in str(mov) else "READY",
                                       "for": (str(ref), str(mov))})

    monkeypatch.setattr(rs.api, "preflight", preflight)
    setup = rs.RunSetup()
    setup.output_dir = tmp_path / "not_yet"
    page = rs.RunPage(setup)
    qtbot.addWidget(page)
    setup.set_inspection("ref", _inspection(Path("/x/ref.cs")))
    setup.set_inspection("mov", _inspection(Path("/x/good_mov.cs")))  # check A starts after the 400 ms debounce
    qtbot.wait(700)  # A is running
    setup.set_inspection("mov", _inspection(Path("/x/bad_mov.cs")))  # new input, whose check would be BLOCKED
    qtbot.wait(400)  # A finishes during the new debounce
    assert setup.preflight is None and page.can_start()[0] is False and not page.start_button.isEnabled()
    qtbot.waitUntil(lambda: setup.preflight is not None, timeout=5000)
    assert setup.preflight["for"] == ("/x/ref.cs", "/x/bad_mov.cs") and setup.readiness == "BLOCKED"
    assert page.can_start() == (False, "Preflight is BLOCKED; see the Inputs page.")


def test_a_stopped_run_names_its_own_staging_folder(qtbot, tmp_path):
    setup = rs.RunSetup()
    page = rs.RunPage(setup)
    qtbot.addWidget(page)
    page._run_output_dir = tmp_path / "run_a"  # as recorded when the run started
    (tmp_path / ".run_a.preparing-1234").mkdir()
    setup.output_dir = tmp_path / "run_b"  # edited on the Parameters page during the run
    page._failed({"code": "worker_stopped", "message": "The cryoROLE worker process was stopped."})
    box = page.findChild(object, "error_worker_stopped")
    assert str(tmp_path / ".run_a.preparing-1234") in box.text()
