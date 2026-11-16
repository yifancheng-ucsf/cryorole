import sys

from PySide6.QtCore import QProcess

from cryorole.gui.process import WorkerProcess


def test_worker_reports_failed_start_and_can_be_reused(qtbot, monkeypatch, tmp_path):
    worker = WorkerProcess()
    failures = []
    worker.failed.connect(failures.append)

    with monkeypatch.context() as patch:
        patch.setattr(sys, "executable", str(tmp_path / "missing-python"))
        with qtbot.waitSignal(worker.stopped, timeout=3000):
            worker.start({"command": "run"})

    assert [error["code"] for error in failures] == ["worker_start_failed"]
    assert not worker.running
    assert worker.process is None

    with qtbot.waitSignal(worker.stopped, timeout=10000):
        worker.start({"command": "unknown"})

    assert [error["code"] for error in failures] == ["worker_start_failed", "internal"]
    assert worker.process is None
    qtbot.waitUntil(lambda: not worker.findChildren(QProcess))
