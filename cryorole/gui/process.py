"""Qt side of the worker child process (:mod:`cryorole.gui.worker`)."""

from __future__ import annotations

import json

from PySide6.QtCore import QObject, QProcess, Signal

from cryorole.gui.childproc import module_command
from cryorole.gui.worker import PREFIX


class WorkerProcess(QObject):
    """Runs one worker command; turns its stdout lines into signals.

    ``cancel()`` asks the worker to stop at the next stage boundary (a cancelled
    run publishes nothing and removes its staging folder). ``kill()`` stops the
    process at once; for a run, its hidden staging folder may then remain.
    """

    event = Signal(dict)  # ProgressEvent fields plus "fraction"
    log = Signal(str, str)  # level, message
    succeeded = Signal(dict)
    failed = Signal(dict)  # code, message, remedy
    cancelled = Signal()
    stopped = Signal()  # the process has exited (after any of the outcome signals)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.process: QProcess | None = None
        self._buffer = b""
        self._outcome_sent = False
        self._stderr = b""

    @property
    def running(self) -> bool:
        return self.process is not None and self.process.state() != QProcess.ProcessState.NotRunning

    def start(self, request: dict) -> None:
        if self.running:
            raise RuntimeError("a worker is already running")
        self._buffer, self._stderr, self._outcome_sent = b"", b"", False
        process = QProcess(self)
        program, *arguments = module_command("cryorole.gui.worker")
        process.setProgram(program)
        process.setArguments(arguments)
        process.readyReadStandardOutput.connect(self._read_stdout)
        process.readyReadStandardError.connect(self._read_stderr)
        process.errorOccurred.connect(self._process_error)
        process.finished.connect(self._finished)
        self.process = process
        process.start()
        process.write((json.dumps(request) + "\n").encode())

    def cancel(self) -> None:
        if self.running:
            self.process.write(b"cancel\n")

    def kill(self) -> None:
        if self.running:
            self.process.kill()

    def _read_stdout(self) -> None:
        self._buffer += bytes(self.process.readAllStandardOutput())
        *lines, self._buffer = self._buffer.split(b"\n")
        for raw in lines:
            line = raw.decode("utf-8", "replace")
            if not line.startswith(PREFIX):
                continue
            try:
                message = json.loads(line[len(PREFIX):])
            except ValueError:
                continue
            kind = message.pop("type", None)
            if kind == "event":
                self.event.emit(message)
            elif kind == "log":
                self.log.emit(message.get("level", "INFO"), message.get("message", ""))
            elif kind == "result":
                self._outcome_sent = True
                self.succeeded.emit(message)
            elif kind == "error":
                self._outcome_sent = True
                self.failed.emit(message)
            elif kind == "cancelled":
                self._outcome_sent = True
                self.cancelled.emit()

    def _read_stderr(self) -> None:
        self._stderr = (self._stderr + bytes(self.process.readAllStandardError()))[-4000:]

    def _process_error(self, error: QProcess.ProcessError) -> None:
        if error != QProcess.ProcessError.FailedToStart:
            return  # other process errors are handled on exit
        process = self.process
        self.process = None
        self._outcome_sent = True
        process.deleteLater()
        try:
            self.failed.emit({
                "code": "worker_start_failed",
                "message": f"Could not start the cryoROLE worker: {process.errorString()}",
                "remedy": "Check that the current Python executable is available and can be run.",
            })
        finally:
            self.stopped.emit()

    def _finished(self, exit_code: int, status) -> None:
        process = self.process
        self._read_stdout()
        try:
            self._report_outcome(exit_code, status)
        finally:
            if self.process is process:
                self.process = None
            process.deleteLater()
            self.stopped.emit()

    def _report_outcome(self, exit_code: int, status) -> None:
        if not self._outcome_sent:  # crashed, killed or out of memory before reporting
            tail = self._stderr.decode("utf-8", "replace").strip().splitlines()[-3:]
            how = "was stopped" if status == QProcess.ExitStatus.CrashExit else f"exited with status {exit_code}"
            self.failed.emit({
                "code": "worker_stopped",
                "message": f"The cryoROLE worker process {how} before finishing." + (
                    " " + " ".join(tail) if tail else ""),
                "remedy": "Check free memory (see the preflight estimate) and try again.",
            })
