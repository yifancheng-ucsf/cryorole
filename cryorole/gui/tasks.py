"""Run short blocking calls (loading a landscape, API calls) off the GUI thread."""

from __future__ import annotations

from typing import Any, Callable

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal


class _Signals(QObject):
    done = Signal(object)
    failed = Signal(object)


class _Task(QRunnable):
    def __init__(self, function: Callable[[], Any]) -> None:
        super().__init__()
        self.function = function
        self.signals = _Signals()

    def run(self) -> None:
        try:
            result = self.function()
        except Exception as exc:  # noqa: BLE001 - every failure is reported to the GUI thread
            outcome = (self.signals.failed, exc)
        else:
            outcome = (self.signals.done, result)
        try:
            outcome[0].emit(outcome[1])
        except RuntimeError:  # the GUI object was deleted (e.g. the window closed) while this ran
            pass


def run_in_background(
    function: Callable[[], Any],
    on_done: Callable[[Any], None],
    on_error: Callable[[BaseException], None],
) -> _Task:
    """Call ``function`` in the thread pool; the callbacks run on the GUI thread.

    The caller keeps the returned task referenced until a callback has run.
    """

    task = _Task(function)
    task.setAutoDelete(False)
    task.signals.done.connect(on_done)
    task.signals.failed.connect(on_error)
    QThreadPool.globalInstance().start(task)
    return task


class LatestWinsRunner:
    """Runs one background call at a time; a newer request replaces any waiting one.

    Results of a call that was overtaken by a newer request are dropped, so the
    GUI only ever shows the answer to the latest request (e.g. a live count while
    a radius slider moves).
    """

    def __init__(self, on_done: Callable[[Any], None], on_error: Callable[[BaseException], None]) -> None:
        self._on_done = on_done
        self._on_error = on_error
        self._pending: Callable[[], Any] | None = None
        self._running: _Task | None = None
        self._discard_running = False

    @property
    def busy(self) -> bool:
        return self._running is not None or self._pending is not None

    def submit(self, function: Callable[[], Any]) -> None:
        self._pending = function
        if self._running is None:
            self._start_next()

    def discard(self) -> None:
        """Drop queued work and suppress any result from the active call."""
        self._pending = None
        self._discard_running = True

    def _start_next(self) -> None:
        function, self._pending = self._pending, None
        self._discard_running = False
        self._running = run_in_background(function, self._done, self._failed)

    def _done(self, result: Any) -> None:
        self._running = None
        if self._pending is not None:
            self._start_next()  # this result is already out of date
        elif not self._discard_running:
            self._on_done(result)

    def _failed(self, exc: BaseException) -> None:
        self._running = None
        if self._pending is not None:
            self._start_next()
        elif not self._discard_running:
            self._on_error(exc)
