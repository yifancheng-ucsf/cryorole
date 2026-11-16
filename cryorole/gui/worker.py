"""Child process that runs `api.run`, `api.canonicalize` or a figure export for the GUI.

Long, memory-heavy operations run here, so a crash or out-of-memory cannot take
the GUI down, and provenance is exactly that of the API/CLI. Figures are saved
here too: matplotlib running in a thread of the GUI process, beside Qt painting,
can crash or hang it (seen with PySide6 on Python 3.10/3.11). Protocol (one
JSON object per line):

* stdin, first line: ``{"command": "run" | "canonicalize" | "save_figures", ...arguments}``;
  later lines: ``cancel`` (sets the run's CancelToken).
* stdout: ``CRYOROLE {"type": ...}`` lines, where type is ``event`` (a
  ProgressEvent), ``log`` (a ``cryorole`` logger record), ``result``,
  ``cancelled`` or ``error`` (``code``, ``message``, ``remedy``).

Started as ``python -m cryorole.gui.worker`` by :mod:`cryorole.gui.process`.
"""

from __future__ import annotations

from dataclasses import asdict
import json
import logging
import sys
import threading

PREFIX = "CRYOROLE "
_lock = threading.Lock()


def emit(message_type: str, /, **payload) -> None:
    with _lock:
        print(PREFIX + json.dumps({"type": message_type, **payload}, default=str), flush=True)


class _LogForwarder(logging.Handler):
    def emit(self, record: logging.LogRecord) -> None:  # noqa: A003 - logging API
        emit("log", level=record.levelname, message=record.getMessage())


def _watch_stdin(token) -> None:
    for line in sys.stdin:
        if line.strip() == "cancel":
            token.cancel()
            return


def main() -> int:
    from cryorole import api
    from cryorole.errors import CancelledError, classify_exception
    from cryorole.logs import LOGGER_NAME

    logger = logging.getLogger(LOGGER_NAME)
    logger.addHandler(_LogForwarder())
    logger.setLevel(logging.INFO)
    logger.propagate = False
    try:
        request = json.loads(sys.stdin.readline())
        command = request.pop("command")
    except (ValueError, KeyError) as exc:
        emit("error", code="internal", message=f"Bad worker request: {exc}", remedy=None)
        return 2
    token = api.CancelToken()
    threading.Thread(target=_watch_stdin, args=(token,), daemon=True).start()

    def progress(event) -> None:
        emit("event", fraction=event.fraction, **asdict(event))

    try:
        if command == "run":
            result = api.run(
                request.pop("ref"), request.pop("mov"), output_dir=request.pop("output_dir"),
                progress=progress, cancel=token, **request,
            )
            emit("result", output_dir=str(result.output_dir), matched_count=result.matched_count)
        elif command == "canonicalize":
            result = api.canonicalize(request.pop("run_dir"), **request)
            emit("result", output_dir=str(result.output_dir))
        elif command == "save_figures":  # request = LandscapeView.figure_request(...)
            view = api.open_landscape(request["run_dir"], space=request["space"], canonical_id=request["canonical_id"])
            result = view.save_figures_from_request(request)
            emit("result", output_dir=str(result.output_dir), report=result.report)
        else:
            emit("error", code="internal", message=f"Unknown worker command {command!r}", remedy=None)
            return 2
    except CancelledError:
        emit("cancelled")
        return 1
    except BaseException as exc:  # noqa: BLE001 - every failure goes back to the GUI as a message
        error = classify_exception(exc)
        emit("error", code=error.code, message=error.message, remedy=error.remedy)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
