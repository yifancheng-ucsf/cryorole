"""cryoROLE logging: services log, the CLI decides how messages are shown.

Services never ``print``. User-relevant warnings go to the ``cryorole`` logger
(``warn_user``); the CLI installs a stderr handler that renders them as
``[cryorole] warning: …``; API callers can attach their own handler.
"""

from __future__ import annotations

import logging
import sys

LOGGER_NAME = "cryorole"
logger = logging.getLogger(LOGGER_NAME)
logger.addHandler(logging.NullHandler())


def warn_user(message: str) -> None:
    """Report a user-relevant warning without printing from service code."""

    logger.warning(message)


def notify_user(message: str) -> None:
    """Report user-relevant progress or information (shown by the CLI unless --quiet)."""

    logger.info(message)


class _CliFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        if record.levelno >= logging.WARNING:
            return f"[cryorole] warning: {record.getMessage()}"
        return record.getMessage()


class _CurrentStderrHandler(logging.StreamHandler):
    """Write to whatever ``sys.stderr`` is at emit time (robust to redirection)."""

    def __init__(self) -> None:
        super().__init__(sys.stderr)

    def emit(self, record: logging.LogRecord) -> None:
        self.stream = sys.stderr
        super().emit(record)


def configure_cli_logging(*, quiet: bool = False) -> None:
    """Install one stderr handler for CLI use (idempotent)."""

    for handler in list(logger.handlers):
        if getattr(handler, "_cryorole_cli", False):
            logger.removeHandler(handler)
    handler = _CurrentStderrHandler()
    handler._cryorole_cli = True  # type: ignore[attr-defined]
    handler.setFormatter(_CliFormatter())
    handler.setLevel(logging.INFO)
    logger.addHandler(handler)
    logger.setLevel(logging.WARNING if quiet else logging.INFO)
    logger.propagate = False
