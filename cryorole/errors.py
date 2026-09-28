"""One user-facing error type for cryoROLE services, the CLI and the API.

``CryoroleError`` carries a stable machine-readable ``code``, a plain-language
``message`` and, when there is one, a ``remedy`` telling the user what to do.
It subclasses ``ValueError`` so existing callers that catch ``ValueError`` keep
working while services migrate.
"""

from __future__ import annotations

from typing import Any, Mapping


class CryoroleError(ValueError):
    """A problem the user can act on."""

    def __init__(
        self,
        code: str,
        message: str,
        remedy: str | None = None,
        *,
        details: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.remedy = remedy
        self.details = dict(details or {})

    def __str__(self) -> str:
        return self.message if not self.remedy else f"{self.message} {self.remedy}"

    def to_dict(self) -> dict[str, Any]:
        return {"code": self.code, "message": self.message, "remedy": self.remedy, "details": self.details}


class CancelledError(CryoroleError):
    """Raised when a caller's cancel token was triggered."""

    def __init__(self, stage: str | None = None) -> None:
        where = f" during {stage}" if stage else ""
        super().__init__(
            "cancelled",
            f"The operation was cancelled{where}; nothing was published.",
            "Start it again when ready; partial work was discarded.",
            details={"stage": stage},
        )


# Stable codes used across services. New codes are added here, never renamed.
ERROR_CODES = {
    "cancelled": "The caller cancelled the operation.",
    "input_not_found": "An input file or directory does not exist.",
    "input_invalid": "An input file exists but cannot be used as given.",
    "output_exists": "The output location exists and --overwrite was not given.",
    "run_dir_unresolved": "No run bundle was given and none could be resolved unambiguously.",
    "id_unresolved": "A canonical frame or selection id was not given and could not be resolved unambiguously.",
    "blocked": "Preflight found a blocking problem.",
    "internal": "An unexpected internal error.",
}


def classify_exception(exc: BaseException) -> CryoroleError:
    """Map a raw exception to a ``CryoroleError`` for display (never loses the message)."""

    if isinstance(exc, CryoroleError):
        return exc
    if isinstance(exc, FileExistsError):
        return CryoroleError(
            "output_exists",
            str(exc),
            None if "--overwrite" in str(exc) else "Use --overwrite to replace it, or choose another output name.",
        )
    if isinstance(exc, (FileNotFoundError, IsADirectoryError, NotADirectoryError)):
        return CryoroleError("input_not_found", str(exc), "Check the path (relative paths are resolved from the current directory).")
    if isinstance(exc, PermissionError):
        return CryoroleError("input_invalid", str(exc), "Check file permissions, or choose a writable output location.")
    if isinstance(exc, OSError):
        return CryoroleError("input_invalid", str(exc), "Check free disk space and that the location is writable.")
    if isinstance(exc, (ValueError, RuntimeError, NotImplementedError)):
        return CryoroleError("input_invalid", str(exc))
    return CryoroleError("internal", f"{type(exc).__name__}: {exc}",
                         "This is a bug. Please report it with the command you ran; set CRYOROLE_DEBUG=1 for a traceback.")
