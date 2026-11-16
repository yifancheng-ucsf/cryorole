"""cryoROLE desktop GUI (PySide6 + pyqtgraph), installed with ``pip install "cryorole[gui]"``.

The GUI is a frontend, like the CLI: it calls ``cryorole.api`` and contains
no science. Importing this package does not import Qt; :func:`launch` checks
the optional dependencies and the display first, so a missing library or an
unusable display becomes a :class:`~cryorole.errors.CryoroleError` instead of
a Qt abort.
"""

from __future__ import annotations

import importlib.util
import os
import sys

from cryorole.errors import CryoroleError
from cryorole.logs import notify_user, warn_user

GUI_DEPENDENCIES = ("PySide6", "pyqtgraph", "OpenGL")


def missing_dependencies() -> list[str]:
    return [name for name in GUI_DEPENDENCIES if importlib.util.find_spec(name) is None]


def check_display_environment(environ: dict[str, str] | None = None, platform: str | None = None) -> None:
    """Fail early, with a remedy, when a Linux session has no display at all."""

    environ = os.environ if environ is None else environ
    platform = sys.platform if platform is None else platform
    if not platform.startswith("linux") or environ.get("QT_QPA_PLATFORM"):
        return
    if not environ.get("DISPLAY") and not environ.get("WAYLAND_DISPLAY"):
        raise CryoroleError(
            "gui_no_display",
            "No desktop display is available (DISPLAY and WAYLAND_DISPLAY are not set).",
            "Start `cryorole gui` from a terminal on the workstation's own desktop, or use a remote desktop.",
        )


def launch(*, no_gl: bool = False, run_dir: str | None = None) -> int:
    """Check dependencies and display, then run the GUI until its window closes."""

    missing = missing_dependencies()
    if missing:
        raise CryoroleError(
            "gui_unavailable",
            f"The GUI needs packages that are not installed: {', '.join(missing)}.",
            'Install them with: pip install "cryorole[gui]"',
        )
    check_display_environment()
    from cryorole.gui.glprobe import probe

    result = probe(check_gl=not no_gl)
    if not result.platform_ok:
        raise CryoroleError("gui_display_unavailable", f"Qt cannot open a window: {result.reason}", result.remedy)
    if result.gl_status == "failed":
        warn_user(f"3D view disabled: {result.reason}. {result.remedy}")
    elif result.gl_status == "software":
        warn_user(f"3D view uses software OpenGL ({result.renderer}); large point counts will be slow. {result.remedy}")
    else:
        notify_user(result.summary())
    from cryorole.gui.app import run_app

    return run_app(result, run_dir=run_dir)
