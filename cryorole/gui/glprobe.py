"""Qt platform and OpenGL probe, run in a child process before the GUI opens.

Qt aborts the whole process when it cannot open a window (for example a
missing ``libxcb-cursor0``) or when a GL widget finds no usable GLX
configuration (for example a display driven by an onboard BMC chip while the
NVIDIA GPUs are only secondary outputs). Python cannot catch either abort, so
the GUI asks a short-lived child process first and decides from its result:

* platform failed: the GUI does not start; the error says what to install;
* OpenGL failed: the GUI starts without the 3D view;
* software OpenGL (llvmpipe, remote desktop): 3D works but is slow.

The child (``python -m cryorole.gui.glprobe``) prints one ``GLPROBE {json}``
line after Qt starts and one after the OpenGL check, so a later abort still
tells the parent how far it got.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
import shutil
import subprocess
import sys

from cryorole.gui.childproc import module_command

PROBE_TIMEOUT_SECONDS = 30.0
_PREFIX = "GLPROBE "
_SOFTWARE_RENDERERS = (
    "llvmpipe",
    "softpipe",
    "swrast",
    "software rasterizer",
    "swiftshader",
    "microsoft basic render",
)
PRIME_OFFLOAD_ENV = {"__NV_PRIME_RENDER_OFFLOAD": "1", "__GLX_VENDOR_LIBRARY_NAME": "nvidia"}


@dataclass(frozen=True)
class ProbeResult:
    """What the child process could do on this display."""

    platform_ok: bool
    gl_status: str  # "ok", "software", "failed" or "skipped"
    platform_name: str | None = None
    vendor: str | None = None
    renderer: str | None = None
    version: str | None = None
    reason: str | None = None
    remedy: str | None = None

    @property
    def gl_available(self) -> bool:
        return self.gl_status in ("ok", "software")

    def summary(self) -> str:
        if not self.platform_ok:
            return f"Qt cannot open a window: {self.reason}"
        if self.gl_status == "skipped":
            return "3D view off (--no-gl)"
        if self.gl_status == "failed":
            return f"3D view off: {self.reason}"
        text = f"3D: {self.renderer} (OpenGL {self.version})"
        return text + " - software rendering, slow at large point counts" if self.gl_status == "software" else text


def probe(*, check_gl: bool = True, timeout: float = PROBE_TIMEOUT_SECONDS) -> ProbeResult:
    """Run the probe child with the current environment and interpret its output."""

    command = module_command("cryorole.gui.glprobe")
    if not check_gl:
        command.append("--platform-only")
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        return ProbeResult(
            platform_ok=False,
            gl_status="failed",
            reason=f"the display check did not finish within {timeout:.0f} s",
            remedy="Check that the desktop session is responsive, or start with --no-gl.",
        )
    return interpret(completed.returncode, completed.stdout, completed.stderr, check_gl=check_gl)


def interpret(
    returncode: int,
    stdout: str,
    stderr: str,
    *,
    check_gl: bool = True,
    nvidia_gpus: int | None = None,
    environ: dict[str, str] | None = None,
    xcb_cursor_missing: bool | None = None,
) -> ProbeResult:
    """Turn the child's output into a :class:`ProbeResult` (pure; unit-tested)."""

    stages = {}
    for line in stdout.splitlines():
        if line.startswith(_PREFIX):
            try:
                payload = json.loads(line[len(_PREFIX):])
            except ValueError:
                continue
            stages[payload.get("stage")] = payload
    platform = stages.get("platform")
    if platform is None:
        reason = (
            _display_connection_error(stderr)
            or _last_message(stderr)
            or f"Qt exited with status {returncode} before opening a window"
        )
        remedy = _platform_remedy(stderr, xcb_cursor_missing)
        return ProbeResult(platform_ok=False, gl_status="failed", reason=reason, remedy=remedy)
    name = platform.get("platform")
    if not check_gl:
        return ProbeResult(platform_ok=True, gl_status="skipped", platform_name=name)
    gl = stages.get("gl")
    if gl is None or gl.get("error"):
        reason = (gl or {}).get("error") or _last_message(stderr) or "the OpenGL check did not finish"
        return ProbeResult(
            platform_ok=True,
            gl_status="failed",
            platform_name=name,
            reason=reason,
            remedy=_gl_remedy(name, nvidia_gpus, environ),
        )
    renderer = gl.get("renderer") or ""
    software = any(token in renderer.lower() for token in _SOFTWARE_RENDERERS)
    return ProbeResult(
        platform_ok=True,
        gl_status="software" if software else "ok",
        platform_name=name,
        vendor=gl.get("vendor"),
        renderer=renderer,
        version=gl.get("version"),
        reason="software OpenGL (no GPU acceleration on this display)" if software else None,
        remedy=_gl_remedy(name, nvidia_gpus, environ) if software else None,
    )


def count_nvidia_gpus() -> int:
    """Number of NVIDIA GPUs listed by ``nvidia-smi -L`` (0 if none or unavailable)."""

    executable = shutil.which("nvidia-smi")
    if executable is None:
        return 0
    try:
        listing = subprocess.run([executable, "-L"], capture_output=True, text=True, timeout=10, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return 0
    return sum(1 for line in listing.stdout.splitlines() if line.startswith("GPU "))


def prime_offload_command() -> str:
    return " ".join(f"{key}={value}" for key, value in PRIME_OFFLOAD_ENV.items()) + " cryorole gui"


def _gl_remedy(platform: str | None, nvidia_gpus: int | None, environ: dict[str, str] | None) -> str:
    environ = os.environ if environ is None else environ
    if platform != "xcb":  # PRIME offload and ssh forwarding are X11 matters
        return "The 2D and 1D views work without OpenGL."
    gpus = count_nvidia_gpus() if nvidia_gpus is None else nvidia_gpus
    offload_set = all(environ.get(key) == value for key, value in PRIME_OFFLOAD_ENV.items())
    if gpus and not offload_set:
        return (
            "This display does not use the NVIDIA GPU for OpenGL (for example a monitor on an onboard chip). "
            f"Start the GUI with: {prime_offload_command()}"
        )
    if environ.get("SSH_CONNECTION") and environ.get("DISPLAY", "").startswith("localhost:"):
        return "The display is forwarded over ssh; run the GUI at the workstation's own desktop for GPU rendering."
    return "The 2D and 1D views work without OpenGL."


def _platform_remedy(stderr: str, xcb_cursor_missing: bool | None) -> str:
    # Qt prints its libxcb-cursor0 hint whenever the xcb plugin fails, so check the
    # display connection first and the library itself before blaming it.
    if _display_connection_error(stderr):
        return "Start `cryorole gui` from a terminal on the workstation's own desktop, or use a remote desktop."
    if xcb_cursor_missing is None:
        xcb_cursor_missing = "xcb-cursor" in stderr.lower() and _library_missing("libxcb-cursor.so.0")
    if xcb_cursor_missing:
        return (
            "Qt 6.5 and later need the system library libxcb-cursor0: `sudo apt install libxcb-cursor0` "
            "(or ask the administrator)."
        )
    return "Run `cryorole gui` from a desktop session. Set QT_DEBUG_PLUGINS=1 to see why Qt failed."


def _display_connection_error(stderr: str) -> str | None:
    for line in stderr.splitlines():
        if "could not connect to display" in line.lower() or "cannot open display" in line.lower():
            return line.strip().removeprefix("qt.qpa.xcb: ").strip()
    return None


def _library_missing(soname: str) -> bool:
    import ctypes

    try:
        ctypes.CDLL(soname)
    except OSError:
        return True
    return False


def _last_message(stderr: str) -> str | None:
    lines = [line.strip() for line in stderr.splitlines() if line.strip()]
    for line in reversed(lines):
        lowered = line.lower()
        if "could not" in lowered or "failed" in lowered or "error" in lowered or "cannot" in lowered:
            return line.removeprefix("WARNING: ")
    return lines[-1] if lines else None


def _emit(payload: dict[str, object]) -> None:
    print(_PREFIX + json.dumps(payload), flush=True)


def _gl_strings() -> dict[str, str]:
    from OpenGL import GL

    values = {}
    for name, code in (("vendor", GL.GL_VENDOR), ("renderer", GL.GL_RENDERER), ("version", GL.GL_VERSION)):
        raw = GL.glGetString(code)
        values[name] = raw.decode("utf-8", "replace") if isinstance(raw, bytes) else str(raw or "")
    return values


FRAME_TIMEOUT_MS = 5_000
# Qt platforms without OpenGL widgets: report at once instead of waiting for a frame.
_NO_GL_PLATFORMS = ("offscreen", "minimal", "vnc", "linuxfb")


def _child_main(argv: list[str]) -> int:
    # The OpenGL check builds what the 3D view builds: a pyqtgraph GLViewWidget,
    # hidden from the screen, rendering one frame. A bare offscreen context is not
    # enough: it can succeed on displays where the real widget aborts.
    from PySide6.QtCore import Qt, QTimer
    from PySide6.QtWidgets import QApplication

    app = QApplication(argv[:1] or ["cryorole-glprobe"])  # Qt aborts here if no platform plugin loads
    _emit({"stage": "platform", "platform": app.platformName()})
    if "--platform-only" in argv:
        return 0
    if app.platformName() in _NO_GL_PLATFORMS:
        _emit({"stage": "gl", "error": f"the Qt platform '{app.platformName()}' has no OpenGL widgets"})
        return 0
    import pyqtgraph.opengl as gl

    result: dict[str, str] = {}

    class _ProbeView(gl.GLViewWidget):
        def paintGL(self, *args, **kwargs):  # noqa: N802 - Qt override
            super().paintGL(*args, **kwargs)
            if not result:
                try:
                    result.update(_gl_strings())
                except Exception as exc:  # noqa: BLE001 - any failure here means "no usable OpenGL"
                    result["error"] = f"OpenGL query failed: {exc}"
                QTimer.singleShot(0, app.quit)

    view = _ProbeView()
    view.setAttribute(Qt.WidgetAttribute.WA_DontShowOnScreen, True)
    view.resize(64, 64)
    view.show()  # a GLX failure aborts the process here or at the first frame
    QTimer.singleShot(FRAME_TIMEOUT_MS, app.quit)
    app.exec()
    if not result:
        _emit({"stage": "gl", "error": f"no OpenGL frame was drawn within {FRAME_TIMEOUT_MS // 1000} s"})
    elif "error" in result:
        _emit({"stage": "gl", "error": result["error"]})
    elif not result.get("renderer"):
        _emit({"stage": "gl", "error": "OpenGL returned no renderer name"})
    else:
        _emit({"stage": "gl", **result})
    return 0


if __name__ == "__main__":
    raise SystemExit(_child_main(sys.argv))
