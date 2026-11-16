"""``cryorole gui`` start-up checks (no Qt needed: the app module is replaced)."""

from __future__ import annotations

import logging
import sys
import types

import pytest

import cryorole.gui as gui
from cryorole.cli.main import main
from cryorole.errors import CryoroleError
from cryorole.gui.glprobe import ProbeResult


@pytest.fixture
def fake_app(monkeypatch):
    calls = []
    module = types.ModuleType("cryorole.gui.app")
    module.run_app = lambda probe, run_dir=None: calls.append(probe) or 0
    monkeypatch.setitem(sys.modules, "cryorole.gui.app", module)
    monkeypatch.setattr(gui, "missing_dependencies", lambda: [])
    monkeypatch.setattr(gui, "check_display_environment", lambda: None)
    return calls


def _set_probe(monkeypatch, result: ProbeResult, seen: list | None = None):
    def fake_probe(*, check_gl=True):
        if seen is not None:
            seen.append(check_gl)
        return result

    monkeypatch.setattr("cryorole.gui.glprobe.probe", fake_probe)


def test_missing_extra_is_a_cryorole_error_with_install_remedy(monkeypatch):
    monkeypatch.setattr(gui, "missing_dependencies", lambda: ["PySide6", "pyqtgraph"])
    with pytest.raises(CryoroleError) as info:
        gui.launch()
    assert info.value.code == "gui_unavailable"
    assert 'pip install "cryorole[gui]"' in info.value.remedy


def test_linux_without_any_display_fails_before_qt_starts():
    with pytest.raises(CryoroleError) as info:
        gui.check_display_environment({}, "linux")
    assert info.value.code == "gui_no_display"


@pytest.mark.parametrize(
    "environ, platform",
    [({"DISPLAY": ":1"}, "linux"), ({"WAYLAND_DISPLAY": "wayland-0"}, "linux"),
     ({"QT_QPA_PLATFORM": "offscreen"}, "linux"), ({}, "darwin")],
)
def test_display_check_accepts_desktop_sessions(environ, platform):
    gui.check_display_environment(environ, platform)


def test_platform_failure_stops_launch(monkeypatch, fake_app):
    _set_probe(monkeypatch, ProbeResult(platform_ok=False, gl_status="failed", reason="no plugin", remedy="install X"))
    with pytest.raises(CryoroleError) as info:
        gui.launch()
    assert info.value.code == "gui_display_unavailable"
    assert info.value.remedy == "install X"
    assert fake_app == []


def test_gl_failure_starts_without_3d_and_warns(monkeypatch, fake_app, caplog):
    result = ProbeResult(platform_ok=True, gl_status="failed", platform_name="xcb",
                         reason="Could not initialize GLX", remedy="use PRIME")
    _set_probe(monkeypatch, result)
    with caplog.at_level(logging.WARNING, logger="cryorole"):
        assert gui.launch() == 0
    assert fake_app == [result]
    assert "3D view disabled: Could not initialize GLX. use PRIME" in caplog.text


def test_no_gl_flag_runs_platform_check_only(monkeypatch, fake_app):
    seen = []
    _set_probe(monkeypatch, ProbeResult(platform_ok=True, gl_status="skipped", platform_name="xcb"), seen)
    assert main(["gui", "--no-gl"]) == 0
    assert seen == [False]


def test_cli_reports_missing_extra_cleanly(monkeypatch, capsys):
    monkeypatch.setattr(gui, "missing_dependencies", lambda: ["PySide6"])
    with pytest.raises(SystemExit) as info:
        main(["gui"])
    assert info.value.code == 2
    err = capsys.readouterr().err
    assert "cryorole: error: The GUI needs packages that are not installed: PySide6." in err
    assert 'pip install "cryorole[gui]"' in err
