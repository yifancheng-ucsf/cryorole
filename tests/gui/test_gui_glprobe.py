"""Interpretation of the Qt/OpenGL probe child output (no Qt needed).

The stderr texts are the ones seen on Linux workstations during the GUI spike:
one without libxcb-cursor0, one with its monitor on the ASPEED BMC chip and no
PRIME offload, and an ssh session without a display.
"""

from __future__ import annotations

import pytest

from cryorole.gui import glprobe
from cryorole.gui.glprobe import interpret

PLATFORM_XCB = 'GLPROBE {"stage": "platform", "platform": "xcb"}\n'
XCB_PLUGIN_FAILURE = (
    "qt.qpa.plugin: From 6.5.0, xcb-cursor0 or libxcb-cursor0 is needed to load the Qt xcb platform plugin.\n"
    'qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.\n'
    "This application failed to start because no Qt platform plugin could be initialized.\n"
)
NO_GLX = (
    "qt.glx: qglx_findConfig: Failed to finding matching FBConfig for QSurfaceFormat(version 2.0)\n"
    "Could not initialize GLX\n"
)
PRIME = {"__NV_PRIME_RENDER_OFFLOAD": "1", "__GLX_VENDOR_LIBRARY_NAME": "nvidia"}


def _gl_line(renderer: str, vendor: str = "NVIDIA Corporation", version: str = "4.6.0") -> str:
    return f'GLPROBE {{"stage": "gl", "vendor": "{vendor}", "renderer": "{renderer}", "version": "{version}"}}\n'


@pytest.fixture(autouse=True)
def _no_real_nvidia_smi(monkeypatch):
    def fail():
        raise AssertionError("nvidia-smi must not be queried in this test")

    monkeypatch.setattr(glprobe, "count_nvidia_gpus", fail)


def test_missing_libxcb_cursor_blocks_the_gui_with_install_remedy():
    result = interpret(-6, "", XCB_PLUGIN_FAILURE, xcb_cursor_missing=True)
    assert not result.platform_ok
    assert "no Qt platform plugin could be initialized" in result.reason
    assert "libxcb-cursor0" in result.remedy


def test_no_display_is_not_blamed_on_libxcb_cursor():
    stderr = "qt.qpa.xcb: could not connect to display :77\n" + XCB_PLUGIN_FAILURE
    result = interpret(-6, "", stderr, xcb_cursor_missing=True)
    assert not result.platform_ok
    assert result.reason == "could not connect to display :77"
    assert "own desktop" in result.remedy
    assert "libxcb-cursor0" not in result.remedy


def test_qt_cursor_hint_is_ignored_when_the_library_is_present():
    result = interpret(-6, "", XCB_PLUGIN_FAILURE, xcb_cursor_missing=False)
    assert not result.platform_ok
    assert "libxcb-cursor0" not in result.remedy


def test_glx_abort_after_platform_start_disables_3d_and_suggests_prime():
    result = interpret(-6, PLATFORM_XCB, NO_GLX, nvidia_gpus=2, environ={})
    assert result.platform_ok
    assert result.gl_status == "failed"
    assert not result.gl_available
    assert result.reason == "Could not initialize GLX"
    assert "__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia cryorole gui" in result.remedy
    assert result.summary() == "3D view off: Could not initialize GLX"


def test_prime_is_not_suggested_again_when_already_set():
    result = interpret(-6, PLATFORM_XCB, NO_GLX, nvidia_gpus=2, environ=dict(PRIME))
    assert result.gl_status == "failed"
    assert "PRIME" not in result.remedy and "__NV" not in result.remedy


def test_software_renderer_keeps_3d_but_warns():
    stdout = PLATFORM_XCB + _gl_line("llvmpipe (LLVM 20.1.2, 256 bits)", vendor="Mesa", version="4.5 Mesa 25.2.8")
    result = interpret(0, stdout, "", nvidia_gpus=2, environ={})
    assert result.gl_status == "software"
    assert result.gl_available
    assert "llvmpipe" in result.summary() and "software rendering" in result.summary()
    assert "__NV_PRIME_RENDER_OFFLOAD" in result.remedy


def test_nvidia_renderer_is_ok_without_remedy():
    stdout = PLATFORM_XCB + _gl_line("NVIDIA GeForce RTX 4090/PCIe/SSE2", version="4.6.0 NVIDIA 595.58.03")
    result = interpret(0, stdout, "", nvidia_gpus=2, environ={})
    assert result.gl_status == "ok"
    assert result.remedy is None
    assert result.summary() == "3D: NVIDIA GeForce RTX 4090/PCIe/SSE2 (OpenGL 4.6.0 NVIDIA 595.58.03)"


def test_offscreen_platform_never_suggests_prime_or_queries_gpus():
    stdout = (
        'GLPROBE {"stage": "platform", "platform": "offscreen"}\n'
        "GLPROBE {\"stage\": \"gl\", \"error\": \"the Qt platform 'offscreen' has no OpenGL widgets\"}\n"
    )
    result = interpret(0, stdout, "")  # count_nvidia_gpus would fail if called
    assert result.gl_status == "failed"
    assert result.remedy == "The 2D and 1D views work without OpenGL."


def test_forwarded_ssh_display_gets_its_own_remedy():
    environ = {"SSH_CONNECTION": "10.0.0.1 1 10.0.0.2 22", "DISPLAY": "localhost:10.0"}
    result = interpret(-6, PLATFORM_XCB, NO_GLX, nvidia_gpus=0, environ=environ)
    assert "forwarded over ssh" in result.remedy


def test_platform_only_check_skips_gl():
    result = interpret(0, PLATFORM_XCB, "", check_gl=False)
    assert result.platform_ok and result.gl_status == "skipped"
    assert result.summary() == "3D view off (--no-gl)"


def test_malformed_probe_lines_are_ignored():
    result = interpret(-6, "GLPROBE {not json\n", "", xcb_cursor_missing=False)
    assert not result.platform_ok
    assert "exited with status -6" in result.reason
