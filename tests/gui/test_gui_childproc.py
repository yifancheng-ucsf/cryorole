"""The GUI's child processes import the installed modules, not files where the GUI was started."""

from __future__ import annotations

import subprocess
import sys

from cryorole.gui.childproc import module_command


def _shadowing_folder(tmp_path):
    (tmp_path / "json.py").write_text("raise ImportError('shadowed json')\n", encoding="utf-8")
    (tmp_path / "cryorole").mkdir()
    (tmp_path / "cryorole" / "__init__.py").write_text("raise ImportError('shadowed cryorole')\n", encoding="utf-8")
    return tmp_path


def _run(command, cwd):
    return subprocess.run(command, input="", capture_output=True, text=True, cwd=cwd, timeout=300, check=False)


def test_worker_command_ignores_modules_in_the_working_directory(tmp_path):
    folder = _shadowing_folder(tmp_path)
    plain = _run([sys.executable, "-m", "cryorole.gui.worker"], folder)
    assert "shadowed" in plain.stderr  # what `python -m` did
    child = _run(module_command("cryorole.gui.worker"), folder)
    assert "shadowed" not in child.stderr, child.stderr


def test_module_command_passes_arguments():
    command = module_command("cryorole.gui.glprobe", "--platform-only")
    assert command[0] == sys.executable and command[-1] == "--platform-only"
