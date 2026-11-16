"""The GUI is a frontend: Qt stays inside cryorole/gui, science stays out of it."""

from __future__ import annotations

import ast
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
GUI_DIR = ROOT / "cryorole" / "gui"
QT_STACK = ("PySide6", "PySide2", "PyQt5", "PyQt6", "pyqtgraph", "OpenGL", "vispy")
# cryoROLE modules the GUI may import. Everything scientific is reached through cryorole.api.
GUI_ALLOWED = ("cryorole.api", "cryorole.errors", "cryorole.logs", "cryorole.gui")


def _imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            if node.module == "cryorole":  # `from cryorole import api` imports cryorole.api
                modules.update(f"cryorole.{alias.name}" for alias in node.names)
            else:
                modules.add(node.module)
    return modules


def _matches(module: str, prefixes) -> bool:
    return any(module == prefix or module.startswith(prefix + ".") for prefix in prefixes)


def test_qt_stack_is_imported_only_inside_the_gui_package():
    offenders = {}
    for path in sorted((ROOT / "cryorole").rglob("*.py")):
        if GUI_DIR in path.parents:
            continue
        found = sorted(module for module in _imports(path) if _matches(module, QT_STACK))
        if found:
            offenders[str(path.relative_to(ROOT))] = found
    assert not offenders


def test_gui_reaches_cryorole_only_through_the_api():
    offenders = {}
    for path in sorted(GUI_DIR.rglob("*.py")):
        found = sorted(
            module
            for module in _imports(path)
            if _matches(module, ("cryorole",)) and not _matches(module, GUI_ALLOWED)
        )
        if found:
            offenders[str(path.relative_to(ROOT))] = found
    assert not offenders


def test_importing_cryorole_gui_and_the_cli_does_not_import_qt():
    code = (
        "import sys, cryorole.gui, cryorole.cli.main\n"
        f"loaded = [m for m in sys.modules if m.split('.')[0] in {QT_STACK!r}]\n"
        "assert not loaded, loaded\n"
    )
    completed = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=ROOT)
    assert completed.returncode == 0, completed.stderr
