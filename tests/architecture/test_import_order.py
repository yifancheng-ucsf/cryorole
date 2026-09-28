"""Entry-point modules import on their own in a fresh interpreter (no import cycles).

``cryorole.io.writers.landscape_store`` used to fail when imported first
(``landscape_store`` -> ``export`` -> ``export.metadata_subset`` ->
``landscape_store``); scripts had to ``import cryorole.cli.main`` first.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
MODULES = (
    "cryorole.io.writers.landscape_store",
    "cryorole.io.writers",
    "cryorole.export",
    "cryorole.export.landscape",
    "cryorole.export.metadata_subset",
    "cryorole.api",
    "cryorole.align",
    "cryorole.canonicalize.service",
    "cryorole.preflight",
    "cryorole.select",
    "cryorole.visualize",
    "cryorole.workflow",
    "cryorole.workflows.run_service",
)


@pytest.mark.parametrize("module", MODULES)
def test_module_imports_first_in_a_fresh_interpreter(module: str) -> None:
    result = subprocess.run([sys.executable, "-c", f"import {module}"], capture_output=True, text=True, cwd=ROOT)
    assert result.returncode == 0, result.stderr
