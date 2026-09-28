"""Regression tests for the CLI/service ownership boundary."""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from cryorole.cli.main import build_parser


ROOT = Path(__file__).resolve().parents[2]
SERVICE_FILES = (
    ROOT / "cryorole" / "workflows" / "run_service.py",
    ROOT / "cryorole" / "canonicalize" / "service.py",
    ROOT / "cryorole" / "visualize" / "service.py",
    ROOT / "cryorole" / "select" / "service.py",
    ROOT / "cryorole" / "preflight" / "service.py",
)


def _imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            modules.add(node.module)
    return modules


def test_domain_services_do_not_import_cli() -> None:
    offenders = {
        str(path.relative_to(ROOT)): sorted(
            module for module in _imports(path) if module == "cryorole.cli" or module.startswith("cryorole.cli.")
        )
        for path in SERVICE_FILES
    }
    assert not {path: modules for path, modules in offenders.items() if modules}


def test_root_main_has_no_scientific_array_stack_imports() -> None:
    imports = _imports(ROOT / "cryorole" / "cli" / "main.py")
    forbidden = {
        module
        for module in imports
        if module == "numpy"
        or module.startswith("numpy.")
        or module == "pandas"
        or module.startswith("pandas.")
        or module == "scipy"
        or module.startswith("scipy.")
    }
    assert not forbidden


@pytest.mark.parametrize(
    "command",
    [
        "align",
        "preflight",
        "run",
        "status",
        "next",
        "guide",
        "explore",
        "visualize",
        "animate",
        "canonical-views",
        "canonicalize",
        "select",
        "export",
        "manifest",
    ],
)
def test_every_public_subcommand_help_renders(command: str) -> None:
    parser = build_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args([command, "--help"])
    assert exc_info.value.code == 0


def test_input_policy_resolver_is_shared_by_run_preflight_and_guide() -> None:
    paths = (
        ROOT / "cryorole" / "workflows" / "run_service.py",
        ROOT / "cryorole" / "preflight" / "service.py",
        ROOT / "cryorole" / "workflow" / "guide.py",
    )
    for path in paths:
        source = path.read_text(encoding="utf-8")
        assert "resolve_input_policies" in source, path


def test_run_quicklook_uses_typed_visualization_service() -> None:
    source = (ROOT / "cryorole" / "workflows" / "run_service.py").read_text(
        encoding="utf-8"
    )
    assert "from cryorole.visualize import QuickLookRequest, write_quicklook" in source
    assert "write_landscape_visualizations" not in source


# Library code below ``cryorole/`` never prints, exits or parses argv: the CLI
# and API frontends own presentation. ``workflows/rotate_landscape.py`` is a
# legacy stand-alone diagnostic script (``scripts/rotate_landscape.py``) whose
# ``main()`` is the only exception; its library function
# ``rotate_landscape_bundle`` is print-free.
_FRONTEND_EXCEPTIONS = {
    ROOT / "cryorole" / "workflows" / "rotate_landscape.py": {"print", "argparse"},
}


def _frontend_calls(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Name) and func.id in {"print", "exit", "quit", "input"}:
                found.add(func.id)
            elif (
                isinstance(func, ast.Attribute)
                and isinstance(func.value, ast.Name)
                and func.value.id == "sys"
                and func.attr == "exit"
            ):
                found.add("sys.exit")
        elif isinstance(node, ast.Import) and any(alias.name == "argparse" for alias in node.names):
            found.add("argparse")
        elif isinstance(node, ast.ImportFrom) and node.module == "argparse":
            found.add("argparse")
    return found


def test_library_code_outside_cli_has_no_print_exit_or_argparse() -> None:
    offenders = {}
    for path in sorted((ROOT / "cryorole").rglob("*.py")):
        if (ROOT / "cryorole" / "cli") in path.parents:
            continue
        found = _frontend_calls(path) - _FRONTEND_EXCEPTIONS.get(path, set())
        if found:
            offenders[str(path.relative_to(ROOT))] = sorted(found)
    assert not offenders


def test_api_package_does_not_import_cli() -> None:
    api_dir = ROOT / "cryorole" / "api"
    offenders = {
        str(path.relative_to(ROOT)): sorted(
            module for module in _imports(path) if module == "cryorole.cli" or module.startswith("cryorole.cli.")
        )
        for path in api_dir.rglob("*.py")
    }
    assert not {path: modules for path, modules in offenders.items() if modules}
