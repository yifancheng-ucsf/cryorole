"""The CLI error boundary, --version, logging and the flag resolver."""

from __future__ import annotations

import pytest

import cryorole.cli.main as cli_main
from cryorole.errors import CancelledError, CryoroleError, classify_exception
from cryorole.logs import configure_cli_logging, notify_user, warn_user
from cryorole.workflow.resolve import (
    resolve_canonical_id,
    resolve_run_dir,
    resolve_selection_id,
)


def _patch_status_handler(monkeypatch, exc: BaseException) -> None:
    def boom(args):
        raise exc

    original = cli_main.build_parser

    def build():
        parser = original()
        parser._subparsers._group_actions[0].choices["status"].set_defaults(handler=boom)
        return parser

    monkeypatch.setattr(cli_main, "build_parser", build)


def test_version_prints_package_version(capsys) -> None:
    from cryorole import __version__

    with pytest.raises(SystemExit) as excinfo:
        cli_main.main(["--version"])
    assert excinfo.value.code == 0
    out = capsys.readouterr().out.strip()
    assert out.startswith("cryorole ")
    assert out.split()[1]  # installed distribution version or cryorole.__version__
    assert __version__


def test_root_help_shows_workflow_examples(capsys) -> None:
    with pytest.raises(SystemExit):
        cli_main.main(["--help"])
    out = capsys.readouterr().out
    assert "Typical workflow:" in out
    assert "cryorole preflight --ref ref.cs --mov mov.cs" in out


def test_expected_errors_print_message_and_remedy_exit_2(monkeypatch, capsys) -> None:
    _patch_status_handler(monkeypatch, FileExistsError("Output directory already exists: out"))
    with pytest.raises(SystemExit) as excinfo:
        cli_main.main(["status", "--run-dir", "x"])
    err = capsys.readouterr().err
    assert excinfo.value.code == 2
    assert "cryorole: error: Output directory already exists: out" in err
    assert "  what to do: Use --overwrite to replace it" in err
    assert "Traceback" not in err


def test_cryorole_error_remedy_is_printed(monkeypatch, capsys) -> None:
    _patch_status_handler(monkeypatch, CryoroleError("input_invalid", "Bad thing.", "Do the other thing."))
    with pytest.raises(SystemExit) as excinfo:
        cli_main.main(["status", "--run-dir", "x"])
    err = capsys.readouterr().err
    assert excinfo.value.code == 2
    assert "cryorole: error: Bad thing.\n  what to do: Do the other thing.\n" in err


def test_unexpected_errors_exit_1_with_bug_report_hint(monkeypatch, capsys) -> None:
    _patch_status_handler(monkeypatch, KeyError("boom"))
    with pytest.raises(SystemExit) as excinfo:
        cli_main.main(["status", "--run-dir", "x"])
    err = capsys.readouterr().err
    assert excinfo.value.code == 1
    assert "KeyError" in err and "CRYOROLE_DEBUG=1" in err


def test_debug_env_reraises(monkeypatch) -> None:
    monkeypatch.setenv("CRYOROLE_DEBUG", "1")
    _patch_status_handler(monkeypatch, KeyError("boom"))
    with pytest.raises(KeyError):
        cli_main.main(["status", "--run-dir", "x"])


def test_keyboard_interrupt_exits_130(monkeypatch, capsys) -> None:
    _patch_status_handler(monkeypatch, KeyboardInterrupt())
    with pytest.raises(SystemExit) as excinfo:
        cli_main.main(["status", "--run-dir", "x"])
    assert excinfo.value.code == 130
    assert "interrupted" in capsys.readouterr().err


def test_classify_exception_codes() -> None:
    assert classify_exception(FileNotFoundError("x")).code == "input_not_found"
    assert classify_exception(ValueError("x")).code == "input_invalid"
    assert classify_exception(KeyError("x")).code == "internal"
    error = CryoroleError("blocked", "m", "r", details={"a": 1})
    assert classify_exception(error) is error
    assert error.to_dict() == {"code": "blocked", "message": "m", "remedy": "r", "details": {"a": 1}}
    assert isinstance(error, ValueError)  # existing ``except ValueError`` callers keep working
    assert CancelledError("density").details == {"stage": "density"}


def test_warn_user_goes_to_stderr_only_after_cli_logging(capsys) -> None:
    warn_user("not shown without a handler")
    assert capsys.readouterr().err == ""
    configure_cli_logging()
    warn_user("parent landscape lacks metadata")
    notify_user("plain notice")
    err = capsys.readouterr().err
    assert "[cryorole] warning: parent landscape lacks metadata" in err
    assert "plain notice" in err
    configure_cli_logging(quiet=True)
    notify_user("hidden when quiet")
    assert "hidden when quiet" not in capsys.readouterr().err


# --- resolver rules ---------------------------------------------------------


def _bundle(path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    (path / "run_manifest.json").write_text("{}", encoding="utf-8")


def test_resolve_run_dir_order(tmp_path) -> None:
    assert resolve_run_dir("given", cwd=tmp_path).resolved_by == "explicit"
    with pytest.raises(CryoroleError) as excinfo:
        resolve_run_dir(None, cwd=tmp_path)
    assert excinfo.value.code == "run_dir_unresolved"
    assert excinfo.value.details["candidates"] == []

    _bundle(tmp_path / "b1")
    _bundle(tmp_path / "b2")
    with pytest.raises(CryoroleError) as excinfo:
        resolve_run_dir(None, cwd=tmp_path)
    assert excinfo.value.details["candidates"] == ["b1", "b2"]  # listed, never picked

    _bundle(tmp_path / "cryorole_outputs")
    resolved = resolve_run_dir(None, cwd=tmp_path)
    assert (resolved.value, resolved.resolved_by) == ("cryorole_outputs", "default_output_dir")

    resolved = resolve_run_dir(None, cwd=tmp_path / "b1")
    assert (resolved.value, resolved.resolved_by) == (".", "current_directory")


def test_resolve_canonical_and_selection_ids(tmp_path) -> None:
    _bundle(tmp_path)
    with pytest.raises(CryoroleError, match="no canonical frame yet"):
        resolve_canonical_id(tmp_path, None)
    with pytest.raises(CryoroleError, match="no selections yet"):
        resolve_selection_id(tmp_path, None)

    (tmp_path / "canonical" / "only").mkdir(parents=True)
    (tmp_path / "canonical" / "only" / "canonical_frame.json").write_text("{}", encoding="utf-8")
    (tmp_path / "canonical" / "empty_dir").mkdir()  # not a frame
    assert resolve_canonical_id(tmp_path, None).value == "only"
    assert resolve_canonical_id(tmp_path, "other").resolved_by == "explicit"

    for name in ("s1", "s2"):
        (tmp_path / "selections" / name).mkdir(parents=True)
        (tmp_path / "selections" / name / "selection.json").write_text("{}", encoding="utf-8")
    with pytest.raises(CryoroleError) as excinfo:
        resolve_selection_id(tmp_path, None)
    assert excinfo.value.details["candidates"] == ["s1", "s2"]


def test_every_help_example_parses() -> None:
    import shlex

    from cryorole.cli.parsers import _COMMAND_EXAMPLES

    parser = cli_main.build_parser()
    for examples in _COMMAND_EXAMPLES.values():
        for line in examples:
            parser.parse_args(shlex.split(line.split("#")[0])[1:])
