"""Public terminology and mode-aware selection guidance."""

import pytest

from cryorole.cli.main import build_parser
from cryorole.select.service import SelectRequest, create_selection
from cryorole.visualize.service import VisualizationRequest, visualize


@pytest.mark.parametrize("option", ["--opacity", "--alpha"])
@pytest.mark.parametrize("value", ["0", "0.4", "1"])
def test_opacity_aliases_resolve_identically(option, value):
    args = build_parser().parse_args(["visualize", "--run-dir", "RUN", option, value])
    assert VisualizationRequest.from_namespace(args).alpha == float(value)


def test_opacity_help_and_conflict(capsys):
    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["visualize", "--help"])
    assert exc.value.code == 0
    help_text = capsys.readouterr().out
    assert "--opacity" in help_text
    assert "--alpha" not in help_text
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["visualize", "--run-dir", "RUN", "--alpha", "0.5", "--opacity", "0.5"])
    assert exc.value.code == 2


def test_selection_help_groups_and_missing_name(capsys):
    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["select", "--help"])
    assert exc.value.code == 0
    help_text = capsys.readouterr().out
    for mode in ("radius", "threshold", "range", "random", "metadata"):
        assert f"{mode} mode" in help_text
        assert f"--mode {mode}" in help_text
    assert "No default" in help_text
    with pytest.raises(SystemExit):
        parser.parse_args(["select", "--run-dir", "RUN"])
    assert "region_01" in capsys.readouterr().err


@pytest.mark.parametrize("mode,options,expected", [
    ("radius", [], "--center"),
    ("radius", ["--center", "0", "0", "0"], "--radius"),
    ("threshold", [], "--sld-min"),
    ("range", [], "--range-bound"),
    ("random", [], "--fraction"),
    ("metadata", [], "--metadata-domain"),
    ("random", ["--fraction", "0.5", "--center-representation", "euler"], "--center-representation"),
    ("threshold", ["--sld-min", "1", "--metric", "so3"], "--metric"),
    ("radius", ["--center", "0", "0", "0", "--radius", "10", "--seed", "0"], "--seed"),
])
def test_selection_request_errors_before_reading_files(tmp_path, mode, options, expected):
    run = tmp_path / "missing_run"
    args = build_parser().parse_args([
        "select", "--run-dir", str(run), "--selection-id", "region_01", "--mode", mode, *options,
    ])
    with pytest.raises(ValueError, match=expected):
        create_selection(SelectRequest.from_namespace(args))
    assert not run.exists()


def test_parser_defaults_do_not_count_as_explicit_mode_options():
    args = build_parser().parse_args([
        "select", "--run-dir", "RUN", "--selection-id", "sample", "--mode", "random", "--fraction", "0.5",
    ])
    request = SelectRequest.from_namespace(args)
    assert "metric" not in request.explicit_options
    assert "center_representation" not in request.explicit_options


@pytest.mark.parametrize("value", [-0.01, 1.01, float("nan"), float("inf")])
def test_opacity_invalid_values_use_public_name(tmp_path, value):
    with pytest.raises(ValueError, match="--opacity"):
        visualize(VisualizationRequest(run_dir=str(tmp_path), alpha=value))
    assert not list(tmp_path.iterdir())


def test_opacity_omission_preserves_style_default():
    args = build_parser().parse_args(["visualize", "--run-dir", "RUN"])
    assert VisualizationRequest.from_namespace(args).alpha is None


def test_next_command_quotes_paths_without_json_escaping(monkeypatch):
    from cryorole.cli.commands import downstream

    value = "C:\\data\\particle's run $1"
    monkeypatch.setattr(downstream.os, "name", "nt")
    assert downstream._quote_argument(value) == "'C:\\data\\particle''s run $1'"
    monkeypatch.setattr(downstream.os, "name", "posix")
    import shlex
    assert shlex.split(downstream._quote_argument(value)) == [value]


@pytest.mark.parametrize("name", [".", "..", "../data", "a/b", "a\\b", "C:run", "a\nb"])
def test_selection_name_cannot_redirect_overwrite(tmp_path, name):
    with pytest.raises(ValueError, match="name, not a path"):
        create_selection(SelectRequest(run_dir=str(tmp_path), selection_id=name,
                                       selection_mode="threshold", sld_min=0, overwrite=True))
    assert not list(tmp_path.iterdir())
