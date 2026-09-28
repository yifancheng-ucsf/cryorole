"""cryorole.api: the typed Python entry points used by notebooks and the future GUI."""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cryorole import api

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run"))
from test_core_productionization import _write_cs  # noqa: E402

N = 300
UIDS = list(range(7000, 7000 + N))


@pytest.fixture()
def inputs(tmp_path):
    rng = np.random.default_rng(11)
    ref_rv = rng.normal(0.0, 0.6, (N, 3))
    mov_rv = (Rotation.from_rotvec(ref_rv) * Rotation.from_rotvec(rng.normal(0.0, 0.15, (N, 3)))).as_rotvec()
    _write_cs(tmp_path / "ref.cs", UIDS, ref_rv)
    _write_cs(tmp_path / "mov.cs", UIDS, mov_rv)
    return tmp_path / "ref.cs", tmp_path / "mov.cs"


def test_api_end_to_end_is_silent_and_reports_progress(inputs, tmp_path, capsys, monkeypatch) -> None:
    ref, mov = inputs
    check = api.preflight(ref, mov, k_neighbors=10)
    assert check.report["readiness"] in {"READY", "READY_WITH_WARNINGS"}

    events: list[api.ProgressEvent] = []
    result = api.run(ref, mov, output_dir=tmp_path / "cryorole_outputs", k_neighbors=10, no_visualize=True,
                     progress=events.append)
    assert result.matched_count == N
    assert {"stage_start", "stage_done"} <= {event.kind for event in events}
    assert all(event.command == "run" for event in events)

    monkeypatch.chdir(tmp_path)  # run_dir omitted from here on: ./cryorole_outputs
    canonical = api.canonicalize(no_visualize=True)
    assert canonical.output_dir.name == "default"
    assert canonical.summary["resolved_by"]["run_dir"]["resolved_by"] == "default_output_dir"

    selected = api.select(selection_id="core", mode="radius", space="canonical", center=(0, 0, 0), radius=180)
    assert selected.selected_counts == (N,)
    summary = json.loads((selected.output_dir / "selection_summary.json").read_text())
    assert summary["resolved_by"]["canonical_id"] == {"value": "default", "resolved_by": "only_candidate"}

    report = api.export(domain="ref")  # the only selection is resolved
    assert report["resolved_by"]["selection_id"]["resolved_by"] == "only_candidate"

    state = api.status()
    assert state["bundle_status"] == "completed"
    assert [item["selection_id"] for item in state["selections"]] == ["core"]
    assert api.next_actions()

    captured = capsys.readouterr()
    assert captured.out == "" and captured.err == ""  # the API never prints


def test_api_cancel_publishes_nothing(inputs, tmp_path) -> None:
    ref, mov = inputs
    token = api.CancelToken()
    seen: list[str] = []

    def progress(event: api.ProgressEvent) -> None:
        seen.append(event.kind)
        if event.kind == "stage_done":
            token.cancel()

    with pytest.raises(api.CancelledError) as excinfo:
        api.run(ref, mov, output_dir=tmp_path / "out", k_neighbors=10, no_visualize=True,
                progress=progress, cancel=token)
    assert excinfo.value.code == "cancelled"
    assert not (tmp_path / "out").exists()
    assert not [p for p in tmp_path.iterdir() if p.name.startswith(".out")]  # no staging or failed bundle left


def test_api_errors_are_cryorole_errors(inputs, tmp_path, monkeypatch) -> None:
    ref, mov = inputs
    with pytest.raises(api.CryoroleError, match="Unknown run option"):
        api.run(ref, mov, output_dir=tmp_path / "x", not_an_option=1)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(api.CryoroleError) as excinfo:
        api.status()
    assert excinfo.value.code == "run_dir_unresolved"
    with pytest.raises(api.CryoroleError, match="canonical_id requires space='canonical'"):
        api.select(tmp_path, selection_id="s", canonical_id="x", center=(0, 0, 0), radius=1)


def test_api_warnings_reach_the_cryorole_logger(caplog) -> None:
    from cryorole.logs import warn_user

    logging.getLogger("cryorole").propagate = True  # restored by the conftest fixture
    with caplog.at_level(logging.WARNING, logger="cryorole"):
        warn_user("something to know")
    assert "something to know" in caplog.text
