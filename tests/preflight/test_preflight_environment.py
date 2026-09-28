"""Preflight machine checks: versions, writable output, existing output, memory."""

from __future__ import annotations

import cryorole.preflight.environment as env


def test_writable_new_output_has_no_findings(tmp_path) -> None:
    result = env.check_environment(tmp_path / "new" / "run")
    assert result["errors"] == [] and result["warnings"] == []
    record = result["environment"]
    assert record["output_writable"] is True
    assert record["output_parent_checked"] == str(tmp_path)
    assert record["packages"]["cryorole"]
    assert record["packages"]["numpy"]
    assert not (tmp_path / "new").exists()  # side-effect free


def test_existing_output_warns_about_overwrite(tmp_path) -> None:
    (tmp_path / "run").mkdir()
    result = env.check_environment(tmp_path / "run")
    assert any(w.startswith("[OUTPUT_EXISTS]") for w in result["warnings"])


def test_unwritable_output_is_an_error(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(env.os, "access", lambda path, mode: False)
    result = env.check_environment(tmp_path / "run")
    assert any(e.startswith("[OUTPUT_NOT_WRITABLE]") for e in result["errors"])


def test_low_memory_warns(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(env, "_available_memory_bytes", lambda: (100 * 2**20, "test"))
    result = env.check_environment(tmp_path / "run", estimated_peak_memory_bytes=500 * 2**20)
    assert any(w.startswith("[LOW_MEMORY]") for w in result["warnings"])


def test_preflight_report_carries_environment(tmp_path) -> None:
    import sys
    from pathlib import Path

    import numpy as np

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run"))
    from test_core_productionization import _write_cs

    from cryorole.preflight.service import PreflightRequest, run_preflight

    uids = list(range(100, 160))
    rng = np.random.default_rng(3)
    _write_cs(tmp_path / "ref.cs", uids, rng.normal(0, 0.5, (60, 3)))
    _write_cs(tmp_path / "mov.cs", uids, rng.normal(0, 0.5, (60, 3)))
    (tmp_path / "out").mkdir()
    report = run_preflight(PreflightRequest(ref=tmp_path / "ref.cs", mov=tmp_path / "mov.cs", output_dir=tmp_path / "out")).report
    assert report["environment"]["output_exists"] is True
    assert report["readiness"] == "READY_WITH_WARNINGS"
    assert any("[OUTPUT_EXISTS]" in w for w in report["warnings"])
