"""Diagnostic persistence and scientific invariance across the public workflow."""

import hashlib
import json

import numpy as np
import pytest

from cryorole.cli.main import main


def write_inputs(tmp_path, group_size):
    values = np.zeros(60, dtype=[("uid", "<u8"), ("alignments3D/pose", "<f8", (3,))])
    values["uid"] = np.arange(60)
    ref = tmp_path / "ref.cs"
    mov = tmp_path / "mov.cs"
    with ref.open("wb") as handle:
        np.save(handle, values)
    values["alignments3D/pose"] = np.random.default_rng(7).normal(0, 0.15, (60, 3))
    values["alignments3D/pose"][:group_size] = [0.1, 0.2, 0.3]
    with mov.open("wb") as handle:
        np.save(handle, values)
    return ref, mov


@pytest.mark.parametrize("group_size,severity", [(3, "info"), (10, "warning")])
def test_diagnostic_report_backend_parity_and_unchanged_sld(tmp_path, group_size, severity):
    ref, mov = write_inputs(tmp_path, group_size)
    reports = []
    arrays = []
    for backend in ("array_native", "dataframe_compat"):
        run = tmp_path / backend
        assert main(["run", "--ref", str(ref), "--mov", str(mov), "--output-dir", str(run),
                     "--run-backend", backend, "--k-neighbors", "15", "--no-visualize"]) == 0
        report = json.loads((run / "reports/density_report.json").read_text())["report"]
        reports.append(report)
        with np.load(run / "data/raw_landscape.npz", allow_pickle=False) as data:
            arrays.append({key: data[key] for key in ("coordinates_analysis", "sld_raw", "sld_display", "particle_key")})
        diagnostic = report["ro_coordinate_diagnostics"]
        assert diagnostic["severity"] == severity
        assert diagnostic["coincident_rows"] == group_size
        assert diagnostic["particle_identity_assessed"] is False
        assert ("RO_COORDINATE_CONCENTRATION" in report["warning_codes"]) == (severity == "warning")
        assert not any("many_near_duplicate" in warning for warning in report["warnings"])
        text = (run / "run_report.md").read_text(encoding="utf-8")
        assert "RO coordinate coincidences" in text
        assert "not duplicate particle identity" in text
        manifest = json.loads((run / "run_manifest.json").read_text())
        assert manifest["reports"]["density_report"]["ro_coordinate_diagnostics"] == diagnostic
        assert manifest["active_policies"]["density_policy"]["ro_coincidence_min_group_size"] == 10
        summary = json.loads((run / "run_summary.json").read_text())
        assert summary["ro_coordinate_diagnostics"] == diagnostic
    assert reports[0]["ro_coordinate_diagnostics"] == reports[1]["ro_coordinate_diagnostics"]
    for name in arrays[0]:
        if name == "particle_key":
            np.testing.assert_array_equal(arrays[0][name], arrays[1][name])
        else:
            np.testing.assert_allclose(arrays[0][name], arrays[1][name], rtol=1e-12, atol=1e-14)


def test_run_canonical_visualize_select_export_preserves_source_and_parent(tmp_path, capsys):
    ref, mov = write_inputs(tmp_path, 3)
    run = tmp_path / "run"
    assert main(["run", "--ref", str(ref), "--mov", str(mov), "--output-dir", str(run),
                 "--no-visualize", "--k-neighbors", "15"]) == 0
    paths = [ref, mov, run / "data/raw_landscape.npz", run / "data/raw_landscape.csv"]
    fingerprints = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    assert main(["canonicalize", "--run-dir", str(run), "--no-visualize"]) == 0
    assert main(["visualize", "--run-dir", str(run), "--space", "canonical",
                 "--opacity", "0.4", "--representation", "rotvec", "--max-points", "10"]) == 0
    assert main(["select", "--run-dir", str(run), "--selection-id", "all_rows",
                 "--mode", "threshold", "--sld-min", "0"]) == 0
    assert main(["export", "--run-dir", str(run), "--selection-id", "all_rows"]) == 0
    for domain, path in (("ref", ref), ("mov", mov)):
        np.testing.assert_array_equal(np.load(path, allow_pickle=False), np.load(
            run / "exports/all_rows" / domain / f"selected_{domain}.cs", allow_pickle=False))
    assert fingerprints == {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    stderr = capsys.readouterr().err
    assert "60 particles" in stderr
    assert "Next: cryorole export" in stderr
