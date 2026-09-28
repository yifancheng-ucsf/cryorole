"""Section 3.2 input-sanity warnings in preflight, run_summary.json and run_report.md."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from cryorole.cli.main import build_parser, run_command
from cryorole.preflight.service import PreflightRequest, run_preflight

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_core_productionization import _write_cs  # noqa: E402

N = 400
UIDS = list(range(5000, 5000 + N))


def _run(ref: Path, mov: Path, out: Path) -> int:
    args = build_parser().parse_args([
        "run", "--ref", str(ref), "--mov", str(mov), "--output-dir", str(out),
        "--k-neighbors", "10", "--no-visualize",
    ])
    return run_command(args)


def _ref_rv() -> np.ndarray:
    return np.random.default_rng(7).normal(0.0, 0.6, (N, 3))


def _perturbed(ref_rv: np.ndarray, sigma_rad: float, *, identical_rows: int = 0) -> np.ndarray:
    from scipy.spatial.transform import Rotation

    rng = np.random.default_rng(8)
    delta = Rotation.from_rotvec(rng.normal(0.0, sigma_rad, (N, 3)))
    mov = (Rotation.from_rotvec(ref_rv) * delta).as_rotvec()
    mov[:identical_rows] = ref_rv[:identical_rows]
    return mov


def _summary(out: Path) -> dict:
    return json.loads((out / "run_summary.json").read_text(encoding="utf-8"))


def test_preflight_blocks_same_file_before_matching(tmp_path) -> None:
    ref = tmp_path / "ref.cs"
    _write_cs(ref, UIDS, _ref_rv())
    result = run_preflight(PreflightRequest(ref=ref, mov=ref))
    assert result.report["readiness"] == "BLOCKED"
    assert result.report["input_sanity"]["same_file"]["code"] == "SAME_INPUT_FILE"
    assert any("same file" in error for error in result.report["errors"])
    assert result.array_preflight is None


def test_preflight_blocks_byte_identical_copy(tmp_path) -> None:
    ref, mov = tmp_path / "ref.cs", tmp_path / "copy.cs"
    _write_cs(ref, UIDS, _ref_rv())
    mov.write_bytes(ref.read_bytes())
    result = run_preflight(PreflightRequest(ref=ref, mov=mov))
    assert result.report["readiness"] == "BLOCKED"
    assert result.report["input_sanity"]["same_file"]["reason"] == "same_sha256"


def test_run_with_same_file_stops_with_the_sanity_explanation(tmp_path) -> None:
    ref = tmp_path / "ref.cs"
    _write_cs(ref, UIDS, _ref_rv())
    with pytest.raises(ValueError, match=r"\[SAME_INPUT_FILE\].*\[IDENTICAL_POSES\].*SLD\) is undefined"):
        _run(ref, ref, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_run_with_mostly_identical_poses_completes_with_strong_warning(tmp_path) -> None:
    ref_rv = _ref_rv()
    _write_cs(tmp_path / "ref.cs", UIDS, ref_rv)
    _write_cs(tmp_path / "mov.cs", UIDS, _perturbed(ref_rv, 0.05, identical_rows=N - 2))

    assert _run(tmp_path / "ref.cs", tmp_path / "mov.cs", tmp_path / "out") == 0

    sanity = _summary(tmp_path / "out")["input_sanity"]
    assert sanity["level"] == "strong_warning"
    assert [f["code"] for f in sanity["findings"]] == ["IDENTICAL_POSES"]
    report = (tmp_path / "out" / "run_report.md").read_text(encoding="utf-8")
    assert "## Input sanity" in report
    assert "**STRONG WARNING** [IDENTICAL_POSES]" in report


def test_run_with_nearly_identical_orientations_warns(tmp_path) -> None:
    ref_rv = _ref_rv()
    _write_cs(tmp_path / "ref.cs", UIDS, ref_rv)
    _write_cs(tmp_path / "mov.cs", UIDS, _perturbed(ref_rv, np.radians(0.25)))

    assert _run(tmp_path / "ref.cs", tmp_path / "mov.cs", tmp_path / "out") == 0

    sanity = _summary(tmp_path / "out")["input_sanity"]
    assert sanity["level"] == "warning"
    assert sanity["findings"][0]["code"] == "NEARLY_IDENTICAL_ORIENTATIONS"
    assert "**WARNING** [NEARLY_IDENTICAL_ORIENTATIONS]" in (tmp_path / "out" / "run_report.md").read_text(
        encoding="utf-8"
    )


def test_normal_run_reports_ro_angle_summary_without_warnings(tmp_path) -> None:
    ref_rv = _ref_rv()
    _write_cs(tmp_path / "ref.cs", UIDS, ref_rv)
    _write_cs(tmp_path / "mov.cs", UIDS, _perturbed(ref_rv, np.radians(8.0)))

    assert _run(tmp_path / "ref.cs", tmp_path / "mov.cs", tmp_path / "out") == 0

    sanity = _summary(tmp_path / "out")["input_sanity"]
    assert sanity["level"] == "none"
    assert sanity["ro_angle_summary"]["n"] == N
    report = (tmp_path / "out" / "run_report.md").read_text(encoding="utf-8")
    assert sanity["ro_angle_summary"]["message"] in report
    assert "No input-sanity warnings" in report


def test_run_with_more_than_k_identical_ros_no_longer_crashes(tmp_path) -> None:
    """P0-3 end to end: 60 identical ROs with k=10 on the default array-native backend."""

    ref_rv = _ref_rv()
    mov = _perturbed(ref_rv, 0.3)
    mov[:60] = ref_rv[:60]
    _write_cs(tmp_path / "ref.cs", UIDS, ref_rv)
    _write_cs(tmp_path / "mov.cs", UIDS, mov)

    assert _run(tmp_path / "ref.cs", tmp_path / "mov.cs", tmp_path / "out") == 0
    density = json.loads((tmp_path / "out" / "reports" / "density_report.json").read_text(encoding="utf-8"))
    assert json.dumps(density).count("n_inf_sld_unfloored") >= 1
    assert "unfloored SLD (`sld_unfloored`) is +inf" in (tmp_path / "out" / "run_report.md").read_text(
        encoding="utf-8"
    )


def test_identical_poses_report_inf_sld_when_ro_product_is_not_exactly_symmetric(tmp_path, monkeypatch) -> None:
    """Regression for the macOS CI failure.

    Some BLAS builds (the macOS NumPy wheels) round R_ref^T @ R_mov differently from Linux,
    so identical poses give ROs ~1e-17 rad from the identity instead of exactly the identity.
    Simulate that asymmetric rounding and require the same report as on Linux.
    """

    import cryorole.workflows.run_pipeline as run_pipeline

    class _AsymmetricRoundingNumpy:
        def __getattr__(self, name):
            return getattr(np, name)

        @staticmethod
        def matmul(a, b, *args, **kwargs):
            out = np.matmul(a, b, *args, **kwargs)
            if isinstance(out, np.ndarray) and out.ndim == 3 and out.shape[1:] == (3, 3):
                out = out.copy()
                out[:, 0, 1] += np.random.default_rng(0).uniform(-1.0, 1.0, len(out)) * 2.2e-17
            return out

    monkeypatch.setattr(run_pipeline, "np", _AsymmetricRoundingNumpy())
    ref_rv = _ref_rv()
    mov = _perturbed(ref_rv, 0.3)
    mov[:60] = ref_rv[:60]
    _write_cs(tmp_path / "ref.cs", UIDS, ref_rv)
    _write_cs(tmp_path / "mov.cs", UIDS, mov)

    assert _run(tmp_path / "ref.cs", tmp_path / "mov.cs", tmp_path / "out") == 0
    density = json.loads((tmp_path / "out" / "reports" / "density_report.json").read_text(encoding="utf-8"))
    assert '"n_inf_sld_unfloored": 60' in json.dumps(density)
    assert "unfloored SLD (`sld_unfloored`) is +inf" in (tmp_path / "out" / "run_report.md").read_text(
        encoding="utf-8"
    )
