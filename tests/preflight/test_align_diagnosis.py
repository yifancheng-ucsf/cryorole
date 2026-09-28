"""A4: preflight suggestions for `cryorole align` (never binding)."""

from __future__ import annotations

import pytest

import shutil
from pathlib import Path

from cryorole.preflight.align_diagnosis import diagnose_star_matching
from cryorole.preflight.service import PreflightRequest, run_preflight

pytestmark = pytest.mark.private_fixtures("relion_recenter")

F = Path(__file__).resolve().parents[1] / "fixtures" / "relion_recenter"


def _star(path: Path, names: list[str], column: str = "_rlnImageName") -> None:
    rows = [f"{name} 10 20 30" for name in names]
    path.write_text("\n".join([
        "data_optics", "loop_", "_rlnOpticsGroup #1", "_rlnImagePixelSize #2", "1 1.0", "",
        "data_particles", "loop_", f"{column} #1", "_rlnAngleRot #2", "_rlnAngleTilt #3", "_rlnAnglePsi #4",
        *rows, "",
    ]), encoding="utf-8")


def test_preflight_suggests_key_pair_for_subtracted_refinement(tmp_path) -> None:
    result = run_preflight(PreflightRequest(ref=F / "job022_run_data_subset.star",
                                            mov=F / "job043_coords_exchange_job042_subset.star"))
    diagnosis = result.report["align_diagnosis"]
    assert result.report["readiness"] == "BLOCKED"  # suggestions never change the verdict
    assert diagnosis["recommended"] == "_rlnImageName = _rlnImageOriginalName"
    assert "--key-pair _rlnImageName=_rlnImageOriginalName" in result.report["recommended_next_command"]
    assert diagnosis["sampling"]["side"] == "ref"


def test_one_sided_sampling_estimates_overlap_without_bias(tmp_path) -> None:
    n = 4000
    _star(tmp_path / "ref.star", [f"{i}@a.mrcs" for i in range(n)])
    _star(tmp_path / "mov.star", [f"{i}@a.mrcs" for i in range(0, 2 * n, 2)])  # half of ref present
    diagnosis = diagnose_star_matching(tmp_path / "ref.star", tmp_path / "mov.star", sample_size=400, seed=1)
    exact = next(c for c in diagnosis["candidates"] if c["name"] == "exact _rlnImageName")
    assert abs(exact["ref_overlap_estimate"] - 0.5) < 0.08
    assert exact["sampled_ref_rows"] == 400


def test_basename_candidate_reports_stack_merge(tmp_path) -> None:
    _star(tmp_path / "ref.star", ["1@/proj/A/stack.mrcs", "1@/proj/B/stack.mrcs"])
    _star(tmp_path / "mov.star", ["1@/other/A/stack.mrcs", "1@/other/B/stack.mrcs"])
    diagnosis = diagnose_star_matching(tmp_path / "ref.star", tmp_path / "mov.star")
    basename = next(c for c in diagnosis["candidates"] if c["name"].endswith("basename"))
    assert basename["stack_merge"] == 2
    assert "merged" in basename["warning"]
    suffix = next(c for c in diagnosis["candidates"] if c["name"].endswith("suffix:2"))
    assert suffix["estimated_matched"] == 2 and not suffix.get("stack_merge")
    assert diagnosis["recommended"] == suffix["name"]


def test_recentered_exact_candidate_from_extraction_note(tmp_path) -> None:
    job = tmp_path / "Extract" / "job055"
    job.mkdir(parents=True)
    shutil.copy(F / "job055_particles_subset.star", job / "particles.star")
    shutil.copy(F / "Extract_job055_note.txt", job / "note.txt")
    diagnosis = diagnose_star_matching(F / "job043_coords_exchange_job042_subset.star", job / "particles.star")
    exact = next(c for c in diagnosis["candidates"] if c["name"].startswith("recentered-exact"))
    assert exact["recenter_shift_px"] == [-9.0, 22.0, -102.0]
    assert exact["verified_fraction"] == 1.0
    assert "--coordinate-match recentered-exact --recenter-shift -9 22 -102" in exact["command"]
    # A name link earlier in the ladder (both files keep the job018 original name) is preferred when unambiguous.
    assert diagnosis["recommended"] == "_rlnImageOriginalName = _rlnImageOriginalName"


def test_stale_subtraction_warning_from_subtract_note(tmp_path) -> None:
    job = tmp_path / "Subtract" / "job041"
    job.mkdir(parents=True)
    shutil.copy(F / "job041_particles_subtracted_subset.star", job / "particles_subtracted.star")
    shutil.copy(F / "Subtract_job041_note.txt", job / "note.txt")
    result = run_preflight(PreflightRequest(ref=F / "job022_run_data_subset.star", mov=job / "particles_subtracted.star"))
    warnings = [w for w in result.report["warnings"] if "STALE_SUBTRACTION_COORDINATES" in w]
    assert len(warnings) == 1
    assert "--fix-subtract-coordinates" in warnings[0]
