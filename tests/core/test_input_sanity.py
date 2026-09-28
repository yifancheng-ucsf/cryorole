"""Section 3.2 input-sanity heuristics."""

from __future__ import annotations

import numpy as np
import pytest

from cryorole.core.input_sanity import (
    InputSanityPolicy,
    assess_ro_angles,
    same_input_file,
    sanity_warning_lines,
    summarize_ro_angles,
)


def test_summary_is_always_reported_with_percentiles() -> None:
    angles = np.radians(np.linspace(0.0, 40.0, 1001))
    summary = summarize_ro_angles(angles)
    assert summary["n"] == 1001
    assert summary["median_deg"] == pytest.approx(20.0)
    assert summary["message"].startswith("RO angle: median 20.0°, 90% 36.0°, 99% 39.6° (n = 1,001)")


def test_j75_j80_like_distribution_triggers_nothing() -> None:
    # Percentiles 1/10/50/90/99 of the J75/J80 example: 1.0/2.9/8.8/20.7/36.2 degrees.
    rng = np.random.default_rng(0)
    degrees = np.interp(rng.uniform(0, 100, 200_000), [0, 1, 10, 50, 90, 99, 100],
                        [0.1, 1.0, 2.9, 8.8, 20.7, 36.2, 60.0])
    report = assess_ro_angles(np.radians(degrees))
    assert report["level"] == "none"
    assert report["findings"] == []
    assert report["ro_angle_summary"]["median_deg"] == pytest.approx(8.8, abs=0.2)


def test_identical_poses_are_a_strong_warning() -> None:
    angles = np.zeros(1000)
    angles[:5] = 0.2  # 99.5 % identical
    report = assess_ro_angles(angles)
    assert report["level"] == "strong_warning"
    assert [f["code"] for f in report["findings"]] == ["IDENTICAL_POSES"]
    assert "identical for 99.5% of particles" in report["findings"][0]["message"]


def test_nearly_identical_orientations_are_a_warning() -> None:
    rng = np.random.default_rng(1)
    angles = np.radians(np.abs(rng.normal(0.0, 0.5, 10_000)))
    report = assess_ro_angles(angles)
    assert report["level"] == "warning"
    finding = report["findings"][0]
    assert finding["code"] == "NEARLY_IDENTICAL_ORIENTATIONS"
    assert "Check the masks" in finding["message"]


@pytest.mark.parametrize("median_deg,p99_deg", [(0.8, 2.5), (1.2, 1.5)])
def test_near_identical_rule_needs_both_conditions(median_deg, p99_deg) -> None:
    rng = np.random.default_rng(2)
    degrees = np.interp(rng.uniform(0, 100, 20_000), [0, 50, 99, 100], [0.0, median_deg, p99_deg, p99_deg])
    assert assess_ro_angles(np.radians(degrees))["level"] == "none"


def test_policy_thresholds_are_recorded_and_validated() -> None:
    report = assess_ro_angles(np.radians([5.0, 6.0]), policy=InputSanityPolicy(near_identical_median_deg=10.0,
                                                                                near_identical_p99_deg=10.0))
    assert report["level"] == "warning"
    assert report["policy"]["heuristic"] is True
    assert report["policy"]["near_identical_median_deg"] == 10.0
    with pytest.raises(ValueError):
        InputSanityPolicy(identical_min_fraction=0.0)


def test_same_file_by_path_or_hash(tmp_path) -> None:
    path = tmp_path / "a.cs"
    path.write_bytes(b"x")
    by_path = same_input_file({}, {}, ref_path=path, mov_path=tmp_path / "." / "a.cs")
    assert by_path["reason"] == "same_resolved_path"
    by_hash = same_input_file(
        {"resolved_path": "/a", "sha256": "f" * 64}, {"resolved_path": "/b", "sha256": "f" * 64}
    )
    assert by_hash["reason"] == "same_sha256"
    assert "identical content" in by_hash["message"]
    assert same_input_file({"resolved_path": "/a", "sha256": "1"}, {"resolved_path": "/b", "sha256": "2"}) is None


def test_same_file_finding_joins_the_report_as_strong_warning() -> None:
    finding = same_input_file({"resolved_path": "/a"}, {"resolved_path": "/a"})
    report = assess_ro_angles(np.radians([10.0, 20.0]), same_file=finding)
    assert report["level"] == "strong_warning"
    assert sanity_warning_lines(report) == (
        "[SAME_INPUT_FILE] The reference and moving inputs are the same file. Select the two different refinements.",
    )
