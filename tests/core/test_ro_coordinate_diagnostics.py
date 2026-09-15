"""Quantized-coordinate concentration is independent of particle identity."""

from dataclasses import replace

import numpy as np
import pytest

from cryorole.core.density import _ro_coordinate_diagnostic, compute_landscape_density_arrays
from cryorole.models.policies import DensityPolicy


def coordinates(group_sizes, total):
    sizes = [*group_sizes, *([1] * (total - sum(group_sizes)))]
    x = np.repeat(np.arange(len(sizes), dtype=float) * 1e-4 + 0.1, sizes)
    return np.column_stack((x, np.full(total, 0.2), np.full(total, 0.3)))


@pytest.mark.parametrize("groups,total,concentrated,severity", [
    ([], 100, 0, "info"),
    ([2] * 100 + [3] * 100, 1000, 0, "info"),
    ([9] * 20, 1000, 0, "info"),
    ([10], 1001, 10, "info"),
    ([10], 1000, 10, "warning"),
    ([99], 10000, 99, "info"),
    ([100], 10001, 100, "warning"),
    ([10] * 10, 20000, 100, "warning"),
    ([], 0, 0, "info"),
])
def test_approved_concentration_boundaries(groups, total, concentrated, severity):
    result = _ro_coordinate_diagnostic(coordinates(groups, total), DensityPolicy())
    assert result["coincident_rows"] == sum(groups)
    assert result["largest_group"] == max(groups, default=0)
    assert result["concentrated_rows"] == concentrated
    assert result["severity"] == severity
    assert result["particle_identity_assessed"] is False
    assert "not duplicate particle identity" in result["message"]
    assert result["policy"] == {"min_group_size": 10, "min_rows": 100, "min_fraction": 0.01}


def test_diagnostic_thresholds_are_explicit_and_disable_preserves_zero_stats():
    policy = replace(DensityPolicy(), ro_coincidence_min_group_size=3)
    result = _ro_coordinate_diagnostic(coordinates([3] * 40, 1000), policy)
    assert result["severity"] == "warning"
    assert result["policy"]["min_group_size"] == 3
    result = _ro_coordinate_diagnostic(coordinates([100], 100),
                                      replace(policy, near_duplicate_coordinate_tolerance_rad=0))
    assert result["coincident_rows"] == 0
    assert result["severity"] == "info"


def test_severity_policy_does_not_change_any_landscape_arrays():
    coords = coordinates([10] * 10, 200)
    keys = np.asarray([f"p{i}" for i in range(len(coords))])
    results = [compute_landscape_density_arrays(keys, coords, np.arange(len(coords)), np.arange(len(coords)), policy=policy) for policy in (
        DensityPolicy(), replace(DensityPolicy(), ro_coincidence_min_group_size=1000),
    )]
    assert results[0][1].ro_coordinate_diagnostics["severity"] == "warning"
    assert results[1][1].ro_coordinate_diagnostics["severity"] == "info"
    for field, value in vars(results[0][0]).items():
        if isinstance(value, np.ndarray):
            np.testing.assert_array_equal(value, getattr(results[1][0], field))


@pytest.mark.parametrize("override", [
    {"ro_coincidence_min_group_size": 1}, {"ro_coincidence_min_group_size": 2.5},
    {"ro_coincidence_min_rows": 0}, {"ro_coincidence_min_rows": True},
    {"ro_coincidence_min_fraction": 0}, {"ro_coincidence_min_fraction": float("nan")},
])
def test_invalid_diagnostic_policies_fail(override):
    with pytest.raises(ValueError, match="ro_coincidence"):
        compute_landscape_density_arrays(np.asarray(["a", "b"]), coordinates([], 2), np.arange(2), np.arange(2),
                                         policy=replace(DensityPolicy(), **override))
