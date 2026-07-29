from __future__ import annotations

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cryorole.animation.transforms import (
    canonical_rotation_to_raw,
    conjugate_active_delta,
    pivoted_scene_delta,
    ro_target_to_active_delta,
)
from cryorole.canonicalize.transforms import apply_canonical_transform


def test_canonical_to_raw_round_trip():
    canonical_transform = Rotation.from_euler("x", 37, degrees=True).as_matrix()
    raw = np.asarray([[0.2, -0.3, 0.4]])
    canonical = apply_canonical_transform(raw, canonical_transform)[0]
    rebuilt = canonical_rotation_to_raw(
        Rotation.from_rotvec(canonical),
        canonical_transform,
    )
    assert np.allclose(rebuilt.as_rotvec(), raw[0], atol=1e-12)


def test_active_delta_direction_and_nonidentity_baseline():
    identity = Rotation.identity()
    target = Rotation.from_euler("z", 30, degrees=True)
    delta = ro_target_to_active_delta(target, identity)
    point = delta.apply([1.0, 0.0, 0.0])
    assert np.allclose(point, [np.cos(np.pi / 6), -np.sin(np.pi / 6), 0], atol=1e-12)

    baseline = Rotation.from_euler("x", 20, degrees=True)
    expected = (target * baseline.inv()).inv()
    assert np.allclose(
        ro_target_to_active_delta(target, baseline).as_matrix(),
        expected.as_matrix(),
    )


def test_raw_to_scene_conjugation():
    raw_to_scene = Rotation.from_euler("x", 90, degrees=True).as_matrix()
    raw_delta = Rotation.from_euler("z", 90, degrees=True).as_matrix()
    scene_delta = conjugate_active_delta(raw_delta, raw_to_scene)
    expected = raw_to_scene @ raw_delta @ raw_to_scene.T
    assert np.allclose(scene_delta, expected)


def test_pivoted_scene_delta_invariants():
    pivot = np.asarray([3.0, -2.0, 7.0])
    rotation = Rotation.from_euler("y", 63, degrees=True).as_matrix()
    transform = pivoted_scene_delta(rotation, pivot)
    pivot_h = np.r_[pivot, 1.0]
    assert np.allclose(transform @ pivot_h, pivot_h)

    point = np.asarray([5.0, 1.0, 9.0, 1.0])
    moved = transform @ point
    assert np.isclose(np.linalg.norm(point[:3] - pivot), np.linalg.norm(moved[:3] - pivot))

    identity = pivoted_scene_delta(np.eye(3), pivot)
    assert np.allclose(identity, np.eye(4))
    assert np.isfinite(transform).all()


@pytest.mark.parametrize(
    "rotation,pivot",
    [
        (np.full((3, 3), np.nan), [0, 0, 0]),
        (np.eye(3), [0, np.inf, 0]),
        (np.eye(2), [0, 0, 0]),
    ],
)
def test_transform_helpers_reject_invalid_values(rotation, pivot):
    with pytest.raises(ValueError):
        pivoted_scene_delta(rotation, pivot)
