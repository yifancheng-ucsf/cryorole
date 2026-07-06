import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
import pytest

from cryorole.normalize.pose_normalizer import PoseNormalizer


def test_cryosparc_native_rotvec_pose_converts_to_active_matrix() -> None:
    rotvec = np.array([0.1, -0.2, 0.3])
    particles = pd.DataFrame(
        {
            "uid": [123],
            "alignments3D/pose": [rotvec],
            "alignments3D/shift": [np.array([1.5, -2.0])],
        }
    )

    pose_table = PoseNormalizer().normalize_cryosparc(particles, domain_name="domain-a")

    expected = Rotation.from_rotvec(rotvec).as_matrix()
    np.testing.assert_allclose(
        pose_table.data.loc[0, "rotation_matrix_active"],
        expected,
    )


def test_cryosparc_shift_and_uid_are_preserved() -> None:
    shift = np.array([2.25, -4.5])
    particles = pd.DataFrame(
        {
            "uid": [456],
            "alignments3D/pose": [np.zeros(3)],
            "alignments3D/shift": [shift],
        }
    )

    pose_table = PoseNormalizer().normalize_cryosparc(particles, domain_name="domain-a")

    assert pose_table.data.loc[0, "uid"] == 456
    np.testing.assert_allclose(pose_table.data.loc[0, "shift_xy"], shift)
    assert pose_table.data.loc[0, "raw_metadata"]["uid"] == 456
    assert pose_table.data.loc[0, "raw_metadata"]["_cryorole_source_row_id"] == 0


def test_cryosparc_precomputed_matrix_without_native_pose_is_rejected() -> None:
    matrix = Rotation.from_euler("ZYX", [10, 20, 30], degrees=True).as_matrix()
    particles = pd.DataFrame(
        {
            "uid": [789],
            "precomputed_matrix": [matrix],
        }
    )

    with pytest.raises(ValueError, match="alignments3D/pose"):
        PoseNormalizer().normalize_cryosparc(particles, domain_name="domain-a")


def test_cryosparc_missing_native_pose_is_rejected() -> None:
    particles = pd.DataFrame({"uid": [1]})

    with pytest.raises(ValueError, match="alignments3D/pose"):
        PoseNormalizer().normalize_cryosparc(particles, domain_name="domain-a")


def test_cryosparc_missing_uid_is_rejected() -> None:
    particles = pd.DataFrame({"alignments3D/pose": [np.zeros(3)]})

    with pytest.raises(ValueError, match="uid"):
        PoseNormalizer().normalize_cryosparc(particles, domain_name="domain-a")


def test_cryosparc_pose_must_be_length_three() -> None:
    particles = pd.DataFrame(
        {
            "uid": [123],
            "alignments3D/pose": [np.array([0.1, 0.2])],
        }
    )

    with pytest.raises(ValueError, match="length-3 rotation vector"):
        PoseNormalizer().normalize_cryosparc(particles, domain_name="domain-a")


def test_cryosparc_shift_must_be_length_two_when_present() -> None:
    particles = pd.DataFrame(
        {
            "uid": [123],
            "alignments3D/pose": [np.zeros(3)],
            "alignments3D/shift": [np.array([1.0, 2.0, 3.0])],
        }
    )

    with pytest.raises(ValueError, match="length-2 shift vector"):
        PoseNormalizer().normalize_cryosparc(particles, domain_name="domain-a")
