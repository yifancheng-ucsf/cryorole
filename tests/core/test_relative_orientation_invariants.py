import numpy as np
import pandas as pd
import pytest
from scipy.spatial.transform import Rotation

from cryorole.core import compute_relative_orientations
from cryorole.core.euler_conventions import resolve_euler_convention
from cryorole.core.representations import derive_rotation_representations
from cryorole.models.match_table import MatchTable
from cryorole.models.pose_table import PoseTable


def relative_orientation_passive(r_ref: np.ndarray, r_mov: np.ndarray) -> np.ndarray:
    return np.linalg.inv(r_ref) @ r_mov


def _pose_table(keys: list[str], matrices: list[np.ndarray], source_rows: list[int]) -> PoseTable:
    return PoseTable(
        pd.DataFrame(
            [
                {
                    "particle_key": key,
                    "domain_name": "domain",
                    "rotation_matrix_active": matrix,
                    "shift_xy": None,
                    "source_type": "synthetic",
                    "source_row_id": source_row,
                    "uid": None,
                    "raw_metadata": {"source_row_id": source_row},
                }
                for key, matrix, source_row in zip(keys, matrices, source_rows)
            ]
        )
    )


def _match_table(keys: list[str], ref_rows: list[int], mov_rows: list[int]) -> MatchTable:
    return MatchTable(
        pd.DataFrame(
            {
                "particle_key": keys,
                "domain_a_row": ref_rows,
                "domain_b_row": mov_rows,
                "match_status": ["matched"] * len(keys),
            }
        )
    )


def test_ro_semantics_are_locked_to_rref_inverse_rmov() -> None:
    r_ref = Rotation.from_euler("ZYZ", [20, 30, 40], degrees=True).as_matrix()
    r_mov = Rotation.from_euler("ZYZ", [-15, 25, 70], degrees=True).as_matrix()
    reference = _pose_table(["p1"], [r_ref], [5])
    moving = _pose_table(["p1"], [r_mov], [9])
    matches = _match_table(["p1"], [5], [9])

    result = compute_relative_orientations(reference, moving, matches)

    np.testing.assert_allclose(result.data.loc[0, "R_ro"], r_ref.T @ r_mov)
    assert result.data.loc[0, "reference_source_row_id"] == 5
    assert result.data.loc[0, "moving_source_row_id"] == 9


def test_ro_swap_inverts_relative_orientation() -> None:
    r_ref = Rotation.from_euler("ZYZ", [20, 30, 40], degrees=True).as_matrix()
    r_mov = Rotation.from_euler("ZYZ", [-15, 25, 70], degrees=True).as_matrix()
    reference = _pose_table(["p1"], [r_ref], [0])
    moving = _pose_table(["p1"], [r_mov], [0])
    matches = _match_table(["p1"], [0], [0])

    ro = compute_relative_orientations(reference, moving, matches).data.loc[0, "R_ro"]
    swapped = compute_relative_orientations(moving, reference, matches).data.loc[0, "R_ro"]

    np.testing.assert_allclose(swapped, np.linalg.inv(ro))


def test_common_global_frame_change_preserves_rotation_angle() -> None:
    r_ref = Rotation.from_euler("ZYZ", [20, 30, 40], degrees=True).as_matrix()
    r_mov = Rotation.from_euler("ZYZ", [-15, 25, 70], degrees=True).as_matrix()
    global_frame = Rotation.from_euler("ZYZ", [75, -10, 5], degrees=True).as_matrix()
    matches = _match_table(["p1"], [0], [0])

    ro = compute_relative_orientations(
        _pose_table(["p1"], [r_ref], [0]),
        _pose_table(["p1"], [r_mov], [0]),
        matches,
    ).data.loc[0, "R_ro"]
    transformed = compute_relative_orientations(
        _pose_table(["p1"], [r_ref @ global_frame], [0]),
        _pose_table(["p1"], [r_mov @ global_frame], [0]),
        matches,
    ).data.loc[0, "R_ro"]

    np.testing.assert_allclose(transformed, global_frame.T @ ro @ global_frame)
    np.testing.assert_allclose(
        Rotation.from_matrix(transformed).magnitude(),
        Rotation.from_matrix(ro).magnitude(),
    )


def test_rotation_representation_roundtrip_stability() -> None:
    ro = Rotation.from_euler("ZYZ", [35, 40, -20], degrees=True).as_matrix()
    derived = derive_rotation_representations(ro)
    rotvec = derived["rotvec_ro"]
    roundtrip = Rotation.from_rotvec(rotvec).as_matrix()

    np.testing.assert_allclose(roundtrip, ro)


def test_euler_convention_helper_maps_public_names_to_scipy_sequences() -> None:
    extrinsic = resolve_euler_convention("extrinsic_zyx", source="cli_default")
    intrinsic = resolve_euler_convention("intrinsic_zyx", source="cli_override")

    assert extrinsic.scipy_euler_sequence == "zyx"
    assert extrinsic.metadata()["euler_convention"] == "extrinsic_zyx"
    assert intrinsic.scipy_euler_sequence == "ZYX"
    assert intrinsic.metadata()["euler_convention_source"] == "cli_override"


def test_ro_computation_derives_outputs_from_internal_active_matrix() -> None:
    r_ref = np.eye(3)
    r_mov = Rotation.from_euler("ZYZ", [10, 20, 30], degrees=True).as_matrix()
    result = compute_relative_orientations(
        _pose_table(["p1"], [r_ref], [0]),
        _pose_table(["p1"], [r_mov], [0]),
        _match_table(["p1"], [0], [0]),
    )

    row = result.data.loc[0]
    np.testing.assert_allclose(row["R_ro"], r_mov)
    np.testing.assert_allclose(row["quat_ro"], Rotation.from_matrix(r_mov).as_quat())
    np.testing.assert_allclose(row["rotvec_ro"], Rotation.from_matrix(r_mov).as_rotvec())
    np.testing.assert_allclose(row["euler_ro"], Rotation.from_matrix(r_mov).as_euler("zyx", degrees=True))
    assert row["angle_ro"] == pytest.approx(np.degrees(Rotation.from_matrix(r_mov).magnitude()))


def test_ro_computation_can_preserve_intrinsic_zyx_compatibility() -> None:
    r_ref = np.eye(3)
    r_mov = Rotation.from_euler("ZYZ", [10, 20, 30], degrees=True).as_matrix()
    result = compute_relative_orientations(
        _pose_table(["p1"], [r_ref], [0]),
        _pose_table(["p1"], [r_mov], [0]),
        _match_table(["p1"], [0], [0]),
        euler_sequence="ZYX",
    )

    row = result.data.loc[0]
    np.testing.assert_allclose(
        row["euler_ro"],
        Rotation.from_matrix(r_mov).as_euler("ZYX", degrees=True),
    )


def test_ro_computation_accepts_only_normalized_matched_inputs() -> None:
    pose = _pose_table(["p1"], [np.eye(3)], [0])
    match = _match_table(["p1"], [0], [0])

    with pytest.raises(TypeError, match="normalized PoseTable"):
        compute_relative_orientations(pd.DataFrame(), pose, match)

    with pytest.raises(TypeError, match="MatchTable"):
        compute_relative_orientations(pose, pose, pd.DataFrame())

    bad_match = MatchTable(
        pd.DataFrame(
            {
                "particle_key": ["p1"],
                "domain_a_row": [0],
                "domain_b_row": [0],
                "match_status": ["unmatched"],
            }
        )
    )
    with pytest.raises(ValueError, match="matched rows only"):
        compute_relative_orientations(pose, pose, bad_match)
