"""P0-6: CryoSPARC rotation vectors with norm above pi.

The fixture holds real poses from the J75/J80 example refinements, including
the rows whose ``alignments3D/pose`` norm exceeds pi (up to ~10 rad). These
tests check convention-independent properties; the convention itself is locked
by ``test_cryosparc_relion_convention_golden.py`` (P0-5).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cryorole.cli.main import build_parser, run_command
from cryorole.io.writers.landscape_store import read_landscape_npz_arrays


FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "j75_j80_out_of_range_poses.npz"


def _write_slim_cs(path: Path, uid: np.ndarray, pose: np.ndarray) -> None:
    dtype = np.dtype([("uid", "<u8"), ("alignments3D/pose", "<f4", (3,))])
    values = np.zeros(len(uid), dtype=dtype)
    values["uid"] = uid
    values["alignments3D/pose"] = pose
    with path.open("wb") as handle:
        np.save(handle, values)


def _run(ref: Path, mov: Path, out: Path, *extra: str) -> None:
    args = build_parser().parse_args([
        "run", "--ref", str(ref), "--mov", str(mov), "--output-dir", str(out),
        "--k-neighbors", "5", "--no-visualize", *extra,
    ])
    assert run_command(args) == 0


@pytest.fixture(scope="module")
def fixture_data():
    data = np.load(FIXTURE)
    return {name: data[name] for name in data.files}


def test_fixture_contains_real_out_of_range_poses(fixture_data) -> None:
    ref_norm = np.linalg.norm(fixture_data["ref_pose"].astype(float), axis=1)
    mov_norm = np.linalg.norm(fixture_data["mov_pose"].astype(float), axis=1)
    assert (ref_norm > np.pi).sum() >= 30
    assert (mov_norm > np.pi).sum() >= 30
    assert max(ref_norm.max(), mov_norm.max()) > 2 * np.pi


def _independent_ro_angle(ref_pose: np.ndarray, mov_pose: np.ndarray) -> np.ndarray:
    """Geodesic angle between the two poses; identical for R^T R' and R R'^T."""

    q_ref = Rotation.from_rotvec(ref_pose.astype(float)).as_quat()
    q_mov = Rotation.from_rotvec(mov_pose.astype(float)).as_quat()
    dot = np.clip(np.abs(np.sum(q_ref * q_mov, axis=1)), 0.0, 1.0)
    return 2.0 * np.arccos(dot)


@pytest.mark.parametrize("backend", ["array_native", "dataframe_compat"])
def test_out_of_range_poses_give_valid_ro_on_both_backends(tmp_path, fixture_data, backend) -> None:
    ref, mov = tmp_path / "ref.cs", tmp_path / "mov.cs"
    _write_slim_cs(ref, fixture_data["uid"], fixture_data["ref_pose"])
    _write_slim_cs(mov, fixture_data["uid"], fixture_data["mov_pose"])

    _run(ref, mov, tmp_path / "out", "--run-backend", backend)

    arrays = read_landscape_npz_arrays(tmp_path / "out" / "data" / "raw_landscape.npz")
    assert arrays.n_points == len(fixture_data["uid"])
    order = arrays.ref_source_row_id
    ro_norm = np.linalg.norm(arrays.coordinates_analysis, axis=1)
    assert np.all(ro_norm <= np.pi + 1e-12)
    expected = _independent_ro_angle(fixture_data["ref_pose"][order], fixture_data["mov_pose"][order])
    assert np.allclose(ro_norm, expected, atol=1e-6)


def test_wrapping_raw_poses_into_the_pi_ball_does_not_change_the_landscape(tmp_path, fixture_data) -> None:
    """A rotation vector and its wrapped equivalent are the same rotation."""

    uid = fixture_data["uid"]
    ref_pose, mov_pose = fixture_data["ref_pose"], fixture_data["mov_pose"]
    wrapped_ref = Rotation.from_rotvec(ref_pose.astype(float)).as_rotvec().astype(np.float32)
    wrapped_mov = Rotation.from_rotvec(mov_pose.astype(float)).as_rotvec().astype(np.float32)
    assert np.linalg.norm(wrapped_ref.astype(float), axis=1).max() <= np.pi + 1e-6

    _write_slim_cs(tmp_path / "ref.cs", uid, ref_pose)
    _write_slim_cs(tmp_path / "mov.cs", uid, mov_pose)
    _write_slim_cs(tmp_path / "ref_wrapped.cs", uid, wrapped_ref)
    _write_slim_cs(tmp_path / "mov_wrapped.cs", uid, wrapped_mov)
    _run(tmp_path / "ref.cs", tmp_path / "mov.cs", tmp_path / "raw")
    _run(tmp_path / "ref_wrapped.cs", tmp_path / "mov_wrapped.cs", tmp_path / "wrapped")

    raw = read_landscape_npz_arrays(tmp_path / "raw" / "data" / "raw_landscape.npz")
    wrapped = read_landscape_npz_arrays(tmp_path / "wrapped" / "data" / "raw_landscape.npz")
    assert raw.particle_key.tolist() == wrapped.particle_key.tolist()
    # float32 re-quantization of the wrapped vectors limits agreement to ~1e-6 rad.
    raw_matrices = Rotation.from_rotvec(raw.coordinates_analysis).as_matrix()
    wrapped_matrices = Rotation.from_rotvec(wrapped.coordinates_analysis).as_matrix()
    assert np.allclose(raw_matrices, wrapped_matrices, atol=2e-6)
