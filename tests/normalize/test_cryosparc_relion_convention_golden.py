"""P0-5: CryoSPARC pose convention locked against pyem csparc2star.py.

The fixture is a deterministic 2,000-particle subset of the J75/J80 example
(uid + ``alignments3D/pose``) and the STAR files that pyem's real
``csparc2star.py`` writes for it (see scripts/make_cryosparc_convention_fixture.py).
RO computed from the ``.cs`` files must equal RO computed from the converted
STAR files through the RELION path.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cryorole.cli.main import build_parser, run_command
from cryorole.io.readers import read_relion_star
from cryorole.io.writers.landscape_store import read_landscape_npz_arrays
from cryorole.models.policies import ConventionPolicy
from cryorole.normalize.conventions import ConventionResolver

pytestmark = pytest.mark.private_fixtures("j75_j80_subset_2000.npz", "j75_subset_2000_pyem.star", "j80_subset_2000_pyem.star")

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
EULER = ["_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi"]


@pytest.fixture(scope="module")
def subset():
    data = np.load(FIXTURES / "j75_j80_subset_2000.npz")
    return {name: data[name] for name in data.files}


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
        "--k-neighbors", "10", "--no-visualize", *extra,
    ])
    assert run_command(args) == 0


def _star_angles(name: str) -> np.ndarray:
    particles = read_relion_star(FIXTURES / name).particles
    return particles[EULER].to_numpy(dtype=float)


@pytest.mark.parametrize("label,key", [("j75", "ref_pose"), ("j80", "mov_pose")])
def test_cryosparc_bridge_matches_pyem_converted_star_per_particle(subset, label, key) -> None:
    cs_active = ConventionResolver(ConventionPolicy.cryosparc_default()).rotvec_to_active_matrices(subset[key])
    star_active = ConventionResolver(ConventionPolicy.relion_default()).euler_to_active_matrices(
        _star_angles(f"{label}_subset_2000_pyem.star")
    )
    # pyem writes angles with 6 decimals (degrees): agreement is ~1e-7.
    assert np.abs(cs_active - star_active).max() < 1e-5


def test_cryosparc_bridge_is_the_transpose_of_scipy_from_rotvec() -> None:
    pose = np.array([[0.1, -0.2, 0.3], [0.0, 0.0, 8.9]])
    active = ConventionResolver(ConventionPolicy.cryosparc_default()).rotvec_to_active_matrices(pose)
    np.testing.assert_allclose(active, np.swapaxes(Rotation.from_rotvec(pose).as_matrix(), 1, 2))


@pytest.mark.parametrize("backend", ["array_native", "dataframe_compat"])
def test_cryosparc_and_pyem_converted_star_yield_identical_relative_orientation(
    tmp_path, subset, backend
) -> None:
    _write_slim_cs(tmp_path / "ref.cs", subset["uid"], subset["ref_pose"])
    _write_slim_cs(tmp_path / "mov.cs", subset["uid"], subset["mov_pose"])
    _run(tmp_path / "ref.cs", tmp_path / "mov.cs", tmp_path / "cs_run", "--run-backend", backend)
    _run(
        FIXTURES / "j75_subset_2000_pyem.star",
        FIXTURES / "j80_subset_2000_pyem.star",
        tmp_path / "star_run",
        "--run-backend",
        backend,
    )

    cs = read_landscape_npz_arrays(tmp_path / "cs_run" / "data" / "raw_landscape.npz")
    star = read_landscape_npz_arrays(tmp_path / "star_run" / "data" / "raw_landscape.npz")
    cs_order = np.argsort(cs.ref_source_row_id)
    star_order = np.argsort(star.ref_source_row_id)
    assert np.array_equal(cs.ref_source_row_id[cs_order], star.ref_source_row_id[star_order])
    cs_ro = Rotation.from_rotvec(cs.coordinates_analysis[cs_order]).as_matrix()
    star_ro = Rotation.from_rotvec(star.coordinates_analysis[star_order]).as_matrix()
    assert np.abs(cs_ro - star_ro).max() < 1e-5

    manifest = json.loads((tmp_path / "cs_run" / "run_manifest.json").read_text(encoding="utf-8"))
    text = json.dumps(manifest)
    assert "Rotation.from_rotvec(pose).as_matrix().T" in text
    assert '"backprojection"' in text


def test_resolver_rejects_wrong_cryosparc_rule() -> None:
    policy = ConventionPolicy(
        source_software="cryosparc",
        source_euler_sequence="rotvec",
        source_semantics="backprojection",
        degrees=False,
        conversion_rule="active_matrix = scipy Rotation.from_rotvec(pose).as_matrix()",
    )
    with pytest.raises(ValueError, match="pyem-consistent bridge"):
        ConventionResolver(policy)
