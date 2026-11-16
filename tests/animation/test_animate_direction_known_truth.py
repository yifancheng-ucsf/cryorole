"""Known-truth regression for the direction in which ``cryorole animate`` moves the map.

Physical setup: a reference body and an asymmetric moving body (point "blobs") in one
box-centred frame; the moving body is turned by a known ``Q`` about an off-centre pivot.
Per-particle poses are derived from that geometry with RELION's own
``Euler_angles2matrix`` (``x_image = A @ x_map``, see ``cryorole.align.relion_geometry``),
not from cryoROLE's normalization. ``cryorole run`` then gives ``RO``, and
``cryorole animate`` (script-only) must place the moving body at ``Q`` about the pivot.

The same geometry was checked end to end on 2026-10-12 with relion_project /
relion_reconstruct (moving-pose reconstruction r = 0.9999 for ``A_ref @ Q``) and with the
generated ChimeraX renderer resampled onto the reference grid.
"""

from __future__ import annotations

import csv
import json

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cryorole.align.relion_geometry import relion_euler_matrices
from cryorole.animation.transforms import (
    RO_DENSITY_DIRECTION_ACTIVE,
    RO_DENSITY_DIRECTION_LEGACY_INVERSE,
    pivoted_scene_delta,
    ro_target_to_active_delta,
)
from cryorole.cli.main import main

APIX = 3.0
BOX = 64
SCENE_OFFSET = np.full(3, BOX // 2 * APIX)  # ChimeraX scene of a map with origin 0
PIVOT = np.array([18.0, -12.0, 9.0])  # box-centred Angstrom, off-centre
Q = Rotation.from_rotvec(np.radians(25.0) * np.array([0.48, 0.80, 0.36]) / np.linalg.norm([0.48, 0.80, 0.36]))
MOVING_CENTRES = np.array(
    [[0, 0, 18], [27, -9, 30], [-21, 15, 39], [6, 30, 48], [-12, -27, 27], [15, 18, 60]],
    dtype=float,
)


def _kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation M and shift t minimizing |M source + t - target|."""

    sc, tc = source.mean(0), target.mean(0)
    u, _, vt = np.linalg.svd((target - tc).T @ (source - sc))
    d = np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt))])
    rotation = u @ d @ vt
    return rotation, tc - rotation @ sc


def _relion_angles(matrix: np.ndarray) -> np.ndarray:
    angles = Rotation.from_matrix(matrix.T).as_euler("ZYZ", degrees=True)
    # Independent check against RELION's own trigonometric Euler_angles2matrix.
    assert np.allclose(relion_euler_matrices(*angles), matrix, atol=1e-9)
    return angles


def _physical_poses(n_particles: int = 200, seed: int = 11):
    """Per-particle truth: the moving body turned by Q * E_i (E_i ~ 2 degree jitter) about the pivot."""

    rng = np.random.default_rng(seed)
    ref_poses = Rotation.random(n_particles, random_state=rng).as_matrix()
    jitter = Rotation.from_rotvec(rng.normal(scale=np.radians(2.0), size=(n_particles, 3)))
    relatives = Q * jitter
    ref_angles, mov_angles, mov_origins = [], [], []
    for a_ref, relative in zip(ref_poses, relatives):
        a_ref = relion_euler_matrices(*_relion_angles(a_ref))
        moved = relative.apply(MOVING_CENTRES - PIVOT) + PIVOT  # truth, reference frame
        image_centres = moved @ a_ref.T  # where the moving blobs land in this image
        a_mov, shift = _kabsch(MOVING_CENTRES, image_centres)
        ref_angles.append(_relion_angles(a_ref))
        mov_angles.append(_relion_angles(a_mov))
        mov_origins.append(-shift[:2])  # RELION origin re-centres the particle
    return np.asarray(ref_angles), np.asarray(mov_angles), np.asarray(mov_origins), relatives


def _write_star(path, angles, origins):
    lines = [
        "",
        "data_optics",
        "",
        "loop_",
        "_rlnOpticsGroupName #1",
        "_rlnOpticsGroup #2",
        "_rlnImagePixelSize #3",
        "_rlnImageSize #4",
        "_rlnImageDimensionality #5",
        f"opticsGroup1 1 {APIX} {BOX} 2",
        "",
        "data_particles",
        "",
        "loop_",
        "_rlnImageName #1",
        "_rlnAngleRot #2",
        "_rlnAngleTilt #3",
        "_rlnAnglePsi #4",
        "_rlnOriginXAngst #5",
        "_rlnOriginYAngst #6",
        "_rlnOpticsGroup #7",
    ]
    for index, (angle, origin) in enumerate(zip(angles, origins)):
        lines.append(
            f"{index + 1:06d}@particles.mrcs {angle[0]:.9f} {angle[1]:.9f} {angle[2]:.9f} "
            f"{origin[0]:.6f} {origin[1]:.6f} 1"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_cs(path, angles, origins):
    matrices = relion_euler_matrices(angles[:, 0], angles[:, 1], angles[:, 2])
    records = np.zeros(
        len(angles),
        dtype=[("uid", "<u8"), ("alignments3D/pose", "<f4", (3,)), ("alignments3D/shift", "<f4", (2,))],
    )
    records["uid"] = 1000 + np.arange(len(angles))
    # CryoSPARC pose = back-projection rotation: pyem expmap(pose) = A, i.e. rotvec(A^T).
    records["alignments3D/pose"] = Rotation.from_matrix(np.swapaxes(matrices, 1, 2)).as_rotvec()
    records["alignments3D/shift"] = origins / APIX
    with path.open("wb") as handle:
        np.save(handle, records)


def _path_csv(path, target: Rotation):
    rv = target.as_rotvec()
    path.write_text(
        "label,rv_x_rad,rv_y_rad,rv_z_rad\nidentity,0,0,0\n"
        f"target,{rv[0]:.12f},{rv[1]:.12f},{rv[2]:.12f}\n",
        encoding="utf-8",
    )


def _final_scene_transform(output_dir) -> np.ndarray:
    with (output_dir / "chimerax" / "frame_transforms.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    transform = np.eye(4)
    transform[:3, :] = [[float(rows[-1][f"m{r}{c}"]) for c in range(4)] for r in range(3)]
    return transform


def _animate(tmp_path, run_dir, name, extra=()):
    session = tmp_path / "scene.cxs"
    session.write_bytes(b"session")
    path_csv = tmp_path / f"{name}_path.csv"
    _path_csv(path_csv, Q)
    output_dir = tmp_path / name
    pivot_scene = PIVOT + SCENE_OFFSET
    code = main(
        [
            "animate", "--run-dir", str(run_dir), "--path-csv", str(path_csv), "--path-space", "rv",
            "--chimerax-session", str(session), "--reference-model-id", "#1", "--moving-model-id", "#2",
            "--pivot", *(str(v) for v in pivot_scene), "--baseline-ro", "identity", "--map-frame", "raw",
            "--frames-per-segment", "3", "--canvas-width", "300", "--canvas-height", "150",
            "--structure-width", "60", "--structure-height", "40", "--render-mode", "script-only",
            "--output-dir", str(output_dir), *extra,
        ]
    )
    assert code == 0
    return output_dir


def _placed(transform: np.ndarray) -> np.ndarray:
    scene = MOVING_CENTRES + SCENE_OFFSET
    return (scene @ transform[:3, :3].T + transform[:3, 3]) - SCENE_OFFSET


@pytest.mark.parametrize("source", ["star", "cs"])
def test_animate_moves_the_moving_map_by_plus_ro(tmp_path, source, monkeypatch):
    monkeypatch.setenv("MPLBACKEND", "Agg")
    ref_angles, mov_angles, mov_origins, relatives = _physical_poses()
    suffix = ".star" if source == "star" else ".cs"
    writer = _write_star if source == "star" else _write_cs
    ref_path, mov_path = tmp_path / f"ref{suffix}", tmp_path / f"mov{suffix}"
    writer(ref_path, ref_angles, np.zeros((len(ref_angles), 2)))
    writer(mov_path, mov_angles, mov_origins)
    run_dir = tmp_path / "run"
    assert main(
        ["run", "--ref", str(ref_path), "--mov", str(mov_path), "--output-dir", str(run_dir), "--no-visualize"]
    ) == 0

    landscape = np.load(run_dir / "data" / "raw_landscape.npz")
    ro = Rotation.from_rotvec(landscape["coordinates_analysis"])
    truth_ro = relatives[np.asarray(landscape["ref_source_row_id"], dtype=int)]
    assert np.degrees((ro.inv() * truth_ro).magnitude()).max() < 1e-3  # RO is the physical turn

    truth = Q.apply(MOVING_CENTRES - PIVOT) + PIVOT
    inverse_truth = Q.inv().apply(MOVING_CENTRES - PIVOT) + PIVOT

    active_dir = _animate(tmp_path, run_dir, "active")
    placed = _placed(_final_scene_transform(active_dir))
    assert np.abs(placed - truth).max() < 1e-4
    assert np.abs(placed - inverse_truth).max() > 5.0
    manifest = json.loads((active_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "9"
    assert manifest["ro_density_direction"] == RO_DENSITY_DIRECTION_ACTIVE
    assert manifest["ro_density_direction_policy"]["formula"] == (
        "delta_active_raw = RO_target @ inverse(RO_baseline)"
    )
    assert manifest["animation_policy"]["ro_density_direction"] == RO_DENSITY_DIRECTION_ACTIVE

    legacy_dir = _animate(tmp_path, run_dir, "legacy", ("--ro-direction", "legacy-inverse"))
    legacy = _placed(_final_scene_transform(legacy_dir))
    assert np.abs(legacy - inverse_truth).max() < 1e-4
    legacy_manifest = json.loads((legacy_dir / "manifest.json").read_text(encoding="utf-8"))
    assert legacy_manifest["ro_density_direction"] == RO_DENSITY_DIRECTION_LEGACY_INVERSE
    assert any("legacy-inverse" in warning for warning in legacy_manifest["warnings"])


def test_nonidentity_baseline_session_reaches_target():
    baseline = Rotation.from_rotvec([0.1, -0.25, 0.3])
    pivot = PIVOT + SCENE_OFFSET
    baseline_place = pivoted_scene_delta(baseline.as_matrix(), pivot)  # session saved at RO_baseline
    delta = ro_target_to_active_delta(Q, baseline, direction=RO_DENSITY_DIRECTION_ACTIVE)
    absolute = pivoted_scene_delta(delta.as_matrix(), pivot) @ baseline_place
    expected = pivoted_scene_delta(Q.as_matrix(), pivot)
    assert np.allclose(absolute, expected, atol=1e-12)
