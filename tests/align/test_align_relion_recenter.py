"""Stage 2: RELION recentring geometry, recentered-exact, chain, subtract correction.

Real-data fixtures (tests/fixtures/relion_recenter/) are deterministic subsets
(first 120 micrographs, 11 columns) of one RELION project:

    Refine3D/job022 (consensus) --relion_particle_subtract --center 9 -22 102--> Subtract/job041
    job043_coords_exchange_job042 (refinement of the subtracted particles, coordinates corrected by hand)
    --relion_preprocess --reextract --recenter -9 22 -102--> Extract/job055
    job052: refinement from an unrelated subtraction (negative control)
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cryorole.align.relion_geometry import (
    predict_reextraction,
    predict_subtraction,
    relion_euler_matrices,
    relion_round,
)
from cryorole.align.star_align import align_star_files
from cryorole.align.subtract_fix import fix_subtract_coordinates
from cryorole.cli.main import align_command, build_parser
from cryorole.io.readers.star_reader import read_relion_star

F = Path(__file__).resolve().parents[1] / "fixtures" / "relion_recenter"
J22 = F / "job022_run_data_subset.star"
J41 = F / "job041_particles_subtracted_subset.star"
J43 = F / "job043_coords_exchange_job042_subset.star"
J55 = F / "job055_particles_subset.star"
J52 = F / "job052_run_data_subset.star"
RECENTER = (-9.0, 22.0, -102.0)
COORDS = ["_rlnCoordinateX", "_rlnCoordinateY"]


def _particles(path: Path):
    return read_relion_star(path).particles


def _write_synthetic_star(path: Path, rows: dict[str, list], *, image_angpix: float, mic_angpix: float) -> None:
    columns = list(rows)
    n = len(rows[columns[0]])
    lines = [
        "data_optics", "", "loop_", "_rlnOpticsGroupName #1", "_rlnOpticsGroup #2",
        "_rlnMicrographOriginalPixelSize #3", "_rlnImagePixelSize #4",
        f"opticsGroup1 1 {mic_angpix:.6f} {image_angpix:.6f}", "", "data_particles", "", "loop_",
        *[f"{c} #{i + 1}" for i, c in enumerate(columns)],
    ]
    for i in range(n):
        lines.append(" ".join(str(rows[c][i]) for c in columns))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------- geometry


def test_relion_matrix_is_relion_euler_angles2matrix() -> None:
    angles = np.random.default_rng(0).uniform(-180, 180, (50, 3))
    ours = relion_euler_matrices(*angles.T)
    # RELION's A equals the transpose of SciPy's intrinsic ZYZ matrix (independent check).
    reference = Rotation.from_euler("ZYZ", angles, degrees=True).as_matrix().transpose(0, 2, 1)
    np.testing.assert_allclose(ours, reference, atol=1e-12)


def test_relion_round_is_half_away_from_zero() -> None:
    np.testing.assert_array_equal(relion_round(np.array([0.5, 1.5, 2.5, -0.5, -1.5, 0.49])), [1, 2, 3, -1, -2, 0])


def _array(path: Path, names):
    return _particles(path)[list(names)].to_numpy(float)


def test_reextraction_prediction_reproduces_job055_exactly() -> None:
    angles = ["_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi"]
    origins = ["_rlnOriginXAngst", "_rlnOriginYAngst"]
    prediction = predict_reextraction(
        coordinates=_array(J43, COORDS),
        origins_angst=_array(J43, origins),
        angles_deg=_array(J43, angles),
        recenter_ref_px=np.array(RECENTER),
        ref_angpix=0.835,
        particle_angpix=0.835,
        micrograph_angpix=0.835,
    )
    np.testing.assert_array_equal(prediction.coordinates, _array(J55, COORDS))
    assert np.abs(prediction.origins_angst - _array(J55, origins)).max() < 2e-3
    np.testing.assert_array_equal(_array(J43, angles), _array(J55, angles))


def test_transposed_matrix_negative_control_fails() -> None:
    angles = _array(J43, ["_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi"])
    matrices = relion_euler_matrices(*angles.T).transpose(0, 2, 1)
    shift = np.einsum("nij,j->ni", matrices[:, :2, :], np.array(RECENTER) * 0.835)
    q = (_array(J43, ["_rlnOriginXAngst", "_rlnOriginYAngst"]) - shift) / 0.835
    agreement = np.all(_array(J43, COORDS) - relion_round(q) == _array(J55, COORDS), axis=1).mean()
    assert agreement < 0.01


def test_subtraction_prediction_reproduces_job041_and_the_job043_correction() -> None:
    prediction = predict_subtraction(
        coordinates=_array(J22, COORDS),
        origins_angst=_array(J22, ["_rlnOriginXAngst", "_rlnOriginYAngst"]),
        angles_deg=_array(J22, ["_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi"]),
        center_model_px=np.array([9.0, -22.0, 102.0]),
        model_angpix=0.835,
        particle_angpix=0.835,
        micrograph_angpix=0.835,
    )
    assert np.abs(prediction.origins_angst - _array(J41, ["_rlnOriginXAngst", "_rlnOriginYAngst"])).max() < 2e-3
    np.testing.assert_array_equal(_array(J41, COORDS), _array(J22, COORDS))  # RELION leaves them stale
    np.testing.assert_array_equal(prediction.corrected_coordinates, _array(J43, COORDS))


# --------------------------------------------------------------------------- recentered-exact


def test_recentered_exact_golden_job043_to_job055(tmp_path) -> None:
    report = align_star_files(
        ref=J43, mov=J55, coordinate_match="recentered-exact", recenter_shift=RECENTER, output_dir=tmp_path / "out"
    )
    match = report["coordinate_match"]
    assert report["strategy"] == "recentered-exact"
    assert match["verified_fraction_of_smaller_input"] == 1.0
    assert match["unverified_count"] == 0
    assert match["max_origin_error_angst"] < 2e-3
    # Duplicate picks that converged to identical geometry cannot be told apart; they are excluded, not paired.
    assert match["indistinguishable_input_duplicates"] == report["ambiguous_ref_row_count"] > 0
    assert report["matched_count"] + report["ambiguous_ref_row_count"] == 1057
    assert any("indistinguishable duplicates" in w for w in report["warnings"])
    ref_rows = _particles(tmp_path / "out" / "aligned_ref.star")
    mov_rows = _particles(tmp_path / "out" / "aligned_mov.star")
    assert ref_rows["_rlnImageOriginalName"].tolist() == mov_rows["_rlnImageOriginalName"].tolist()
    table = (tmp_path / "out" / "match_table.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "coord_verified" in table and "origin_error_A" in table


def test_recentered_exact_angle_mismatch_disqualifies_the_pair(tmp_path) -> None:
    text = J55.read_text(encoding="utf-8").splitlines()
    first_row = next(i for i, line in enumerate(text) if line.startswith("000001@"))
    fields = text[first_row].split()
    fields[4] = f"{float(fields[4]) + 0.01:.6f}"  # _rlnAngleRot
    text[first_row] = " ".join(fields)
    mov = tmp_path / "job055_modified.star"
    mov.write_text("\n".join(text) + "\n", encoding="utf-8")

    report = align_star_files(
        ref=J43, mov=mov, coordinate_match="recentered-exact", recenter_shift=RECENTER, output_dir=tmp_path / "out"
    )

    assert report["coordinate_match"]["unverified_reasons"]["angle_mismatch"] == 1
    assert "000001@job055.mrcs" not in [r for r in _particles(tmp_path / "out" / "aligned_mov.star")["_rlnImageName"]]
    assert _particles(tmp_path / "out" / "unverified_mov.star")["_rlnImageName"].tolist() == ["000001@job055.mrcs"]


@pytest.mark.parametrize("shift", [(9.0, -22.0, 102.0), (-9.0, 22.0, 102.0)])
def test_recentered_exact_wrong_vector_writes_nothing(tmp_path, shift) -> None:
    with pytest.raises(ValueError, match="not an extraction input/output pair"):
        align_star_files(ref=J43, mov=J55, coordinate_match="recentered-exact", recenter_shift=shift,
                         output_dir=tmp_path / "out")
    assert not (tmp_path / "out" / "aligned_ref.star").exists()


def test_recentered_exact_negative_control_job052(tmp_path) -> None:
    with pytest.raises(ValueError, match="not an extraction input/output pair"):
        align_star_files(ref=J52, mov=J55, coordinate_match="recentered-exact", recenter_shift=RECENTER,
                         output_dir=tmp_path / "out")


def test_recentered_exact_uses_micrograph_pixel_rounding_and_ref_angpix(tmp_path) -> None:
    rng = np.random.default_rng(5)
    n = 40
    a_part, a_mic, a_ref = 1.67, 0.835, 1.2
    coords = np.column_stack([rng.integers(200, 3000, n), rng.integers(200, 3000, n)]).astype(float)
    origins = rng.uniform(-10, 10, (n, 2))
    angles = rng.uniform(-180, 180, (n, 3))
    angles[:, 1] = rng.uniform(0, 180, n)
    shift = np.array([5.0, -7.0, 31.0])
    out = predict_reextraction(coordinates=coords, origins_angst=origins, angles_deg=angles, recenter_ref_px=shift,
                               ref_angpix=a_ref, particle_angpix=a_part, micrograph_angpix=a_mic)
    base = {"_rlnMicrographName": ["m1"] * n, "_rlnAngleRot": [f"{a:.6f}" for a in angles[:, 0]],
            "_rlnAngleTilt": [f"{a:.6f}" for a in angles[:, 1]], "_rlnAnglePsi": [f"{a:.6f}" for a in angles[:, 2]]}
    angles = np.round(angles, 6)
    out = predict_reextraction(coordinates=coords, origins_angst=origins, angles_deg=angles, recenter_ref_px=shift,
                               ref_angpix=a_ref, particle_angpix=a_part, micrograph_angpix=a_mic)
    _write_synthetic_star(tmp_path / "in.star", {
        "_rlnImageName": [f"{i + 1}@in.mrcs" for i in range(n)], **base,
        "_rlnCoordinateX": [f"{v:.6f}" for v in coords[:, 0]], "_rlnCoordinateY": [f"{v:.6f}" for v in coords[:, 1]],
        "_rlnOriginXAngst": [f"{v:.6f}" for v in origins[:, 0]], "_rlnOriginYAngst": [f"{v:.6f}" for v in origins[:, 1]],
    }, image_angpix=a_part, mic_angpix=a_mic)
    order = rng.permutation(n)
    _write_synthetic_star(tmp_path / "out.star", {
        "_rlnImageName": [f"{i + 1}@out.mrcs" for i in range(n)],
        **{k: [v[i] for i in order] for k, v in base.items()},
        "_rlnCoordinateX": [f"{out.coordinates[i, 0]:.6f}" for i in order],
        "_rlnCoordinateY": [f"{out.coordinates[i, 1]:.6f}" for i in order],
        "_rlnOriginXAngst": [f"{out.origins_angst[i, 0]:.6f}" for i in order],
        "_rlnOriginYAngst": [f"{out.origins_angst[i, 1]:.6f}" for i in order],
    }, image_angpix=a_part, mic_angpix=a_mic)

    report = align_star_files(ref=tmp_path / "in.star", mov=tmp_path / "out.star", coordinate_match="recentered-exact",
                              recenter_shift=shift, ref_angpix=a_ref, output_dir=tmp_path / "ok")
    assert report["matched_count"] == n
    with pytest.raises(ValueError, match="not an extraction input/output pair"):
        align_star_files(ref=tmp_path / "in.star", mov=tmp_path / "out.star", coordinate_match="recentered-exact",
                         recenter_shift=shift, output_dir=tmp_path / "no_ref_angpix")


# --------------------------------------------------------------------------- chain and keys


def test_via_extraction_chain_pairs_only_rows_verified_on_every_link(tmp_path) -> None:
    report = align_star_files(ref=J22, mov=J55, via_extraction=[J43, J55], recenter_shift=RECENTER,
                              output_dir=tmp_path / "chain")
    links = report["chain"]["links"]
    assert links[0]["key"] == "_rlnImageName=_rlnImageOriginalName"
    assert links[2]["key"] == "_rlnImageName"
    assert report["chain"]["rows_lost_per_link"]["input_to_output"] == links[1]["ambiguous_input_rows"]
    ref_rows = _particles(tmp_path / "chain" / "aligned_ref.star")
    mov_rows = _particles(tmp_path / "chain" / "aligned_mov.star")
    # Ground truth for this project: job055 keeps the job022 image name in _rlnImageOriginalName.
    assert ref_rows["_rlnImageName"].tolist() == mov_rows["_rlnImageOriginalName"].tolist()
    assert len(ref_rows) == report["matched_count"] == links[1]["verified_count"]


def test_key_pair_links_consensus_to_subtracted_refinement(tmp_path) -> None:
    report = align_star_files(ref=J22, mov=J43, key_pairs=["_rlnImageName=_rlnImageOriginalName"],
                              output_dir=tmp_path / "kp")
    assert report["strategy"] == "key-pair"
    assert report["matched_count"] == 1057
    assert report["key_columns"] == ["_rlnImageName=_rlnImageOriginalName"]
    with pytest.raises(ValueError, match="REF_COLUMN=MOV_COLUMN"):
        align_star_files(ref=J22, mov=J43, key_pairs=["_rlnImageName"], output_dir=tmp_path / "bad")


def test_coordinate_unchanged_matches_subtraction_without_moving_coordinates(tmp_path) -> None:
    report = align_star_files(ref=J22, mov=J41, coordinate_match="unchanged", output_dir=tmp_path / "u")
    assert report["strategy"] == "coordinate-unchanged"
    assert report["matched_count"] == 1057


def test_float_tolerance_neighbour_is_ambiguous(tmp_path) -> None:
    def star(path, rows):
        _write_synthetic_star(path, {"_rlnMicrographName": [r[0] for r in rows],
                                     "_rlnCoordinateX": [r[1] for r in rows], "_rlnCoordinateY": [r[2] for r in rows]},
                              image_angpix=1.0, mic_angpix=1.0)

    star(tmp_path / "ref.star", [("m", 10.0, 10.0), ("m", 50.0, 50.0)])
    star(tmp_path / "mov.star", [("m", 10.02, 10.0), ("m", 10.6, 10.0), ("m", 50.0, 50.0)])
    report = align_star_files(ref=tmp_path / "ref.star", mov=tmp_path / "mov.star",
                              key_columns=["_rlnMicrographName", "_rlnCoordinateX", "_rlnCoordinateY"],
                              float_tolerances={"_rlnCoordinateX": 0.5, "_rlnCoordinateY": 0.5},
                              output_dir=tmp_path / "out")
    assert report["matched_count"] == 1
    assert report["ambiguous_ref_row_count"] == 1
    assert len(_particles(tmp_path / "out" / "ambiguous_ref.star")) == 1


def test_default_output_location_report_hashes_and_next_command(tmp_path) -> None:
    ref = tmp_path / "project" / "ref.star"
    ref.parent.mkdir()
    shutil.copy(J22, ref)
    report = align_star_files(ref=ref, mov=J43, key_pairs=["_rlnImageName=_rlnImageOriginalName"])
    out = ref.parent / "cryorole_alignments" / "default"
    assert Path(report["output_dir"]) == out
    assert set(report["sha256"]) == {"ref_source", "mov_source", "aligned_ref", "aligned_mov", "match_table"}
    assert report["next_command"].startswith("cryorole run --ref ")
    assert str((out / "aligned_ref.star").resolve()) in report["next_command"]
    assert report["next_command"].endswith("--row-aligned")
    assert json.loads((out / "align_report.json").read_text(encoding="utf-8"))["schema_version"] == "2"


# --------------------------------------------------------------------------- subtract correction


def _subtract_project(tmp_path: Path, *, note: bool = True) -> Path:
    project = tmp_path / "project"
    (project / "Subtract" / "job041").mkdir(parents=True)
    (project / "Refine3D" / "job022").mkdir(parents=True)
    shutil.copy(J41, project / "Subtract" / "job041" / "particles_subtracted.star")
    shutil.copy(J22, project / "Refine3D" / "job022" / "run_data.star")
    if note:
        shutil.copy(F / "Subtract_job041_note.txt", project / "Subtract" / "job041" / "note.txt")
    return project


def test_fix_subtract_coordinates_golden_matches_job043(tmp_path) -> None:
    project = _subtract_project(tmp_path)
    before = (project / "Subtract" / "job041" / "particles_subtracted.star").read_bytes()

    report = fix_subtract_coordinates(project / "Subtract" / "job041", log=lambda *_: None)

    corrected = _particles(Path(report["outputs"]["corrected_star"]))
    np.testing.assert_array_equal(corrected[COORDS].to_numpy(float), _array(J43, COORDS))
    original = _particles(J41)
    assert corrected.drop(columns=COORDS).equals(original.drop(columns=COORDS))
    assert report["verification"]["verified_fraction"] == 1.0
    assert report["verification"]["stale_coordinate_rows"] == 1057
    assert "note.txt" in report["sources"]["center"]
    assert Path(report["outputs"]["corrected_star"]).parent == (project / "cryorole_alignments" / "fix_subtract_job041").resolve()
    assert (project / "Subtract" / "job041" / "particles_subtracted.star").read_bytes() == before


def test_fix_subtract_coordinates_apply_to_downstream_refinement(tmp_path) -> None:
    project = _subtract_project(tmp_path)
    report = fix_subtract_coordinates(project / "Subtract" / "job041", apply_to=J43, output_dir=tmp_path / "o",
                                      log=lambda *_: None)
    out = _particles(Path(report["outputs"]["corrected_star"]))
    assert out.equals(_particles(J43))  # job043 is exactly job042 with these corrected coordinates
    assert report["apply_to"]["corrected_rows"] == 1057


def test_fix_subtract_coordinates_wrong_parameters_write_nothing(tmp_path) -> None:
    project = _subtract_project(tmp_path)
    with pytest.raises(ValueError, match="Nothing was written"):
        fix_subtract_coordinates(project / "Subtract" / "job041", center=[-9, 22, -102], output_dir=tmp_path / "o1",
                                 log=lambda *_: None)
    with pytest.raises(ValueError, match="Nothing was written"):
        fix_subtract_coordinates(project / "Subtract" / "job041", model_angpix=1.67, output_dir=tmp_path / "o2",
                                 log=lambda *_: None)
    assert not (tmp_path / "o1").exists() and not (tmp_path / "o2").exists()


def test_fix_subtract_coordinates_without_note_needs_explicit_values(tmp_path) -> None:
    project = _subtract_project(tmp_path, note=False)
    with pytest.raises(ValueError, match="--center"):
        fix_subtract_coordinates(project / "Subtract" / "job041", log=lambda *_: None)
    report = fix_subtract_coordinates(
        project / "Subtract" / "job041" / "particles_subtracted.star",
        center=[9, -22, 102],
        subtract_input=project / "Refine3D" / "job022" / "run_data.star",
        output_dir=tmp_path / "o",
        log=lambda *_: None,
    )
    assert report["sources"]["center"] == "--center"


def test_fix_subtract_rounds_in_particle_pixels_when_binned(tmp_path) -> None:
    rng = np.random.default_rng(3)
    n = 60
    a_part, a_mic = 1.67, 0.835
    coords = np.column_stack([rng.integers(300, 3000, n), rng.integers(300, 3000, n)]).astype(float)
    origins = np.round(rng.uniform(-12, 12, (n, 2)), 6)
    angles = np.round(np.column_stack([rng.uniform(-180, 180, n), rng.uniform(0, 180, n), rng.uniform(-180, 180, n)]), 6)
    center = np.array([7.0, -11.0, 40.0])
    sub = predict_subtraction(coordinates=coords, origins_angst=origins, angles_deg=angles, center_model_px=center,
                              model_angpix=a_part, particle_angpix=a_part, micrograph_angpix=a_mic)
    common = {"_rlnMicrographName": ["m1"] * n, "_rlnCoordinateX": [f"{v:.6f}" for v in coords[:, 0]],
              "_rlnCoordinateY": [f"{v:.6f}" for v in coords[:, 1]], "_rlnAngleRot": [f"{v:.6f}" for v in angles[:, 0]],
              "_rlnAngleTilt": [f"{v:.6f}" for v in angles[:, 1]], "_rlnAnglePsi": [f"{v:.6f}" for v in angles[:, 2]]}
    _write_synthetic_star(tmp_path / "input.star", {"_rlnImageName": [f"{i + 1}@p.mrcs" for i in range(n)], **common,
                          "_rlnOriginXAngst": [f"{v:.6f}" for v in origins[:, 0]],
                          "_rlnOriginYAngst": [f"{v:.6f}" for v in origins[:, 1]]}, image_angpix=a_part, mic_angpix=a_mic)
    _write_synthetic_star(tmp_path / "subtracted.star", {
        "_rlnImageName": [f"{i + 1}@s.mrcs" for i in range(n)], **common,
        "_rlnOriginXAngst": [f"{v:.6f}" for v in sub.origins_angst[:, 0]],
        "_rlnOriginYAngst": [f"{v:.6f}" for v in sub.origins_angst[:, 1]],
        "_rlnImageOriginalName": [f"{i + 1}@p.mrcs" for i in range(n)]}, image_angpix=a_part, mic_angpix=a_mic)

    report = fix_subtract_coordinates(tmp_path / "subtracted.star", center=center, subtract_input=tmp_path / "input.star",
                                      output_dir=tmp_path / "o", log=lambda *_: None)
    corrected = _particles(Path(report["outputs"]["corrected_star"]))[COORDS].to_numpy(float)
    moved = coords - corrected
    assert np.all(np.mod(moved, 2) == 0)  # whole particle pixels = even micrograph pixels at 2x binning
    np.testing.assert_array_equal(corrected, sub.corrected_coordinates)


def test_align_cli_fix_subtract_and_next_command(tmp_path, capsys) -> None:
    project = _subtract_project(tmp_path)
    args = build_parser().parse_args(["align", "--fix-subtract-coordinates", str(project / "Subtract" / "job041")])
    assert align_command(args) == 0
    out = capsys.readouterr()
    assert out.out.strip().endswith("particles_subtracted_coords_corrected.star")
    assert "verification: 1057/1057" in out.err

    args = build_parser().parse_args(["align", "--ref", str(J22), "--mov", str(J43), "--key-pair",
                                      "_rlnImageName=_rlnImageOriginalName", "--output-dir", str(tmp_path / "a")])
    assert align_command(args) == 0
    assert "Next: cryorole run --ref" in capsys.readouterr().err
    with pytest.raises(ValueError, match="only used with --fix-subtract-coordinates"):
        align_command(build_parser().parse_args(["align", "--ref", str(J22), "--mov", str(J43), "--center", "1", "2", "3"]))
