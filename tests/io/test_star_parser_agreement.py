"""P0-1 / P0-2: one STAR structure parser for run, selection, and export paths."""

from __future__ import annotations

from pathlib import Path

import pytest

from cryorole.export.star_subset import write_relion_star_subset
from cryorole.io.readers.star_reader import (
    locate_particle_loop,
    read_relion_star,
    read_relion_star_particle_columns,
)

POSE = ("_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi")


def _relion31_spa() -> str:
    return "\n".join([
        "# version 30001",
        "data_optics",
        "loop_",
        "_rlnOpticsGroup #1",
        "_rlnVoltage #2",
        "1 300",
        "2 300",
        "",
        "# version 30001",
        "data_particles",
        "",
        "loop_",
        "_rlnImageName #1",
        "_rlnAngleRot #2",
        "_rlnAngleTilt #3",
        "_rlnAnglePsi #4",
        "_rlnOpticsGroup #5",
        "1@a.mrcs 10 20 30 1",
        "2@a.mrcs 11 21 31 1",
        "# interleaved comment",
        "3@a.mrcs 12 22 32 2",
        "",
    ])


def _relion4_tomo_trailing_tomograms() -> str:
    return "\n".join([
        "data_optics",
        "loop_",
        "_rlnOpticsGroup #1",
        "_rlnTomoTiltSeriesPixelSize #2",
        "1 1.35",
        "",
        "data_particles",
        "loop_",
        "_rlnTomoName #1",
        "_rlnTomoParticleName #2",
        "_rlnAngleRot #3",
        "_rlnAngleTilt #4",
        "_rlnAnglePsi #5",
        "tomo1 tomo1/1 10 20 30",
        "tomo1 tomo1/2 11 21 31",
        "tomo2 tomo2/3 12 22 32",
        "",
        "data_tomograms",
        "loop_",
        "_rlnTomoName #1",
        "_rlnVoltage #2",
        "tomo1 300",
        "tomo2 300",
        "tomo3 300",
        "tomo4 300",
        "tomo5 300",
        "",
    ])


def _relion5_general_block() -> str:
    return "\n".join([
        "data_general",
        "",
        "_rlnTomoSubTomosAre2DStacks 1",
        "",
        "data_optics",
        "loop_",
        "_rlnOpticsGroup #1",
        "1",
        "",
        "data_particles",
        "loop_",
        "_rlnTomoParticleName #1",
        "_rlnAngleRot #2",
        "_rlnAngleTilt #3",
        "_rlnAnglePsi #4",
        "t/1 10 20 30",
        "t/2 11 21 31",
        "t/3 12 22 32",
        "",
    ])


def _optics_last() -> str:
    return "\n".join([
        "data_particles",
        "loop_",
        "_rlnImageName #1",
        "_rlnAngleRot #2",
        "_rlnAngleTilt #3",
        "_rlnAnglePsi #4",
        "1@a.mrcs 10 20 30",
        "2@a.mrcs 11 21 31",
        "3@a.mrcs 12 22 32",
        "",
        "data_optics",
        "loop_",
        "_rlnOpticsGroup #1",
        "_rlnVoltage #2",
        "1 300",
    ])


REALISTIC = {
    "relion31_spa": _relion31_spa,
    "relion4_tomo": _relion4_tomo_trailing_tomograms,
    "relion5_general": _relion5_general_block,
    "optics_last": _optics_last,
}


@pytest.mark.parametrize("name", sorted(REALISTIC))
def test_star_readers_agree_on_row_ids_for_realistic_fixtures(tmp_path: Path, name: str) -> None:
    path = tmp_path / f"{name}.star"
    path.write_text(REALISTIC[name](), encoding="utf-8")

    full = read_relion_star(path)
    streamed = read_relion_star_particle_columns(path, columns=POSE)
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    loop = locate_particle_loop(lines)

    assert full.report.row_count == streamed.report.row_count == loop.row_count == 3
    assert list(full.particles.index) == list(streamed.particles.index) == [0, 1, 2]
    for column in POSE:
        assert full.particles[column].tolist() == streamed.particles[column].tolist()
    assert full.particles["_rlnAngleRot"].tolist() == ["10", "11", "12"]
    assert tuple(streamed.report.columns) == tuple(full.particles.columns)


def _two_particle_blocks() -> str:
    return "\n".join([
        "data_particles",
        "loop_",
        "_rlnImageName #1",
        "_rlnAngleRot #2",
        "1@a.mrcs 10",
        "",
        "data_particles",
        "loop_",
        "_rlnImageName #1",
        "_rlnAngleRot #2",
        "2@a.mrcs 11",
        "3@a.mrcs 12",
    ])


def _two_loops_in_particle_block() -> str:
    return "\n".join([
        "data_particles",
        "loop_",
        "_rlnImageName #1",
        "_rlnAngleRot #2",
        "1@a.mrcs 10",
        "loop_",
        "_rlnImageName #1",
        "_rlnAngleRot #2",
        "2@a.mrcs 11",
        "3@a.mrcs 12",
    ])


@pytest.mark.parametrize("builder", [_two_particle_blocks, _two_loops_in_particle_block])
def test_every_star_path_rejects_multiple_particle_loops(tmp_path: Path, builder) -> None:
    path = tmp_path / "dup.star"
    path.write_text(builder(), encoding="utf-8")

    with pytest.raises(ValueError, match="multiple data_particles"):
        read_relion_star(path)
    with pytest.raises(ValueError, match="multiple data_particles"):
        read_relion_star_particle_columns(path, columns=("_rlnAngleRot",))
    with pytest.raises(ValueError, match="multiple data_particles"):
        write_relion_star_subset(path, tmp_path / "out.star", [0], row_id_field="ref_source_row_id")


def test_full_reader_keeps_optics_when_loop_follows_without_block_change(tmp_path: Path) -> None:
    """Regression: the old reader reset buffers on ``loop_`` without flushing."""

    path = tmp_path / "optics.star"
    path.write_text(_relion31_spa(), encoding="utf-8")

    data = read_relion_star(path)

    assert list(data.optics["_rlnOpticsGroup"]) == ["1", "2"]
    assert data.particles["_rlnImageName"].tolist() == ["1@a.mrcs", "2@a.mrcs", "3@a.mrcs"]


def test_star_subset_export_selects_data_particles_loop_in_tomo_star_without_image_name(
    tmp_path: Path,
) -> None:
    source = tmp_path / "tomo.star"
    source.write_text(_relion4_tomo_trailing_tomograms(), encoding="utf-8")
    before = source.read_text(encoding="utf-8")
    output = tmp_path / "selected.star"

    result = write_relion_star_subset(
        source,
        output,
        [0, 2],
        row_id_field="ref_source_row_id",
        expected_source_row_count=3,
        required_headers=POSE,
    )

    exported = read_relion_star(output)
    assert result.selected_row_count == 2
    assert result.source_row_count == 3
    assert exported.particles["_rlnTomoParticleName"].tolist() == ["tomo1/1", "tomo2/3"]
    text = output.read_text(encoding="utf-8")
    # The trailing tomogram table is copied verbatim, not subset.
    for tomo in ("tomo1 300", "tomo2 300", "tomo3 300", "tomo4 300", "tomo5 300"):
        assert tomo in text
    assert source.read_text(encoding="utf-8") == before


def test_star_subset_export_fails_loudly_on_recorded_row_count_mismatch(tmp_path: Path) -> None:
    source = tmp_path / "tomo.star"
    source.write_text(_relion4_tomo_trailing_tomograms(), encoding="utf-8")

    with pytest.raises(ValueError, match="has 3 rows but the run recorded 5; refusing to export"):
        write_relion_star_subset(
            source,
            tmp_path / "out.star",
            [0],
            row_id_field="ref_source_row_id",
            expected_source_row_count=5,
        )
    assert not (tmp_path / "out.star").exists()


def test_star_subset_export_requires_pose_columns_in_resolved_loop(tmp_path: Path) -> None:
    source = tmp_path / "noposes.star"
    source.write_text(
        "\n".join(["data_particles", "loop_", "_rlnImageName #1", "1@a.mrcs", "2@a.mrcs"]),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="missing required columns"):
        write_relion_star_subset(
            source, tmp_path / "out.star", [0], row_id_field="ref_source_row_id", required_headers=POSE
        )


def test_star_without_particle_block_has_no_particle_table(tmp_path: Path) -> None:
    source = tmp_path / "legacy.star"
    source.write_text(
        "\n".join(["data_", "loop_", "_rlnImageName #1", "_rlnAngleRot #2", "1@a.mrcs 10"]),
        encoding="utf-8",
    )
    assert read_relion_star(source).report.row_count == 0
    assert read_relion_star_particle_columns(source, columns=("_rlnAngleRot",)).report.row_count == 0
    with pytest.raises(ValueError, match="no 'data_particles' block"):
        write_relion_star_subset(source, tmp_path / "out.star", [0], row_id_field="ref_source_row_id")


def test_star_subset_last_row_without_trailing_newline_is_not_glued(tmp_path: Path) -> None:
    source = tmp_path / "nonl.star"
    source.write_text(
        "\n".join([
            "data_particles", "loop_", "_rlnImageName #1", "_rlnAngleRot #2",
            "1@a.mrcs 10", "2@a.mrcs 11", "3@a.mrcs 12",
        ]),
        encoding="utf-8",
    )
    output = tmp_path / "out.star"
    write_relion_star_subset(source, output, [2, 0], row_id_field="ref_source_row_id")
    assert read_relion_star(output).particles["_rlnImageName"].tolist() == ["3@a.mrcs", "1@a.mrcs"]
