from pathlib import Path

from cryorole.io.readers.star_reader import read_relion_star


def test_star_reader_parses_raw_tables_without_interpretation(tmp_path: Path) -> None:
    star_path = tmp_path / "particles.star"
    star_path.write_text(
        "\n".join(
            [
                "data_optics",
                "",
                "loop_",
                "_rlnOpticsGroup #1",
                "_rlnVoltage #2",
                "1 300",
                "",
                "data_particles",
                "",
                "loop_",
                "_rlnImageName #1",
                "_rlnParticleName #2",
                "_rlnAngleRot #3",
                "_rlnAngleTilt #4",
                "_rlnAnglePsi #5",
                "_rlnOriginXAngst #6",
                "_rlnOriginYAngst #7",
                "000001@stack.mrcs p001 10 20 30 1.5 -2.0",
                "000002@stack.mrcs p002 40 50 60 0.0 3.0",
            ]
        ),
        encoding="utf-8",
    )

    data = read_relion_star(star_path)

    assert data.report.source_type == "relion"
    assert data.report.row_count == 2
    assert data.particles.index.name == "source_row_id"
    assert list(data.particles.index) == [0, 1]
    assert list(data.optics.columns) == ["_rlnOpticsGroup", "_rlnVoltage"]
    assert list(data.particles.columns) == [
        "_rlnImageName",
        "_rlnParticleName",
        "_rlnAngleRot",
        "_rlnAngleTilt",
        "_rlnAnglePsi",
        "_rlnOriginXAngst",
        "_rlnOriginYAngst",
    ]
    assert data.particles.loc[0, "_rlnAngleRot"] == "10"
    assert "_rotation_matrix_active" not in data.particles.columns
