from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from cryorole.align.star_align import align_star_files
from cryorole.cli.main import align_command, build_parser


@pytest.fixture(autouse=True)
def _work_in_tmp_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)


def _write_star(
    path: Path,
    rows: list[dict[str, object]],
    *,
    columns: list[str] | None = None,
    block: str = "data_particles",
    include_optics: bool = True,
) -> None:
    columns = columns or [
        "_rlnImageName",
        "_rlnTomoParticleName",
        "_rlnMicrographName",
        "_rlnCoordinateX",
        "_rlnCoordinateY",
        "_rlnDefocusU",
        "_rlnAngleRot",
    ]
    lines: list[str] = []
    if include_optics:
        lines.extend(
            [
                "data_optics",
                "",
                "loop_",
                "_rlnOpticsGroup #1",
                "_rlnOpticsGroupName #2",
                "1 opticsGroup1",
                "",
            ]
        )
    lines.extend([block, "", "loop_"])
    for index, column in enumerate(columns, start=1):
        lines.append(f"{column} #{index}")
    for row in rows:
        lines.append(" ".join(str(row.get(column, row.get(column.lstrip('_'), ""))) for column in columns))
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def _particle_rows(path: Path) -> list[str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    rows: list[str] = []
    in_particles = False
    in_loop = False
    for line in lines:
        stripped = line.strip()
        if stripped == "data_particles":
            in_particles = True
            continue
        if in_particles and stripped == "loop_":
            in_loop = True
            continue
        if in_particles and in_loop:
            if not stripped:
                continue
            if stripped.startswith("_"):
                continue
            if stripped.startswith("data_") or stripped == "loop_":
                break
            rows.append(stripped)
    return rows


def _basic_rows(prefix: str = "img") -> list[dict[str, object]]:
    return [
        {
            "_rlnImageName": f"{idx:06d}@/path/{prefix}_{idx}.mrcs",
            "_rlnTomoParticleName": f"particle_{idx}",
            "_rlnMicrographName": f"micro_{idx // 2}",
            "_rlnCoordinateX": float(idx * 10),
            "_rlnCoordinateY": float(idx * 20),
            "_rlnDefocusU": 1000.0 + idx,
            "_rlnAngleRot": idx,
        }
        for idx in range(4)
    ]


def test_align_help_shows_compact_public_params(capsys) -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["align", "--help"])

    help_text = capsys.readouterr().out
    for option in (
        "--ref",
        "--mov",
        "--align-id",
        "--key",
        "--float-tol",
        "--path-mode",
        "--duplicate-policy",
        "--overwrite",
    ):
        assert option in help_text
    assert "run" in help_text
    assert "--row-aligned" in help_text
    assert "--output" not in help_text
    assert "occurrence" not in help_text
    assert " error" not in help_text


def test_align_requires_ref_and_mov() -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["align"])


def test_non_star_input_fails_clearly(tmp_path) -> None:
    ref = tmp_path / "ref.cs"
    mov = tmp_path / "mov.star"
    ref.write_text("x", encoding="utf-8")
    _write_star(mov, _basic_rows())

    with pytest.raises(ValueError, match="only supports STAR"):
        align_star_files(ref=ref, mov=mov)


def test_auto_uses_tomo_particle_name_and_reorders_mov(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    ref_rows = _basic_rows()
    mov_rows = [ref_rows[2], ref_rows[0], ref_rows[1], ref_rows[3]]
    _write_star(ref, ref_rows)
    _write_star(mov, mov_rows)

    report = align_star_files(ref=ref, mov=mov)

    output_dir = tmp_path / "alignments" / "default"
    assert report["key_columns"] == ["_rlnTomoParticleName"]
    assert (output_dir / "aligned_ref.star").exists()
    assert (output_dir / "aligned_mov.star").exists()
    assert "data_optics" in (output_dir / "aligned_ref.star").read_text(encoding="utf-8")
    assert _particle_rows(output_dir / "aligned_ref.star") == [
        " ".join(str(row[column]) for column in [
            "_rlnImageName",
            "_rlnTomoParticleName",
            "_rlnMicrographName",
            "_rlnCoordinateX",
            "_rlnCoordinateY",
            "_rlnDefocusU",
            "_rlnAngleRot",
        ])
        for row in ref_rows
    ]
    assert [line.split()[1] for line in _particle_rows(output_dir / "aligned_mov.star")] == [
        "particle_0",
        "particle_1",
        "particle_2",
        "particle_3",
    ]
    table = pd.read_csv(output_dir / "match_table.csv")
    assert table.columns.tolist() == [
        "particle_key",
        "match_key",
        "ref_source_row_id",
        "mov_source_row_id",
        "ref_output_row_id",
        "mov_output_row_id",
        "match_status",
    ]
    assert table["ref_source_row_id"].tolist() == [0, 1, 2, 3]
    assert table["mov_source_row_id"].tolist() == [1, 2, 0, 3]


def test_auto_falls_back_to_image_name(tmp_path) -> None:
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    columns = ["_rlnImageName", "_rlnAngleRot"]
    rows = [
        {"_rlnImageName": "1@a.mrcs", "_rlnAngleRot": 1},
        {"_rlnImageName": "2@a.mrcs", "_rlnAngleRot": 2},
    ]
    _write_star(ref, rows, columns=columns)
    _write_star(mov, list(reversed(rows)), columns=columns)

    report = align_star_files(ref=ref, mov=mov)

    assert report["key_columns"] == ["_rlnImageName"]
    assert report["matched_count"] == 2


def test_auto_missing_safe_key_fails_with_explicit_key_hint(tmp_path) -> None:
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    columns = ["_rlnCoordinateX", "_rlnCoordinateY"]
    _write_star(ref, [{"_rlnCoordinateX": 1, "_rlnCoordinateY": 2}], columns=columns)
    _write_star(mov, [{"_rlnCoordinateX": 1, "_rlnCoordinateY": 2}], columns=columns)

    with pytest.raises(ValueError, match="--key _rlnMicrographName _rlnCoordinateX _rlnCoordinateY"):
        align_star_files(ref=ref, mov=mov)


def test_image_name_zero_overlap_fails_without_aligned_outputs(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    _write_star(ref, [{"_rlnImageName": "1@a.mrcs"}], columns=["_rlnImageName"])
    _write_star(mov, [{"_rlnImageName": "2@b.mrcs"}], columns=["_rlnImageName"])

    with pytest.raises(ValueError, match="zero overlap"):
        align_star_files(ref=ref, mov=mov)

    assert not (tmp_path / "alignments" / "default" / "aligned_ref.star").exists()


def test_low_overlap_warns_but_writes_outputs(tmp_path) -> None:
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    ref_rows = [{"_rlnImageName": f"{idx}@a.mrcs"} for idx in range(5)]
    mov_rows = [{"_rlnImageName": "0@a.mrcs"}, {"_rlnImageName": "x@b.mrcs"}]
    _write_star(ref, ref_rows, columns=["_rlnImageName"])
    _write_star(mov, mov_rows, columns=["_rlnImageName"])

    report = align_star_files(ref=ref, mov=mov)

    assert report["matched_count"] == 1
    assert report["low_overlap_warning"] is True
    assert any("Low overlap" in warning for warning in report["warnings"])
    assert (Path("alignments") / "default" / "aligned_ref.star").exists()


def test_explicit_coordinate_key_with_float_tolerance(tmp_path) -> None:
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    columns = ["_rlnMicrographName", "_rlnCoordinateX", "_rlnCoordinateY"]
    _write_star(
        ref,
        [{"_rlnMicrographName": "m1", "_rlnCoordinateX": 10.0, "_rlnCoordinateY": 20.0}],
        columns=columns,
    )
    _write_star(
        mov,
        [{"_rlnMicrographName": "m1", "_rlnCoordinateX": 10.04, "_rlnCoordinateY": 20.03}],
        columns=columns,
    )

    report = align_star_files(
        ref=ref,
        mov=mov,
        key_columns=columns,
        float_tolerances={"_rlnCoordinateX": 0.1, "_rlnCoordinateY": 0.1},
    )

    assert report["matched_count"] == 1
    assert any("user-selected tracking key" in warning for warning in report["warnings"])


def test_float_tol_invalid_and_nonnumeric_fail(tmp_path) -> None:
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "align",
                "--ref",
                "ref.star",
                "--mov",
                "mov.star",
                "--float-tol",
                "_rlnCoordinateX",
            ]
        )

    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    columns = ["_rlnImageName"]
    _write_star(ref, [{"_rlnImageName": "a"}], columns=columns)
    _write_star(mov, [{"_rlnImageName": "a"}], columns=columns)
    with pytest.raises(ValueError, match="not numeric"):
        align_star_files(
            ref=ref,
            mov=mov,
            key_columns=["_rlnImageName"],
            float_tolerances={"_rlnImageName": 0.1},
        )


def test_explicit_defocus_key_works_but_warns(tmp_path) -> None:
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    columns = ["_rlnDefocusU", "_rlnDefocusV"]
    rows = [{"_rlnDefocusU": 1000.0, "_rlnDefocusV": 1100.0}]
    _write_star(ref, rows, columns=columns)
    _write_star(mov, rows, columns=columns)

    report = align_star_files(ref=ref, mov=mov, key_columns=columns)

    assert report["matched_count"] == 1
    assert any("user-selected tracking key" in warning for warning in report["warnings"])


def test_missing_key_column_lists_available_columns(tmp_path) -> None:
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    _write_star(ref, [{"_rlnImageName": "a"}], columns=["_rlnImageName"])
    _write_star(mov, [{"_rlnImageName": "a"}], columns=["_rlnImageName"])

    with pytest.raises(ValueError, match="Available columns"):
        align_star_files(ref=ref, mov=mov, key_columns=["_rlnMissing"])


def test_path_modes_exact_basename_and_suffix(tmp_path) -> None:
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    columns = ["_rlnImageName"]
    _write_star(ref, [{"_rlnImageName": "1@/a/b/stack.mrcs"}], columns=columns)
    _write_star(mov, [{"_rlnImageName": "1@/x/y/stack.mrcs"}], columns=columns)
    with pytest.raises(ValueError, match="zero overlap"):
        align_star_files(ref=ref, mov=mov, key_columns=columns, path_mode="exact")

    basename_report = align_star_files(
        ref=ref,
        mov=mov,
        key_columns=columns,
        path_mode="basename",
        align_id="basename",
    )
    assert basename_report["matched_count"] == 1

    _write_star(mov, [{"_rlnImageName": "1@/x/b/stack.mrcs"}], columns=columns)
    suffix_report = align_star_files(
        ref=ref,
        mov=mov,
        key_columns=columns,
        path_mode="suffix:2",
        align_id="suffix",
    )
    assert suffix_report["matched_count"] == 1


def test_invalid_path_mode_fails(tmp_path) -> None:
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    _write_star(ref, [{"_rlnImageName": "a"}], columns=["_rlnImageName"])
    _write_star(mov, [{"_rlnImageName": "a"}], columns=["_rlnImageName"])

    with pytest.raises(ValueError, match="--path-mode"):
        align_star_files(ref=ref, mov=mov, path_mode="suffix:0")


def test_duplicate_policy_exclude_is_default(tmp_path) -> None:
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    columns = ["_rlnImageName", "_rlnAngleRot"]
    ref_rows = [
        {"_rlnImageName": "dup@a.mrcs", "_rlnAngleRot": 1},
        {"_rlnImageName": "dup@a.mrcs", "_rlnAngleRot": 2},
        {"_rlnImageName": "keep@a.mrcs", "_rlnAngleRot": 3},
    ]
    mov_rows = [
        {"_rlnImageName": "dup@a.mrcs", "_rlnAngleRot": 4},
        {"_rlnImageName": "keep@a.mrcs", "_rlnAngleRot": 5},
    ]
    _write_star(ref, ref_rows, columns=columns)
    _write_star(mov, mov_rows, columns=columns)

    report = align_star_files(ref=ref, mov=mov)

    output_dir = Path("alignments") / "default"
    assert report["duplicate_policy"] == "exclude"
    assert report["matched_count"] == 1
    assert report["duplicate_ref_row_count"] == 2
    assert report["duplicate_mov_row_count"] == 1
    assert len(_particle_rows(output_dir / "duplicate_ref.star")) == 2
    assert len(_particle_rows(output_dir / "duplicate_mov.star")) == 1
    assert any("Duplicate" in warning for warning in report["warnings"])


def test_duplicate_policy_first_keeps_first_occurrence(tmp_path) -> None:
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    columns = ["_rlnImageName", "_rlnAngleRot"]
    rows = [
        {"_rlnImageName": "dup@a.mrcs", "_rlnAngleRot": 1},
        {"_rlnImageName": "dup@a.mrcs", "_rlnAngleRot": 2},
    ]
    _write_star(ref, rows, columns=columns)
    _write_star(mov, rows, columns=columns)

    report = align_star_files(ref=ref, mov=mov, duplicate_policy="first")

    assert report["matched_count"] == 1
    assert report["duplicate_ref_row_count"] == 1
    assert report["duplicate_mov_row_count"] == 1
    assert any("first" in warning for warning in report["warnings"])


def test_no_public_occurrence_or_error_duplicate_policy() -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            ["align", "--ref", "ref.star", "--mov", "mov.star", "--duplicate-policy", "occurrence"]
        )
    with pytest.raises(SystemExit):
        parser.parse_args(
            ["align", "--ref", "ref.star", "--mov", "mov.star", "--duplicate-policy", "error"]
        )


def test_output_dir_exists_requires_overwrite(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    rows = [{"_rlnImageName": "a"}]
    _write_star(ref, rows, columns=["_rlnImageName"])
    _write_star(mov, rows, columns=["_rlnImageName"])
    (tmp_path / "alignments" / "default").mkdir(parents=True)

    with pytest.raises(FileExistsError, match="Align output directory already exists"):
        align_star_files(ref=ref, mov=mov)

    report = align_star_files(ref=ref, mov=mov, overwrite=True)
    assert report["matched_count"] == 1


def test_align_command_writes_outputs_and_preserves_sources(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    ref = tmp_path / "ref.star"
    mov = tmp_path / "mov.star"
    rows = _basic_rows()
    _write_star(ref, rows)
    _write_star(mov, rows)
    ref_before = ref.read_text(encoding="utf-8")
    mov_before = mov.read_text(encoding="utf-8")
    parser = build_parser()
    args = parser.parse_args(["align", "--ref", str(ref), "--mov", str(mov)])

    assert align_command(args) == 0

    output_dir = tmp_path / "alignments" / "default"
    report = json.loads((output_dir / "align_report.json").read_text(encoding="utf-8"))
    for key in (
        "artifact_type",
        "schema_version",
        "status",
        "ref_path",
        "mov_path",
        "align_id",
        "output_dir",
        "output_paths",
        "key_policy",
        "key_columns",
        "key_source",
        "path_mode",
        "float_tolerances",
        "duplicate_policy",
        "ref_row_count",
        "mov_row_count",
        "matched_count",
        "ref_only_count",
        "mov_only_count",
        "duplicate_ref_row_count",
        "duplicate_mov_row_count",
        "duplicate_ref_key_count",
        "duplicate_mov_key_count",
        "ref_overlap_fraction",
        "mov_overlap_fraction",
        "low_overlap_warning",
        "warnings",
    ):
        assert key in report
    assert ref.read_text(encoding="utf-8") == ref_before
    assert mov.read_text(encoding="utf-8") == mov_before
