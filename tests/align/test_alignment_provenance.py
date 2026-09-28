"""A2: verified align lineage in `run --row-aligned`; never a requirement for running."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from cryorole.align.star_align import align_star_files
from cryorole.cli.main import build_parser, run_command

pytestmark = pytest.mark.private_fixtures("relion_recenter")

F = Path(__file__).resolve().parents[1] / "fixtures" / "relion_recenter"


def _aligned(tmp_path: Path) -> tuple[Path, Path, Path]:
    project = tmp_path / "project"
    project.mkdir()
    ref = project / "consensus.star"
    mov = project / "body.star"
    shutil.copy(F / "job022_run_data_subset.star", ref)
    shutil.copy(F / "job043_coords_exchange_job042_subset.star", mov)
    report = align_star_files(ref=ref, mov=mov, key_pairs=["_rlnImageName=_rlnImageOriginalName"])
    out = Path(report["output_dir"])
    return out / "aligned_ref.star", out / "aligned_mov.star", out


def _run(ref: Path, mov: Path, out: Path) -> dict:
    args = build_parser().parse_args([
        "run", "--ref", str(ref), "--mov", str(mov), "--output-dir", str(out), "--row-aligned",
        "--k-neighbors", "10", "--no-visualize",
    ])
    assert run_command(args) == 0
    return json.loads((out / "run_summary.json").read_text(encoding="utf-8"))


def test_verified_lineage_is_attached_with_original_status(tmp_path) -> None:
    ref, mov, align_dir = _aligned(tmp_path)
    summary = _run(ref, mov, tmp_path / "run")
    provenance = summary["alignment_provenance"]
    assert provenance["attached"] is True
    assert provenance["lineage"]["strategy"] == "key-pair"
    assert provenance["original_files"] == {
        "ref": {"path": str((tmp_path / "project" / "consensus.star").resolve()), "status": "verified"},
        "mov": {"path": str((tmp_path / "project" / "body.star").resolve()), "status": "verified"},
    }
    manifest = json.loads((tmp_path / "run" / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["results"]["alignment_provenance"]["attached"] is True
    assert "Alignment lineage: `key-pair`" in (tmp_path / "run" / "run_report.md").read_text(encoding="utf-8")


def test_lineage_attaches_when_originals_moved_or_changed(tmp_path) -> None:
    ref, mov, _ = _aligned(tmp_path)
    (tmp_path / "project" / "consensus.star").rename(tmp_path / "moved.star")
    with (tmp_path / "project" / "body.star").open("a", encoding="utf-8") as handle:
        handle.write("\n# edited\n")
    provenance = _run(ref, mov, tmp_path / "run")["alignment_provenance"]
    assert provenance["attached"] is True
    assert provenance["original_files"]["ref"]["status"] == "unavailable"
    assert provenance["original_files"]["mov"]["status"] == "mismatch"


@pytest.mark.parametrize("damage", ["report_json", "match_table", "aligned_file"])
def test_unverifiable_report_warns_but_run_continues(tmp_path, damage, capsys) -> None:
    ref, mov, align_dir = _aligned(tmp_path)
    if damage == "report_json":
        (align_dir / "align_report.json").write_text("{not json", encoding="utf-8")
    elif damage == "match_table":
        with (align_dir / "match_table.csv").open("a", encoding="utf-8") as handle:
            handle.write("x,x,0,0,0,0,matched\n")
    else:
        with mov.open("a", encoding="utf-8") as handle:
            handle.write("\n# edited after align\n")
    provenance = _run(ref, mov, tmp_path / "run")["alignment_provenance"]
    assert provenance["attached"] is False
    assert provenance["reason"]
    assert "alignment provenance not attached" in capsys.readouterr().err


def test_row_aligned_without_any_report_runs_and_records_nothing(tmp_path) -> None:
    ref = tmp_path / "a.star"
    mov = tmp_path / "b.star"
    shutil.copy(F / "job022_run_data_subset.star", ref)
    shutil.copy(F / "job043_coords_exchange_job042_subset.star", mov)
    summary = _run(ref, mov, tmp_path / "run")
    assert summary["alignment_provenance"] is None
