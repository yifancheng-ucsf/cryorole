"""CryoSPARC result-group (.csg) export: mirror the source job's group so Import Result Group works.

The layout follows a real CryoSPARC v5 refinement (the J75 example): the job's
particle file holds alignments3D and ctf, a passthrough file holds blob, and the
job's .csg names both. A real subset of J75 exported this way imported in CryoSPARC.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cryorole import api
from cryorole.cli.main import main
from cryorole.export.cs_group import format_cs_group, read_cs_group

N = 240
PARTICLE_CSG = """created: 2026-09-09 21:21:10.415270+00:00
group:
  description: All particles that were processed, including alignments
  name: particles
  title: All particles
  type: particle
results:
  alignments3D:
    metafile: '>{job}_004_particles.cs'
    num_items: {n}
    type: particle.alignments3D
  blob:
    metafile: '>{job}_passthrough_particles.cs'
    num_items: {n}
    type: particle.blob
  ctf:
    metafile: '>{job}_004_particles.cs'
    num_items: {n}
    type: particle.ctf
version: v5.0.6
"""


def _particles(uids, poses):
    dtype = np.dtype([("uid", "<u8"), ("blob/path", "S40"), ("blob/idx", "<u4"),
                      ("ctf/df1_A", "<f4"), ("alignments3D/pose", "<f4", (3,))])
    values = np.zeros(len(uids), dtype=dtype)
    values["uid"] = uids
    values["blob/path"] = [f"J1/imported/stack_{u % 7}.mrcs".encode() for u in uids]
    values["blob/idx"] = np.asarray(uids) % 100
    values["ctf/df1_A"] = 10000.0 + np.asarray(uids, dtype=float)
    values["alignments3D/pose"] = poses
    return values


def _save(path: Path, array: np.ndarray) -> None:
    with path.open("wb") as handle:
        np.save(handle, array)


def _job(project: Path, job: str, uids, poses, *, passthrough_order=None, csg=True) -> Path:
    folder = project / job
    folder.mkdir(parents=True)
    main_file = folder / f"{job}_004_particles.cs"
    particles = _particles(uids, poses)
    _save(main_file, particles)
    passthrough = particles[["uid", "blob/path", "blob/idx"]]
    if passthrough_order is not None:
        passthrough = passthrough[passthrough_order]
    _save(folder / f"{job}_passthrough_particles.cs", np.array(passthrough.tolist(), dtype=passthrough.dtype))
    if csg:
        (folder / f"{job}_particles.csg").write_text(PARTICLE_CSG.format(job=job, n=len(uids)), encoding="utf-8")
    return main_file


@pytest.fixture
def project(tmp_path):
    return build_project(tmp_path)


def build_project(tmp_path: Path) -> dict[str, Path]:
    """A CryoSPARC project with two jobs (J5 with a .csg, J6 without), a run on them and a selection."""

    root = tmp_path / "CS-project"
    root.mkdir()
    (root / "project.json").write_text(json.dumps({"uid": "P9", "project_dir": str(root)}), encoding="utf-8")
    rng = np.random.default_rng(3)
    uids = rng.choice(np.arange(10**6, 10**6 + 10 * N), N, replace=False).astype(np.uint64)
    ref = rng.normal(0.0, 0.5, (N, 3))
    mov = (Rotation.from_rotvec(ref) * Rotation.from_rotvec(rng.normal(0.0, 0.1, (N, 3)))).as_rotvec()
    ref_file = _job(root, "J5", uids, ref, passthrough_order=rng.permutation(N))
    mov_file = _job(root, "J6", uids, mov, csg=False)
    run = tmp_path / "run"
    api.run(ref_file, mov_file, output_dir=run, k_neighbors=10, no_visualize=True)
    api.select(run, selection_id="dense", mode="threshold", sld_min=1.0)
    return {"root": root, "run": run, "ref": ref_file, "mov": mov_file}


def test_export_mirrors_the_source_jobs_group_for_import(project, tmp_path):
    report = api.export(project["run"], selection_id="dense", domain="ref")
    domain = json.loads(Path(report["domain_reports"]["ref"]).read_text(encoding="utf-8"))
    group = domain["cryosparc_group"]
    assert group["resolved_by"] == "only_candidate"
    assert group["source_group"] == str(project["ref"].parent / "J5_particles.csg")
    out = Path(domain["output_path"]).parent
    main_subset = np.load(out / "selected_ref.cs")
    passthrough = np.load(out / "selected_ref_J5_passthrough_particles.cs")
    assert np.array_equal(passthrough["uid"], main_subset["uid"])  # same order despite the shuffled source
    assert np.array_equal(passthrough["blob/path"], main_subset["blob/path"])
    written = read_cs_group(out / "selected_ref.csg")
    assert {slot: written.metafile_name(slot) for slot in written.results} == {
        "alignments3D": "selected_ref.cs", "ctf": "selected_ref.cs",
        "blob": "selected_ref_J5_passthrough_particles.cs",
    }
    assert {entry["num_items"] for entry in written.results.values()} == {len(main_subset)}
    assert written.version == "v5.0.6" and written.group["type"] == "particle"
    assert "dense" in written.group["title"]
    steps = "\n".join(report["next_steps"]["ref"])
    destination = project["root"] / "exports" / "groups" / "cryorole_J5_dense_ref"
    assert f"mkdir -p {destination}" in steps
    assert str(destination / "selected_ref.csg") in steps
    assert (Path(report["output_dir"]) / "NEXT_STEPS.txt").is_file()


def test_explicit_group_elsewhere_and_missing_group_are_reported(project, tmp_path):
    elsewhere = tmp_path / "downloaded"
    elsewhere.mkdir()
    shutil.copy(project["ref"].parent / "J5_particles.csg", elsewhere / "J5_particles.csg")
    shutil.copy(project["ref"].parent / "J5_passthrough_particles.cs", elsewhere / "J5_passthrough_particles.cs")
    report = api.export(project["run"], selection_id="dense", domain="both", output_dir=tmp_path / "out",
                        cryosparc_group_ref=elsewhere / "J5_particles.csg")
    ref = json.loads(Path(report["domain_reports"]["ref"]).read_text(encoding="utf-8"))
    mov = json.loads(Path(report["domain_reports"]["mov"]).read_text(encoding="utf-8"))
    assert ref["cryosparc_group"]["resolved_by"] == "explicit"
    assert ref["cryosparc_group"]["project_dir"] is None  # not inside a project: generic copy advice
    assert "inside your CryoSPARC project directory" in report["next_steps"]["ref"][0]
    assert mov["cryosparc_group"] == {"resolved_by": "none", "candidates": []}
    assert not (tmp_path / "out" / "mov" / "selected_mov.csg").exists()
    assert any("--cryosparc-group-mov" in warning for warning in report["warnings"])


def test_inconsistent_groups_fail_loudly(project, tmp_path):
    folder = project["ref"].parent
    csg = folder / "J5_particles.csg"
    good = PARTICLE_CSG.format(job="J5", n=N)
    csg.write_text(good.replace(f"num_items: {N}", "num_items: 999"), encoding="utf-8")
    with pytest.raises(ValueError, match="not the same file"):
        api.export(project["run"], selection_id="dense", domain="ref", output_dir=tmp_path / "a")
    passthrough = np.load(folder / "J5_passthrough_particles.cs")
    _save(folder / "J5_passthrough_particles.cs", passthrough[1:])
    csg.write_text(good, encoding="utf-8")
    with pytest.raises(ValueError, match="its .csg says"):
        api.export(project["run"], selection_id="dense", domain="ref", output_dir=tmp_path / "b")
    _save(folder / "J5_passthrough_particles.cs", passthrough)
    shutil.copy(csg, folder / "J5_other.csg")
    with pytest.raises(ValueError, match="Several CryoSPARC particle groups"):
        api.export(project["run"], selection_id="dense", domain="ref", output_dir=tmp_path / "c")


def test_missing_uid_in_a_group_file_is_an_error(project, tmp_path):
    folder = project["ref"].parent
    api.select(project["run"], selection_id="everything", mode="random", fraction=1.0, seed=1)
    passthrough = np.load(folder / "J5_passthrough_particles.cs")
    passthrough["uid"][0] = 5  # a uid that no selected particle has
    _save(folder / "J5_passthrough_particles.cs", passthrough)
    with pytest.raises(ValueError, match="1 selected uids are missing"):
        api.export(project["run"], selection_id="everything", domain="ref", output_dir=tmp_path / "x")


def test_preview_reports_sources_and_groups_without_writing(project):
    before = sorted(p.name for p in project["run"].iterdir())
    preview = api.export_preview(project["run"])
    assert preview["ref"]["source_type"] == "cryosparc_cs"
    assert preview["ref"]["cryosparc_group"].endswith("J5_particles.csg")
    assert preview["mov"]["cryosparc_group"] is None
    assert "--cryosparc-group-mov" in preview["mov"]["cryosparc_group_note"]
    assert sorted(p.name for p in project["run"].iterdir()) == before


def test_cli_prints_import_steps_on_stderr(project, tmp_path, capsys):
    assert main(["export", "--run-dir", str(project["run"]), "--selection-id", "dense", "--domain", "ref",
                 "--output-dir", str(tmp_path / "cli")]) == 0
    captured = capsys.readouterr()
    assert captured.out.strip() == str(tmp_path / "cli")
    assert "[ref] next steps:" in captured.err and "Import Result Group" in captured.err


def test_csg_round_trip_and_unsupported_content(tmp_path):
    text = PARTICLE_CSG.format(job="J75", n=1162304)
    path = tmp_path / "J75_particles.csg"
    path.write_text(text, encoding="utf-8")
    group = read_cs_group(path)
    assert format_cs_group(group=group.group, results=group.results, version=group.version,
                           created=group.created) == text
    path.write_text(text.replace("type: particle\n", "type: [particle]\n", 1), encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported .csg value"):
        read_cs_group(path)
    path.write_text(text.replace("'>J75_passthrough_particles.cs'", "'/abs/J75_passthrough_particles.cs'"),
                    encoding="utf-8")
    with pytest.raises(ValueError, match="not a file next to the .csg"):
        read_cs_group(path)


# Regression coverage retained from the October GUI review.

WRAPPED = "  description: A stack of imported particles. May or may not contain data, ctfs, pick\n    locations, etc.\n"


def test_csg_values_wrapped_over_lines_are_read(tmp_path):
    """CryoSPARC's YAML writer folds long values (e.g. Import Particles descriptions) onto lines
    indented deeper than their key; quoted values continue until their closing quote."""
    text = (PARTICLE_CSG.format(job="J7", n=12)
            .replace("  description: All particles that were processed, including alignments\n", WRAPPED)
            .replace("  title: All particles\n", "  title: 'Imported particles: a title quoted\n    over two lines'\n"))
    path = tmp_path / "J7_imported_particles.csg"
    path.write_text(text, encoding="utf-8")
    group = read_cs_group(path)
    assert group.group["description"] == (
        "A stack of imported particles. May or may not contain data, ctfs, pick locations, etc."
    )
    assert group.group["title"] == "Imported particles: a title quoted over two lines"
    assert group.results["blob"]["num_items"] == 12 and group.version == "v5.0.6"
    path.write_text(text.replace("results:\n", "results:\n    stray line\n"), encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported .csg line"):  # not the continuation of a value
        read_cs_group(path)
    # long plain values are written folded the same way, so CryoSPARC's file round-trips byte for byte
    plain = PARTICLE_CSG.format(job="J7", n=12).replace(
        "  description: All particles that were processed, including alignments\n", WRAPPED)
    path.write_text(plain, encoding="utf-8")
    again = read_cs_group(path)
    assert format_cs_group(group=again.group, results=again.results, version=again.version,
                           created=again.created) == plain


def test_wrapped_group_is_found_and_mirrored(project, tmp_path):
    csg = project["ref"].parent / "J5_particles.csg"
    csg.write_text(csg.read_text(encoding="utf-8").replace(
        "  description: All particles that were processed, including alignments\n", WRAPPED), encoding="utf-8")
    chosen = api.inspect_input(project["ref"].parent).chosen
    assert chosen.group == csg and chosen.reason.startswith("particle group J5_particles.csg")
    report = api.export(project["run"], selection_id="dense", domain="ref", output_dir=tmp_path / "out")
    domain = json.loads(Path(report["domain_reports"]["ref"]).read_text(encoding="utf-8"))
    assert domain["cryosparc_group"]["resolved_by"] == "only_candidate"
    assert (tmp_path / "out" / "ref" / "selected_ref.csg").is_file()


def test_group_problems_stop_the_export_before_it_writes(project, tmp_path):
    folder = project["ref"].parent
    shutil.copy(folder / "J5_particles.csg", folder / "J5_other.csg")  # two candidates: never guess
    target = tmp_path / "export_try"
    with pytest.raises(ValueError, match="Several CryoSPARC particle groups"):
        api.export(project["run"], selection_id="dense", domain="ref", output_dir=target)
    assert not target.exists()
    report = api.export(project["run"], selection_id="dense", domain="ref", output_dir=target,  # the remedy works
                        cryosparc_group_ref=folder / "J5_particles.csg")
    assert Path(report["output_dir"]) == target and (target / "ref" / "selected_ref.csg").is_file()
    (folder / "J5_other.csg").unlink()
    (folder / "J5_passthrough_particles.cs").rename(folder / "archived.cs")
    with pytest.raises(ValueError, match="missing next to its .csg.*--no-cryosparc-group"):
        api.export(project["run"], selection_id="dense", domain="ref", output_dir=tmp_path / "b")
    assert not (tmp_path / "b").exists()
    preview = api.export_preview(project["run"])["ref"]
    assert preview["cryosparc_group_usable"] is False and "would stop before writing" in preview["cryosparc_group_note"]


def test_no_cryosparc_group_writes_the_cs_subset_alone(project, tmp_path):
    folder = project["ref"].parent
    (folder / "J5_passthrough_particles.cs").rename(folder / "archived.cs")
    report = api.export(project["run"], selection_id="dense", domain="ref", output_dir=tmp_path / "a",
                        mirror_cryosparc_group=False)
    domain = json.loads(Path(report["domain_reports"]["ref"]).read_text(encoding="utf-8"))
    assert domain["cryosparc_group"] == {"resolved_by": "disabled", "candidates": []}
    assert (tmp_path / "a" / "ref" / "selected_ref.cs").is_file()
    assert not (tmp_path / "a" / "ref" / "selected_ref.csg").exists() and not report["warnings"]
    assert main(["export", "--run-dir", str(project["run"]), "--selection-id", "dense", "--domain", "ref",
                 "--output-dir", str(tmp_path / "cli"), "--no-cryosparc-group"]) == 0
    with pytest.raises(ValueError, match="cannot be combined"):
        api.export(project["run"], selection_id="dense", domain="ref", output_dir=tmp_path / "c",
                   mirror_cryosparc_group=False, cryosparc_group_ref=folder / "J5_particles.csg")


def test_late_failure_removes_only_a_folder_the_export_created(project, tmp_path):
    folder = project["ref"].parent
    api.select(project["run"], selection_id="everything", mode="random", fraction=1.0, seed=1)
    passthrough = np.load(folder / "J5_passthrough_particles.cs")
    passthrough["uid"][0] = 5  # found only while writing: a selected uid is missing from the passthrough
    _save(folder / "J5_passthrough_particles.cs", passthrough)
    with pytest.raises(ValueError, match="selected uids are missing"):
        api.export(project["run"], selection_id="everything", domain="ref", output_dir=tmp_path / "new")
    assert not (tmp_path / "new").exists()
    existing = tmp_path / "existing"
    existing.mkdir()
    (existing / "keep.txt").write_text("user file", encoding="utf-8")
    with pytest.raises(ValueError, match="selected uids are missing"):
        api.export(project["run"], selection_id="everything", domain="ref", output_dir=existing, overwrite=True)
    assert (existing / "keep.txt").read_text(encoding="utf-8") == "user file"
