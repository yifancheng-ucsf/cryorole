"""Particle-file detection in CryoSPARC and RELION job folders (GUI Inputs page)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cryorole import api
from cryorole.errors import CryoroleError

CSG = """created: 2026-09-09 21:21:10+00:00
group:
  description: ''
  name: {name}
  title: {title}
  type: particle
results:
{slots}version: v5.0.6
"""


def _cs(path: Path, n: int) -> Path:
    array = np.zeros(n, dtype=[("uid", "<u8"), ("alignments3D/pose", "<f4", (3,))])
    array["uid"] = np.arange(n)
    with path.open("wb") as handle:
        np.save(handle, array)
    return path


def _csg(path: Path, name: str, slots: dict[str, tuple[str, int]], title: str = "All particles") -> None:
    body = "".join(
        f"  {slot}:\n    metafile: '>{file}'\n    num_items: {n}\n    type: particle.{slot}\n"
        for slot, (file, n) in slots.items()
    )
    path.write_text(CSG.format(name=name, title=title, slots=body), encoding="utf-8")


def test_cryosparc_job_folder_uses_its_particle_group(tmp_path):
    job = tmp_path / "J75"
    job.mkdir()
    _cs(job / "J75_004_particles.cs", 12)
    _cs(job / "J75_passthrough_particles.cs", 12)
    _cs(job / "J75_002_particles.cs", 12)  # an earlier iteration: not named by the group
    _csg(job / "J75_particles.csg", "particles", {"alignments3D": ("J75_004_particles.cs", 12),
                                                  "blob": ("J75_passthrough_particles.cs", 12)})
    _csg(job / "J75_volume.csg", "volume", {"map": ("J75_004_volume_map.cs", 1)})  # not a particle group
    result = api.inspect_input(job)
    assert result.chosen.path == job / "J75_004_particles.cs"
    assert result.chosen.software == "cryosparc" and result.chosen.particle_count == 12
    assert result.chosen.group == job / "J75_particles.csg"
    assert len(result.candidates) == 1


def test_several_particle_groups_are_listed_not_guessed(tmp_path):
    job = tmp_path / "J80"
    job.mkdir()
    for k in (0, 1):
        _cs(job / f"J80_class_{k:02d}_particles.cs", 5)
        _csg(job / f"J80_particles_class_{k}.csg", f"particles_class_{k}",
             {"alignments3D": (f"J80_class_{k:02d}_particles.cs", 5)}, title=f"Class {k}")
    result = api.inspect_input(job)
    assert result.chosen is None and len(result.candidates) == 2
    assert "choose one" in result.note


def test_cryosparc_folder_without_groups_never_guesses_an_iteration(tmp_path):
    job = tmp_path / "J9"
    job.mkdir()
    _cs(job / "J9_002_particles.cs", 4)
    _cs(job / "J9_passthrough_particles.cs", 4)
    single = api.inspect_input(job)
    assert single.chosen.path.name == "J9_002_particles.cs" and single.chosen.group is None
    _cs(job / "J9_004_particles.cs", 4)
    several = api.inspect_input(job)
    assert several.chosen is None
    assert [c.path.name for c in several.candidates] == ["J9_002_particles.cs", "J9_004_particles.cs"]


def test_relion_job_folder_prefers_run_data_star(tmp_path):
    job = tmp_path / "Refine3D" / "job022"
    job.mkdir(parents=True)
    for name in ("run_it010_data.star", "run_it012_data.star"):
        (job / name).write_text("data_\n", encoding="utf-8")
    several = api.inspect_input(job)
    assert several.chosen is None and len(several.candidates) == 2
    (job / "run_data.star").write_text("data_\n", encoding="utf-8")
    final = api.inspect_input(job)
    assert final.chosen.path == job / "run_data.star" and final.chosen.software == "relion"
    extract = tmp_path / "Extract" / "job042"
    extract.mkdir(parents=True)
    (extract / "particles.star").write_text("data_\n", encoding="utf-8")
    assert api.inspect_input(extract).chosen.path.name == "particles.star"


def test_files_and_bad_inputs(tmp_path):
    cs = _cs(tmp_path / "ref.cs", 3)
    chosen = api.inspect_input(cs).chosen
    assert chosen.path == cs and chosen.reason == "the file you chose" and chosen.particle_count == 3
    (tmp_path / "notes.txt").write_text("x", encoding="utf-8")
    with pytest.raises(CryoroleError) as wrong:
        api.inspect_input(tmp_path / "notes.txt")
    assert wrong.value.code == "input_invalid"
    with pytest.raises(CryoroleError) as missing:
        api.inspect_input(tmp_path / "nothing")
    assert missing.value.code == "input_not_found"
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(CryoroleError, match="No CryoSPARC or RELION particle output"):
        api.inspect_input(empty)


def test_folder_with_cryosparc_and_relion_outputs_is_listed_not_guessed(tmp_path):
    job = tmp_path / "mixed"
    job.mkdir()
    _cs(job / "J9_004_particles.cs", 8)
    _csg(job / "J9_particles.csg", "particles", {"alignments3D": ("J9_004_particles.cs", 8)})
    (job / "run_data.star").write_text("data_particles\n", encoding="utf-8")
    result = api.inspect_input(job)
    assert result.chosen is None and "both CryoSPARC and RELION" in result.note
    assert {candidate.software for candidate in result.candidates} == {"cryosparc", "relion"}
