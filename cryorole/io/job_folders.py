"""Find the particle metadata file in a CryoSPARC or RELION job folder (for the GUI Inputs page).

A file is taken as given. In a folder, a file is chosen only when there is
exactly one candidate (D4: locations are filled in only when unambiguous; the
choice and the reason are reported, and the alternatives are listed):

* CryoSPARC job folder: the job's particle result groups (``*.csg``) that have
  an ``alignments3D`` slot; the chosen file is that slot's ``.cs``. Without any
  ``.csg`` (older exports), the ``J*_particles.cs`` files that are not
  passthrough files are the candidates; iteration numbers are never guessed.
* RELION job folder: ``run_data.star`` (the job's result); otherwise the
  ``*_data.star`` files and ``particles.star`` are the candidates.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cryorole.errors import CryoroleError
from cryorole.export.cs_group import read_cs_group

CRYOSPARC_SUFFIX = ".cs"
RELION_SUFFIX = ".star"


@dataclass(frozen=True)
class InputCandidate:
    path: Path
    software: str  # "cryosparc" or "relion"
    reason: str
    job: str | None = None
    group: Path | None = None  # CryoSPARC result group naming this file
    particle_count: int | None = None


@dataclass(frozen=True)
class InputInspection:
    given: Path
    chosen: InputCandidate | None
    candidates: tuple[InputCandidate, ...]
    note: str | None = None


def inspect_input(path: str | Path) -> InputInspection:
    """Which particle metadata file to use for ``path`` (a .cs/.star file or a job folder)."""

    given = Path(path).expanduser()
    if not given.exists():
        raise CryoroleError("input_not_found", f"Input not found: {given}", "Check the path.")
    if given.is_file():
        candidate = _file_candidate(given, reason="the file you chose")
        return InputInspection(given=given, chosen=candidate, candidates=(candidate,))
    candidates = _cryosparc_candidates(given) + _relion_candidates(given)
    if not candidates:
        raise CryoroleError(
            "input_invalid",
            f"No CryoSPARC or RELION particle output found in {given}.",
            "Choose a refinement job folder, or its particle .cs / .star file.",
        )
    if len({candidate.software for candidate in candidates}) > 1:  # never prefer one program's files
        return InputInspection(
            given=given,
            chosen=None,
            candidates=tuple(candidates),
            note=f"{given.name} holds both CryoSPARC and RELION particle outputs; choose one.",
        )
    preferred = [candidate for candidate in candidates if candidate.reason.startswith("RELION job result")]
    if preferred:
        return InputInspection(given=given, chosen=preferred[0], candidates=tuple(candidates))
    if len(candidates) == 1:
        return InputInspection(given=given, chosen=candidates[0], candidates=tuple(candidates))
    return InputInspection(
        given=given,
        chosen=None,
        candidates=tuple(candidates),
        note=f"{len(candidates)} particle outputs in {given.name}; choose one.",
    )


def _file_candidate(path: Path, *, reason: str) -> InputCandidate:
    suffix = path.suffix.lower()
    if suffix == CRYOSPARC_SUFFIX:
        groups = [group for group in _particle_groups(path.parent) if path.name in group.metafile_names]
        return InputCandidate(
            path=path,
            software="cryosparc",
            reason=reason,
            job=path.parent.name,
            group=groups[0].path if len(groups) == 1 else None,
            particle_count=_cs_rows(path),
        )
    if suffix == RELION_SUFFIX:
        return InputCandidate(path=path, software="relion", reason=reason, job=path.parent.name)
    raise CryoroleError(
        "input_invalid",
        f"{path.name} is not a CryoSPARC .cs or RELION .star file.",
        "Choose a particle .cs or .star file, or a job folder.",
    )


def _particle_groups(folder: Path):
    groups = []
    for csg in sorted(folder.glob("*.csg")):
        try:
            group = read_cs_group(csg)
        except (ValueError, OSError):
            continue
        if group.group.get("type") == "particle" and "alignments3D" in group.results:
            groups.append(group)
    return groups


def _cryosparc_candidates(folder: Path) -> list[InputCandidate]:
    candidates = []
    for group in _particle_groups(folder):
        name = group.metafile_name("alignments3D")
        path = folder / name
        if path.is_file():
            candidates.append(InputCandidate(
                path=path,
                software="cryosparc",
                reason=f"particle group {group.path.name} ({group.group.get('title', group.group.get('name', ''))})",
                job=folder.name,
                group=group.path,
                particle_count=int(group.results["alignments3D"]["num_items"]),
            ))
    if candidates:
        return candidates
    for path in sorted(folder.glob("*_particles.cs")):
        if "passthrough" not in path.name:
            candidates.append(InputCandidate(
                path=path, software="cryosparc", reason="particle file (no .csg in the folder)",
                job=folder.name, particle_count=_cs_rows(path),
            ))
    return candidates


def _relion_candidates(folder: Path) -> list[InputCandidate]:
    result = folder / "run_data.star"
    if result.is_file():
        return [InputCandidate(path=result, software="relion", reason="RELION job result run_data.star",
                               job=folder.name)]
    names = sorted({*folder.glob("*_data.star"), *folder.glob("particles*.star")})
    return [InputCandidate(path=path, software="relion", reason="particle STAR file", job=folder.name)
            for path in names if path.is_file()]


def _cs_rows(path: Path) -> int | None:
    try:
        return int(len(np.load(path, mmap_mode="r", allow_pickle=False)))
    except (OSError, ValueError):
        return None
