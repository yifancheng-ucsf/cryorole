"""Read recentring parameters recorded by RELION job directories.

RELION writes the executed command line to ``<job>/note.txt`` and (3.1+) the
job options to ``<job>/job.star``. cryoROLE only *reads* these files to
suggest or fill in parameters; every value is reported with its source and is
verified against the particle data before it is used for anything.
"""

from __future__ import annotations

import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

EXTRACT_PROGRAMS = ("relion_preprocess", "relion_preprocess_mpi")
SUBTRACT_PROGRAMS = ("relion_particle_subtract", "relion_particle_subtract_mpi")


@dataclass(frozen=True)
class RelionJobCommand:
    """The last command recorded in a RELION ``note.txt``."""

    source: str
    program: str
    options: dict[str, Any] = field(default_factory=dict)

    def get(self, name: str, default: Any = None) -> Any:
        return self.options.get(name, default)

    def vector(self, prefix: str) -> tuple[float, float, float] | None:
        keys = (f"{prefix}_x", f"{prefix}_y", f"{prefix}_z")
        if not all(key in self.options for key in keys):
            return None
        return tuple(float(self.options[key]) for key in keys)  # type: ignore[return-value]


def _tokens_of_command(text: str) -> list[str]:
    cleaned = text.replace("`which ", "").replace("`", " ")
    try:
        return shlex.split(cleaned)
    except ValueError:
        return cleaned.split()


def parse_note_txt(path: str | Path) -> RelionJobCommand | None:
    """Parse the last executed RELION command in ``note.txt`` (``None`` if absent)."""

    note = Path(path)
    if not note.is_file():
        return None
    commands: list[str] = []
    for line in note.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if "relion_" in stripped and not stripped.startswith("++++"):
            commands.append(stripped)
    if not commands:
        return None
    tokens = _tokens_of_command(commands[-1])
    program = next((token for token in tokens if Path(token).name.startswith("relion_")), "")
    options: dict[str, Any] = {}
    index = tokens.index(program) + 1 if program in tokens else 0
    while index < len(tokens):
        token = tokens[index]
        if token.startswith("--"):
            name = token[2:]
            if index + 1 < len(tokens) and not tokens[index + 1].startswith("--"):
                options[name] = tokens[index + 1]
                index += 2
                continue
            options[name] = True
        index += 1
    return RelionJobCommand(source=str(note), program=Path(program).name, options=options)


def find_job_command(job_dir: str | Path) -> RelionJobCommand | None:
    """Return the recorded command of a RELION job directory, if any."""

    return parse_note_txt(Path(job_dir) / "note.txt")


def project_dir_for(job_dir: str | Path) -> Path:
    """RELION project directory: the parent of ``<JobType>/jobNNN``."""

    return Path(job_dir).resolve().parent.parent


def resolve_project_path(value: str, job_dir: str | Path) -> Path:
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    return project_dir_for(job_dir) / candidate


def extraction_parameters(command: RelionJobCommand) -> dict[str, Any] | None:
    """``--recenter`` vector and ``--ref_angpix`` of a re-extraction, or ``None``."""

    if command.program not in EXTRACT_PROGRAMS or not command.get("recenter"):
        return None
    vector = command.vector("recenter")
    if vector is None:
        return None
    ref_angpix = command.get("ref_angpix")
    return {
        "recenter_shift_px": list(vector),
        "ref_angpix": float(ref_angpix) if ref_angpix not in (None, True) and float(ref_angpix) > 0 else None,
        "input_star": command.get("reextract_data_star"),
        "output_star": command.get("part_star"),
        "source": command.source,
    }


def subtraction_parameters(command: RelionJobCommand) -> dict[str, Any] | None:
    """``--center_x/y/z`` and inputs of a particle subtraction, or ``None``."""

    if command.program not in SUBTRACT_PROGRAMS:
        return None
    vector = command.vector("center")
    return {
        "center_px": list(vector) if vector is not None and all(abs(v) < 9999 for v in vector) else None,
        "recenter_on_mask": bool(command.get("recenter_on_mask", False)),
        "optimiser": command.get("i"),
        "output_dir": command.get("o"),
        "new_box": command.get("new_box"),
        "source": command.source,
    }
