"""Resolve omitted ``--run-dir``, ``--canonical-id`` and ``--selection-id``.

A value is filled in only when exactly one candidate exists; otherwise the
caller gets a ``CryoroleError`` that lists the candidates. Every resolution
records how it was made (``resolved_by``) so reports can show it:

``explicit``
    the user gave the value.
``current_directory``
    ``--run-dir`` omitted and the current directory is a run bundle.
``default_output_dir``
    ``--run-dir`` omitted and ``./cryorole_outputs`` is a run bundle (the
    default ``run --output-dir``).
``only_candidate``
    ``--canonical-id`` / ``--selection-id`` omitted and the bundle has exactly
    one canonical frame / selection.

Nothing here reads scientific data or writes files.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from cryorole.errors import CryoroleError

DEFAULT_RUN_OUTPUT_DIR = "cryorole_outputs"
RUN_BUNDLE_MARKER = "run_manifest.json"


@dataclass(frozen=True)
class Resolved:
    """A resolved option value and how it was obtained."""

    option: str
    value: str
    resolved_by: str

    @property
    def implicit(self) -> bool:
        return self.resolved_by != "explicit"

    def notice(self) -> str:
        how = {
            "current_directory": "the current directory is a run bundle",
            "default_output_dir": f"./{DEFAULT_RUN_OUTPUT_DIR} is the default run output",
            "only_candidate": "it is the only one in this run bundle",
        }.get(self.resolved_by, self.resolved_by)
        return f"using {self.option} {self.value} ({how})"

    def record(self) -> dict[str, str]:
        return {"value": self.value, "resolved_by": self.resolved_by}


def is_run_bundle(path: str | Path) -> bool:
    return (Path(path) / RUN_BUNDLE_MARKER).is_file()


def resolve_run_dir(explicit: str | Path | None, *, cwd: str | Path | None = None) -> Resolved:
    """Return the run bundle to use; never guesses between several bundles."""

    if explicit:
        return Resolved("--run-dir", str(explicit), "explicit")
    here = Path(cwd) if cwd is not None else Path.cwd()
    if is_run_bundle(here):
        return Resolved("--run-dir", ".", "current_directory")
    default = here / DEFAULT_RUN_OUTPUT_DIR
    if is_run_bundle(default):
        return Resolved("--run-dir", DEFAULT_RUN_OUTPUT_DIR, "default_output_dir")
    candidates = sorted(child.name for child in here.iterdir() if child.is_dir() and is_run_bundle(child)) if here.is_dir() else []
    if candidates:
        listed = ", ".join(candidates[:10]) + (" …" if len(candidates) > 10 else "")
        raise CryoroleError(
            "run_dir_unresolved",
            f"No --run-dir given and the current directory is not a run bundle; bundles here: {listed}.",
            f"Add --run-dir {candidates[0]} (or another bundle).",
            details={"candidates": candidates},
        )
    raise CryoroleError(
        "run_dir_unresolved",
        f"No --run-dir given, and neither the current directory nor ./{DEFAULT_RUN_OUTPUT_DIR} is a run bundle.",
        "Add --run-dir RUN, or create one with `cryorole run --ref REF --mov MOV --output-dir RUN`.",
        details={"candidates": []},
    )


def canonical_frame_ids(run_dir: str | Path) -> list[str]:
    root = Path(run_dir) / "canonical"
    if not root.is_dir():
        return []
    return sorted(
        child.name
        for child in root.iterdir()
        if child.is_dir() and ((child / "canonical_frame.json").is_file() or (child / "canonical_landscape.npz").is_file())
    )


def selection_ids(run_dir: str | Path) -> list[str]:
    root = Path(run_dir) / "selections"
    if not root.is_dir():
        return []
    return sorted(child.name for child in root.iterdir() if child.is_dir() and (child / "selection.json").is_file())


def resolve_canonical_id(run_dir: str | Path, explicit: str | None) -> Resolved:
    """Return the canonical frame to read; never guesses between several frames."""

    if explicit:
        return Resolved("--canonical-id", str(explicit), "explicit")
    frames = canonical_frame_ids(run_dir)
    if len(frames) == 1:
        return Resolved("--canonical-id", frames[0], "only_candidate")
    if not frames:
        raise CryoroleError(
            "id_unresolved",
            f"The run bundle {run_dir} has no canonical frame yet.",
            f"Run `cryorole canonicalize --run-dir {run_dir}` first, or use --space raw.",
            details={"candidates": []},
        )
    raise CryoroleError(
        "id_unresolved",
        f"The run bundle has {len(frames)} canonical frames: {', '.join(frames)}.",
        f"Add --canonical-id {frames[0]} (or another of them).",
        details={"candidates": frames},
    )


def resolve_selection_id(run_dir: str | Path, explicit: str | None) -> Resolved:
    """Return the selection to read (export / visualize); never guesses between several."""

    if explicit:
        return Resolved("--selection-id", str(explicit), "explicit")
    found = selection_ids(run_dir)
    if len(found) == 1:
        return Resolved("--selection-id", found[0], "only_candidate")
    if not found:
        raise CryoroleError(
            "id_unresolved",
            f"The run bundle {run_dir} has no selections yet.",
            f"Create one with `cryorole select --run-dir {run_dir} --selection-id NAME ...` first.",
            details={"candidates": []},
        )
    raise CryoroleError(
        "id_unresolved",
        f"The run bundle has {len(found)} selections: {', '.join(found[:10])}{' …' if len(found) > 10 else ''}.",
        f"Add --selection-id {found[0]} (or another of them).",
        details={"candidates": found},
    )
