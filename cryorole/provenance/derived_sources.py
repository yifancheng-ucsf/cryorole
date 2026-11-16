"""Input paths consumed by derived canonical/display artifacts."""
from __future__ import annotations

import json
from pathlib import Path

from cryorole.provenance.output_paths import protect_input_directory, protect_input_file


def protect_derived_directory(output: str | Path, inputs) -> None:
    """Check an in-place writer's existing directories and file aliases."""
    sources = tuple(Path(value) for value in inputs if value is not None)
    protect_input_directory(output, sources)
    for entry in Path(output).rglob("*"):
        if entry.is_dir():
            # rglob does not descend through directory symlinks; inspect their
            # resolved destination before a nested writer can follow them.
            protect_input_directory(entry, sources)
        elif entry.is_file() or entry.is_symlink():
            for source in sources:
                if source.is_file():
                    protect_input_file(entry, source)


def frame_input_paths(value: str | Path | None) -> tuple[Path, ...]:
    """Retain lexical frame names and any referenced transform payload.

    Geometry/schema validation remains with the frame reader. Directory frames
    use the diagnostic reader's NPZ-then-JSON preference.
    """
    if value is None:
        return ()
    current = Path(value)
    paths: list[Path] = []
    seen: set[Path] = set()
    while True:
        resolved = current.resolve()
        if resolved in seen:
            raise ValueError("Canonical frame references contain a cycle")
        seen.add(resolved)
        paths.append(current)
        if current.is_dir():
            candidate = next((current / name for name in ("canonical_frame.npz", "canonical_frame.json")
                              if (current / name).is_file()), None)
        elif current.suffix.lower() == ".json" and current.is_file():
            with current.open(encoding="utf-8") as stream:
                payload = json.load(stream)
            reference = payload.get("canonical_frame_npz") if isinstance(payload, dict) else None
            candidate = None
            if isinstance(payload, dict) and "canonical_transform" not in payload and reference:
                candidate = Path(str(reference))
                if not candidate.is_absolute():
                    candidate = current.parent / candidate
        else:
            candidate = None
        if candidate is None:
            return tuple(paths)
        current = candidate
