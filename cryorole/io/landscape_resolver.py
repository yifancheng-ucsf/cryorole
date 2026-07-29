"""Dependency-light landscape artifact path resolution."""

from __future__ import annotations

from pathlib import Path


def resolve_landscape_path(
    path_or_run_dir: str | Path,
    *,
    space: str = "raw",
    canonical_id: str | None = None,
) -> Path:
    """Resolve run-dir landscape defaults or return an explicit artifact path."""

    path = Path(path_or_run_dir)
    if path.is_file():
        return path
    if not path.is_dir():
        raise ValueError(f"Landscape path or run directory does not exist: {path}")
    if space not in {"raw", "canonical"}:
        raise ValueError("space must be 'raw' or 'canonical'")
    if space == "raw":
        return _first_existing(
            (
                path / "data" / "raw_landscape.npz",
                path / "data" / "raw_landscape.csv",
                path / "debug" / "landscape_debug.json",
                path / "landscape.json",
            )
        )
    canonical = canonical_id or "default"
    return _first_existing(
        (
            path / "canonical" / canonical / "canonical_landscape.npz",
            path / "canonical" / canonical / "canonical_landscape.csv",
            path / "canonical_landscape.npz",
            path / "canonical_landscape.csv",
            path / "debug" / f"canonical_{canonical}_landscape_debug.json",
            path / "canonical_landscape.json",
        )
    )


def _first_existing(candidates: tuple[Path, ...]) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    formatted = ", ".join(str(candidate) for candidate in candidates)
    raise ValueError(f"No supported landscape artifact found; checked: {formatted}")
