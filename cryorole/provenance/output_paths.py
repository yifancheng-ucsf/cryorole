"""Keep output publication separate from original input paths and aliases."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable


def _input_aliases(source: str | Path) -> set[Path]:
    """Keep names reached through parent links before resolving the final link."""
    path = Path(source).expanduser().absolute()
    aliases = {path, path.resolve()}
    # A path may enter the output through one parent link and leave through
    # another. Full resolve alone loses the intermediate input entry names.
    for parent in path.parents:
        aliases.add(parent.resolve() / path.relative_to(parent))
    return aliases


def protect_input_directory(output_dir: str | Path, inputs: Iterable[str | Path | None]) -> None:
    """Refuse replacing a directory containing an input name or its resolved target."""
    roots = {Path(output_dir).expanduser().absolute(), Path(output_dir).expanduser().resolve()}
    for source in inputs:
        if source is None:
            continue
        candidates = _input_aliases(source)
        if any(candidate.is_relative_to(root) for candidate in candidates for root in roots):
            raise ValueError(f'Output directory contains an input source: {source}. Choose a separate output directory.')


def protect_input_file(output_path: str | Path, input_path: str | Path) -> None:
    """Refuse a source file, symlink alias or hard-link alias as its own output."""
    output, source = Path(output_path).expanduser(), Path(input_path).expanduser()
    same = output.resolve() == source.resolve()
    if not same and output.exists() and source.exists():
        same = output.samefile(source)
    if same:
        raise ValueError('Output file would overwrite its input source; choose a separate output path.')
