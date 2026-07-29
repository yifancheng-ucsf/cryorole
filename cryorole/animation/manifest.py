"""Animation-bundle manifest persistence."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


ANIMATION_MANIFEST_SCHEMA_VERSION = "5"


def write_animation_manifest(payload: dict[str, Any], path: str | Path) -> Path:
    """Atomically replace the manifest inside one animation output directory."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default),
        encoding="utf-8",
    )
    temporary.replace(output)
    return output


def _json_default(value: object):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")
