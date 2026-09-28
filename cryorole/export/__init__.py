"""Export helpers for cryoROLE derived artifacts.

Names are loaded lazily (PEP 562) so that ``cryorole.export.landscape`` can be
imported by ``cryorole.io.writers.landscape_store`` without importing the
rest of the package first; this removes the old
``landscape_store`` <-> ``export.metadata_subset`` import cycle.
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING, Any

_EXPORTS: dict[str, str] = {
    "read_landscape_json": "cryorole.export.landscape",
    "read_landscape_json_metadata": "cryorole.export.landscape",
    "write_canonical_landscape_csv": "cryorole.export.landscape",
    "write_json_artifact": "cryorole.export.landscape",
    "write_landscape_json": "cryorole.export.landscape",
    "write_raw_landscape_csv": "cryorole.export.landscape",
    "write_report_json": "cryorole.export.landscape",
    "build_run_manifest_payload": "cryorole.export.manifest",
    "write_run_manifest": "cryorole.export.manifest",
    "export_selection_metadata_subset": "cryorole.export.metadata_subset",
    "export_selection": "cryorole.export.selection_export",
    "read_selection_json": "cryorole.export.selection_export",
    "write_landscape_visualizations": "cryorole.export.visualization",
    "read_landscape": "cryorole.io.writers.landscape_store",
    "read_landscape_csv": "cryorole.io.writers.landscape_store",
    "read_landscape_metadata": "cryorole.io.writers.landscape_store",
    "read_landscape_npz": "cryorole.io.writers.landscape_store",
    "read_landscape_npz_arrays": "cryorole.io.writers.landscape_store",
    "landscape_from_arrays": "cryorole.io.writers.landscape_store",
    "resolve_landscape_path": "cryorole.io.writers.landscape_store",
    "write_canonical_landscape_csv_from_arrays": "cryorole.io.writers.landscape_store",
    "write_canonical_landscape_csv_from_npz": "cryorole.io.writers.landscape_store",
    "write_landscape_npz": "cryorole.io.writers.landscape_store",
    "write_landscape_npz_arrays": "cryorole.io.writers.landscape_store",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str) -> Any:
    module_name = _EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(importlib.import_module(module_name), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(_EXPORTS))


if TYPE_CHECKING:  # pragma: no cover - static analysers see the real names
    from cryorole.export.landscape import (  # noqa: F401
        read_landscape_json,
        read_landscape_json_metadata,
        write_canonical_landscape_csv,
        write_json_artifact,
        write_landscape_json,
        write_raw_landscape_csv,
        write_report_json,
    )
    from cryorole.export.manifest import (  # noqa: F401
        build_run_manifest_payload,
        write_run_manifest,
    )
    from cryorole.export.metadata_subset import (  # noqa: F401
        export_selection_metadata_subset,
    )
    from cryorole.export.selection_export import (  # noqa: F401
        export_selection,
        read_selection_json,
    )
    from cryorole.export.visualization import (  # noqa: F401
        write_landscape_visualizations,
    )
    from cryorole.io.writers.landscape_store import (  # noqa: F401
        read_landscape,
        read_landscape_csv,
        read_landscape_metadata,
        read_landscape_npz,
        read_landscape_npz_arrays,
        landscape_from_arrays,
        resolve_landscape_path,
        write_canonical_landscape_csv_from_arrays,
        write_canonical_landscape_csv_from_npz,
        write_landscape_npz,
        write_landscape_npz_arrays,
    )
