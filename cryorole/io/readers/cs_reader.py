"""CryoSPARC .cs reader.

The reader loads raw structured arrays into a DataFrame without assigning
scientific interpretation to pose fields.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from cryorole.io.import_report import ImportReport


@dataclass(frozen=True)
class CryoSparcData:
    """Raw CryoSPARC .cs table."""

    path: str
    particles: pd.DataFrame
    report: ImportReport


def read_cryosparc_cs(path: str | Path) -> CryoSparcData:
    """Parse a CryoSPARC .cs file into a raw particle table."""

    cs_path = Path(path)
    if not cs_path.is_file():
        raise FileNotFoundError(f"CryoSPARC .cs file does not exist: {cs_path}")

    array = np.load(cs_path, allow_pickle=False)
    if array.dtype.names:
        columns = {}
        for name in array.dtype.names:
            values = array[name]
            columns[name] = list(values) if getattr(values, "ndim", 1) > 1 else values
        particles = pd.DataFrame(columns)
    else:
        particles = pd.DataFrame(array)

    report = ImportReport(
        path=str(cs_path),
        source_type="cryosparc",
        row_count=len(particles),
        columns=tuple(str(col) for col in particles.columns),
    )
    return CryoSparcData(path=str(cs_path), particles=particles, report=report)
