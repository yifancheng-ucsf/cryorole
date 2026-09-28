"""Regenerate the P0-5 CryoSPARC/RELION convention golden fixture.

Inputs: ``tests/fixtures/j75_j80_subset_2000.npz`` (uid + alignments3D/pose
for a deterministic 2,000-particle subset of the J75/J80 example,
``numpy.random.default_rng(20260924)``).

Outputs: ``tests/fixtures/j75_subset_2000_pyem.star`` and
``j80_subset_2000_pyem.star``, produced by pyem's real ``csparc2star.py``
(https://github.com/asarnow/pyem, commit 22f28768d1129c41293143b70c588358aed3f64a).
The fixture .cs files only need the fields csparc2star.py requires; placeholder
blob/CTF values do not affect the angles.

Usage::

    PYTHONPATH=/path/to/pyem python scripts/make_cryosparc_convention_fixture.py /path/to/pyem
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


def main(pyem_dir: str) -> None:
    data = np.load(FIXTURES / "j75_j80_subset_2000.npz")
    dtype = np.dtype([
        ("uid", "<u8"), ("blob/path", "S64"), ("blob/idx", "<u4"), ("blob/shape", "<u4", (2,)),
        ("blob/psize_A", "<f4"), ("alignments3D/pose", "<f4", (3,)), ("alignments3D/shift", "<f4", (2,)),
        ("alignments3D/psize_A", "<f4"), ("alignments3D/class_posterior", "<f4"),
        ("alignments3D/class", "<u4"), ("alignments3D/split", "<u4"),
    ])
    with tempfile.TemporaryDirectory() as tmp:
        for label, key in (("j75", "ref_pose"), ("j80", "mov_pose")):
            values = np.zeros(len(data["uid"]), dtype=dtype)
            values["uid"] = data["uid"]
            values["alignments3D/pose"] = data[key]
            values["blob/path"] = b"J1/extract/particles.mrc"
            values["blob/idx"] = np.arange(len(values))
            values["blob/shape"] = 256
            values["blob/psize_A"] = 1.0
            values["alignments3D/psize_A"] = 1.0
            values["alignments3D/class_posterior"] = 1.0
            cs_path = Path(tmp) / f"{label}.cs"
            with cs_path.open("wb") as handle:
                np.save(handle, values)
            subprocess.run(
                [sys.executable, str(Path(pyem_dir) / "pyem" / "cli" / "csparc2star.py"), str(cs_path),
                 str(FIXTURES / f"{label}_subset_2000_pyem.star")],
                check=True,
            )


if __name__ == "__main__":
    main(sys.argv[1])
