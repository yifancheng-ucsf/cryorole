"""Optional real-ChimeraX validation.

Set CRYOROLE_CHIMERAX_BIN and CRYOROLE_CHIMERAX_TEST_SESSION to a prepared
session containing asymmetric reference/moving models #1 and #2. The moving
model should be visibly asymmetric so a +90-degree one-axis target detects
direction and multiplication-order errors.
"""

from __future__ import annotations

import csv
import json
import os

import numpy as np
import pytest

from cryorole.cli.main import main


@pytest.mark.chimerax
def test_real_chimerax_known_rotation_and_pivot(tmp_path):
    executable = os.environ.get("CRYOROLE_CHIMERAX_BIN")
    session = os.environ.get("CRYOROLE_CHIMERAX_TEST_SESSION")
    if not executable or not session:
        pytest.skip(
            "Set CRYOROLE_CHIMERAX_BIN and CRYOROLE_CHIMERAX_TEST_SESSION "
            "to run the asymmetric real-ChimeraX validation"
        )
    run_dir = tmp_path / "run"
    (run_dir / "data").mkdir(parents=True)
    np.savez(
        run_dir / "data" / "raw_landscape.npz",
        artifact_type=np.asarray("raw_landscape"),
        schema_version=np.asarray("1"),
        particle_key=np.asarray(["a", "b"]),
        coordinates_analysis=np.asarray([[0, 0, 0], [0, 0, np.pi / 2]]),
        sld_display=np.asarray([1.0, 2.0]),
        sld_display_is_outlier=np.asarray([False, False]),
    )
    (run_dir / "run_summary.json").write_text(
        json.dumps({"euler_convention": "extrinsic_zyx"}),
        encoding="utf-8",
    )
    waypoints = tmp_path / "waypoints.csv"
    with waypoints.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["label", "rv_x_rad", "rv_y_rad", "rv_z_rad"])
        writer.writerow(["identity", 0, 0, 0])
        writer.writerow(["z90", 0, 0, np.pi / 2])
    output = tmp_path / "animation"
    assert (
        main(
            [
                "animate",
                "--run-dir",
                str(run_dir),
                "--path-csv",
                str(waypoints),
                "--path-space",
                "rv",
                "--chimerax-session",
                session,
                "--reference-model-id",
                "#1",
                "--moving-model-id",
                "#2",
                "--pivot",
                "0",
                "0",
                "0",
                "--baseline-ro",
                "identity",
                "--map-frame",
                "raw",
                "--frames-per-segment",
                "2",
                "--render-mode",
                "execute",
                "--chimerax-bin",
                executable,
                "--no-encode",
                "--output-dir",
                str(output),
            ]
        )
        == 0
    )
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "composite_rendered"
    assert manifest["structure_validation"]["frame_count"] == 2
    assert manifest["compositor"]["frame_count"] == 2
