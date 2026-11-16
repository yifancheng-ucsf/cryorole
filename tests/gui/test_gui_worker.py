"""The GUI worker child process (no Qt): protocol, results, errors, cancel."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run"))
from test_core_productionization import _write_cs  # noqa: E402

N = 300


def _inputs(tmp_path):
    rng = np.random.default_rng(1)
    uids = list(range(10, 10 + N))
    ref = rng.normal(0.0, 0.5, (N, 3))
    mov = (Rotation.from_rotvec(ref) * Rotation.from_rotvec(rng.normal(0.0, 0.1, (N, 3)))).as_rotvec()
    _write_cs(tmp_path / "ref.cs", uids, ref)
    _write_cs(tmp_path / "mov.cs", uids, mov)
    return tmp_path / "ref.cs", tmp_path / "mov.cs"


def _worker(request: dict, *, extra_stdin: str = "") -> list[dict]:
    completed = subprocess.run(
        [sys.executable, "-m", "cryorole.gui.worker"], input=json.dumps(request) + "\n" + extra_stdin,
        capture_output=True, text=True, timeout=300, check=False,
    )
    return [json.loads(line[len("CRYOROLE "):]) for line in completed.stdout.splitlines()
            if line.startswith("CRYOROLE ")]


def test_run_streams_stages_then_the_result(tmp_path):
    ref, mov = _inputs(tmp_path)
    messages = _worker({"command": "run", "ref": str(ref), "mov": str(mov), "output_dir": str(tmp_path / "run"),
                        "k_neighbors": 10, "no_visualize": True})
    starts = [m for m in messages if m["type"] == "event" and m["kind"] == "stage_start"]
    assert starts and starts[0]["stage_number"] == 1 and starts[0]["stage_total"] >= 2
    assert messages[-1] == {"type": "result", "output_dir": str((tmp_path / "run").resolve()), "matched_count": N}


def test_errors_come_back_as_messages(tmp_path):
    ref, mov = _inputs(tmp_path)
    (tmp_path / "run").mkdir()
    messages = _worker({"command": "run", "ref": str(ref), "mov": str(mov), "output_dir": str(tmp_path / "run")})
    assert messages[-1]["type"] == "error" and messages[-1]["code"] == "output_exists"


def test_cancel_publishes_nothing_and_leaves_no_staging_folder(tmp_path):
    ref, mov = _inputs(tmp_path)
    messages = _worker({"command": "run", "ref": str(ref), "mov": str(mov), "output_dir": str(tmp_path / "run"),
                        "k_neighbors": 10, "no_visualize": True}, extra_stdin="cancel\n")
    assert messages[-1]["type"] == "cancelled"
    assert not (tmp_path / "run").exists()
    assert not list(tmp_path.glob(".run.preparing-*"))


def test_canonicalize_command(tmp_path):
    ref, mov = _inputs(tmp_path)
    _worker({"command": "run", "ref": str(ref), "mov": str(mov), "output_dir": str(tmp_path / "run"),
             "k_neighbors": 10, "no_visualize": True})
    messages = _worker({"command": "canonicalize", "run_dir": str(tmp_path / "run"), "canonical_id": "fit30",
                        "fit_top_fraction": 0.3, "no_visualize": True})
    assert messages[-1]["type"] == "result" and messages[-1]["output_dir"].endswith("canonical/fit30")
