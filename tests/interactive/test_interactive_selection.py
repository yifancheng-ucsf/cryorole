from __future__ import annotations

import json
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

from cryorole.cli.main import build_parser, select_command
from cryorole.interactive import ExploreSession, create_explore_server
from cryorole.models.landscape import Landscape
from cryorole.models.policies import SelectionPolicy
from cryorole.select import select_particles


def _bundle(tmp_path: Path, n: int = 20) -> Path:
    run = tmp_path / "run"
    (run / "data").mkdir(parents=True)
    coordinates = np.zeros((n, 3), dtype=float)
    coordinates[:, 2] = np.linspace(0, 0.5, n)
    np.savez_compressed(
        run / "data" / "raw_landscape.npz",
        artifact_type=np.asarray("raw_landscape"), schema_version=np.asarray("1"),
        particle_key=np.asarray([f"p{i}" for i in range(n)]),
        coordinates_analysis=coordinates, coordinates_display=coordinates,
        sld_unfloored=np.ones(n), sld_raw=np.linspace(1, 5, n),
        sld_display=np.linspace(1, 5, n), sld_display_is_outlier=np.zeros(n, bool),
        sld_was_floored=np.zeros(n, bool), sld_local_k_mean=np.ones(n),
        sld_effective_local_k_mean=np.ones(n), sld_distance_floor=np.full(n, 1e-4),
        ref_source_row_id=np.arange(n), mov_source_row_id=np.arange(n),
    )
    manifest = {"run_id": "run-ux", "bundle_transaction": {"state": "completed"}}
    (run / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    (run / "run_summary.json").write_text(
        json.dumps({"run_id": "run-ux", "euler_convention": "extrinsic_zyx", "scipy_euler_sequence": "zyx"}),
        encoding="utf-8",
    )
    (run / ".cryorole_bundle_complete").write_text("{}", encoding="utf-8")
    return run


def _add_canonical(run: Path, n: int = 20) -> None:
    target = run / "canonical" / "default"
    target.mkdir(parents=True)
    raw = np.zeros((n, 3), dtype=float)
    canonical = np.zeros((n, 3), dtype=float)
    canonical[:, 0] = np.linspace(0, 0.5, n)
    np.savez_compressed(
        target / "canonical_landscape.npz",
        artifact_type=np.asarray("canonical_landscape"), schema_version=np.asarray("1"),
        particle_key=np.asarray([f"p{i}" for i in range(n)]),
        coordinates_analysis=raw, coordinates_display=canonical,
        coordinates_canonical=canonical, canonical_transform=np.eye(3),
        sld_unfloored=np.ones(n), sld_raw=np.ones(n), sld_display=np.ones(n),
        sld_display_is_outlier=np.zeros(n, bool), sld_was_floored=np.zeros(n, bool),
        sld_local_k_mean=np.ones(n), sld_effective_local_k_mean=np.ones(n),
        sld_distance_floor=np.full(n, 1e-4), ref_source_row_id=np.arange(n),
        mov_source_row_id=np.arange(n),
    )


def test_draft_exact_evaluation_matches_cli_core_and_downsample_is_display_only(tmp_path) -> None:
    run = _bundle(tmp_path)
    session = ExploreSession(run, max_display_points=5)
    assert session.session_payload()["displayed_point_count"] == 5
    assert session.session_payload()["full_candidate_count"] == 20
    assert not (run / "selections").exists()

    draft = session.evaluate(center=[0, 0, 0.2], representation="rotvec", radius_deg=5)
    coordinates = session.coordinates
    landscape = Landscape(data=pd.DataFrame({
        "particle_key": session.arrays.particle_key,
        "coordinates_analysis": list(coordinates),
        "coordinates_display": list(coordinates),
        "sld_unfloored": np.ones(len(coordinates)),
        "sld_raw": np.ones(len(coordinates)),
        "sld_display": np.ones(len(coordinates)),
        "sld_was_floored": np.zeros(len(coordinates), dtype=bool),
        "sld_local_k_mean": np.ones(len(coordinates)),
        "sld_effective_local_k_mean": np.ones(len(coordinates)),
        "sld_distance_floor": np.full(len(coordinates), 1e-4),
    }))
    expected = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="radius_around_center", center_input=(0, 0, 0.2),
            center_input_representation="rotvec", evaluation_space="analysis",
            metric="so3_geodesic", radius=5, radius_unit="degrees", selection_id="x",
        ),
    )
    assert draft["selected_count"] == expected.selected_count
    assert draft["exact"] is True
    assert not (run / "selections").exists()


def test_confirm_requires_id_writes_standard_artifact_and_refuses_overwrite(tmp_path) -> None:
    run = _bundle(tmp_path)
    session = ExploreSession(run, max_display_points=5)
    session.evaluate(center=[0, 0, 0.2], representation="rotvec", radius_deg=5)

    result = session.confirm("state_1")

    selection_path = run / "selections" / "state_1" / "selection.json"
    payload = json.loads(selection_path.read_text(encoding="utf-8"))
    assert result["selection_id"] == "state_1"
    assert payload["parent_run_id"] == "run-ux"
    assert payload["parent_landscape_metadata"]["sha256"]
    assert payload["interaction_provenance"]["displayed_point_count"] == 5
    try:
        session.confirm("state_1")
    except FileExistsError:
        pass
    else:
        raise AssertionError("interactive confirm must not overwrite a selection")


def test_interactive_and_cli_share_standard_selection_artifact_membership(tmp_path) -> None:
    run = _bundle(tmp_path)
    session = ExploreSession(run, max_display_points=5)
    session.evaluate(center=[0, 0, 0.2], representation="rotvec", radius_deg=5)
    session.confirm("interactive")
    args = build_parser().parse_args(
        [
            "select", "--run-dir", str(run), "--selection-id", "cli",
            "--center-representation", "rotvec", "--center", "0", "0", "0.2",
            "--radius", "5", "--metric", "so3",
        ]
    )
    select_command(args)

    interactive_dir = run / "selections" / "interactive"
    cli_dir = run / "selections" / "cli"
    for filename in ("selected_particle_keys.csv", "selection.csv", "selected_landscape_rows.csv"):
        assert (interactive_dir / filename).read_text(encoding="utf-8") == (
            cli_dir / filename
        ).read_text(encoding="utf-8")


def test_canonical_space_and_identity_changes_are_validated(tmp_path) -> None:
    run = _bundle(tmp_path)
    _add_canonical(run)
    session = ExploreSession(run, space="canonical", max_display_points=5)
    draft = session.evaluate(
        center=[0, 0, 0], representation="rotvec", radius_deg=2,
        expected_run_id="run-ux", expected_landscape_sha256=session.landscape_sha256,
    )
    assert draft["space"] == "canonical"
    try:
        session.evaluate(center=[0, 0, 0], representation="rotvec", radius_deg=2, expected_run_id="wrong")
    except ValueError as exc:
        assert "run_id mismatch" in str(exc)
    else:
        raise AssertionError("run mismatch must be rejected")
    with session.landscape_path.open("ab") as handle:
        handle.write(b"changed")
    try:
        session.confirm("changed_landscape")
    except ValueError as exc:
        assert "landscape changed" in str(exc)
    else:
        raise AssertionError("changed landscape must be rejected before confirm")


def test_local_server_binds_loopback_rejects_traversal_and_shuts_down(tmp_path) -> None:
    session = ExploreSession(_bundle(tmp_path), max_display_points=5)
    server = create_explore_server(session, port=0)
    thread = server.start_in_thread()
    try:
        assert server.host == "127.0.0.1"
        with urlopen(server.url + "api/health", timeout=3) as response:
            assert json.load(response)["status"] == "ok"
        with urlopen(server.url, timeout=3) as response:
            html = response.read().decode("utf-8")
        assert "https://" not in html
        with urlopen(server.url + "assets/app.js", timeout=3) as response:
            javascript = response.read().decode("utf-8")
        assert "http://" not in javascript
        assert "https://" not in javascript
        assert "center_rv" in javascript
        assert "center_euler" in javascript
        assert "viridis" not in javascript  # palette is selected by session policy
        try:
            urlopen(server.url + "../run_manifest.json", timeout=3)
        except HTTPError as exc:
            assert exc.code == 404
        else:
            raise AssertionError("path traversal must be rejected")
        request = Request(
            server.url + "api/shutdown",
            data=json.dumps({"session_token": session.session_token}).encode(),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urlopen(request, timeout=3) as response:
            assert json.load(response)["status"] == "shutting_down"
        thread.join(timeout=3)
        assert not thread.is_alive()
    finally:
        server.shutdown()


def test_http_evaluate_and_confirm_protocol_uses_exact_python_service(tmp_path) -> None:
    run = _bundle(tmp_path)
    session = ExploreSession(run, max_display_points=4, colormap="rainbow_r")
    server = create_explore_server(session, port=0)
    server.start_in_thread()
    identity = {
        "session_token": session.session_token,
        "run_id": session.run_id,
        "landscape_sha256": session.landscape_sha256,
    }
    try:
        evaluate = Request(
            server.url + "api/evaluate",
            data=json.dumps({
                **identity, "center": [0, 0, 0.2],
                "representation": "rotvec", "radius_deg": 5,
            }).encode(),
            headers={"Content-Type": "application/json"}, method="POST",
        )
        with urlopen(evaluate, timeout=3) as response:
            draft = json.load(response)
        assert draft["exact"] is True
        assert draft["full_candidate_count"] == 20
        assert len(draft["center_rv"]) == len(draft["center_euler"]) == 3

        confirm = Request(
            server.url + "api/confirm",
            data=json.dumps({**identity, "selection_id": "from_http"}).encode(),
            headers={"Content-Type": "application/json"}, method="POST",
        )
        with urlopen(confirm, timeout=3) as response:
            result = json.load(response)
        assert result["selected_count"] == draft["selected_count"]
        assert (run / "selections" / "from_http" / "selection.json").is_file()
    finally:
        server.shutdown()


def test_empty_display_filter_does_not_change_full_evaluation(tmp_path) -> None:
    run = _bundle(tmp_path)
    session = ExploreSession(run, display_threshold=1000, max_display_points=5)
    assert session.session_payload()["displayed_point_count"] == 0
    draft = session.evaluate(center=[0, 0, 0.2], representation="rotvec", radius_deg=5)
    assert draft["full_candidate_count"] == 20
    assert draft["selected_count"] > 0
