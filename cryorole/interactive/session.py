"""Array-native draft/evaluate/confirm session for offline exploration."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shlex
from typing import Any
from uuid import uuid4

import numpy as np
from scipy.spatial.transform import Rotation

from cryorole.export import read_landscape_npz_arrays, read_selection_json
from cryorole.io.landscape_resolver import resolve_landscape_path
from cryorole.models.policies import SelectionPolicy
from cryorole.models.selection import Selection
from cryorole.provenance import stream_sha256
from cryorole.run_bundle import validate_completed_run_bundle
from cryorole.select import (
    SelectedRowProvenance,
    SelectionArtifactRequest,
    center_to_rotvec,
    evaluate_radius_mask,
    write_selection_artifact,
)


_SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class ExploreSession:
    """One bounded display sample backed by one full NPZ array load."""

    def __init__(
        self,
        run_dir: str | Path,
        *,
        space: str = "raw",
        canonical_id: str = "default",
        selection_id: str | None = None,
        max_display_points: int = 50_000,
        display_threshold: float | None = None,
        display_top_fraction: float | None = None,
        colormap: str = "viridis",
    ) -> None:
        if max_display_points < 1:
            raise ValueError("max_display_points must be positive")
        if display_threshold is not None and display_top_fraction is not None:
            raise ValueError("display threshold and top fraction are mutually exclusive")
        if display_top_fraction is not None and not 0 < display_top_fraction <= 1:
            raise ValueError("display_top_fraction must be > 0 and <= 1")
        self.run_dir = Path(run_dir).expanduser().resolve()
        validate_completed_run_bundle(self.run_dir)
        self.space = space
        self.canonical_id = canonical_id
        self.landscape_path = resolve_landscape_path(
            self.run_dir, space=space, canonical_id=canonical_id
        ).resolve()
        if self.landscape_path.suffix.lower() != ".npz":
            raise ValueError("interactive explore requires an array-native NPZ landscape")
        self.landscape_sha256 = stream_sha256(self.landscape_path)
        self.arrays = read_landscape_npz_arrays(self.landscape_path)
        self.coordinates = self._coordinates_for_space()
        self.run_id, self.euler_sequence, self.euler_convention = self._run_metadata()
        self.session_token = uuid4().hex
        self.session_id = f"explore-{uuid4().hex[:12]}"
        self.colormap = colormap
        self.display_filter = self._display_filter(display_threshold, display_top_fraction)
        eligible = self._display_filter_indices(display_threshold, display_top_fraction)
        self.display_indices = _deterministic_sample(eligible, max_display_points)
        self.overlay_selection_id = selection_id
        self.overlay_keys = self._selection_overlay(selection_id)
        self._draft: dict[str, Any] | None = None

    def session_payload(self) -> dict[str, Any]:
        indices = self.display_indices
        coordinates = self.coordinates[indices]
        euler = Rotation.from_rotvec(coordinates).as_euler(
            self.euler_sequence, degrees=True
        )
        keys = self.arrays.particle_key[indices]
        return {
            "artifact_type": "cryorole_explore_session",
            "schema_version": "1.0",
            "session_id": self.session_id,
            "session_token": self.session_token,
            "run_id": self.run_id,
            "run_dir": str(self.run_dir),
            "landscape_path": str(self.landscape_path),
            "landscape_sha256": self.landscape_sha256,
            "space": self.space,
            "canonical_id": self.canonical_id if self.space == "canonical" else None,
            "euler_convention": self.euler_convention,
            "scipy_euler_sequence": self.euler_sequence,
            "sld_field": "sld_display",
            "metric": "so3_geodesic",
            "radius_unit": "degrees",
            "full_candidate_count": self.arrays.n_points,
            "displayed_point_count": len(indices),
            "sampling_policy": "deterministic_even_index_after_display_filter",
            "display_filter": self.display_filter,
            "selection_evaluation": "exact_full_parent_landscape",
            "colormap": self.colormap,
            "points": {
                "display_index": indices.tolist(),
                "particle_key": keys.tolist(),
                "rv": coordinates.tolist(),
                "euler": euler.tolist(),
                "sld": self.arrays.sld_display[indices].tolist(),
                "sld_display_is_outlier": self.arrays.sld_display_is_outlier[indices].tolist(),
                "overlay_selected": [str(key) in self.overlay_keys for key in keys],
            },
        }

    def evaluate(
        self,
        *,
        center: list[float] | tuple[float, float, float],
        representation: str,
        radius_deg: float,
        expected_run_id: str | None = None,
        expected_landscape_sha256: str | None = None,
    ) -> dict[str, Any]:
        self._validate_identity(expected_run_id, expected_landscape_sha256, rehash=False)
        center_rv = center_to_rotvec(
            center,
            representation=representation,
            scipy_euler_sequence=self.euler_sequence,
            euler_degrees=True,
        )
        selected, distances = evaluate_radius_mask(
            self.coordinates,
            center_rv,
            radius=radius_deg,
            radius_unit="degrees",
            metric="so3_geodesic",
        )
        display_selected = selected[self.display_indices]
        center_euler = Rotation.from_rotvec(center_rv).as_euler(
            self.euler_sequence, degrees=True
        )
        self._draft = {
            "artifact_type": "cryorole_selection_draft",
            "schema_version": "1.0",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "run_id": self.run_id,
            "landscape_sha256": self.landscape_sha256,
            "space": self.space,
            "canonical_id": self.canonical_id if self.space == "canonical" else None,
            "center_input": [float(value) for value in center],
            "center_input_representation": representation,
            "center_rv": center_rv.tolist(),
            "center_euler": center_euler.tolist(),
            "euler_convention": self.euler_convention,
            "radius": float(radius_deg),
            "radius_unit": "degrees",
            "metric": "so3_geodesic",
            "full_candidate_count": self.arrays.n_points,
            "selected_count": int(np.count_nonzero(selected)),
            "selected_fraction": float(np.mean(selected)),
            "exact": True,
            "displayed_point_count": len(self.display_indices),
            "display_selected": display_selected.tolist(),
            "display_local_neighborhood": (
                distances[self.display_indices] <= np.deg2rad(float(radius_deg) * 2)
            ).tolist(),
            "sampling_policy": "deterministic_even_index_after_display_filter",
            "display_filter": self.display_filter,
            "command_template": self._selection_command(center_rv, radius_deg, "SELECTION_ID"),
        }
        self._draft["_selected_mask"] = selected
        return {key: value for key, value in self._draft.items() if not key.startswith("_")}

    def confirm(
        self,
        selection_id: str,
        *,
        expected_run_id: str | None = None,
        expected_landscape_sha256: str | None = None,
    ) -> dict[str, Any]:
        if self._draft is None:
            raise ValueError("Evaluate a draft before confirming a scientific selection")
        if not _SAFE_ID.fullmatch(selection_id):
            raise ValueError("selection_id must use only letters, numbers, '.', '_' or '-'")
        self._validate_identity(expected_run_id, expected_landscape_sha256, rehash=True)
        output_dir = self.run_dir / "selections" / selection_id
        if output_dir.exists():
            raise FileExistsError(f"Selection already exists and interactive confirm never overwrites: {output_dir}")
        selected_mask = np.asarray(self._draft["_selected_mask"], dtype=bool)
        selected_keys = tuple(self.arrays.particle_key[selected_mask].tolist())
        interaction = {
            "source": "cryorole_explore",
            "session_id": self.session_id,
            "displayed_point_count": len(self.display_indices),
            "full_candidate_count": self.arrays.n_points,
            "sampling_policy": "deterministic_even_index_after_display_filter",
            "display_filter": self.display_filter,
            "selection_evaluation": "exact_full_parent_landscape",
            "overlay_selection_id": self.overlay_selection_id,
        }
        parent_metadata = {
            "path": str(self.landscape_path),
            "sha256": self.landscape_sha256,
            "size_bytes": self.landscape_path.stat().st_size,
            "row_count": self.arrays.n_points,
            "space": self.space,
            "canonical_id": self.canonical_id if self.space == "canonical" else None,
        }
        policy = SelectionPolicy(
            selection_mode="radius_around_center",
            center_input=tuple(self._draft["center_input"]),
            center_input_representation=str(self._draft["center_input_representation"]),
            center_input_space="evaluation",
            center_euler_convention=self.euler_convention,
            center_scipy_euler_sequence=self.euler_sequence,
            center_degrees=True,
            evaluation_space="canonical" if self.space == "canonical" else "analysis",
            metric="so3_geodesic",
            radius=float(self._draft["radius"]),
            radius_unit="degrees",
            top_fraction=None,
            parent_landscape_id=(self.canonical_id if self.space == "canonical" else "raw"),
            parent_landscape_metadata=parent_metadata,
            selection_id=selection_id,
        )
        selection = Selection(
            selection_id=selection_id,
            parent_landscape_id=policy.parent_landscape_id,
            parent_run_id=self.run_id,
            created_at=datetime.now(timezone.utc).isoformat(),
            interaction_provenance=interaction,
            parent_landscape_metadata=parent_metadata,
            selection_mode="radius_around_center",
            selection_basis=f"coordinates_{'canonical' if self.space == 'canonical' else 'analysis'}",
            metric="so3_geodesic",
            center_input=policy.center_input,
            center_input_representation=policy.center_input_representation,
            center_input_space="evaluation",
            center_euler_sequence="zyx",
            center_euler_convention=self.euler_convention,
            center_scipy_euler_sequence=self.euler_sequence,
            center_degrees=True,
            center_evaluated=tuple(self._draft["center_rv"]),
            evaluation_space=policy.evaluation_space,
            resolved_evaluation_space=policy.evaluation_space,
            radius=policy.radius,
            radius_unit="degrees",
            density_support_field="sld_raw",
            top_fraction=None,
            selected_particle_keys=selected_keys,
            selected_count=len(selected_keys),
            total_count=self.arrays.n_points,
            active_policy=policy,
        )
        summary = {
            "space": self.space, "canonical_id": parent_metadata["canonical_id"],
            "center_input": selection.center_input,
            "center_input_representation": selection.center_input_representation,
            "center_evaluated_rv": selection.center_evaluated,
            "center_euler": self._draft["center_euler"],
            "euler_convention": self.euler_convention,
            "radius": selection.radius, "radius_unit": selection.radius_unit,
            "metric": selection.metric, "full_candidate_count": selection.total_count,
            "selected_count": selection.selected_count,
            "interaction_provenance": interaction,
            "exact_cli_command": self._selection_command(np.asarray(selection.center_evaluated), float(selection.radius), selection_id),
            "next_commands": {
                "visualize": f"cryorole visualize --run-dir {shlex.quote(str(self.run_dir))} --selection-id {shlex.quote(selection_id)}",
                "export": f"cryorole export --run-dir {shlex.quote(str(self.run_dir))} --selection-id {shlex.quote(selection_id)}",
            },
        }
        selected_indices = np.flatnonzero(selected_mask)
        write_selection_artifact(
            SelectionArtifactRequest(
                selection=selection,
                output_dir=output_dir,
                overwrite=False,
                selected_rows=SelectedRowProvenance(
                    particle_key=self.arrays.particle_key[selected_indices].tolist(),
                    ref_source_row_id=(
                        self.arrays.ref_source_row_id[selected_indices].tolist()
                        if self.arrays.ref_source_row_id is not None
                        else None
                    ),
                    mov_source_row_id=(
                        self.arrays.mov_source_row_id[selected_indices].tolist()
                        if self.arrays.mov_source_row_id is not None
                        else None
                    ),
                ),
                summary=summary,
            )
        )
        return json.loads(
            (output_dir / "selection_summary.json").read_text(encoding="utf-8")
        )

    def _coordinates_for_space(self) -> np.ndarray:
        if self.space == "raw":
            return self.arrays.coordinates_analysis
        if self.space == "canonical":
            if self.arrays.coordinates_canonical is None:
                raise ValueError("canonical explore requires coordinates_canonical")
            return self.arrays.coordinates_canonical
        raise ValueError("space must be raw or canonical")

    def _run_metadata(self) -> tuple[str, str, str]:
        manifest = json.loads((self.run_dir / "run_manifest.json").read_text(encoding="utf-8"))
        summary_path = self.run_dir / "run_summary.json"
        summary = json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.is_file() else {}
        run_id = manifest.get("run_id") or summary.get("run_id")
        if not run_id:
            raise ValueError("interactive explore requires a run_id")
        return str(run_id), str(summary.get("scipy_euler_sequence", "zyx")), str(summary.get("euler_convention", "extrinsic_zyx"))

    def _display_filter(self, threshold, top_fraction) -> dict[str, object]:
        if threshold is not None:
            return {"mode": "threshold", "field": "sld_display", "threshold": float(threshold)}
        if top_fraction is not None:
            return {"mode": "top_fraction", "field": "sld_display", "top_fraction": float(top_fraction)}
        return {"mode": "all", "field": "sld_display"}

    def _display_filter_indices(self, threshold, top_fraction) -> np.ndarray:
        indices = np.arange(self.arrays.n_points, dtype=np.int64)
        if threshold is not None:
            return indices[self.arrays.sld_display >= float(threshold)]
        if top_fraction is not None:
            count = max(1, int(np.ceil(self.arrays.n_points * float(top_fraction))))
            order = np.lexsort((indices, -self.arrays.sld_display))
            return np.sort(order[:count]).astype(np.int64)
        return indices

    def _selection_overlay(self, selection_id: str | None) -> set[str]:
        if selection_id is None:
            return set()
        selection = read_selection_json(self.run_dir / "selections" / selection_id / "selection.json")
        if selection.parent_run_id and selection.parent_run_id != self.run_id:
            raise ValueError("overlay selection belongs to a different run_id")
        return {str(key) for key in selection.selected_particle_keys}

    def _validate_identity(self, expected_run_id, expected_hash, *, rehash: bool) -> None:
        if expected_run_id is not None and expected_run_id != self.run_id:
            raise ValueError("run_id mismatch")
        if expected_hash is not None and expected_hash != self.landscape_sha256:
            raise ValueError("landscape identity mismatch")
        current_manifest = json.loads((self.run_dir / "run_manifest.json").read_text(encoding="utf-8"))
        if current_manifest.get("run_id") != self.run_id:
            raise ValueError("run_id changed during interactive session")
        if rehash and stream_sha256(self.landscape_path) != self.landscape_sha256:
            raise ValueError("landscape changed during interactive session")

    def _selection_command(self, center_rv: np.ndarray, radius_deg: float, selection_id: str) -> str:
        parts = [
            "cryorole", "select", "--run-dir", str(self.run_dir),
            "--selection-id", selection_id, "--space", self.space,
        ]
        if self.space == "canonical":
            parts.extend(["--canonical-id", self.canonical_id])
        parts.extend([
            "--center-representation", "rotvec", "--center",
            *(f"{value:.12g}" for value in center_rv), "--radius", f"{radius_deg:.12g}",
            "--metric", "so3",
        ])
        return shlex.join(parts)


def _deterministic_sample(indices: np.ndarray, maximum: int) -> np.ndarray:
    values = np.asarray(indices, dtype=np.int64)
    if len(values) <= maximum:
        return values
    positions = np.linspace(0, len(values) - 1, num=maximum, dtype=np.int64)
    return values[np.unique(positions)]
