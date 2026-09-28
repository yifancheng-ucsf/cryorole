"""Explicit policy objects for interpretation-changing choices."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Sequence


DEFAULT_MIN_MATCH_OVERLAP = 0.50


@dataclass(frozen=True)
class InputPolicy:
    """Controls source type inference and strictness for raw input parsing."""

    source_type: str | None = None
    strict: bool = True
    required_fields: tuple[str, ...] = ()


@dataclass(frozen=True)
class ConventionPolicy:
    """Controls conversion from source pose conventions to internal matrices."""

    source_software: str
    source_euler_sequence: str
    source_semantics: str = "passive"
    internal_semantics: str = "active"
    degrees: bool = True
    conversion_rule: str = (
        "active_matrix = scipy Rotation.from_euler('ZYZ', angles).as_matrix().T"
    )

    @classmethod
    def relion_default(cls) -> "ConventionPolicy":
        """Return the required RELION convention policy."""

        return cls(source_software="relion", source_euler_sequence="ZYZ")

    @classmethod
    def cryosparc_default(cls) -> "ConventionPolicy":
        """Return the required CryoSPARC ``alignments3D/pose`` convention policy.

        CryoSPARC stores an axis-angle vector (radians) describing the rotation
        applied during back-projection (image to volume); RELION's Rot/Tilt/Psi
        describe the opposite direction (reference volume to particle image).
        The two matrices are therefore transposes of each other. cryoROLE uses
        the pyem ``csparc2star.py`` mapping, so a ``.cs`` input gives the same
        internal active matrix as its pyem-converted STAR through the RELION
        bridge.
        """

        return cls(
            source_software="cryosparc",
            source_euler_sequence="rotvec",
            source_semantics="backprojection",
            internal_semantics="active",
            degrees=False,
            conversion_rule="active_matrix = scipy Rotation.from_rotvec(pose).as_matrix().T",
        )


@dataclass(frozen=True)
class IdentityPolicy:
    """Controls particle identity resolution and ambiguity behavior."""

    identity_mode: str
    identity_columns: tuple[str, ...] = ()
    column_normalization_rules: Mapping[str, str] = field(default_factory=dict)
    duplicate_policy: str = "fail"
    overlap_threshold: float = 0.0
    failure_mode: str = "fail"
    mapping_file: str | None = None

    @classmethod
    def cryosparc_uid(cls) -> "IdentityPolicy":
        """Return the CryoSPARC uid identity policy."""

        return cls(identity_mode="cryosparc_uid", identity_columns=("uid",))

    @classmethod
    def relion_image_name(cls) -> "IdentityPolicy":
        """Return the public RELION particle identity policy.

        The resolver prefers tomo particle names when present and falls back to
        image names for single-particle STAR files.
        """

        return cls(identity_mode="relion_image_name", identity_columns=("_rlnImageName",))

    @classmethod
    def relion_user_columns(
        cls,
        columns: Sequence[str],
        *,
        duplicate_policy: str = "fail",
        overlap_threshold: float = 0.0,
        failure_mode: str = "fail",
    ) -> "IdentityPolicy":
        """Return a RELION user-column identity policy."""

        return cls(
            identity_mode="relion_user_columns",
            identity_columns=tuple(columns),
            duplicate_policy=duplicate_policy,
            overlap_threshold=overlap_threshold,
            failure_mode=failure_mode,
        )


@dataclass(frozen=True)
class MatchPolicy:
    """Controls joining of resolved particle identities."""

    join_type: str = "inner"
    duplicate_handling: str = "fail"
    low_overlap_behavior: str = "fail"
    overlap_threshold: float = DEFAULT_MIN_MATCH_OVERLAP
    low_overlap_allowed: bool = False


@dataclass(frozen=True)
class DensityPolicy:
    """Controls analysis-space density computation and display normalization."""

    density_metric: str = "sld"
    sld_metric: str = "rotvec_euclidean"
    coordinate_source: str = "rotvec_ro"
    k_neighbors: int = 50
    distance_floor_mode: str = "relative_to_global_local_k_mean"
    distance_floor_fraction: float = 1e-4
    duplicate_policy: str = "report_only"
    display_normalization_mode: str = "identity"
    display_outlier_mode: str = "tail_jump"
    tail_search_fraction: float = 0.01
    tail_jump_factor: float = 5.0
    max_display_outlier_fraction: float = 0.002
    near_identity_ro_tolerance_rad: float = 1e-8
    near_duplicate_coordinate_tolerance_rad: float = 1e-8
    ro_coincidence_min_group_size: int = 10
    ro_coincidence_min_rows: int = 100
    ro_coincidence_min_fraction: float = 0.01


@dataclass(frozen=True)
class RepresentationPolicy:
    """Controls derived representation output conventions."""

    euler_convention: str = "extrinsic_zyx"
    scipy_euler_sequence: str = "zyx"
    euler_degrees: bool = True


@dataclass(frozen=True)
class CanonicalizationPolicy:
    """Controls policy-driven canonical frame alignment of landscapes."""

    method: str = "motion_axis_pca"
    axis_method: str = "pca_svd"
    axis_assignment: str = "pc123_to_alpha_beta_gamma"
    fit_subset: str = "top_fraction"
    fit_top_fraction: float = 0.40
    density_support_field: str = "sld_raw"
    sign_rule: str = "density_weighted_skewness"
    positive_side: str = "low_density_skew"
    sign_weight_field: str = "sld_raw"
    sign_ambiguity_threshold: float = 1e-8
    handedness_rule: str = "right_handed"
    origin_policy: str = "preserve"


@dataclass(frozen=True)
class AnimationPolicy:
    """Controls non-destructive Phase 1-4 offline animation export."""

    coordinate_set: str = "raw"
    path_space: str = "rv"
    euler_convention: str = "auto"
    interpolation: str = "so3_quaternion_slerp"
    frames_per_segment: int = 30
    fps: float = 30.0
    hold_frames: int = 0
    reverse: bool = False
    ping_pong: bool = False
    sld_field: str = "sld_display"
    map_frame: str = "raw"
    matrix_vector_convention: str = "column"
    quaternion_storage_order: str = "wxyz"
    quaternion_internal_order: str = "xyzw"
    render_mode: str = "script-only"
    composite_width: int = 1920
    composite_height: int = 1080
    composite_layout: str = "stacked"
    composite_background_color: str = "#ffffff"
    encode_movie: bool = True
    crf: int = 18
    movie_name: str = "animation.mp4"


@dataclass(frozen=True)
class SelectionPolicy:
    """Controls policy-driven particle selection from an existing Landscape."""

    selection_mode: str = "top_fraction_by_density"
    density_support_field: str = "sld_raw"
    density_artifact_policy: str = "include_all"
    top_fraction: float = 0.40
    random_fraction: float | None = None
    random_seed: int | None = None
    metadata_domain: str | None = None
    metadata_column: str | None = None
    metadata_values: tuple[str, ...] = ()
    metadata_source_details: Mapping[str, object] = field(default_factory=dict)
    metadata_source_file: str | None = None
    metadata_source_row_id_field: str | None = None
    split_by_metadata: bool = False
    threshold: float | None = None
    sld_min: float | None = None
    sld_max: float | None = None
    threshold_operator: str = ">="
    tie_break_rule: str = "stable_input_order"
    center_input: tuple[float, ...] | None = None
    center_input_representation: str = "rotvec"
    center_input_space: str = "evaluation"
    center_euler_sequence: str = "zyx"
    center_euler_convention: str = "extrinsic_zyx"
    center_scipy_euler_sequence: str = "zyx"
    center_degrees: bool = True
    evaluation_space: str = "canonical_if_available"
    metric: str = "so3_geodesic"
    radius: float | None = None
    radius_unit: str | None = None
    range_coordinate_source: str = "canonical_if_available"
    range_representation: str = "euler"
    range_euler_sequence: str = "zyx"
    range_euler_convention: str = "extrinsic_zyx"
    range_scipy_euler_sequence: str = "zyx"
    range_degrees: bool = True
    range_bounds: Mapping[str, tuple[float | None, float | None] | None] = field(
        default_factory=dict
    )
    parent_landscape_id: str | None = None
    parent_landscape_metadata: Mapping[str, object] = field(default_factory=dict)
    selection_id: str | None = None


@dataclass(frozen=True)
class SelectionExportPolicy:
    """Controls non-destructive export of Selection metadata and particle keys."""

    output_dir: str | Path
    selected_keys_filename: str = "selected_particle_keys.csv"
    selection_json_filename: str = "selection.json"
    export_report_filename: str = "export_report.json"
    overwrite: bool = False
    include_selected_particle_keys_in_json: bool = True
    export_reconstructed_maps: bool = False
    invoke_external_tools: bool = False


@dataclass(frozen=True)
class SelectionMetadataExportPolicy:
    """Controls source-format-preserving metadata subset export for a Selection."""

    output_dir: str | Path | None = None
    overwrite: bool = False
    domain: str = "both"
    format: str = "auto"
    run_dir: str | Path | None = None
    relocated_ref: str | Path | None = None
    relocated_mov: str | Path | None = None
    allow_unverified_source: bool = False


@dataclass(frozen=True)
class RunManifestPolicy:
    """Controls non-destructive writing of run provenance manifests."""

    output_path: str | Path
    overwrite: bool = False
    manifest_type: str = "cryorole_run_manifest"
    schema_version: str = "1"
    workflow_name: str | None = None
    command: str | None = None
    compute_file_hashes: bool = False
    hash_algorithm: str = "sha256"
    include_selected_particle_keys: bool = False
