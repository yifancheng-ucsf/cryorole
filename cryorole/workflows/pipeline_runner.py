"""Thin Phase 1 pipeline shell.

This shell supports read -> normalize -> identity -> match, plus explicit
helpers for RO, density, canonicalization, selection, selection export, and
manifest writing from already normalized results. It does not compute GUI or
reconstruction.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from cryorole.canonicalize import canonicalize_landscape
from cryorole.core import compute_landscape_density, compute_relative_orientations
from cryorole.export import export_selection, write_run_manifest
from cryorole.io.readers import read_cryosparc_cs, read_relion_star
from cryorole.match import ResolvedIdentity, match_resolved_identities, resolve_identity
from cryorole.models.export_report import ExportReport
from cryorole.models.landscape import Landscape
from cryorole.models.manifest_report import ManifestReport
from cryorole.models.match_table import MatchReport, MatchTable
from cryorole.models.policies import (
    CanonicalizationPolicy,
    ConventionPolicy,
    DensityPolicy,
    IdentityPolicy,
    MatchPolicy,
    RunManifestPolicy,
    SelectionExportPolicy,
    SelectionPolicy,
)
from cryorole.models.pose_table import PoseTable
from cryorole.models.ro_result import ROResult
from cryorole.models.selection import Selection
from cryorole.normalize import ConventionResolver, PoseNormalizer
from cryorole.select import select_particles


@dataclass
class Phase1PipelineResult:
    """Result returned by the Phase 1 pipeline shell."""

    pose_a: PoseTable
    pose_b: PoseTable
    identity_a: ResolvedIdentity
    identity_b: ResolvedIdentity
    match_table: MatchTable
    match_report: MatchReport


@dataclass
class Phase2ROResult:
    """Explicit Phase 2 RO computation result."""

    phase1: Phase1PipelineResult
    ro_result: ROResult


@dataclass
class Phase3DensityResult:
    """Density result computed from an RO result without canonicalization."""

    phase2: Phase2ROResult
    landscape: Landscape


@dataclass
class Phase4CanonicalizationResult:
    """Canonicalized landscape result computed from a density result."""

    phase3: Phase3DensityResult
    landscape: Landscape


@dataclass
class Phase5SelectionResult:
    """Selection result computed from an existing Landscape."""

    landscape: Landscape
    selection: Selection


@dataclass
class Phase6SelectionExportResult:
    """Export report for a Selection."""

    selection: Selection
    export_report: ExportReport


@dataclass
class Phase7ManifestResult:
    """Report for a written manifest."""

    manifest_report: ManifestReport


class PipelineRunner:
    """Run only the Phase 1 read/normalize/match portion of cryoROLE."""

    def __init__(self) -> None:
        self.normalizer = PoseNormalizer()

    def _read_and_normalize(
        self,
        path: str,
        *,
        domain_name: str,
        convention_policy: ConventionPolicy | None,
    ):
        suffix = Path(path).suffix.lower()
        if suffix == ".star":
            raw = read_relion_star(path)
            resolver = ConventionResolver(convention_policy or ConventionPolicy.relion_default())
            pose = self.normalizer.normalize_relion(
                raw.particles,
                domain_name=domain_name,
                convention_resolver=resolver,
            )
            return raw.particles, "relion", pose
        if suffix == ".cs":
            raw = read_cryosparc_cs(path)
            pose = self.normalizer.normalize_cryosparc(raw.particles, domain_name=domain_name)
            return raw.particles, "cryosparc", pose
        raise ValueError(f"Unsupported input suffix for Phase 1: {suffix}")

    def run_phase1(
        self,
        path_a: str,
        path_b: str,
        *,
        domain_a_name: str,
        domain_b_name: str,
        identity_policy_a: IdentityPolicy | None = None,
        identity_policy_b: IdentityPolicy | None = None,
        convention_policy_a: ConventionPolicy | None = None,
        convention_policy_b: ConventionPolicy | None = None,
        match_policy: MatchPolicy | None = None,
    ) -> Phase1PipelineResult:
        """Run the Phase 1 shell and return matched, normalized inputs."""

        raw_a, source_a, pose_a = self._read_and_normalize(
            path_a,
            domain_name=domain_a_name,
            convention_policy=convention_policy_a,
        )
        raw_b, source_b, pose_b = self._read_and_normalize(
            path_b,
            domain_name=domain_b_name,
            convention_policy=convention_policy_b,
        )
        self._validate_row_aligned_policies(
            raw_a,
            raw_b,
            source_a=source_a,
            source_b=source_b,
            identity_policy_a=identity_policy_a,
            identity_policy_b=identity_policy_b,
        )
        identity_a = resolve_identity(raw_a, source_type=source_a, policy=identity_policy_a)
        identity_b = resolve_identity(raw_b, source_type=source_b, policy=identity_policy_b)
        match_table, match_report = match_resolved_identities(
            identity_a,
            identity_b,
            policy=match_policy,
        )
        return Phase1PipelineResult(
            pose_a=pose_a.with_particle_keys(identity_a.keys),
            pose_b=pose_b.with_particle_keys(identity_b.keys),
            identity_a=identity_a,
            identity_b=identity_b,
            match_table=match_table,
            match_report=match_report,
        )

    def _validate_row_aligned_policies(
        self,
        raw_a,
        raw_b,
        *,
        source_a: str,
        source_b: str,
        identity_policy_a: IdentityPolicy | None,
        identity_policy_b: IdentityPolicy | None,
    ) -> None:
        row_aligned_a = (
            identity_policy_a is not None and identity_policy_a.identity_mode == "row_aligned"
        )
        row_aligned_b = (
            identity_policy_b is not None and identity_policy_b.identity_mode == "row_aligned"
        )
        if not row_aligned_a and not row_aligned_b:
            return
        if not row_aligned_a or not row_aligned_b:
            raise ValueError("row_aligned identity mode must be requested for both domains")
        if len(raw_a) != len(raw_b):
            raise ValueError(
                "row_aligned identity mode requires equal reference and moving row counts; "
                f"got {len(raw_a)} and {len(raw_b)}. Use cryorole align or manually "
                "pre-align metadata before running RO computation."
            )

    def compute_ro_for_phase1_result(
        self,
        phase1: Phase1PipelineResult,
        *,
        euler_sequence: str | None = None,
        degrees: bool = True,
    ) -> Phase2ROResult:
        """Compute RO from normalized matched Phase 1 outputs only."""

        ro_result = compute_relative_orientations(
            phase1.pose_a,
            phase1.pose_b,
            phase1.match_table,
            euler_sequence=euler_sequence,
            degrees=degrees,
        )
        return Phase2ROResult(phase1=phase1, ro_result=ro_result)

    def compute_density_for_ro_result(
        self,
        phase2: Phase2ROResult,
        *,
        density_policy: DensityPolicy | None = None,
    ) -> Phase3DensityResult:
        """Compute density from RO-derived analysis coordinates only."""

        landscape = compute_landscape_density(
            phase2.ro_result,
            policy=density_policy,
        )
        return Phase3DensityResult(phase2=phase2, landscape=landscape)

    def canonicalize_density_result(
        self,
        phase3: Phase3DensityResult,
        *,
        canonicalization_policy: CanonicalizationPolicy | None = None,
    ) -> Phase4CanonicalizationResult:
        """Canonicalize an existing density-bearing Landscape only."""

        landscape = canonicalize_landscape(
            phase3.landscape,
            policy=canonicalization_policy,
        )
        return Phase4CanonicalizationResult(phase3=phase3, landscape=landscape)

    def select_from_landscape(
        self,
        landscape: Landscape,
        *,
        selection_policy: SelectionPolicy | None = None,
    ) -> Phase5SelectionResult:
        """Create a non-destructive Selection from an existing Landscape."""

        selection = select_particles(landscape, policy=selection_policy)
        return Phase5SelectionResult(landscape=landscape, selection=selection)

    def export_selection_result(
        self,
        selection: Selection,
        *,
        export_policy: SelectionExportPolicy,
    ) -> Phase6SelectionExportResult:
        """Export a Selection without source mutation or reconstruction."""

        export_report = export_selection(selection, policy=export_policy)
        return Phase6SelectionExportResult(selection=selection, export_report=export_report)

    def write_manifest(
        self,
        *,
        manifest_policy: RunManifestPolicy,
        **manifest_sections,
    ) -> Phase7ManifestResult:
        """Write a JSON manifest from available partial workflow sections."""

        manifest_report = write_run_manifest(
            policy=manifest_policy,
            **manifest_sections,
        )
        return Phase7ManifestResult(manifest_report=manifest_report)
