# AGENTS.md

## Purpose

This repository contains **cryoROLE 2.0**, a policy-driven cryo-EM relative-orientation analysis platform.

The core workflow is:

```text
input -> normalize -> match -> compute -> canonicalize -> visualize -> select -> export -> manifest
```

The primary engineering goal is to make this workflow stable, auditable, scalable, reproducible, and safe for large cryo-EM particle sets. Do not add new analysis variants before preserving the scientific and artifact contracts below.

---

## Required read order

For any behavior-changing work, read documents in this order:

1. `AGENTS.md` — non-negotiable guardrails and current priority.
2. `docs/architecture.md` — stable command, artifact, and policy contract.
3. `docs/production_run_plan.md` — current implementation plan for production-scale `cryorole run` work.
4. `docs/roadmap.md` — phase priorities and release direction.

If a task changes persistence, CLI behavior, artifact layout, visualization, selection, export, or memory behavior, update the relevant document together with code and tests.

---

## Current active priority

The active engineering priority is:

```text
Make cryorole run production-scale.
```

This means reducing memory pressure in the `run` path by moving the production backend toward array-native data structures, vectorized RO computation, batched SLD density calculation, and chunked raw CSV export.

Canonicalization, visualization, selection, and export already have a working command-layer contract. Do not destabilize them while refactoring the run backend.

---

## Non-negotiable scientific invariants

1. Do not change the scientific definition of relative orientation:

```text
RO = R_ref^-1 R_mov
```

2. Do not silently change the passive-frame interpretation currently used by cryoROLE.
3. Keep source convention conversion centralized in the normalization layer.
4. Do not scatter `.T`, `.inv()`, active/passive conversion, or Euler convention logic across downstream modules.
5. Internal rotation truth is a `3 x 3` active rotation matrix. Rotation vectors, quaternions, and Euler angles are derived representations.
6. Euler angle output convention is a derived-representation policy. Public cryoROLE RO/RV Euler output uses extrinsic fixed-axis ZYX; intrinsic ZYX is legacy/internal compatibility only.
7. Do not conflate analysis-space quantities with display-space quantities.
8. Do not write display-only reparameterizations back into STAR/CS metadata.
9. Do not overwrite original input files.
10. Any behavior that changes interpretation must be controlled by an explicit policy and recorded in reports/manifests.

---

## Input and matching guardrails

### RELION `.star`

- RELION `Rot/Tilt/Psi` must be parsed as intrinsic `ZYZ` Euler angles.
- SciPy parsing must use uppercase `"ZYZ"`.
- RELION passive-to-active conversion must remain centralized and tested.
- Public `cryorole run` may default STAR identity matching to `_rlnTomoParticleName` first, then `_rlnImageName` / `rlnImageName`, only when the resolved column is present and safe.
- Never assume row-order correspondence between two STAR files unless the user explicitly requests row-aligned mode.
- Never use Euler angles, origin shifts, or other refinement-result columns as default particle identity keys.
- If RELION identity ambiguity remains unresolved, fail before RO computation.

### CryoSPARC `.cs`

- Treat `.cs` as a native input format.
- Default identity key is `uid`.
- Native pose field is `alignments3D/pose`.
- Native pose encoding is rotation vector / axis-angle.
- Internal active matrices must be derived from `Rotation.from_rotvec(...).as_matrix()`.
- Preserve vector-valued `.cs` fields without flattening them incorrectly.

### Matching

- Particle matching must be explicit and auditable.
- `--row-aligned` is an explicit user assertion for any supported input type; when set, `run` must not key-match, reorder, or drop rows, and must fail if ref/mov row counts differ.
- Key-based matching that reorders rows or drops unmatched particles must emit a prominent warning and record counts in reports/manifests.
- Duplicate-key and low-overlap behavior must be policy-controlled and reported.
- Match tables must preserve source-row provenance for later selection/export backtracking.

---

## Artifact and command guardrails

### Run

- `cryorole run` must create a run bundle.
- Public `cryorole run` should use explicit input arguments: `--ref REF_METADATA --mov MOV_METADATA`.
- Run metadata domain labels default to `ref` and `mov`, with optional explicit overrides.
- `data/raw_landscape.npz` is the production machine-readable raw landscape.
- `data/raw_landscape.csv` is the user-facing raw table.
- `cryorole run` must derive user-facing Euler columns with the recorded Euler convention policy. Public output uses extrinsic fixed-axis ZYX.
- Raw CSV Euler column names may remain simplified as `raw_ea_zyx_*`; reports/manifests own the intrinsic/extrinsic interpretation.
- Full landscape JSON is debug-only and must be opt-in.
- Production run backend should be array-native; DataFrame/object-column paths are compatibility/debug paths.
- Raw visualization is a default convenience and should write both all-particle and `sld_raw > 1.5` preview views using legacy rainbow style. Display `vmax` is `min(displayed max SLD, 100)` unless disabled with `--no-visualize`.
- Translation distance is future report-only diagnostic work. Do not make translation distance a default `run` analysis coordinate or let it affect RO, SLD, canonicalization, selection, or export by default.

### Canonicalize

- `cryorole canonicalize` derives a canonical coordinate frame and must not overwrite raw artifacts.
- Canonicalization should use array-native backend for NPZ/run-dir inputs.
- Public canonicalization should default to `cryorole canonicalize --run-dir RUN`, writing under `canonical/<canonical_id>/`.
- `--fit-top-fraction` and its compact alias `--fit-top` control the high-`sld_raw` subset used for canonical fitting and must be recorded.
- Public canonicalization uses density-weighted skewness for axis sign disambiguation. The default positive side is `low_density_skew`; public `--positive-side low|high` may choose the direction.
- Canonicalization should write a reusable canonical frame artifact and may apply one from another run without refitting.
- Canonicalization writes display-only quick-look previews by default unless disabled with `--no-visualize`: all particles, `sld_raw > 1.5`, and the top-`sld_raw` fraction used for frame fitting.

### Visualize

- Visualization is display-only.
- Public visualization uses the source landscape's recorded Euler convention; do not expose Euler convention overrides in the compact public command.
- Display filters, axis ranges, downsampling, and style defaults must not alter raw/canonical landscapes or selections.
- Public visualization should default to legacy rainbow style, PNG output, fixed `sld_display` coloring, optional display colormap choice, and equal display units within each figure so plots do not distort.
- Public visualization may write basic 1D distributions for displayed Euler/rotvec coordinates; these use the same display-filtered rows and remain display-only.
- Public visualization outputs should stay under the run bundle `visualizations/` directory, using a user-provided `visual_id` or `default`.
- Visualization may consume a selected-derived landscape when explicitly requested; this is still display-only and must not mutate the parent raw/canonical landscape or the selection.
- `sld_display` is display-only; do not let it drive scientific defaults.
- Do not use fixed-percentile SLD clipping as the default display-scaling strategy. Default `run` quick-look figures may cap color display with `vmax = min(displayed max SLD, 100)` without changing stored density values.
- Tail-jump display outliers must remain in raw/canonical landscapes and plots. Do not overwrite `sld_raw`; do not rewrite outlier `sld_display` values merely to clip colors.

### Select

- Selection is a first-class scientific decision object.
- Radius-based orientation selection defaults to SO(3) geodesic distance.
- Public radius selection uses `--radius` for degrees or `--radius-rad` for radians.
- Do not use naive Euler-space Euclidean distance for scientific radius selection.
- Euler center input or Euler range selection must record and obey the explicit Euler convention policy.
- Top-density selection defaults to including all rows, including display-only SLD outliers; excluding display outliers must be an explicit selection policy and recorded.
- Random-fraction selection is allowed only as an explicit selection mode with a recorded fraction, seed, candidate count, and selected count.
- Source-metadata selection is allowed only against the run-time ref/mov source metadata, using recorded `ref_source_row_id` or `mov_source_row_id`; the metadata domain, column, selected values, and source-row policy must be recorded.
- Public metadata value selection may accept comma-separated values such as `--metadata-value 1,3`.
- Metadata split selection may create one selection per unique metadata value, but each child selection must remain a standard selection artifact and preserve export backtracking.
- Selection may optionally write a selected-derived landscape for visualization. Recomputed SLD on that subset must be an explicit policy and must preserve parent SLD fields separately.
- Public `--overwrite` may replace artifacts for the same `selection_id` only; it must not modify raw/canonical landscapes or unrelated selections.

### Export

- The primary user-facing command is direct `cryorole export`.
- `cryorole export selection` may remain only as a compatibility alias.
- STAR/CS metadata subset export must consume an explicit selection and use recorded `ref_source_row_id` / `mov_source_row_id` backtracking.
- Do not reselect or rematch particles during export.
- Do not modify source STAR/CS files.
- First-pass metadata subset export is a pure subset operation; do not write canonical/display coordinates as physical poses.

---

## Production-scale run guardrails

During production-scale `run` work:

1. Preserve RO, convention, matching, SLD, canonicalization, selection, and export semantics.
2. Prefer compact arrays over per-row Python objects.
3. Avoid storing full source rows as per-particle dictionaries in production paths.
4. Compute RO in vectorized form, using transpose for rotation-matrix inverse:

```python
R_ro = np.matmul(np.swapaxes(R_ref, 1, 2), R_mov)
```

5. SLD/kNN density must support batched query to avoid materializing huge `n x k` arrays at once.
6. Raw CSV writing should be post-computation and chunked/streaming where practical.
7. Preserve the existing SLD distance-floor stabilization, and add tail-jump SLD display-outlier detection without changing `sld_raw`.
8. Density reports should warn about large near-duplicate or near-identity RO concentrations and should count display-only SLD outliers.
9. Memory profiling and benchmarks should report peak RSS, wall time, row counts, and artifact sizes.

Detailed implementation plan: `docs/production_run_plan.md`.

---

## Things Codex must not do

- Do not change RO definition.
- Do not change passive/active interpretation silently.
- Do not move convention logic into plotting, CLI, export, GUI, or downstream helpers.
- Do not silently fall back to intrinsic ZYX for RO/RV Euler output, visualization, or Euler-based selection when the recorded landscape policy is extrinsic fixed-axis ZYX.
- Do not invent additional default RELION particle identity rules beyond the documented `_rlnTomoParticleName`, `_rlnImageName`, and `rlnImageName` run defaults.
- Do not use row index as a matching key unless row-aligned mode is explicitly requested, row counts match, and the policy is reported.
- Do not restore full `landscape.json` as the production persistence contract.
- Do not treat display-filtered rows as scientific selections.
- Do not treat random downsampling, visualization downsampling, or display filtering as equivalent operations; only explicit random selection creates a scientific selection artifact.
- Do not implement metadata selection by accepting an arbitrary new STAR file by default; use the run-time ref/mov source metadata and recorded source-row provenance. External annotation files require a future explicit join-key policy.
- Do not infer whether metadata selection should use ref or mov when both domains may contain the requested column; require an explicit metadata domain.
- Do not overwrite raw/canonical landscapes when writing a selected-derived landscape.
- Do not silently replace parent `sld_raw` with recomputed subset SLD; record recomputation policy and preserve parent SLD provenance.
- Do not use `sld_display` as the default scientific selection/canonical fitting field.
- Do not silently exclude display-only SLD outliers from scientific selection; require an explicit density-artifact policy such as excluding display outliers.
- Do not conflate coordinate/shift deltas or defocus proxies with RO, and do not use translation diagnostics to drive scientific defaults without an explicit recorded policy.
- Do not silently map PC1 / primary motion to the tertiary displayed coordinate.
- Do not claim density-weighted skewness defines an absolute biological clockwise/counterclockwise direction.
- Do not implement STAR/CS export by reselecting, rematching, or rewriting source poses.
- Do not remove reports/manifests or downgrade test coverage.

---

## Testing discipline

Before changing behavior, add or update tests for the relevant invariant group.

Minimum expectation for non-trivial changes:

1. Unit tests for new policy/model/helper behavior.
2. CLI tests for new or changed command arguments.
3. Persistence tests for artifact shape, paths, and row counts.
4. Regression tests protecting RO, convention handling, SLD, canonical axis mapping, selection semantics, and export backtracking when touched.
5. For production-scale run work, include memory/benchmark-oriented tests or scripts when practical; do not put very large benchmarks in the normal fast test suite.

Every implementation summary should list changed files and the tests run.

---

## Working style

- Prefer the smallest safe refactor that advances the current priority.
- Keep scientific core, CLI, visualization, selection, export, and persistence responsibilities separated.
- Record all interpretation-changing policies in reports/manifests.
- Keep documents short enough to stay authoritative: update `architecture.md` for stable contracts, `production_run_plan.md` for current sprint details, and `roadmap.md` for phase priority.
