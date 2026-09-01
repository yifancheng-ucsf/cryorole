# cryoROLE 2.0 Roadmap

**Status:** Living roadmap
**Current code state:** Architecture-complete beta with Core Productionization,
Workflow UX, and the typed CLI/service boundary refactor implemented.

---

## 1. Current position

cryoROLE 2.0 now has a complete command-layer workflow:

```text
preflight -> run -> status/next -> [canonicalize] -> explore/visualize -> select -> export
```

Implemented or substantially established:

```text
Run Bundle layout
raw landscape NPZ/CSV persistence
array-native canonicalization for NPZ/run-dir inputs
canonicalize --fit-top-fraction
display-only visualization separation
visualization-compatible rotated-landscape diagnostic script
legacy-compatible visualization style
SO(3)-aware selection
source-row-based STAR/CS metadata subset export
transactional run bundles with rollback-safe whole-bundle overwrite
streamed source identity and export-time SHA-256 verification
50% public matching safety gate with explicit audited override
array-native match-first run backend, batched SLD, and chunked raw CSV
shared side-effect-free preflight and run --dry-run
artifact-derived status/next and conservative guide
localhost-only offline interactive explorer and explicit Confirm boundary
policy/report/manifest-oriented architecture
offline animation trajectory/transform core (Phase 1)
offline three-panel landscape animation (Phase 2)
ChimeraX script-only and optional execute export (Phase 3)
animation display-policy parity and Linux offscreen completion checks
deterministic composite frames and validated H.264 MP4 export (Phase 4)
resolution-scaled projection marker/text and tighter panel spacing
optional validated two- and three-session structure views
repeatable rigid reference/moving model groups with per-model validation
shared resolved input policy across preflight/run/guide
content-level SourceIdentityGuard at formal run consumption
typed run/canonicalize/visualize/select services
standard Selection artifact writer shared by CLI and interactive Confirm
array-first NPZ visualization filtering and mode-minimal NPZ selection
```

Animation follow-up work:

```text
Phase 5 examples, publication validation, and public release notes
```

The Phase 4 presentation refinement is implemented: the default composite is
structure-first stacked, the three-projection strip is compact, and its
optional overlay is coordinate-only.

The camera-only `cryorole canonical-views` helper is implemented for exporting
deterministic canonical +X/+Y/+Z views from a raw or explicitly mapped
ChimeraX session without moving either density model. Optional
`--save-sessions` writes the same three views as derived camera-only `.cxs`
sessions while preserving the source session and model transforms.

The second presentation refinement is implemented: projection gaps are tighter,
the marker and coordinate line scale for the final composite, and an optional
secondary ChimeraX session provides synchronized dual views after independent
completion and initial-model-transform parity checks.

The third presentation refinement is implemented: dual-view stacked output may
crop a fixed horizontal fraction from both structure sources and anchor the
views inward. The animation colorbar now displays `SLD` without renaming the
underlying `sld_display` artifact field.

The fourth presentation refinement is implemented: an optional tertiary
ChimeraX session extends the existing single/dual workflow to three equal
stacked structure views with independent completion checks and shared baseline
transform validation. Fixed horizontal/vertical multi-view crop controls and a
three-view 58/42 structure/landscape split reduce unused canvas space while
preserving the existing single/dual layout.

Rigid model groups are implemented for three-body and similar sessions:
`--reference-model-id` and `--moving-model-id` are repeatable, all moving
models share one absolute delta from their own baselines, and completion plus
multi-view parity checks cover every declared model.

The implemented Phase 1–4 contract is documented
in `docs/animation_export.md`. Animation remains downstream and does not replace
the current production-hardening priority.

The bounded Phase 3 stabilization is complete:

1. reuses `visualize` filtering, legacy style, color scaling, point order,
   projection order, and Euler viewport policy for animation backgrounds;
2. makes `--threshold` consistent with `visualize`, retaining
   `--sld-threshold` only as a compatibility alias;
3. uses ChimeraX `--offscreen` for Linux headless rendering;
4. requires an explicit renderer-completion artifact in addition to process and
   PNG validation.

The next priority is not new science. It is release validation, larger 500k/1M
capacity checks, and continued gradual removal of DataFrame compatibility paths
where equivalence coverage already exists.

---

## 2. Release framing

Suggested version framing:

```text
2.0-beta.1   architecture-complete, command workflow functional
2.0-beta.2   production-scale run backend and export hardening
2.0-rc1      documentation, examples, benchmarks, paper validation outputs
2.0.0        stable public release
```

This is a planning frame, not a strict versioning requirement.

---

## 3. Completed milestone: Core Productionization

Goal:

```text
Make cryorole run robust for hundreds of thousands to millions of particles.
```

Core tasks:

1. Make the public `run` command compact: `cryorole run --ref REF --mov MOV`, with optional `--output-dir RUN_DIR` for naming the run bundle.
2. Default run domain labels to `ref` and `mov`, with explicit overrides when needed.
3. Default matching to `.cs` `uid` and STAR `_rlnTomoParticleName` / `_rlnImageName`; warn on reordering/dropped particles and fail unsafe matches with an align/manual-prealignment hint.
4. Add/strengthen array-native pose and matched-pose data models.
5. Vectorize RO computation.
6. Add batched kNN SLD computation.
7. Add opt-in `--sld-metric so3_geodesic`; keep `rotvec_euclidean` as the default and record the resolved metric.
8. Write raw NPZ directly from arrays.
9. Write raw CSV as chunked post-export.
10. Keep default raw visualization on as seven flat quick-look PNGs: Euler/RV triptychs for all particles, `sld_raw >= 1`, and the top 40%, plus a full-data log-SLD histogram. Use legacy rainbow display and `vmax = min(full-landscape max SLD, 100)`; `--no-visualize` must leave downstream commands usable and still write `run_report.md`.
11. Keep `cryorole run --help` compact by hiding visualization and developer/debug controls from normal help.
12. Add memory/timing reporting and synthetic benchmark scripts for both SLD metrics.
13. Assign stable `run_id`, record streamed source hashes, verify sources before
    export, support verified relocation, and gate legacy unverified export.
14. Publish run bundles transactionally and replace overwrite targets only after
    the new bundle is complete.

Detailed plan: `docs/production_run_plan.md`.

Exit criteria:

```text
full fast test suite passes
small deterministic array-native run matches compatibility path
default RV-Euclidean SLD is unchanged; opt-in SO(3) SLD matches an exact reference
default run CLI is compact and documented
CS uid and STAR image-name matching are tested
RELION tomo `_rlnTomoParticleName` matching is tested
reordered/dropped matches warn and are reported
default run visualization writes exactly seven flat PNGs and no extra quick-look artifacts
run help shows only public run controls
100k synthetic benchmarks compare both SLD metrics with memory profiles
canonicalize/select/export work from array-native run bundle
fault injection leaves no partial formal run and preserves old overwrite targets
source relocation, replacement, hash mismatch, legacy override, and run-id mismatch are tested
```

Implementation status:

Local sign-off on 2026-08-20 passed the complete suite (`546 passed, 2
skipped`). A synthetic 100k full-run comparison measured 184.84 MiB peak RSS
for `array_native` and 932.67 MiB for `dataframe_compat`; detailed parameters
and limitations are recorded in `docs/production_run_plan.md`.

```text
complete: source identity schema and export verification
complete: zero/low-overlap matching safety and reporting
complete: PoseArrays/MatchedPoseArrays/ROArrays/LandscapeArrays production path
complete: vectorized RO and batched RV/SO(3) SLD
complete: direct NPZ and chunked raw CSV
complete: RunBundleWriter staging, validation, commit, failure report, rollback
complete: auto -> array_native; explicit dataframe_compat reference backend
complete: canonical display-outlier preservation and canonical positive-side default parity
complete: full fast suite and recorded 100k production benchmark
```

---

## 3.1 Completed milestone: Workflow UX

Goal:

```text
Make the production workflow safe and approachable without hiding scientific decisions.
```

Implemented:

1. `cryorole preflight` and `cryorole run --dry-run` share production run's
   minimal-field, array-native pose/matching validation and create no bundle.
2. Schema-1.0 readiness reports use `READY`, `READY_WITH_WARNINGS`, or
   `BLOCKED`, with source identities, conventions, matching diagnostics,
   resource formulas/assumptions, and exact commands.
3. Production run consumes the same result and rejects an input whose content
   changed after inspection.
4. `status` derives integrity and progress from actual artifacts; `next`
   distinguishes required, recommended, and optional actions.
5. `guide` wraps preflight/status/next conservatively and requires an explicit
   run execution choice. It never makes a scientific selection or export.
6. `explore` serves packaged offline assets on `127.0.0.1`, renders linked
   EA/RV projections, and clearly separates draft evaluation from Confirm.
7. Display filtering/downsampling affects only presentation. Exact counts and
   Confirm use the full parent landscape and shared Python SO(3) evaluator.
8. Confirm writes the standard Selection with run/landscape hash, coordinate,
   radius, convention, count, timestamp, and UI provenance, and refuses
   overwrite or stale sessions.

Recorded 100k workflow benchmark (2026-08-21, synthetic local dataset,
`max_display_points=50,000`):

```text
preflight: 0.222 s wall, 134.37 MiB peak RSS, 100,000 matched
explore load: 0.052 s
exact full-landscape radius evaluation: 0.060 s
explore subprocess: 111.88 MiB peak RSS
raw CSV required: no
```

These are local planning measurements, not hardware-independent guarantees.
The milestone adds no optional Python dependency: the explorer server is
standard-library based and browser assets are packaged with the project.

---

## 4. Next sprint: public canonicalize cleanup

Goal:

```text
Make cryorole canonicalize a compact public command that writes reusable frames.
```

Tasks:

1. Keep the public command simple: `cryorole canonicalize --run-dir RUN`.
2. Write outputs to `RUN/canonical/<canonical_id>/`; default `canonical_id` is `default`.
3. Use density-weighted skewness as the public sign rule and default `--positive-side low`.
4. Keep `--fit-top-fraction` and add compact alias `--fit-top`.
5. Use public extrinsic fixed-axis ZYX Euler output; do not expose Euler convention choice in normal help.
6. Write default display-only canonical quick-look previews unless `--no-visualize` is set: `all_particles`, `filter_particles_by_sld_gt_1p5`, and `filter_particles_by_top_sld_XXpct` using the effective fit fraction.
7. Write `canonical_frame.json` and `canonical_frame.npz`.
8. Add `--use-frame FRAME` to apply an existing frame and skip fitting.
9. Add focused tests for compact help, default paths, low positive-side default, frame writing, frame reuse, and `--no-visualize`.

Non-goals:

```text
new canonicalization algorithms
translation or SE(3) frame fitting
frame registry/database
changing RO, SLD, selection, or export semantics
```

---

## 5. Implemented: public visualize cleanup

Goal:

```text
Make cryorole visualize a compact display-only command with predictable run-bundle output.
```

Tasks:

1. Keep the public command centered on `cryorole visualize --run-dir RUN`.
2. Remove public `--output-dir`; write under `RUN/visualizations/` using `--visual-id` or `default`.
3. Default to exactly two PNG three-view projections, `sld_display >= 1`, legacy rainbow style, fixed `sld_display` coloring, and equal display units within each figure.
4. Use one composable `--view 2d[,1d,3d]` control. Keep row filters (`--range`, `--top-fraction`, `--sld-threshold`, `--all`) separate from viewport policy (`--axis-limit`).
5. Remove public Euler convention/radian/sequence controls; use the source landscape's recorded Euler convention.
6. Make 1D opt-in, flat, full-filtered-data output with auto-binned percentage histograms; make coordinate KDE an explicit option.
7. Write selection visualizations under `visualizations/selections/<selection_id>/`, including parent-landscape and selected-landscape modes.
8. Provide opt-in self-contained offline interactive 3D and explicit static 3D. Interactive inspection cannot create or confirm a Selection.
9. Record all display filters, sampling policies, viewport policy, colormap policy, per-view counts, 1D statistics, 3D mode, selection provenance, and generated files in `visualization_report.json`.
10. Cover compact help, exact default outputs, view combinations, output paths, selection paths, range/viewport separation, formats, 1D full-data behavior, offline 3D, and display-only safety with focused tests.

Non-goals:

```text
new scientific selection behavior
changing run/canonicalize/select/export semantics
Euler convention overrides in public visualize
1D peak finding or multi-run overlay
```

---

## 6. Following sprint: public select cleanup

Goal:

```text
Make cryorole select a compact scientific-selection command with predictable run-bundle output.
```

Tasks:

1. Center the public command on `cryorole select --run-dir RUN --selection-id ID`; remove public `--output` but keep `--overwrite` for rerunning the same selection id.
2. Make radius-around-center the default mode with compact `--center/-c` and `--radius/-r`; default center input is Euler degrees and default metric is SO(3).
3. Keep only necessary radius advanced controls: `--center-representation`, `--radius-rad`, and `--metric`.
4. Remove public center/range Euler override clutter such as center input space, Euler convention, Euler sequence, and radians flags.
5. Keep compact modes for threshold, range, random, and metadata selection.
6. Threshold mode should support `--sld-min` and `--sld-max`.
7. Range mode should use only `--range-bound AXIS:LOWER:UPPER` plus the selected `--space`.
8. Random mode should use `--fraction F` and optional `--seed`.
9. Metadata mode should use `--metadata-domain`, `--metadata-column`, `--metadata-value VALUE[,VALUE...]`, and `--split-by-value`.
10. Keep optional selected-derived landscape writing with `--write-selected-landscape`; default inherits parent SLD and `--recompute-sld` is explicit.
11. Make select help concise but descriptive enough to explain `--write-selected-landscape`, `--recompute-sld`, `--overwrite`, and each mode's required parameters.
12. Add tests for compact help, default radius behavior, aliases, threshold min/max, range simplification, random reproducibility, metadata multi-value selection, selected-landscape persistence, overwrite behavior, and export backtracking.

Non-goals:

```text
new density definitions
display downsampling as scientific selection
arbitrary external STAR annotation files without an explicit join-key policy
overwriting raw/canonical landscapes
changing export to depend on recomputed SLD
density artifact or evaluation-space policy controls in public select
```

---

## 7. Following sprint: export hardening

Goal:

```text
Make public export simple, auditable, and reconstruction-ready in real
RELION/CryoSPARC workflows.
```

Tasks:

1. Center public usage on `cryorole export --run-dir RUN --selection-id ID`.
2. Default output to `RUN/exports/<selection_id>/`; do not require `--output-dir`.
3. Keep `--domain both` and `--format auto` as the public defaults.
4. Treat `--selection PATH/to/selection.json` and `--output-dir` as advanced compatibility inputs.
5. Fail clearly if `--selection` and `--selection-id` are both provided.
6. Keep `cryorole export selection` only as a compatibility alias; docs should prefer direct `cryorole export`.
7. Enrich `export_report.json` with selection path, run directory, source files, selected counts, source row-id stats, resolved domains, resolved formats, output directory, output files, overwrite policy, row-count checks, and warnings.
8. Add STAR particle-loop and optics-table diagnostics.
9. Add tests with RELION-style optics tables and multiple loops.
10. Add tests with realistic CryoSPARC structured arrays and vector fields.
11. Improve errors for missing selection inputs and missing source-row provenance.
12. Document RELION and CryoSPARC export workflows.

Non-goals:

```text
canonical transform metadata export
rewriting source poses with canonical/display coordinates
reselecting or rematching during export
making export depend on selected-derived or recomputed SLD
complex export registries or profiles
```

---

## 8. Following sprint: public align cleanup

Goal:

```text
Make cryorole align a compact STAR pre-alignment command for producing
row-aligned inputs before run.
```

Tasks:

1. Center public usage on `cryorole align --ref REF.star --mov MOV.star`.
2. Default output to `alignments/<align_id>/`; default `align_id` is `default`.
3. Write `aligned_ref.star`, `aligned_mov.star`, `match_table.csv`, `align_report.json`, and diagnostic STAR files for ref-only, mov-only, duplicate-ref, and duplicate-mov rows.
4. Auto key selection should use `_rlnTomoParticleName` first, then `_rlnImageName` / `rlnImageName`.
5. If auto image-name matching has nonzero overlap below 50%, continue with a prominent warning; if overlap is zero, fail and ask for explicit `--key`.
6. Allow explicit `--key` columns, including defocus columns, but record warnings when user-selected keys are not default identity keys.
7. Support `--float-tol COL=TOL` for numeric key columns and `--path-mode exact|basename|suffix:N` for path-like keys.
8. Default duplicate handling to `--duplicate-policy exclude`: exclude all rows whose key is duplicated in either input and write duplicate diagnostics.
9. Support `--duplicate-policy first`: keep the first occurrence and exclude later redundant rows with warnings.
10. Preserve optics and non-particle STAR blocks from each source file; only filter/reorder the target particle loop.
11. Ensure aligned outputs have equal row counts and can be passed to `cryorole run --row-aligned`.
12. Add focused tests for tomo-particle-name matching, image-name fallback, low-overlap warning, zero-overlap failure, explicit coordinate/defocus keys with float tolerance, path normalization, duplicate exclude/first policies, optics preservation, and output reports.

Non-goals:

```text
RO or SLD computation
selection or export behavior
CS alignment in the first public slice
row-index matching without explicit row-aligned assertion
using defocus, coordinates, Euler angles, shifts, or CTF fields as default identity keys
complex external annotation joins
```

---

## 9. Planned sprint: offline animation export

Goal:

```text
Generate reproducible offline movies that synchronize a trajectory in three
landscape projections with a validated rigid-body ChimeraX rendering.
```

This sprint starts only after the active production-scale `run` work and the
command/artifact cleanup needed for stable input bundles.

Tasks:

1. Implement a ChimeraX-independent trajectory core for EA/RV waypoint CSV
   files, quaternion sign continuity, SLERP, deterministic frame counts, and
   per-frame quaternion/RV/EA output.
2. Resolve raw and canonical landscape artifacts through recorded policy.
   Canonical waypoints must map back through the recorded canonical frame to
   raw RO before physical rendering.
3. Centralize and validate the baseline-relative passive-RO to active-density
   transform. Never accumulate transforms frame by frame.
4. Keep landscape coordinate set, baseline RO, ChimeraX scene basis, physical
   raw-to-scene transform, and pivot coordinate frame separate and auditable.
5. Render three static landscape backgrounds once, then add synchronized
   markers, optional trails, labels, and axis-limit expansion per frame.
6. Generate ChimeraX scripts from a preconfigured session without requiring
   ChimeraX in normal CI.
7. Add optional execute mode with exact model-group validation, fixed
   stationary transforms, absolute moving-model transforms, captured logs, and
   frame-count validation.
8. Preserve deterministic composition and validated
   H.264/yuv420p/faststart MP4 encoding.
9. Preserve a staged manifest that distinguishes `scripts_ready`,
   `structure_rendered`, `composite_rendered`, `movie_encoded`, and `failed`.
10. Add an INO80 P1-P6 example using placeholder paths and session setup
    instructions; do not commit large MRC or `.cxs` files.

Implementation phases:

```text
Phase 1  trajectory, canonical-to-raw mapping, and transform core
Phase 2  production-scale landscape frame renderer
Phase 3  ChimeraX script export and optional execution
Phase 4  compositor, encoder, logs, and staged manifest
Phase 5  examples, public CLI documentation, and validation report
```

Exit criteria:

```text
EA and RV paths use one SO(3) trajectory truth
raw/canonical paths resolve to the correct raw physical RO
reference transform remains unchanged
moving density rotates in the validated direction about the declared pivot
all moving transforms are absolute relative to the saved baseline
normal tests require neither ChimeraX nor FFmpeg
optional ChimeraX test uses a small asymmetric synthetic density
script-only artifacts are complete and honestly report scripts_ready
execute mode stops before encoding on render/frame validation failure
composite, trajectory, and encoded frame counts agree
source run, landscapes, selections, metadata, and maps are unchanged
manifest records every convention, frame assertion, stage, and output
documentation distinguishes rigid-body rendering from reconstruction
```

Non-goals:

```text
real-time GUI or click-to-drive interaction
automatic trajectory, axis, or pivot discovery
translation, SE(3), density morphing, or per-frame reconstruction
selection creation or reconstruction submission
raw-MRC scene construction in the first slice
ambiguous Euler, loop, or map-frame fallbacks
```

Detailed plan: `docs/animation_export.md`.

---

## 10. Documentation sprint

Goal:

```text
Make external users able to run the complete workflow without reading source code.
```

Docs to add:

```text
docs/quick_start.md
docs/cli_reference.md
docs/output_files.md
docs/relion_workflow.md
docs/cryosparc_workflow.md
docs/plotting_external_csv_tools.md
docs/faq.md
docs/animation_export.md (promote from design plan to implemented user guide)
```

Concepts that must be explained clearly:

```text
public run defaults: --ref, --mov, --output-dir, default domains, default matching, --row-aligned
default run quick-look: three Euler/RV triptych pairs plus full-data log-SLD PNG
default canonical previews: all_particles, filter_particles_by_sld_gt_1p5, and fit-top support preview
when to run cryorole align or manually pre-align metadata
align auto keys, explicit keys, duplicate policy, and low-overlap warnings
raw landscape vs canonical landscape
rotation vector vs Euler display
sld_raw vs sld_display
SLD floor stabilization vs tail-jump display-outlier detection
fit-top-fraction vs visualize top-fraction
canonical frame artifacts and --use-frame
visualization filter vs scientific selection
axis limits vs range selection
ref-domain export vs mov-domain export
source-row export backtracking
why canonicalization quick-look figures are display-only
why animation trajectories use SO(3) rather than Euler-linear interpolation
canonical waypoint back-mapping, baseline RO, ChimeraX scene basis, and pivot frame
script-only versus execute animation stages
rigid-body rendering versus independently reconstructed density
```

---

## 11. Paper-support sprint

Goal:

```text
Generate reproducible summary tables and validation outputs for manuscript figures and Methods.
```

Suggested scripts/reports:

```text
selection count table generator
canonicalization summary table
SLD distribution summary
state sampling table
run benchmark table
ref/mov swap validation summary
global frame reorientation validation summary
subsampling stability summary
```

Target datasets:

```text
hFASN validation
INO80-hexasome motion corridor
70S ribosome translocation landscape
```

---

## 12. Future 2.1+ ideas

These are valuable but should not block cryoROLE 2.0 stabilization:

```text
SE(3) translation landscape extension
translation diagnostics before full SE(3): tomo XYZ distance, SPA XY shift distance, explicit defocus-as-Z proxy policy
MPI or distributed density computation
full workflow GUI beyond the offline display-only 3D viewer
interactive ChimeraX integration beyond offline script export
canonical map/metadata transform export
advanced clustering/state annotation
optional Parquet backend
```

Sequence recommendation:

```text
single-node array-native backend first
then batch-level parallelism
then MPI/distributed support if still needed
```

---

## 13. Development rules for roadmap items

1. Preserve RO definition and convention semantics.
2. Keep interpretation-changing behavior policy-driven and reported.
3. Do not introduce hidden defaults that alter scientific meaning.
4. Keep source metadata immutable.
5. Keep display operations display-only.
6. Add tests with each behavior-changing implementation.
7. Keep documentation modular: avoid growing `AGENTS.md` and `architecture.md` with sprint-level details.
8. Keep animation interpolation in SO(3), map canonical waypoints back to raw RO
   before physical rendering, and record baseline/scene/pivot frames explicitly.
