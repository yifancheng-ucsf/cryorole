# cryoROLE 2.0 Roadmap

**Status:** Living roadmap  
**Current code state:** Architecture-complete beta with active production-scale run work.

---

## 1. Current position

cryoROLE 2.0 now has a complete command-layer workflow:

```text
run -> canonicalize -> visualize -> select -> export
```

Implemented or substantially established:

```text
Run Bundle layout
raw landscape NPZ/CSV persistence
array-native canonicalization for NPZ/run-dir inputs
canonicalize --fit-top-fraction
display-only visualization separation
legacy-compatible visualization style
SO(3)-aware selection
source-row-based STAR/CS metadata subset export
policy/report/manifest-oriented architecture
```

The next priority is not new science. The next priority is production hardening.

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

## 3. Active sprint: production-scale run

Goal:

```text
Make cryorole run robust for hundreds of thousands to millions of particles.
```

Core tasks:

1. Make the public `run` command compact: `cryorole run --ref REF --mov MOV`.
2. Default run domain labels to `ref` and `mov`, with explicit overrides when needed.
3. Default matching to `.cs` `uid` and STAR `_rlnTomoParticleName` / `_rlnImageName`; warn on reordering/dropped particles and fail unsafe matches with an align/manual-prealignment hint.
4. Add/strengthen array-native pose and matched-pose data models.
5. Vectorize RO computation.
6. Add batched kNN SLD computation.
7. Write raw NPZ directly from arrays.
8. Write raw CSV as chunked post-export.
9. Keep default raw visualization on, using legacy rainbow all-particle and `sld_raw > 1.5` preview views with display `vmax = min(max displayed SLD, 100)`; `--no-visualize` must leave downstream commands usable.
10. Keep `cryorole run --help` compact by hiding visualization and developer/debug controls from normal help.
11. Add memory/timing reporting and synthetic benchmark scripts.

Detailed plan: `docs/production_run_plan.md`.

Exit criteria:

```text
full fast test suite passes
small deterministic array-native run matches compatibility path
default run CLI is compact and documented
CS uid and STAR image-name matching are tested
RELION tomo `_rlnTomoParticleName` matching is tested
reordered/dropped matches warn and are reported
default run visualization writes all-particle and SLD-threshold 2D/3D preview outputs
run help shows only public run controls
100k synthetic benchmark completes with memory profile
canonicalize/select/export work from array-native run bundle
```

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

## 5. Following sprint: public visualize cleanup

Goal:

```text
Make cryorole visualize a compact display-only command with predictable run-bundle output.
```

Tasks:

1. Keep the public command centered on `cryorole visualize --run-dir RUN`.
2. Remove public `--output-dir`; write under `RUN/visualizations/` using `--visual-id` or `default`.
3. Default to legacy rainbow style, PNG output, fixed `sld_display` coloring, and equal display units within each figure.
4. Keep compact display controls: `--space`, `--canonical-id`, `--selection-id`, `--use-selected-landscape`, `--visual-id`, `--representation`, `--colormap`, `--range`, `--top-fraction`, `--threshold`, `--formats`, `--vmin`, `--vmax`, `--bins`, `--hist-mode`, `--kde-bandwidth`, `--xlim`, `--ylim`, and `--overwrite`.
5. Remove public Euler convention/radian/sequence controls; use the source landscape's recorded Euler convention.
6. Write basic 1D distributions by default from the same displayed rows as 2D/3D views: histogram plus KDE curve only.
7. Write selection visualizations under `visualizations/selections/<selection_id>/`, including parent-landscape and selected-landscape modes.
8. Record all display filters, range viewport policy, colormap policy, 1D distribution policy, selection provenance, and generated files in `visualization_report.json`.
9. Add focused tests for compact help, output paths, selection paths, range behavior, colormap policy, formats, 1D outputs, and display-only safety.

Non-goals:

```text
new plotting backends
interactive viewer
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

## 9. Documentation sprint

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
```

Concepts that must be explained clearly:

```text
public run defaults: --ref, --mov, default domains, default matching, --row-aligned
default run previews: all_particles and filter_particles_by_sld_gt_1p5
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
```

---

## 10. Paper-support sprint

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

## 11. Future 2.1+ ideas

These are valuable but should not block cryoROLE 2.0 stabilization:

```text
SE(3) translation landscape extension
translation diagnostics before full SE(3): tomo XYZ distance, SPA XY shift distance, explicit defocus-as-Z proxy policy
MPI or distributed density computation
interactive viewer / GUI
ChimeraX integration
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

## 12. Development rules for roadmap items

1. Preserve RO definition and convention semantics.
2. Keep interpretation-changing behavior policy-driven and reported.
3. Do not introduce hidden defaults that alter scientific meaning.
4. Keep source metadata immutable.
5. Keep display operations display-only.
6. Add tests with each behavior-changing implementation.
7. Keep documentation modular: avoid growing `AGENTS.md` and `architecture.md` with sprint-level details.
