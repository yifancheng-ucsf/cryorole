# Public `cryorole run` Plan

**Status:** Active implementation plan  
**Goal:** Make `cryorole run` a clean public entry point that is safe, auditable, and production-scale without changing scientific semantics.

---

## 1. Scope

`cryorole run` is the fact-generation command:

```text
input metadata
  -> normalize poses
  -> resolve/match particles
  -> compute RO
  -> compute SLD
  -> write raw NPZ/CSV
  -> write reports/manifest
  -> write default quick-look previews unless disabled
```

It must preserve:

```text
RO = R_ref^-1 R_mov
```

Internal rotation truth remains an active `3 x 3` rotation matrix. Source convention handling stays centralized in normalization.

`run` is not the place for complex alignment, outlier curation, advanced visualization design, or scientific selection. Those belong in future or separate commands:

```text
cryorole align      future complex pre-alignment / match-table workflow
cryorole outlier    future explicit outlier diagnostics and handling
cryorole visualize  advanced display controls
cryorole select     scientific selections
```

Translation distance diagnostics are also out of scope for this public `run` slice. Future work may report tomo XYZ distance and single-particle XY shift distance, but it must not change RO, SLD, selection, or export defaults.

---

## 2. Public CLI

Default public command:

```bash
cryorole run --ref REF_METADATA --mov MOV_METADATA
```

Public controls should stay small:

```text
--ref-domain NAME     default: ref
--mov-domain NAME     default: mov
--row-aligned         user asserts row N matches row N
--no-visualize        skip default quick-look previews
```

Default domain labels are `ref` and `mov`. Overrides are labels/provenance only; they must not change scientific interpretation.

---

## 3. Matching Behavior

### 3.1 Default key-based matching

When `--row-aligned` is not set:

1. CryoSPARC `.cs` inputs match by `uid`.
2. STAR inputs match by `_rlnTomoParticleName` when present, otherwise `_rlnImageName` / `rlnImageName`.
3. Matching must preserve `ref_source_row_id` and `mov_source_row_id`.
4. If matching reorders rows or drops unmatched particles, `run` emits a prominent warning and records counts.
5. If matching cannot be resolved safely, `run` fails before RO computation and suggests `cryorole align` or manual pre-alignment.

Reports must record:

```text
match_key
ref_row_count
mov_row_count
matched_row_count
dropped_ref_only_count
dropped_mov_only_count
matched_rows_reordered
warnings
```

### 3.2 Explicit row-aligned mode

`--row-aligned` is a user assertion. It applies to STAR and CS inputs.

In this mode, `run` must:

1. Check only that ref/mov row counts match.
2. Pair row `N` with row `N`.
3. Not key-match.
4. Not reorder rows.
5. Not drop particles.
6. Fail clearly if row counts differ.
7. Record the row-aligned policy in reports/manifests.

---

## 4. Default Outputs

Required run bundle artifacts:

```text
run_manifest.json
run_summary.json
data/raw_landscape.npz
data/raw_landscape.csv
data/match_table.csv
reports/import_ref_report.json
reports/import_mov_report.json
reports/identity_ref_report.json
reports/identity_mov_report.json
reports/match_report.json
reports/density_report.json
```

`data/raw_landscape.npz` is the production machine-readable raw landscape.  
`data/raw_landscape.csv` is the user-facing raw table.  
Full landscape JSON is debug-only and opt-in.

Raw landscape rows must include source-row provenance for export:

```text
particle_key
ref_source_row_id
mov_source_row_id
coordinates_analysis
sld_raw and SLD diagnostics
```

---

## 5. Default Quick-Look Previews

Unless `--no-visualize` is set, `run` writes display-only previews under:

```text
visualizations/raw_default/
  visualization_report.json
  all_particles/
  filter_particles_by_sld_gt_1p5/
```

Each preview group should include the basic 2D projections and 3D views:

```text
euler_alpha_beta.png
euler_beta_gamma.png
euler_alpha_gamma.png
euler_3view_projection.png
rotvec_xy.png
rotvec_yz.png
rotvec_xz.png
rotvec_3view_projection.png
landscape_3d_euler.png
landscape_3d_rotvec.png
display_table.csv
```

Default style:

```text
colormap: legacy rainbow
low SLD: red
high SLD: blue
preview cutoff: sld_raw > 1.5
display vmax cap: 100
axis units: consistent across 2D and 3D
percentile clipping: none
```

The cap is an upper bound, not a forced value.

The `filter_particles_by_sld_gt_1p5` group is a preview, not a selection. Display `vmax` resolves to `min(max displayed SLD, 100)`. Display filtering must not alter raw landscape rows, `sld_raw`, `sld_display`, reports, selection behavior, canonicalization, or export.

Directory labels should be filesystem-safe. Future examples:

```text
filter_particles_by_sld_gt_2p0
filter_particles_by_top_sld_40pct
```

Exact thresholds and fractions belong in `visualization_report.json`.

---

## 6. SLD and Outlier Diagnostics

`run` computes `sld_raw` using the production SLD formula and existing distance-floor stabilization.

Extreme high-SLD detection may remain as an internal/reporting diagnostic:

```text
sld_display_is_outlier
n_sld_display_outliers
fraction_sld_display_outliers
display_outlier_threshold_sld
largest_tail_jump_ratio
max_sld_raw
max_over_display_vmax
near-identity / near-duplicate warnings
```

This diagnostic must not:

1. Change `sld_raw`.
2. Rewrite `sld_display` merely to clip colors.
3. Drop rows.
4. Exclude particles from scientific selections by default.
5. Become public `run` CLI tuning surface.

Future explicit outlier workflows should live outside the default `run` command.

---

## 7. Developer Controls

These controls are useful for implementation, benchmarking, and debugging, but should be hidden from normal `cryorole run --help`:

```text
--run-backend auto|array_native|dataframe_compat
--quiet
--verbose
--profile-time
--profile-memory
--raw-csv-chunk-size N
--density-query-batch-size N
--no-raw-csv
--write-debug-json
```

Do not expose SLD tail-jump tuning as public `run` options. Keep those as internal/reporting policy until a dedicated outlier or visualization workflow needs them.

---

## 8. Production Backend Checklist

Implementation should move toward array-native production paths:

1. Read source metadata without copying full source rows into per-particle dictionaries.
2. Normalize source poses into compact arrays.
3. Produce matched index/source-row arrays.
4. Compute RO in vectorized form:

```python
R_ro = np.matmul(np.swapaxes(R_ref, 1, 2), R_mov)
```

5. Compute SLD with batched kNN query.
6. Write `raw_landscape.npz` from arrays.
7. Write `raw_landscape.csv` as a post-computation chunked export.
8. Keep preview visualization display-only and skippable.
9. Write timing/memory reports only when requested.
10. Keep progress on stderr and final run directory on stdout.

---

## 9. Acceptance Tests

Minimum tests for this plan:

```text
CLI:
  run requires --ref and --mov
  default domains are ref/mov and are recorded
  explicit domain labels are recorded
  --no-visualize skips preview artifacts

Matching:
  CS defaults to uid matching
  STAR defaults to tomo-particle-name matching when present, otherwise image-name matching
  key matching reorders safely and warns
  key matching with dropped particles warns and records counts
  unsafe STAR matching fails before RO computation
  --row-aligned works for STAR and CS when row counts match
  --row-aligned fails on row-count mismatch

Science:
  RO matches the reference implementation
  RELION/CryoSPARC convention tests still pass
  vectorized RO equals compatibility path on deterministic cases
  batched SLD equals unbatched reference on small datasets

Persistence:
  raw NPZ and CSV row counts match
  match_table preserves source-row provenance
  run_summary records matching policy and preview policy
  density_report records SLD diagnostics without changing sld_raw

Visualization:
  default run writes all_particles previews
  default run writes filter_particles_by_sld_gt_1p5 previews
  preview report records cutoff 1.5, max displayed SLD, vmax cap 100, and resolved vmax
  preview filters do not create selections

Production:
  progress goes to stderr
  final run directory remains stdout
  --profile-time writes timing report
  --profile-memory writes memory report
  synthetic 100k benchmark completes outside the fast test suite
```

---

## 10. First Implementation Task

Start with the smallest behavior-changing slice:

```text
Implement the compact public run behavior and reports without changing RO math.
```

Deliverables:

```text
--ref / --mov required
default ref/mov domain labels
--row-aligned for STAR and CS with row-count preflight
default CS uid matching
default STAR image-name matching
match warning/report fields
default preview policy fields
--no-visualize still skips previews
focused CLI and matching tests
```

After that, continue with array-native pose extraction, vectorized RO, batched SLD, chunked CSV, and benchmark reporting.
