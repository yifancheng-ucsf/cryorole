# Public `cryorole run` Plan

**Status:** Core Productionization implemented; typed run boundary and downstream array-native preparation implemented
**Goal:** Make `cryorole run` a clean public entry point that is safe, auditable, and production-scale without changing default scientific semantics.

---

## 1. Scope

`cryorole run` is the fact-generation command:

```text
input metadata
  -> extract minimal identity and pose fields
  -> resolve/match particles
  -> normalize matched poses only
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

`run` is not the place for complex alignment, outlier curation, advanced visualization design, or scientific selection. Those belong in separate commands or later explicit workflows:

```text
cryorole align      pre-run STAR alignment / match-table workflow
cryorole outlier    future explicit outlier diagnostics and handling
cryorole visualize  advanced display controls
cryorole select     scientific selections
cryorole animate    Phase 1-4 offline trajectory/landscape/ChimeraX/video export
```

Translation distance diagnostics are also out of scope for this public `run` slice. Future work may report tomo XYZ distance and single-particle XY shift distance, but it must not change RO, SLD, selection, or export defaults.

Animation work must consume completed run/canonical artifacts as a downstream
display/export workflow. It must not add trajectory, ChimeraX, frame-rendering,
or video-encoding responsibilities to `cryorole run`.
Animation and `visualize` must share display-policy helpers instead of adding
animation-specific filtering or color semantics.
Animation Phase 4 compositor/MP4 work and its stacked-layout presentation
refinement are implemented.
The downstream `canonical-views` helper may render canonical-axis camera views
from completed canonical frame artifacts; it adds no responsibility to `run`.

The production command now crosses a typed `RunRequest -> execute_run ->
RunResult` service boundary. The service retains transactional
`RunBundleWriter` ownership, uses the shared resolved input-policy service, and
performs content-level source identity verification before consuming the
preflight result. Test/fake sources require an explicit execution context and
cannot change production policy through hidden argparse attributes.

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
--output-dir RUN      default: cryorole_outputs
--row-aligned         user asserts row N matches row N
--allow-low-overlap   explicit override below the public 50% overlap threshold
--sld-metric METRIC   rotvec_euclidean|so3_geodesic; default: rotvec_euclidean
--no-visualize        skip default quick-look previews
```

Default domain labels are `ref` and `mov`. Overrides are labels/provenance only; they must not change scientific interpretation.
`--output-dir` names the run bundle created by `run`; downstream commands use that path through `--run-dir`.

---

## 3. Matching Behavior

### 3.1 Default key-based matching

When `--row-aligned` is not set:

1. CryoSPARC `.cs` inputs match by `uid`.
2. STAR inputs match by `_rlnTomoParticleName` when present, otherwise `_rlnImageName` / `rlnImageName`.
3. Matching must preserve `ref_source_row_id` and `mov_source_row_id`.
4. If matching reorders rows or drops unmatched particles, `run` emits a prominent warning and records counts.
5. Zero matches fail before pose normalization. The centralized public overlap
   threshold is 50%; lower overlap fails unless `--allow-low-overlap` is
   explicit and recorded.
6. If matching cannot be resolved safely, `run` fails before RO computation and suggests `cryorole align --key ...` when an explicit STAR key is needed, or manual pre-alignment.

Reports must record:

```text
match_key
ref_row_count
mov_row_count
matched_row_count
dropped_ref_only_count
dropped_mov_only_count
matched_rows_reordered
ref_coverage
mov_coverage
overlap_smaller_input
low_overlap_allowed
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
run_report.md
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

New bundles additionally contain `bundle_state.json` and
`.cryorole_bundle_complete`. Manifest schema 3.0 records the unique `run_id`,
streamed SHA-256 source identities, and completed transactional state. Failed
runs are never published at the requested output path.

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
visualizations/quicklook/
  all_euler_3view_projection.png
  all_rotvec_3view_projection.png
  sld_ge_1_euler_3view_projection.png
  sld_ge_1_rotvec_3view_projection.png
  top_40pct_euler_3view_projection.png
  top_40pct_rotvec_3view_projection.png
  sld_log_distribution.png
```

The directory contains exactly these seven flat PNG files. Projection pairs
show all particles, `sld_raw >= 1`, and the top 40% by `sld_raw`; cutoff ties
are retained and each range uses deterministic display-only sampling. The
log-SLD histogram uses the complete landscape and marks SLD 1, SLD 100, P99,
and the tail-jump threshold only when detected. Run does not write individual
projections, 3D figures, display tables, subdirectories, alternate formats, or
other distributions.

Default style:

```text
colormap: legacy rainbow
low SLD: red
high SLD: blue
display filter: none
display vmax cap: 100
axis units: consistent across the three panels
percentile clipping: none
```

The cap is an upper bound, not a forced value. Display `vmax` resolves to
`min(full-landscape max SLD, 100)`. The quick-look figures are not selections and
must not alter raw landscape rows, `sld_raw`, `sld_display`, reports, selection
behavior, canonicalization, or export. The quick-look policy, resolved color
bound, range counts, cutoff, and histogram diagnostics are recorded in the run
summary and manifest. More complex display
policies belong to explicit `cryorole visualize` workflows.

The bundle-root `run_report.md` concisely explains inputs, matching, policies,
output paths, quick-look meaning, SLD warnings, and next commands. It is written
even when quick-look rendering is disabled.

---

## 6. SLD and Outlier Diagnostics

`run` computes `sld_raw` using the production SLD formula and existing distance-floor stabilization.

SLD metric plan:

```text
rotvec_euclidean (default): existing Euclidean kNN in RO rotation-vector space
so3_geodesic (opt-in):      cryorole run ... --sld-metric so3_geodesic
```

The SO(3) path converts RO to unit quaternions, performs batched `cKDTree`
candidate lookup over `q` and `-q`, expands the query until `k` unique physical
neighbors are available, and computes exact distances as
`2 * arccos(|q_i dot q_j|)` before averaging.
The default path and existing SLD values must remain unchanged. Reports and
manifests record the requested/resolved metric; selected-landscape SLD
recomputation inherits it.

Implementation order:

1. Add policy, CLI validation, and provenance while protecting the default path.
2. Add a small exact SO(3) reference implementation and boundary tests.
3. Add the batched quaternion `cKDTree` path and verify reference equivalence.
4. Benchmark both metrics at 100k rows before production sign-off.

Expected average complexity remains `O(N log N + Nk)` with `O(N)` index memory;
the SO(3) constant factor is measured, and no `N x N` distance matrix is allowed.

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

Default warnings distinguish `HIGH_SLD_PRESENT` for any `sld_raw > 100`,
`DISCONTINUOUS_SLD_TAIL` for a qualifying tail jump, and
`DISTANCE_FLOOR_APPLIED`. The fixed-threshold warning communicates display
saturation and high dynamic range, not a definitive scientific outlier label.

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

Implemented production path:

1. Read source metadata without copying full source rows into per-particle dictionaries.
2. Normalize source poses into compact arrays.
3. Produce matched index/source-row arrays.
4. Compute RO in vectorized form:

```python
R_ro = np.matmul(np.swapaxes(R_ref, 1, 2), R_mov)
```

5. Compute SLD with batched metric-specific kNN query; keep RV Euclidean as the default.
6. Write `raw_landscape.npz` from arrays.
7. Write `raw_landscape.csv` as a post-computation chunked export.
8. Keep preview visualization display-only and skippable.
9. Write timing/memory reports only when requested.
10. Keep progress on stderr and final run directory on stdout.
11. Resolve `auto` to `array_native`; retain `dataframe_compat` only as an
    explicit reference backend.
12. Publish through `RunBundleWriter`: sibling staging, stage history, artifact
    validation, completion marker, atomic whole-bundle replacement, and rollback.
13. Record `run_id` and ref/mov identity (`original_path`, resolved absolute
    path, type, size, mtime, streamed SHA-256, row count). Verify identity before
    source-row export; relocation and legacy override are explicit policies.

### 8.1 Implemented array models

```text
PoseArrays -> MatchedPoseArrays -> ROArrays -> LandscapeArrays
```

RELION intrinsic-ZYZ passive-to-active conversion remains centralized and now
has a vectorized batch entry point. CryoSPARC pose conversion is batch
`Rotation.from_rotvec`. RO uses `swapaxes(ref, 1, 2) @ mov`. Raw NPZ is written
directly from arrays, and raw CSV derives Euler chunks while streaming rows.

### 8.2 Benchmark contract and memory budget

`benchmarks/benchmark_full_run.py` runs each backend in a separate subprocess
and reports wall time, peak RSS, row count, total artifact bytes, and per-file
artifact sizes. The normal sign-off size is 100k; 500k and 1M are explicit slow
runs. The initial Windows 100k array-native target is peak RSS below 1 GiB.
Query and CSV batch sizes are recorded, and focused tests verify RV query calls
do not exceed the configured batch size.

2026-08-20 local sign-off measurement (Windows 10 19045, Python 3.11.9,
synthetic native CS input, 100,000 matched rows, k=50, visualization disabled):

| Backend / query batch | Wall time | Peak RSS | Bundle bytes |
| --- | ---: | ---: | ---: |
| `array_native` / 25,000 | 3.33 s | 184.84 MiB | 53,277,982 |
| `dataframe_compat` / 25,000 | 43.86 s | 932.67 MiB | 53,275,642 |
| `array_native` / 5,000 | 3.21 s | 179.93 MiB | 53,277,981 |
| `array_native` / 39,000 | 3.23 s | 206.18 MiB | 53,277,982 |

The array-native result is below the 1 GiB 100k target. The 5,000 versus
39,000 query-batch measurements also show that the configured batch bound has
a measurable peak-RSS effect while producing the same scientific artifacts.
These are single local subprocess measurements, not cross-machine performance
guarantees. The 500k and 1M opt-in capacity runs were not part of this sign-off.

2026-08-21 service-boundary sign-off on the same Windows/Python environment
used independent 100k subprocesses. Array-native run (k=50, 25k query/CSV
batches, visualization disabled) measured 3.29 s, 191.55 MiB peak RSS, and
53,277,984 bundle bytes. The downstream benchmark measured:

| Stage | Wall time | Peak RSS | Artifact bytes |
| --- | ---: | ---: | ---: |
| preflight | 0.231 s | 135.80 MiB | report only |
| visualize (100k 2D / 50k 3D) | 12.21 s | 256.61 MiB | 11,916,553 |
| radius select (1,606 selected) | 0.266 s | 144.49 MiB | 103,729 |

These are local single-run measurements, not performance guarantees. They
confirm that downstream selection used all 100k parent rows and visualization
kept independent 2D and 3D sampling limits.

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
  omitted --sld-metric preserves current RV-Euclidean SLD exactly
  batched RV and SO(3) SLD equal their unbatched references on small datasets
  SO(3) neighbor/distance tests cover q/-q equivalence and the rotvec pi boundary

Persistence:
  raw NPZ and CSV row counts match
  match_table preserves source-row provenance
  run_summary records matching policy and preview policy
  density_report and run_manifest record the requested/resolved SLD metric
  density_report records SLD diagnostics without changing sld_raw

Visualization:
  default run writes exactly six range-specific triptychs and one full-data log-SLD PNG
  default run quick-look directory contains no tables, reports, 3D, subdirectories, or extra files
  run_summary records range counts, tie policy, max SLD, vmax cap 100, and resolved vmax
  run_report.md explains core artifacts, quick-look meaning, diagnostics, and next commands
  quick-look rendering does not create selections

Production:
  progress goes to stderr
  final run directory remains stdout
  --profile-time writes timing report
  --profile-memory writes memory report
synthetic 100k benchmarks compare both SLD metrics outside the fast test suite
full-run 100k benchmark compares array_native and dataframe_compat wall/RSS/artifacts
source relocation/hash mismatch/legacy override and selection/run-id mismatch fail safely
fault injection never publishes partial output and overwrite failure preserves the old bundle
```

---

## 10. Milestone completion state

The Core Productionization implementation is complete in code. Local sign-off
on 2026-08-20 passed the complete suite (`546 passed, 2 skipped`) and the
recorded 100k full-run benchmark above. The skipped tests require external
ChimeraX/FFmpeg integration. 500k/1M remain opt-in slow capacity checks rather
than normal CI.

### Workflow UX integration note

The subsequent Workflow UX milestone did not introduce a second run backend.
`cryorole preflight`, `cryorole run --dry-run`, and production `run` share the
array-native minimal-field pose/matching preflight service. CryoSPARC uses a
memory-mapped structured array and RELION streams the particle loop while
retaining only pose fields and the applicable identity candidates. Production run
rechecks source identity before consuming its result, then retains the same
transactional writer, vectorized RO, batched SLD, NPZ, and chunked CSV path.
The interactive explorer also reads NPZ directly and its exact radius
evaluation shares the Python selection evaluator; display sampling is never a
candidate restriction. These additions do not change RO, convention, matching,
SLD, canonicalization, selection, export, or memory contracts above.

Downstream NPZ visualization and selection now avoid the legacy full
object-coordinate Landscape expansion as their first step. Visualization
filters compact arrays before materializing plotting rows, keeps independent
2D/3D counts, and uses all filtered rows only when 1D statistics are explicitly
requested. It does not write a default display table. Selection builds only the evaluator fields required by its mode;
full selected rows are materialized only for an explicitly requested
selected-derived landscape. Legacy JSON/CSV/DataFrame paths remain compatibility
paths.
