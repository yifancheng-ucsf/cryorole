# cryoROLE 2.0 Architecture Contract

**Status:** Stable contract with active production-scale `run` work in progress  
**Audience:** cryoROLE developers, Codex contributors, advanced users  
**Companion documents:**

- `AGENTS.md` — short non-negotiable guardrails.
- `docs/production_run_plan.md` — current production-scale `run` implementation plan.
- `docs/roadmap.md` — release and development priorities.

---

## 1. Executive summary

cryoROLE 2.0 is a policy-driven cryo-EM relative-orientation analysis platform. It compares two independently refined near-rigid domains particle-by-particle and produces a reproducible relative-orientation landscape.

The architecture is organized around a stable **Run Bundle** and user-facing command layers:

```text
align = prepare row-aligned metadata when default run matching is not enough
run = compute raw facts
canonicalize = derive a canonical coordinate frame
visualize = render display-only figures and display tables
select = create auditable particle selections
export = bridge selections or transforms back to external tools
```

The long-term persistence contract is:

```text
NPZ = machine-readable numeric arrays
CSV = user-facing flat tables
JSON = reports, summaries, manifests, and optional debug artifacts
```

Full object-record landscape JSON is not a production-scale persistence format.

---

## 2. Core scientific contract

### 2.1 Relative orientation

cryoROLE preserves the relative-orientation definition:

```text
RO = R_ref^-1 R_mov
```

`R_ref` and `R_mov` are normalized source poses represented internally as active `3 x 3` rotation matrices.

### 2.2 Source conventions

RELION and CryoSPARC metadata conventions are normalized before downstream computation.

RELION rules:

- `Rot/Tilt/Psi` are intrinsic `ZYZ` Euler angles.
- SciPy parsing uses uppercase `"ZYZ"`.
- Passive-to-active conversion is centralized and tested.
- Public `cryorole run` defaults STAR particle identity to `_rlnTomoParticleName` when present, otherwise `_rlnImageName` / `rlnImageName`, when the resolved column is unambiguous.
- If default STAR matching is unavailable or unsafe, `run` must fail before RO computation unless the user explicitly requests row-aligned mode or another supported identity policy.

CryoSPARC rules:

- `.cs` is a native input format.
- Default identity key is `uid`.
- Native pose field is `alignments3D/pose`.
- Pose encoding is rotation vector / axis-angle.
- Internal active matrices are derived with `Rotation.from_rotvec(...).as_matrix()`.

### 2.3 Representation truth

Internal rotation truth is a `3 x 3` active rotation matrix. Rotation vector, quaternion, and Euler coordinates are derived representations.

Rotation vectors are the primary compact analysis coordinate. ZYX Euler angles in degrees are a display/interpretation coordinate.

Public user-facing Euler convention for RO/RV-derived output is extrinsic fixed-axis ZYX. Intrinsic ZYX is legacy/internal compatibility only. Reports, manifests, visualization reports, and selection artifacts must record the resolved Euler convention whenever Euler coordinates are derived or interpreted.

CSV column names intentionally remain compact, for example `raw_ea_zyx_alpha_deg`. The intrinsic/extrinsic interpretation is owned by recorded policy metadata, not by expanded CSV column names.

### 2.4 Analysis versus display

Analysis-space values and display-space values must remain separate.

| Concept | Analysis role | Display role |
|---|---|---|
| rotation vector | statistics, density, canonicalization, geodesic reasoning | may be plotted directly |
| Euler angles | derived coordinate | user-facing labels/figures |
| `sld_raw` | scientific density/threshold/canonical fitting | can be color-mapped |
| `sld_display` | none by default | display-only density/color value, default equal to `sld_raw` |
| `sld_display_is_outlier` | optional explicit artifact-exclusion flag for selection | marks tail-jump display outliers |
| axis limits | none | figure viewport only |
| display filters | none | displayed rows only |
| selection radius | scientific decision | reported, not display-only |

Display reparameterizations must not be written back into external metadata as physical poses.

Tail-jump display-outlier detection must not alter `sld_raw`, raw/canonical rows, source metadata, or canonicalization inputs. It only changes automatic display color scaling and must be recorded as display policy. Display outliers remain plotted; they simply do not determine the automatic colorbar maximum.

---

## 3. Run Bundle contract

A standard run directory is:

```text
<run_dir>/
  run_manifest.json
  run_summary.json

  data/
    raw_landscape.npz
    raw_landscape.csv
    match_table.csv

  reports/
    import_ref_report.json
    import_mov_report.json
    identity_ref_report.json
    identity_mov_report.json
    match_report.json
    density_report.json

  visualizations/
    raw_default/
      visualization_report.json
      all_particles/
        display_table.csv
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
      filter_particles_by_sld_gt_1p5/
        display_table.csv
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
    canonical/
      <canonical_id>/
        visualization_report.json
        all_particles/
          display_table.csv
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
        filter_particles_by_sld_gt_1p5/
          display_table.csv
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
        filter_particles_by_top_sld_40pct/
          display_table.csv
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

  canonical/
    <canonical_id>/
      canonical_landscape.npz
      canonical_landscape.csv
      canonicalization_report.json
      canonicalize_summary.json
      canonical_frame.json
      canonical_frame.npz
  selections/
  exports/
  debug/
```

Required production artifacts:

| Artifact | Role |
|---|---|
| `run_manifest.json` | top-level provenance and artifact index |
| `run_summary.json` | human-readable run summary |
| `data/raw_landscape.npz` | machine-readable raw landscape |
| `data/raw_landscape.csv` | user-facing raw landscape table |
| `data/match_table.csv` | matched-particle provenance |
| `reports/density_report.json` | SLD diagnostics |
| `reports/match_report.json` | matching diagnostics |

`debug/landscape_debug.json` may exist only when explicitly requested.

---

## 4. Align contract

`cryorole align` prepares metadata files that are safe to use with
`cryorole run --row-aligned`. It is a pre-run matching and filtering utility,
not an RO, SLD, selection, or export command.

First public target:

```bash
cryorole align --ref REF.star --mov MOV.star
```

Default output:

```text
alignments/<align_id>/
  aligned_ref.star
  aligned_mov.star
  match_table.csv
  align_report.json
  ref_only.star
  mov_only.star
  duplicate_ref.star
  duplicate_mov.star
```

Small public controls:

```text
--align-id ID                  default: default
--key COL [COL ...]            explicit STAR identity key columns
--float-tol COL=TOL            quantize a numeric key column; may repeat
--path-mode exact|basename|suffix:N
--duplicate-policy exclude|first   default: exclude
--overwrite
```

Rules:

1. Initial public `align` support is STAR-focused.
2. Original STAR inputs must never be modified.
3. Only the particle loop is filtered/reordered; optics and other STAR blocks are preserved from each input.
4. `aligned_mov.star` must be reordered to match `aligned_ref.star` row by row.
5. Auto key selection is conservative: use `_rlnTomoParticleName` when present in both files; otherwise use `_rlnImageName` / `rlnImageName` when present in both files.
6. If no safe auto key exists, fail clearly and ask the user to provide `--key`.
7. Explicit `--key` may use any STAR columns, including defocus columns, but non-identity columns such as defocus must be recorded with a warning because they are user-chosen tracking keys, not cryoROLE defaults.
8. Numeric key columns should use explicit `--float-tol`; this is especially important for coordinate or defocus keys.
9. If the resolved key has zero overlap, fail clearly. If overlap is below 50% but nonzero, continue and emit a prominent warning.
10. Duplicate keys are detected by building the normalized key for every row and counting repeated keys within each input.
11. The default duplicate policy is `exclude`: any key duplicated in either input is excluded from aligned outputs and written to duplicate STAR files.
12. `--duplicate-policy first` keeps the first occurrence for each duplicated key and excludes later redundant rows; it must warn and report excluded duplicate counts.
13. `match_table.csv` must record source-row provenance and aligned-row order for downstream auditing.
14. `align_report.json` must record inputs, key policy, normalization policy, duplicate policy, row counts, matched count, overlap fractions, dropped counts, duplicate counts, warnings, and output paths.

Recommended use after alignment:

```bash
cryorole run \
  --ref alignments/default/aligned_ref.star \
  --mov alignments/default/aligned_mov.star \
  --row-aligned
```

---

## 5. Production `run` contract

`cryorole run` is the fact layer. It reads inputs, normalizes conventions, resolves identities, matches particles, computes RO, computes SLD, writes raw artifacts, and emits run provenance.

Public default command:

```bash
cryorole run --ref REF_METADATA --mov MOV_METADATA
```

Small public controls:

```text
--ref-domain NAME     default: ref
--mov-domain NAME     default: mov
--row-aligned         explicit row-order assumption; requires equal row counts
--no-visualize        skip raw quick-look visualizations
```

`--row-aligned` is a user assertion that the two metadata files are already particle-aligned. It applies to STAR and CS inputs. In this mode `run` must only check that row counts match, then pair row `N` with row `N`; it must not key-match, reorder rows, or drop particles.

Without `--row-aligned`, default matching is key-based:

1. CryoSPARC `.cs` inputs use `uid`.
2. STAR inputs use `_rlnTomoParticleName` when present, otherwise `_rlnImageName` / `rlnImageName`.
3. Matching reports must record the resolved key, input row counts, matched count, dropped ref-only/mov-only counts, and whether rows were reordered.
4. If matching succeeds but reorders rows or drops particles, `run` must emit a prominent warning and continue only with matched pairs.
5. If matching cannot be resolved safely, `run` must fail clearly and suggest `cryorole align` or manual pre-alignment.

Production-scale `run` must move toward an array-native backend:

```text
source metadata
  -> normalized pose arrays
  -> matched index arrays
  -> vectorized RO arrays
  -> batched SLD arrays
  -> raw_landscape.npz
  -> chunked raw_landscape.csv
  -> default raw visualization unless --no-visualize
```

Required production behavior:

1. Raw landscape is the source of truth and must not be overwritten by downstream commands.
2. `raw_landscape.npz` is written from compact arrays.
3. `raw_landscape.csv` is a post-computation user-facing export and should be chunked/streaming.
4. RO computation should use vectorized matrix multiplication:

```python
R_ro = np.matmul(np.swapaxes(R_ref, 1, 2), R_mov)
```

5. `cryorole run` records the Euler convention used for user-facing Euler columns. Public output uses extrinsic fixed-axis ZYX.
6. SLD density must support batched kNN query to avoid materializing large `n x k` intermediate arrays.
7. Raw visualization is generated by default as display-only quick-look output and must be disableable with `--no-visualize`.
8. Full source-row metadata should not be copied into per-particle Python objects in the production path.

Default raw quick-look visualization should use the legacy style: rainbow color mapping with low SLD rendered red and high SLD rendered blue, common axis units across 2D and 3D figures, and no percentile clipping. It writes two preview groups by default:

```text
all_particles                  all matched rows
filter_particles_by_sld_gt_1p5 rows with sld_raw > 1.5
```

These are display previews, not selections. The default run color scale uses `vmax = min(max displayed SLD, 100)`. Rows above `vmax` remain plotted with saturated high-density color, and stored `sld_raw` / `sld_display` values are not changed. This style is a `run` convenience; refined display policies belong in `cryorole visualize`.

Filtered visualization directories should use filesystem-safe labels such as `filter_particles_by_sld_gt_2p0` or `filter_particles_by_top_sld_40pct`; exact thresholds and fractions belong in `visualization_report.json`.

Translation distance is not part of the current default `run` analysis coordinate and must not affect RO, SLD, canonicalization, selection, or export by default.

Future translation diagnostics may add report/raw columns only under an explicit policy:

```text
tomo STAR:    XYZ distance from centered coordinate fields when available
single-particle STAR: XY shift distance from origin fields, with pixel-size conversion when needed
single-particle defocus-as-Z: future explicit proxy policy only
```

These diagnostics must be recorded as metadata/report fields, not interpreted as relative orientation.

Detailed plan and acceptance criteria are in `docs/production_run_plan.md`.

---

## 6. Landscape schema

### 6.1 Raw landscape NPZ

Required arrays:

```text
particle_key                  string[n]
ref_source_row_id              int64[n]
mov_source_row_id              int64[n]
coordinates_analysis           float64[n, 3]   # raw RO rotation vector, radians
sld_unfloored                  float64[n]
sld_raw                        float64[n]
sld_display                    float64[n]
sld_display_is_outlier         bool[n]
sld_was_floored                bool[n]
sld_local_k_mean               float64[n]
sld_effective_local_k_mean     float64[n]
sld_distance_floor             float64[n]
```

Optional arrays may include `ro_angle_rad`, `matched_ref_index`, `matched_mov_index`, or other explicitly documented numeric arrays.

Compatibility readers may tolerate older transitional artifacts containing `sld_display_was_clipped`, but the stable field name for new artifacts is `sld_display_is_outlier`. New writers should not use percentile-clipped `sld_display` as the persistence contract.

### 6.2 Raw landscape CSV

Recommended columns:

```text
particle_key
ref_source_row_id
mov_source_row_id
raw_rv_x_rad
raw_rv_y_rad
raw_rv_z_rad
raw_angle_deg
raw_ea_zyx_alpha_deg
raw_ea_zyx_beta_deg
raw_ea_zyx_gamma_deg
sld_unfloored
sld_raw
sld_display
sld_display_is_outlier
sld_was_floored
sld_local_k_mean
sld_effective_local_k_mean
sld_distance_floor
```

Column names should include representation and units. The `zyx` label names the axis order only; the intrinsic/extrinsic convention is recorded in reports/manifests. Backward-compatible aliases may exist but are not the long-term contract.

### 6.3 Canonical landscape

Canonical artifacts live under:

```text
canonical/<canonical_id>/
  canonical_landscape.npz
  canonical_landscape.csv
  canonicalization_report.json
  canonicalize_summary.json
  canonical_frame.json
  canonical_frame.npz
```

Canonical CSV preserves raw coordinates and adds canonical coordinates:

```text
raw_rv_x_rad
raw_rv_y_rad
raw_rv_z_rad
raw_ea_zyx_alpha_deg
raw_ea_zyx_beta_deg
raw_ea_zyx_gamma_deg
canonical_rv_x_rad
canonical_rv_y_rad
canonical_rv_z_rad
canonical_ea_zyx_alpha_deg
canonical_ea_zyx_beta_deg
canonical_ea_zyx_gamma_deg
```

Canonical CSV Euler columns use the public cryoROLE convention: extrinsic fixed-axis ZYX.

---

## 7. SLD density contract

cryoROLE 2.0 uses **SLD** as the public density metric:

```text
SLD = kNN-scaled Local Density
```

Default formula:

```text
local_k_mean_i = mean distance from p_i to its k nearest neighbors
global_local_k_mean = mean(local_k_mean_i over all points)
distance_floor = distance_floor_fraction * global_local_k_mean
effective_local_k_mean_i = max(local_k_mean_i, distance_floor)
sld_unfloored_i = global_local_k_mean / local_k_mean_i
sld_raw_i = global_local_k_mean / effective_local_k_mean_i
```

Field semantics:

| Field | Meaning |
|---|---|
| `sld_unfloored` | diagnostic ratio without distance floor |
| `sld_raw` | production density with floor stabilization |
| `sld_display` | display-only density/color value, default equal to `sld_raw` |
| `sld_display_is_outlier` | whether tail-jump display policy marked this row as an automatic color-scale outlier |
| `sld_was_floored` | whether floor was applied |
| `sld_local_k_mean` | original local kNN mean distance |
| `sld_effective_local_k_mean` | floor-stabilized distance |
| `sld_distance_floor` | floor value used |

Do not apply log transforms to `sld_raw`. Display transforms must write to display-only fields or display tables.

The SLD distance floor and display-outlier detection are separate mechanisms:

```text
distance floor = numerical stabilization during SLD computation
tail-jump display-outlier detection = high-SLD diagnostic after `sld_raw` exists
```

The distance floor prevents zero or near-zero kNN distances from producing infinite or unbounded `sld_raw`. Tail-jump display-outlier detection marks extreme high-density tails for reports and optional display policy. Display-outlier detection must not change `sld_raw`, `sld_display`, `sld_unfloored`, `sld_was_floored`, or `sld_local_k_mean`.

Default display-outlier policy:

```text
display_outlier_detection_mode = tail_jump
tail_search_fraction = 0.01
tail_jump_factor = 5.0
max_display_outlier_fraction = 0.002
```

Tail-jump detection sorts finite `sld_raw` values ascending, searches only the highest `tail_search_fraction` of values, and looks for an adjacent jump where:

```text
sld_after_jump / sld_before_jump >= tail_jump_factor
n_after_jump / n_total <= max_display_outlier_fraction
```

If a qualifying jump is found, rows above the jump are marked `sld_display_is_outlier = True`. If no qualifying jump is found, no display outliers are marked. Default `run` quick-look figures still use their fixed display cap policy rather than deleting or rewriting outlier rows.

Reports must record the display-outlier mode, parameters, detected threshold, outlier count, outlier fraction, display cap when used, and enough high-density diagnostics to distinguish true dense states from near-duplicate artifacts. Suggested diagnostics include `max_sld_raw`, `p99_sld_raw`, `p99_5_sld_raw`, `max_over_display_vmax`, `largest_tail_jump_ratio`, `display_outlier_threshold_sld`, `n_sld_display_outliers`, `fraction_sld_display_outliers`, and near-duplicate/near-identity RO warnings when applicable.

---

## 8. Canonicalization contract

Canonicalization derives a coordinate frame for interpretation and comparison. It does not change raw RO facts.

Default command:

```bash
cryorole canonicalize --run-dir RUN
```

Small public controls:

```text
--canonical-id ID       default: default
--fit-top F             alias for --fit-top-fraction; default: 0.40
--fit-top-fraction F    high-sld_raw fraction used to fit the frame
--positive-side low|high
--use-frame FRAME       apply an existing canonical frame instead of fitting
--no-visualize          skip canonical quick-look figures
```

Rules:

1. Canonicalization must read from a run bundle and write under `canonical/<canonical_id>/`.
2. It must not overwrite raw artifacts.
3. `canonical_landscape.npz` is the machine-readable canonical landscape; `canonical_landscape.csv` is the user-facing table.
4. Public canonicalization uses extrinsic fixed-axis ZYX for derived Euler columns.
5. `--fit-top-fraction F` / `--fit-top F` controls the high-`sld_raw` subset used for frame fitting. Valid range: `0 < F <= 1`.
6. Public axis sign disambiguation uses density-weighted skewness. The default positive side is `low_density_skew`; public `--positive-side low|high` maps to `low_density_skew` / `high_density_skew`.
7. Canonicalization writes display-only quick-look previews by default under `visualizations/canonical/<canonical_id>/`; `--no-visualize` disables this.
8. Canonicalization reports must record backend, fit subset, fit top fraction, axis assignment, sign rule, positive-side convention, handedness, transform direction, warnings, frame artifact paths, visualization status, and output paths.

Canonical quick-look previews are:

```text
all_particles
filter_particles_by_sld_gt_1p5
filter_particles_by_top_sld_XXpct
```

`filter_particles_by_top_sld_XXpct` uses the effective `fit_top_fraction` recorded by the canonicalization policy. For default fitting this is `filter_particles_by_top_sld_40pct`; if the user fits with `--fit-top 0.50`, it becomes `filter_particles_by_top_sld_50pct`. When `--use-frame` is used, the preview uses the fit fraction recorded in the imported frame metadata. These are previews, not selections.

Default axis assignment:

```text
PC1 -> alpha / primary motion
PC2 -> beta  / secondary motion
PC3 -> gamma / tertiary motion
```

The implementation may map these through canonical RV axes before deriving ZYX Euler display coordinates, but the final user-facing effect must be tested: largest spread appears along alpha, second along beta, third along gamma.

Default sign disambiguation is density-weighted skewness using `sld_raw`:

```text
positive_side = low_density_skew
```

Density-weighted skewness defines a reproducible positive side; it does not define an absolute biological clockwise/counterclockwise direction.

Each canonicalization must write a reusable frame:

```text
canonical_frame.json    audited metadata and paths
canonical_frame.npz     canonical_transform numeric array
```

Minimum frame fields:

```text
canonical_transform
transform_direction = canonical_rv = raw_rv @ canonical_transform
coordinate_space = rotvec_ro_radians
axis_assignment
sign_rule = density_weighted_skewness
positive_side
fit_top_fraction
density_support_field = sld_raw
euler_convention = extrinsic_zyx
source_run_dir
source_canonical_id
schema_version
```

When `--use-frame FRAME` is provided, canonicalize must validate the frame, apply its `canonical_transform`, record frame provenance, and skip PCA fitting/sign selection. Fitting controls such as `--fit-top` and `--positive-side` should conflict with `--use-frame` rather than silently changing the imported frame.

---

## 9. Visualization contract

Visualization is display-only. It must not alter raw/canonical landscapes or create selections.

Default visualization writes:

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
visualization_report.json
distributions_1d/
  euler_alpha_distribution.png
  euler_beta_distribution.png
  euler_gamma_distribution.png
  rotvec_x_distribution.png
  rotvec_y_distribution.png
  rotvec_z_distribution.png
  distribution_1d_report.json
  distribution_1d_stats.csv
```

Public command:

```bash
cryorole visualize --run-dir RUN
```

Small public controls:

```text
--space raw|canonical          default: raw
--canonical-id ID              default: default
--selection-id ID
--use-selected-landscape
--visual-id ID                 default: default
--representation both|euler|rotvec   default: both
--colormap NAME                default: rainbow_r
--range AXIS:LOWER:UPPER       display-only row/range filter; may repeat
--top-fraction F               display top fraction by sld_display
--threshold T                  display rows above density threshold
--formats png[,svg,pdf]        default: png
--vmin V
--vmax V
--bins N                       1D histogram bins; default: 72
--hist-mode count|percent      1D histogram mode; default: count
--kde-bandwidth VALUE          default: scott; scott, silverman, or positive float
--xlim MIN:MAX                 1D x-axis display limit
--ylim MIN:MAX                 1D y-axis display limit
--overwrite
```

Rules:

1. Users do not choose `--output-dir`; visualization output stays under `RUN/visualizations/`.
2. Default style is legacy rainbow: low density is red, high density is blue, and display axes use equal units within each 2D/3D figure so plots do not distort.
3. Public `visualize` colors by `sld_display`. Users may change the display colormap with `--colormap`, but the color field is not a public workflow choice.
4. Default format is PNG; SVG/PDF are opt-in through `--formats`.
5. `--top-fraction` and `--threshold` are display filters, not selections, and must be mutually exclusive.
6. `--range` is a display-only filter and viewport policy. Euler ranges filter the displayed row set and set Euler axis limits for matching axes. Rotvec plots show the same filtered rows in rotvec units rather than converting Euler bounds into fake rotvec bounds. Rotvec ranges likewise set rotvec axis limits only.
7. Public visualization does not expose Euler convention controls. Euler plots use the source landscape's recorded convention, which for public run/canonicalize outputs is extrinsic fixed-axis ZYX.
8. Basic 1D distributions are written by default from the same displayed rows as the 2D/3D views. They include only a histogram and KDE curve; no mean line, peak finding, overlay comparison, or selection is created.
9. 1D outputs follow `--representation`: `euler` writes alpha/beta/gamma, `rotvec` writes x/y/z, and `both` writes all six. `--formats` applies to 1D figures too.

Default output paths:

```text
visualizations/raw/<visual_id>/
visualizations/canonical/<canonical_id>/<visual_id>/
```

Selection-scoped visualization consumes an existing selection with:

```bash
cryorole visualize --run-dir RUN --selection-id SELECTION_ID
```

Selection-scoped visualization is a display-only row subset. It must read particle keys from `selections/<selection_id>/selected_particle_keys.csv` when present, fall back to `selections/<selection_id>/selection.json`, filter the requested raw or canonical landscape by `particle_key`, and then run the normal visualization pipeline. It must not create, modify, or reinterpret selections.

Default selection-scoped output paths:

```text
visualizations/selections/<selection_id>/parent_raw/<visual_id>/
visualizations/selections/<selection_id>/parent_canonical/<canonical_id>/<visual_id>/
```

Selection-scoped visualization reports must record `selection_id`, `selected_count`, and `total_count`, plus enough provenance to identify the selection source artifact and the row count before and after selection filtering. Missing or inconsistent selection keys must be reported clearly and must not be silently dropped.

Visualization may also consume a selected-derived landscape when one exists and the user explicitly requests it, for example:

```bash
cryorole visualize --run-dir RUN --selection-id SELECTION_ID --use-selected-landscape
```

In this mode, visualization reads `selections/<selection_id>/selected_landscape/landscape.npz` or equivalent selected-landscape artifacts instead of filtering the parent raw/canonical landscape at display time. It writes under:

```text
visualizations/selections/<selection_id>/selected_landscape/<visual_id>/
```

Reports must record `selected_landscape_used = true`, the selected-landscape path, the parent landscape recorded by the selected landscape, and whether SLD values were inherited or recomputed. This remains display-only; visualization must not create, modify, or reinterpret selections or parent landscapes.

Visualization reports must record source landscape, coordinate source, representation, inherited Euler convention, fixed display color field, colormap, display density field, display filters, top fraction, threshold, ranges, format list, color scale, 1D distribution policy, selection provenance when used, and generated files.

When visualization uses automatic SLD color scaling, it may exclude rows marked `sld_display_is_outlier` from automatic colorbar scaling by default. Outlier points must remain plotted and present in display tables; values above the chosen `vmax` may render with the colormap's saturated top color. Reports must make display-outlier handling visible via fields such as:

```text
color_scale_mode = tail_jump_exclude_outliers
tail_search_fraction = 0.01
tail_jump_factor = 5.0
max_display_outlier_fraction = 0.002
display_color_vmax
n_sld_display_outliers
fraction_sld_display_outliers
```

---

## 10. Selection contract

Selection is a first-class scientific decision object, distinct from display filtering.

Standard output:

```text
selections/<selection_id>/
  selection.json
  selection.csv
  selected_particle_keys.csv
  selected_landscape_rows.csv
  selection_summary.json
  selected_landscape/              # optional
    landscape.npz
    landscape.csv
    landscape_report.json
```

Rules:

1. Public selection is centered on `cryorole select --run-dir RUN --selection-id ID`.
2. Users do not choose `--output`; outputs always live under `selections/<selection_id>/`.
3. `--overwrite` may be used to replace an existing `selections/<selection_id>/` output directory. It must not modify raw/canonical landscapes or other selections, and the overwrite policy must be recorded.
4. Default mode is radius selection around a center. The compact public form is:

```bash
cryorole select \
  --run-dir RUN \
  --selection-id center_sel \
  --space canonical \
  --center 0 0 0 \
  --radius 15
```

`--center` has alias `-c`; `--radius` has alias `-r`. By default, center input is Euler degrees, the center representation is `euler`, and radius selection uses SO(3) geodesic distance in degrees. `--radius-rad` may be used instead of `--radius`.

Small public controls:

```text
--run-dir RUN                       existing run bundle; required
--selection-id ID                   output id under RUN/selections/; required
--space raw|canonical               parent landscape space; default: raw
--canonical-id ID                   canonical landscape id when --space canonical; default: default
--mode radius|threshold|range|random|metadata
--center / -c A B C                 radius center; default representation is Euler degrees
--radius / -r DEG                   radius in degrees
--radius-rad RAD                    radius in radians; mutually exclusive with --radius
--center-representation euler|rotvec   default: euler
--metric so3|rotvec                 radius metric; default: so3
--sld-min VALUE                     lower sld_raw bound for threshold mode
--sld-max VALUE                     upper sld_raw bound for threshold mode
--range-bound AXIS:LOWER:UPPER      coordinate range bound; may repeat
--fraction F                        random selection fraction; 0 < F <= 1
--seed SEED                         optional random seed
--metadata-domain ref|mov           source metadata domain for metadata mode
--metadata-column COLUMN            source metadata column for metadata mode
--metadata-value VALUE[,VALUE...]   one or more values to include
--split-by-value                    write one child selection per metadata value
--write-selected-landscape          also write a selected-derived landscape for visualization
--recompute-sld                     recompute SLD within the selected-derived landscape
--overwrite                         replace existing artifacts for the same selection id
```

Public select does not expose center Euler override flags. Euler center input uses the parent landscape's recorded Euler convention.

5. Naive Euler-space Euclidean distance must not be the default scientific selection metric.
6. Scientific density thresholding uses `sld_raw`, not `sld_display`. The compact threshold mode supports lower and/or upper bounds:

```bash
cryorole select \
  --run-dir RUN \
  --selection-id sld_window \
  --mode threshold \
  --sld-min 1.5 \
  --sld-max 10
```

7. Range selection is an explicit coordinate policy and must be reported. The compact public form is:

```bash
cryorole select \
  --run-dir RUN \
  --selection-id alpha_window \
  --mode range \
  --space canonical \
  --range-bound alpha:-20:20
```

Public range selection does not expose coordinate-source, representation, Euler sequence, or radians override flags. Euler range interpretation uses the parent landscape's recorded Euler convention.

8. Random selection is an explicit scientific selection artifact, not visualization downsampling. The compact public form uses `--mode random --fraction F`; `--seed` may control reproducibility. It must record fraction, seed, candidate count, selected count, and parent landscape row count. Valid fractions are `0 < F <= 1`.
9. Top-density and threshold selections include all rows by default, including rows with `sld_display_is_outlier == True`. Public compact select does not expose density artifact or evaluation-space policy controls.
10. Source-metadata selection is an explicit `metadata` mode after a run. It must use the run-time ref/mov source metadata recorded by the run bundle and the landscape's `ref_source_row_id` or `mov_source_row_id`; it must not rematch particles or assume identity from a newly supplied STAR file.
11. Metadata selection must require `--metadata-domain ref|mov` and `--metadata-column COLUMN`. If both domains contain a requested column, cryoROLE must not guess which one is intended.
12. `--metadata-value VALUE[,VALUE...]` selects one or more values. For example `--metadata-value 1,3` includes rows whose metadata value equals `1` or `3`.
13. Metadata selection by value must record `metadata_domain`, `metadata_source_file`, `metadata_column`, `metadata_values`, `source_row_id_field`, candidate count, selected count, and missing/invalid metadata counts.
14. Metadata split selection may create one child selection per unique metadata value with `--split-by-value`. Each child selection must be a normal selection artifact with its own `selection.json`, `selected_particle_keys.csv`, optional selected-derived landscape, and export backtracking.

Recommended compact selection CLI:

```bash
cryorole select \
  --run-dir RUN \
  --selection-id ref_class3 \
  --mode metadata \
  --metadata-domain ref \
  --metadata-column rlnClassNumber \
  --metadata-value 3

cryorole select \
  --run-dir RUN \
  --selection-id mov_classes_1_3 \
  --mode metadata \
  --metadata-domain mov \
  --metadata-column rlnClassNumber \
  --metadata-value 1,3

cryorole select \
  --run-dir RUN \
  --mode metadata \
  --metadata-domain ref \
  --metadata-column rlnOpticsGroup \
  --split-by-value
```

`rlnOpticsGroup` selection should use the particle-loop `_rlnOpticsGroup` value. The optics table must remain source metadata preserved by export; optics-table row numbers are not particle identities.

External annotation STAR files are out of scope for the default metadata-selection mode. A future external-annotation mode must require an explicit join key, duplicate policy, overlap reporting, and provenance report; it must not silently row-align an arbitrary file.

### 10.1 Selected-derived landscape

Selection may optionally write a selected-derived landscape for direct visualization and downstream inspection:

```bash
cryorole select \
  --run-dir RUN \
  --space raw \
  --selection-id random_10pct \
  --mode random \
  --fraction 0.10 \
  --seed 123 \
  --write-selected-landscape \
  --recompute-sld
```

`--write-selected-landscape` writes an additional landscape containing only selected rows, so `cryorole visualize --selection-id ID --use-selected-landscape` can inspect the selection directly. By default this derived landscape keeps the parent SLD fields unchanged. `--recompute-sld` is an explicit opt-in that recomputes density inside the selected subset and preserves parent density fields under parent-prefixed columns.

The selected-derived landscape lives under:

```text
selections/<selection_id>/selected_landscape/
  landscape.npz
  landscape.csv
  landscape_report.json
```

Rules:

1. A selected-derived landscape must not overwrite raw or canonical landscape artifacts.
2. The selected-derived landscape row set must exactly match the selected particle keys and preserve `particle_key`, `ref_source_row_id`, and `mov_source_row_id`.
3. By default, selected-derived landscapes inherit the parent landscape's SLD fields and record `density_source = parent_landscape`.
4. `--recompute-sld` recomputes SLD only within the selected subset and records `density_source = recomputed_on_selection`.
5. When SLD is recomputed, parent density values must be preserved under explicit parent-prefixed fields such as `parent_sld_raw`, `parent_sld_display`, and `parent_sld_display_is_outlier`. Recomputed `sld_raw` then means density within the selected subset, not density in the parent landscape.
6. Recomputed SLD must preserve the SLD distance-floor stabilization and tail-jump display-outlier policy. It must not change RO, source conventions, source-row provenance, canonical transforms, or export behavior.
7. Small selections whose row count cannot support the requested kNN density must fail clearly or record an explicit effective-k policy; they must not silently emit misleading density fields.
8. `landscape_report.json` must record the parent landscape path, selected count, coordinate space, selected-landscape schema, density source, SLD recomputation parameters when applicable, and output paths.

---

## 11. Export contract

Export bridges explicit cryoROLE results back to external tools. Source inputs are immutable.

Primary public command:

```bash
cryorole export \
  --run-dir RUN \
  --selection-id SELECTION_ID \
  --domain both \
  --format auto
```

Public export is ID-driven, like `select` and `visualize`. `--selection-id` resolves
`RUN/selections/<selection_id>/selection.json` and writes to
`RUN/exports/<selection_id>/` by default.

`--selection PATH/to/selection.json` may remain as an advanced compatibility input
for direct artifact export. It should not be the documented common path. If
`--selection` and `--selection-id` are both provided, export should fail clearly
rather than silently choosing one.

`--output-dir` may remain as an advanced compatibility override. The default
public workflow should not require it. If omitted with `--run-dir`, output stays
inside the run bundle under `exports/<selection_id>/`.

`cryorole export selection` may remain as a compatibility alias, but
documentation should prefer direct `cryorole export`.

### 11.1 Selection metadata subset export

Default behavior is source-format-preserving subset export:

```text
ref input .star -> ref/selected_ref.star
ref input .cs   -> ref/selected_ref.cs
mov input .star -> mov/selected_mov.star
mov input .cs   -> mov/selected_mov.cs
```

Mixed ref/mov source formats are valid.

Default output:

```text
exports/<selection_id>/
  selected_particle_keys.txt
  selected_landscape_rows.csv
  export_report.json

  ref/
    selected_ref.star or selected_ref.cs
    export_ref_report.json

  mov/
    selected_mov.star or selected_mov.cs
    export_mov_report.json
```

Rules:

1. Export consumes an explicit selection.
2. Export must not reselect or rematch particles.
3. `--domain ref` uses `ref_source_row_id`.
4. `--domain mov` uses `mov_source_row_id`.
5. `--domain both` applies both independently.
6. `--format auto` resolves source format independently for each domain.
7. STAR export preserves source columns, optics data when present, and source pose/origin values.
8. CS export preserves structured-array dtype, fields, vector-valued fields, and source pose/shift values.
9. Export must not write canonical/display coordinates as physical poses.
10. Export must fail clearly if required `ref_source_row_id` / `mov_source_row_id` provenance is missing.
11. Export reports must record selection path, run directory, source files, selected counts, row-id provenance, domain choice, requested and resolved formats, output directory, output files, overwrite policy, and warnings.
12. User-facing errors should point to the compact path first: provide `--run-dir RUN --selection-id ID`, or use `--selection PATH/to/selection.json` for advanced direct artifact export.

### 11.2 Future transform export

Canonical transform export is a separate advanced mode. Only true global frame rotations may be exported as physical metadata/map transforms. Display-only recentering, RV-space recentering, axis limits, or visualization coordinates are non-exportable by default.

---

## 12. Policy system

Any behavior that changes interpretation must be represented as an explicit policy and recorded.

Core policies include:

```text
ConventionPolicy
IdentityPolicy
MatchPolicy
DensityPolicy
RepresentationPolicy
CanonicalizationPolicy
VisualizationPolicy
SelectionPolicy
SourceMetadataSelectionPolicy
SelectionExportPolicy
RawMetadataPolicy
```

Euler convention policy is part of derived representation policy. Public output uses `extrinsic_zyx`; any legacy/internal compatibility mapping must remain centralized so SciPy uppercase/lowercase details do not leak across plotting, selection, or workflow orchestration code.

Density and visualization policies must keep SLD floor stabilization separate from tail-jump display-outlier detection. Selection policy may reference display-outlier rows only as an explicit artifact-exclusion flag; it must continue to use `sld_raw` for density ranking or thresholding.

Source metadata selection policy is a selection policy, not a matching policy. It filters already matched landscape rows by looking up recorded ref/mov source-row metadata from the run-time inputs.

Raw metadata copying can dominate memory. Production mode should not copy full source rows into per-particle Python objects. Full raw metadata is a debug option, not a default.

---

## 13. Code organization

Responsibilities should stay separated:

```text
io/            readers and artifact stores/writers
normalize/     convention resolution and schema normalization
match/         identity resolution and matching
core/          RO, representations, distances, density
canonicalize/  frame alignment logic
select/        geodesic selection and backtracking
visualize/     display logic
export/        derived artifact writing and manifests
models/        typed data objects and policies
workflows/     orchestration
cli/           frontend only
gui/           optional frontend only
```

Do not put convention logic in CLI, plotting, GUI, or export code. Do not put frontend logic in core modules.

---

## 14. Testing contract

Behavior-changing changes must include tests for relevant invariants.

Required test categories:

```text
convention tests
matching tests
RO/math tests
SLD density tests
persistence/run-bundle tests
production-scale run tests/benchmarks
canonicalization tests
visualization tests
selection tests
export tests
CLI regression tests
```

Production-scale tests should protect:

- array-native run backend shape and row counts;
- vectorized RO equivalence to reference implementation;
- batched SLD equivalence to non-batched SLD on small datasets;
- chunked raw CSV row count and schema;
- `--no-visualize` behavior;
- memory profile report generation;
- benchmark scripts for 100k/1M synthetic datasets without putting very large runs in the fast suite.

---

## 15. Documentation ownership

Avoid duplicating long implementation details across multiple files.

- `AGENTS.md` owns short guardrails and active priority.
- `docs/architecture.md` owns stable contracts.
- `docs/production_run_plan.md` owns current production-scale run implementation details.
- `docs/roadmap.md` owns priorities and release sequencing.
- Future `docs/cli_reference.md` should own user-facing command help and examples.
