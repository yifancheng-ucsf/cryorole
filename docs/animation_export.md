# cryoROLE Offline Animation Export Plan

**Status:** Phase 1–4 and presentation refinements, including fixed inward dual-view crop, implemented; real ChimeraX/FFmpeg validation remains environment-dependent
**Command:** `cryorole animate`
**Priority:** downstream feature; it must not interrupt the active production-scale `run` work
**Audience:** cryoROLE developers, validation users, and contributors preparing publication movies

---

## 1. Purpose

`cryorole animate` is an offline display/export command. The implemented
Phase 1–4 path connects
an orientation trajectory in a cryoROLE landscape to a rigid-body rendering of
two density domains:

```text
waypoint CSV
  -> SO(3) trajectory
  -> three synchronized Euler projections
  -> fixed reference-domain density
  -> moving-domain rigid-body rotation about a declared pivot
  -> optional ChimeraX structure frames
  -> deterministic composite frames
  -> optional validated H.264 MP4
```

It is not a new analysis method. It consumes an existing run bundle and writes
a separate animation bundle. It must not modify the run bundle, landscapes,
SLD values, selections, source metadata, or reconstructed maps.

The scientific interpretation must be stated in every user-facing guide:

> The continuously transformed moving-domain density is a rigid-body rendering of the composite-map density, not an independently reconstructed density at every trajectory position.

---

## 2. Design assessment

The engineering brief is compatible with the cryoROLE architecture when the
following boundaries are explicit:

1. Interpolation belongs to SO(3), not Euler coordinate space.
2. Euler angles are derived values used for waypoint input, labels, and the
   three projection markers.
3. Canonical coordinates are not automatically raw physical RO values. A
   canonical waypoint must be mapped through the recorded canonical frame back
   to raw RO before a density transform is derived.
4. Landscape coordinate space and ChimeraX scene/map basis are separate
   concepts. Equality of the words `raw` or `canonical` is not sufficient
   evidence that a session and a trajectory share a physical frame.
5. The saved moving-model transform represents a declared baseline RO. The
   baseline must be recorded and used for every frame.
6. `script-only` can prepare trajectories, landscape frames, and ChimeraX
   scripts without ChimeraX. It cannot promise composite frames or MP4 until
   matching structure frames exist.
7. Animation filtering is display-only. It must never create a selection.

These constraints are implementation gates, not optional refinements.

---

## 3. Scope

### 3.1 Implemented Phase 1–4 slice

The first slice should support:

- an existing raw or canonical run landscape;
- EA or RV waypoint CSV input;
- quaternion SLERP along a waypoint path;
- a preconfigured ChimeraX `.cxs` session;
- explicit reference and moving model IDs;
- an explicit pivot in ChimeraX scene coordinates;
- three synchronized Euler projection panels;
- script generation without a ChimeraX dependency;
- optional local ChimeraX execution;
- deterministic aspect-preserving composition;
- optional FFmpeg encoding and FFprobe validation;
- an auditable manifest and logs.

### 3.2 Non-goals

The first slice does not include:

- a real-time ChimeraX or web GUI;
- click-to-drive landscape interaction;
- principal-curve or density-ridge extraction;
- automatic trajectory, axis, or pivot discovery;
- moving-domain translation or SE(3);
- density morphing;
- per-frame reconstruction or map resampling;
- automatic reconstruction submission;
- creation of scientific selections;
- inference of a unique 3D orientation from three 2D projections;
- direct raw-MRC scene construction;
- cryoDRGN or e2gmm integration.

---

## 4. Implemented CLI

Representative command:

```bash
cryorole animate \
  --run-dir RUN_DIR \
  --coordinate-set canonical \
  --canonical-id default \
  --path-csv trajectory_waypoints.csv \
  --path-space ea \
  --chimerax-session ino80_animation_scene.cxs \
  --reference-model-id "#1" \
  --moving-model-id "#2" \
  --pivot 120.4 98.6 75.2 \
  --baseline-ro first-waypoint \
  --map-frame canonical \
  --frames-per-segment 45 \
  --fps 30 \
  --output-dir animation_ino80
```

Required inputs:

```text
--run-dir PATH
--path-csv PATH
--path-space ea|rv
--chimerax-session PATH
--reference-model-id MODEL_ID
--moving-model-id MODEL_ID
--pivot X,Y,Z
--baseline-ro identity|first-waypoint|ea:A,B,G|rv:X,Y,Z
--map-frame raw|canonical|explicit
--output-dir PATH
```

Coordinate and baseline controls:

```text
--coordinate-set raw|canonical        default: raw
--canonical-id ID                     default: default
--euler-convention auto|extrinsic_zyx|intrinsic_zyx
--map-frame-transform PATH            required with --map-frame explicit
```

Trajectory controls:

```text
--frames-per-segment INT              default: 30; must be >= 2
--fps FLOAT                           default: 30
--hold-frames INT                     default: 0
--reverse
--ping-pong
```

Display controls:

```text
--threshold FLOAT                     same semantics as visualize
--sld-threshold FLOAT                 compatibility alias
--top-fraction FLOAT
--colormap NAME                       default: rainbow_r
--vmin FLOAT
--vmax FLOAT
--range AXIS:LOWER:UPPER
--point-size FLOAT                    default: 1
--landscape-alpha FLOAT               default: 1
--trail-frames INT                    default: 0
--show-current-values
--axis-limits PATH
```

Structure-rendering controls:

```text
--render-mode script-only|execute     default: script-only
--chimerax-bin PATH
--canvas-width INT                    default: 1800
--canvas-height INT                   default: 600
--structure-width INT
--structure-height INT
--overwrite
```

Composition and encoding controls:

```text
--composite-width INT                 default: 1920
--composite-height INT                default: 1080
--layout stacked|side_by_side         default: stacked
--background-color COLOR              default: #ffffff
--no-encode
--ffmpeg-bin PATH
--ffprobe-bin PATH
--crf INT                             default: 18
--movie-name NAME.mp4                 default: animation.mp4
```

Optional dual-view input:

```text
--secondary-chimerax-session PATH
```

Its presence enables two synchronized structure views. The current
`--chimerax-session` remains the primary session and single-view default.

Implemented dual-view crop control:

```text
--dual-structure-horizontal-crop FRACTION   default: 0
```

Rules:

1. `--threshold`/`--sld-threshold` and `--top-fraction` are mutually exclusive.
2. The first public slice uses the source landscape's `sld_display` field for
   background coloring and display filtering. It does not expose an arbitrary
   `--sld-column`.
3. Filtering and style resolution reuse the same public helpers as
   `cryorole visualize`: threshold means `sld_display >= threshold`,
   top-fraction selection is stable, tail-jump outliers do not set automatic
   `vmax`, points are drawn by ascending density, and the default style is
   `rainbow_r` with equal aspect and a bottom horizontal colorbar.
4. `--euler-convention auto` resolves from the selected landscape's recorded
   policy. Column names are not sufficient for inference.
5. An explicit Euler convention is a declaration to validate against recorded
   provenance, not permission to reinterpret an existing landscape.
6. `--chimerax-bin` is required only for `--render-mode execute`.
7. A closed-loop trajectory option is deferred until its meaning is specified.
   MP4 playback looping and adding a geodesic segment from the final waypoint
   back to the first waypoint are different policies and must not share an
   ambiguous `--loop` flag.
8. `--reverse` reverses waypoint order before interpolation.
9. `--ping-pong` appends the reverse traversal without duplicating the turning
   frame. The returned first waypoint is a real final frame, not an adjacent
   duplicate.
10. `--map-frame` describes the physical basis of the prepared session. It is
    not an alias for `--coordinate-set`, and frame mismatch is not bypassed by
    a warning-only flag.
11. Execute mode always writes validated composite frames. Without
    `--no-encode`, explicit FFmpeg and FFprobe paths are required.

### 4.1 Baseline precondition

`--baseline-ro` declares the RO represented by the moving model's saved initial
`scene_position`. It is not a display offset.

- `identity` is valid only when the saved session represents identity RO.
- `first-waypoint` is valid when the saved session represents the first
  waypoint orientation.
- explicit EA/RV baselines must use the same convention and coordinate
  resolution rules as waypoints.

The implementation cannot infer this relationship from a `.cxs` file. The
manifest must record it as a user assertion. The safe first public contract
therefore requires an explicit `--baseline-ro`; it does not silently default to
identity.

### 4.2 Canonical-axis ChimeraX views

`cryorole canonical-views` is a camera-only companion command:

```bash
cryorole canonical-views \
  --run-dir RUN \
  --canonical-id default \
  --chimerax-session scene.cxs \
  --map-frame raw \
  --render-mode execute \
  --chimerax-bin /path/to/chimerax \
  --output-dir canonical_views
```

It accepts `--map-frame raw|explicit`; explicit mode requires
`--map-frame-transform`. Script-only is the default. Width and height default
to 900 pixels.

---

## 5. Output bundle

The command writes a new, immutable-by-default animation bundle:

```text
<output_dir>/
  manifest.json
  trajectory.csv
  logs/
    cryorole_animate.log
    chimerax_render.log        # only when execute is attempted
    chimerax_render_status.json
    ffmpeg_encode.log          # only when encoding is attempted
    ffprobe_validate.log       # only when validation is attempted
  chimerax/
    render_structure.py
    run_render.cxc
    frame_transforms.csv
  frames/
    landscape/
      frame_000000.png
      ...
    structure/
      frame_000000.png
      ...
    composite/
      frame_000000.png
      ...
  animation.mp4                # only after successful encoding and validation
```

Absent artifacts must match the manifest state:

| State | Required artifacts |
|---|---|
| `scripts_ready` | manifest, trajectory, landscape frames, ChimeraX scripts |
| `structure_rendered` | all above plus validated structure frames |
| `composite_rendered` | all above plus validated composite frames |
| `movie_encoded` | all above plus validated H.264 MP4 |
| `failed` | manifest and available logs with the failed stage recorded |

`script-only` ends at `scripts_ready`. Execute with `--no-encode` ends at
`composite_rendered`; default execute ends at `movie_encoded` after FFprobe
validation. No empty structure, composite, or movie artifact is fabricated.

Existing output directories require `--overwrite`. Overwrite applies only to
the exact animation output directory and must never remove or mutate the
source run bundle.

---

## 6. Source landscape and coordinate resolution

### 6.1 Artifact resolution

For `--coordinate-set raw`, read:

```text
RUN/data/raw_landscape.npz
RUN/data/raw_landscape.csv
```

For `--coordinate-set canonical`, read:

```text
RUN/canonical/<canonical_id>/canonical_landscape.npz
RUN/canonical/<canonical_id>/canonical_landscape.csv
RUN/canonical/<canonical_id>/canonical_frame.json
RUN/canonical/<canonical_id>/canonical_frame.npz
```

Use the NPZ/report artifacts for policy and machine-readable coordinate
resolution. Read only the necessary columns from CSV when it is used for the
static background:

```text
selected alpha/beta/gamma columns
sld_display
sld_display_is_outlier when present
```

Do not load full source particle metadata or debug landscape JSON.

### 6.2 Euler convention

EA waypoints and displayed EA coordinates use the selected landscape's recorded
Euler convention. Public cryoROLE output is normally:

```text
extrinsic fixed-axis ZYX in degrees
```

Legacy intrinsic ZYX input is accepted only when provenance explicitly records
that policy or the user declaration matches an auditable legacy artifact. A
missing or conflicting convention is an error.

### 6.3 Canonical waypoint back-mapping

The recorded canonical frame currently defines:

```text
canonical_rv = raw_rv @ canonical_transform
```

Before physical rendering, a canonical waypoint must be mapped back:

```text
raw_rv = canonical_rv @ inverse(canonical_transform)
raw_ro = Rotation.from_rotvec(raw_rv)
```

The implementation must validate shape, finiteness, invertibility, transform
direction, handedness, and the frame artifact's coordinate-space declaration.
It must not treat canonical Euler angles as raw Euler angles or apply the
canonical RV matrix directly to Euler triples.

### 6.4 Landscape space versus scene basis

These fields must remain separate:

```text
waypoint_coordinate_set       raw or canonical landscape coordinates
resolved_ro_space             raw cryoROLE RO used for physical delta
session_scene_basis           physical basis used by ChimeraX scene positions
pivot_coordinate_frame        ChimeraX scene coordinates
baseline_ro                   RO represented by the saved moving-model pose
```

If the ChimeraX maps were globally reoriented, the animation needs an explicit,
audited physical `raw_to_scene` rotation so the raw-frame active delta can be
conjugated into the session scene basis. A canonical RV frame must not be
assumed to be that physical map transform without validation.

`--map-frame raw` resolves `raw_to_scene` to identity. `--map-frame canonical`
may resolve it from the selected canonical frame only after the implementation
validates that the artifact is a proper physical change-of-basis, confirms its
direction, and finds the explicit `physical_change_of_basis: true` audit
declaration. Ordinary canonical RV frames do not carry this declaration and
therefore fail this mode. `--map-frame explicit` reads a dedicated audited
transform from `--map-frame-transform`.

For a column-vector transform `S` satisfying
`x_scene = S @ x_raw`, convert the active delta by conjugation:

```text
delta_active_scene = S @ delta_active_raw @ inverse(S)
```

The first execute-capable implementation must fail when it cannot resolve this
mapping. A warning-only frame mismatch escape hatch is not part of the safe
public default.

### 6.5 Canonical-axis camera views

For `canonical_rv = raw_rv @ C`, the canonical axes expressed in a scene with
column-vector raw-to-scene rotation `S` are the columns of:

```text
axes_scene = S @ C
```

The three deterministic positive-axis views are:

| File | Toward viewer | Screen right | Screen up |
|---|---|---|---|
| `canonical_x_plus.png` | `+Xc` | `+Yc` | `+Zc` |
| `canonical_y_plus.png` | `+Yc` | `+Xc` | `-Zc` |
| `canonical_z_plus.png` | `+Zc` | `+Xc` | `+Yc` |

The renderer changes only the camera, computes every view from the saved
camera baseline, verifies all model scene transforms are unchanged, restores
the camera, and never overwrites the source session.
The implementation follows ChimeraX `Place` camera coordinates: `+X` is
screen-right, `+Y` is screen-up, `+Z` points toward the viewer, and the camera
looks along `-Z`.

Output is:

```text
<output_dir>/
  canonical_views.json
  chimerax/
    render_canonical_views.py
    run_render.cxc
  logs/
    chimerax_render.log              # execute attempted
    chimerax_render_status.json
  views/                             # execute success only
    canonical_x_plus.png
    canonical_y_plus.png
    canonical_z_plus.png
```

The manifest records `C`, `S`, raw/scene axis vectors, screen bases, stage
status, frame dimensions, completion checks, paths, warnings, and errors.
These are camera views of the original composite density, not transformed or
resampled canonical maps.

---

## 7. Waypoint CSV contract

### 7.1 EA waypoints

Required columns:

```csv
label,alpha_deg,beta_deg,gamma_deg
P1,-50.0,0.0,0.0
P2,-30.0,1.2,-0.5
P3,-10.0,1.8,-0.8
```

### 7.2 RV waypoints

Required columns:

```csv
label,rv_x_rad,rv_y_rad,rv_z_rad
P1,0.00,0.00,-0.87
P2,0.01,0.02,-0.52
P3,0.02,0.03,-0.17
```

RV units are always radians.

### 7.3 Optional columns and validation

Optional columns:

```text
hold_frames
segment_frames
annotation
```

Rules:

1. At least two waypoints are required.
2. Labels are non-empty and unique.
3. Coordinates are finite.
4. `hold_frames` is a non-negative integer.
5. `segment_frames` is an integer of at least 2 and describes the segment from
   the current waypoint to the next. The last value is ignored and reported.
6. Missing `segment_frames` uses `--frames-per-segment`.
7. Missing per-waypoint `hold_frames` uses the CLI `--hold-frames` value. A hold
   count means extra frames after the arrival frame.
8. EA and RV schemas must not be mixed.
9. Unknown columns may be preserved as metadata but must not change trajectory
   interpretation.

---

## 8. SO(3) trajectory contract

Waypoint parsing must use the existing centralized cryoROLE representation
helpers:

```text
EA waypoint -> recorded Euler policy -> Rotation
RV waypoint -> radians -> Rotation
```

Adjacent quaternion signs are normalized before interpolation:

```python
if dot(q_i, q_next) < 0:
    q_next = -q_next
```

Each segment uses quaternion SLERP on the shortest SO(3) path. Euler-linear
fallback is forbidden.

For a segment with `m = segment_frames`, generate `m` samples including both
endpoints. At a segment join, keep one copy of the shared waypoint. Without
holds or ping-pong, for `S` segments:

```text
frame_count = sum(segment_frames_s) - (S - 1)
```

`hold_frames = h` adds `h` extra copies after the waypoint arrival frame.

The trajectory table contains:

```text
frame_index
time_sec
segment_index
from_label
to_label
segment_fraction
quat_w
quat_x
quat_y
quat_z
rv_x_rad
rv_y_rad
rv_z_rad
alpha_deg
beta_deg
gamma_deg
is_waypoint
waypoint_label
```

The public CSV order is explicitly `wxyz`. If the internal library uses SciPy's
`xyzw`, conversion must occur at one documented boundary. The manifest records:

```text
quaternion_storage_order = wxyz
quaternion_internal_order = xyzw
```

All quaternion values must be normalized. EA values in `trajectory.csv` are
derived from the same per-frame `Rotation` used for RV and physical rendering.

---

## 9. RO-to-density transform contract

cryoROLE's RO definition remains:

```text
RO = R_ref^-1 R_mov
```

The initial animation bridge defines the relative passive change from the
declared baseline to a target:

```text
delta_passive = RO_target @ inverse(RO_baseline)
delta_active = inverse(delta_passive)
```

This formula must live in one shared helper and nowhere else:

```python
ro_target_to_active_delta(target_ro, baseline_ro)
```

The helper's docstring and tests must state the passive/active convention,
matrix action, multiplication order, and expected direction. Before execute
mode is signed off, a non-symmetric synthetic object and known single-axis
rotations must validate the formula against both cryoROLE normalization and
ChimeraX rendering. Visual plausibility is not validation.

No frame may accumulate a transform from the prior frame.

---

## 10. Pivot and absolute scene transform

The pivot is specified in ChimeraX scene coordinates, normally in angstroms:

```text
--pivot X,Y,Z
```

For a point in the same scene basis:

```text
x_target = pivot + delta_active @ (x_baseline - pivot)
```

The conceptual homogeneous transform is:

```text
T_frame =
  Translate(pivot)
  @ Rotate(delta_active_in_scene_basis)
  @ Translate(-pivot)
  @ T_baseline
```

The actual `Place` composition order must be verified against the ChimeraX API
and a reproducible 90-degree synthetic test. The reference model's initial and
final scene transforms must be identical.

---

## 11. ChimeraX rendering

The first version uses a preconfigured session. Users prepare model placement,
contours, colors, transparency, camera, background, clipping, and lighting
before saving the `.cxs` file.

Generated files:

```text
chimerax/render_structure.py
chimerax/run_render.cxc
logs/chimerax_render_status.json
```

The generated renderer must:

1. open the specified session;
2. resolve the exact reference and moving model IDs;
3. reject identical or ambiguous model matches;
4. preserve both initial scene transforms;
5. read the resolved per-frame RO/transform data;
6. compute every moving-model pose from the saved baseline;
7. export fixed-size PNG frames;
8. verify the reference transform did not change;
9. emit a structured log;
10. close without modifying the source session.

`execute` uses an argument list with `subprocess.run` and captures
stdout/stderr. Linux headless mode invokes:

```text
chimerax --offscreen --exit --nocolor --script run_render.cxc
```

Other platforms use their supported noninteractive mode. The generated
renderer records `success` in `chimerax_render_status.json` only after the final
frame, reference-transform check, and moving-transform restoration; exceptions
write a `failed` status with traceback. Execute succeeds only when the process,
completion artifact, and frame index/count/size checks all succeed. A zero
ChimeraX return code does not override a logged Python exception or missing
completion artifact.

Normal CI must not require ChimeraX. Optional tests use a small asymmetric
synthetic density and the `chimerax` marker.

### 11.1 Dual structure views

**Implemented.**

Dual view uses one trajectory and identical absolute model transforms in two
sessions that differ only in camera or display styling. Both sessions must
resolve the same model IDs and have matching initial reference and moving
scene transforms; baseline RO, map frame, and pivot remain shared. A mismatch
fails before composition.

Each session writes and validates its own completion artifact and contiguous
PNG sequence. Composition begins only after both succeed. Dual view supports
`stacked` layout only and writes:

```text
chimerax/primary/
chimerax/secondary/
frames/structure/primary/
frames/structure/secondary/
logs/chimerax_render_primary.log
logs/chimerax_render_primary_status.json
logs/chimerax_render_secondary.log
logs/chimerax_render_secondary_status.json
```

`script-only` writes both script sets but does not claim transform parity.
Source sessions are never overwritten. Single view retains its original paths.

---

## 12. Landscape frame rendering

Panel order matches public `cryorole visualize`:

```text
alpha-beta
beta-gamma
alpha-gamma
```

Every marker is derived from the same per-frame `Rotation`:

```text
alpha-beta marker = (alpha, beta)
beta-gamma marker = (beta, gamma)
alpha-gamma marker = (alpha, gamma)
```

The static background must use the same resolved display rows and style as an
equivalent `cryorole visualize` command:

```text
fixed field: sld_display
threshold: >=
top fraction: stable selection and source-row order
default colormap: rainbow_r
default point size/alpha: 1/1
color scale: tail-jump-aware visualize policy
draw order: ascending density
axis viewport: visualize legacy Euler limits
colorbar: horizontal bottom
visible colorbar label: SLD
aspect: equal
```

The visible label `SLD` is presentation-only. The artifact field,
filtering input, color-scale provenance, and manifest field remain
`sld_display`.

For production-scale landscapes:

1. read the background data once;
2. apply display-only filtering once;
3. render each static projection background once;
4. add only the marker, optional fading trail, and optional coordinate line per
   frame;
5. never rebuild a million-point scatter or reread CSV for every frame.

Default Euler limits follow the visualize legacy viewport. `--range` applies
the same display-only row filter and viewport semantics as visualize. Explicit
limits or ranges must contain every trajectory point or fail clearly; markers
must never be silently clipped.

Display filters do not change the trajectory and do not create selection
artifacts.

### 12.1 Compact projection presentation

**Implemented.**

- keep the fixed order `alpha-beta`, `beta-gamma`, `alpha-gamma`;
- use three equal-width panels with no more than 2% of the landscape width
  between adjacent panels;
- use one centered horizontal colorbar shared by all three panels;
- render the dynamic text block once at the upper-left of the landscape strip
  as `α=... β=... γ=...`;
- do not render the waypoint label or free-form `annotation` in the default
  projection frame; preserve them in trajectory artifacts;
- keep equal aspect, complete labels, and fixed geometry across all frames.

### 12.2 Projection legibility refinement

**Implemented.**

- reduce adjacent panel gap from 2% to at most 1.25% of landscape width
  without removing labels or changing the visualize viewport;
- keep background `--point-size` independent from the larger trajectory
  marker;
- enlarge the marker to a resolution-scaled high-contrast center and outline,
  targeting a 12–14 px outer diameter at `1920 x 1080`;
- render the single coordinate line in a fixed strip header at 16–18 px at
  `1920 x 1080`, with a light backing for contrast;
- resolve marker, text, and panel geometry once and reuse it for every frame.

---

## 13. Phase 4 composition and encoding

**Implemented.** The compositor uses complete landscape and structure PNGs
without inspecting or cropping Matplotlib subplot boundaries. Each input is
placed with contain/letterbox geometry in a fixed recorded region of a
`1920 x 1080` RGB canvas. Padding, gap, background, destination rectangles,
and Lanczos resizing are recorded in the manifest.

### 13.1 Stacked layout

**Implemented.** The compositor retains `side_by_side` for compatibility and
uses a structure-first `stacked` layout by default:

- structure frame centered in the upper region, using about 64% of canvas
  height;
- compact three-projection strip in the lower region, using about 36%;
- small fixed outer padding and one fixed gap between the two regions;
- aspect-preserving placement with geometry resolved once and reused for every
  frame;
- no per-frame content-aware cropping or scaling.

Exact rectangles, padding, gap, proportions, and layout name are recorded in
the manifest.

### 13.2 Dual-view stacked layout

**Implemented.**

When a secondary session is supplied, the upper structure region is split
into equal primary and secondary rectangles with at most 1.5% canvas-width
gap. Both use fixed aspect-preserving contain geometry. The landscape remains
full-width in the lower region. No titles, content-aware cropping, per-frame
zoom, or dynamic geometry are added by default.

Before composition, both input sequences must be contiguous, non-empty,
valid PNGs with stable dimensions and counts equal to `trajectory.csv`.
Composite frames are validated again before encoding.

### 13.3 Implemented fixed inward dual-view crop

Dual view may opt into:

```text
--dual-structure-horizontal-crop FRACTION
```

`FRACTION` is finite, defaults to `0`, and must satisfy `0 <= F < 0.5`.
For `F > 0`, crop `F` from both horizontal sides of each validated structure
frame, keep the full height, then place primary right-aligned and secondary
left-aligned with the resolved dual-view gap. Crop rectangles and destinations
are resolved once and reused for all frames. Per-side crop pixels use
deterministic half-up rounding of `source_width * FRACTION`.

This is an explicit presentation crop, not automatic content detection.
Source structure frames remain unchanged; there is no per-frame crop, stretch,
camera change, or model-transform change. Single view and `F = 0` retain
existing behavior.

MP4 policy:

```text
codec: H.264
pixel format: yuv420p
frame rate: constant
CRF: 18
faststart: enabled
```

FFmpeg and FFprobe run with argument lists and `shell=False`. FFprobe validates
H.264, `yuv420p`, dimensions, fps, and frame count. When frame count is not
reported, `duration * fps` is accepted only within a recorded 0.5-frame
tolerance. `--no-encode` retains composite PNGs without creating an MP4.
Landscape, structure, and composite frames are never automatically removed.

---

## 14. Manifest contract

`manifest.json` records at least:

The implemented presentation-refinement contract uses animation manifest
schema version `4`.

```text
schema_version
status
command
cryorole_version
git_commit when available
creation_timestamp
run_dir
source_landscape
canonical_frame when used
waypoint_coordinate_set
resolved_ro_space
euler_convention
path_csv
path_space
quaternion orders
baseline_ro declaration
session_scene_basis
raw_to_scene transform provenance when used
pivot and pivot coordinate frame
ChimeraX session path
reference and moving model IDs
trajectory policies
fps and frame count
display filter and resolved SLD field
resolved visualize display policy and style
projection rectangles, gap, annotation policy, and shared colorbar policy
resolved trajectory marker and current-value text sizes
total and displayed landscape rows
axis-limit policy and expansions
canvas and structure frame sizes
trajectory/landscape/script/structure stage statuses
ChimeraX return code
compositor layout/version, regions, fractions, destinations, and frame validation
structure view count, session paths, per-view destinations and validation
dual-view crop fraction, source crop rectangles, inward alignment, and gap
FFmpeg arguments and return code
FFprobe payload, validation method, and return code
output paths
warnings and errors
resolved CLI arguments
```

Do not hash large MRC or `.cxs` inputs by default. A future explicit input-hash
policy may add this behavior.

---

## 15. Failure policy

Fail clearly for:

- fewer than two waypoints;
- missing or mixed EA/RV columns;
- non-finite values or duplicate labels;
- invalid segment/hold frame counts;
- missing raw/canonical artifacts;
- missing, ambiguous, or conflicting Euler provenance;
- invalid canonical frame or failed canonical-to-raw mapping;
- unresolved session scene basis or raw-to-scene transform;
- an invalid baseline declaration;
- missing session or model IDs;
- identical reference and moving model IDs;
- mismatched initial model transforms across requested dual-view sessions;
- invalid pivot;
- missing executables required by the selected mode;
- non-finite rotation or scene transforms;
- an empty displayed landscape;
- clipped trajectory under explicit axis limits;
- frame index, count, or dimension mismatches;
- a non-zero ChimeraX return code;
- a logged ChimeraX/Python rendering exception, missing completion artifact, or
  completion artifact that does not report success;
- an existing output directory without `--overwrite`;

No failed stage may silently continue into a downstream stage. Failure leaves
the manifest and logs needed for diagnosis.

---

## 16. Implementation phases

### Phase 1: trajectory and transform core

**Implemented.**

- waypoint schemas and validation;
- EA/RV parsing through centralized helpers;
- quaternion sign continuity and SLERP;
- frame-count and hold behavior;
- `trajectory.csv`;
- canonical-to-raw waypoint mapping;
- RO-to-active-delta helper;
- pivot/absolute transform math;
- unit and synthetic transform tests;
- no ChimeraX dependency.

### Phase 2: landscape animation

**Implemented.**

- artifact and policy resolution;
- display-only SLD filtering;
- static backgrounds;
- dynamic markers, trail, annotations, and axis limits;
- landscape frames;
- production-scale background benchmark.

### Phase 3: ChimeraX script export

**Implemented.** Normal CI validates generated scripts and subprocess/frame
handling without ChimeraX. Real ChimeraX scene composition still requires the
optional local integration validation.

- session and model validation;
- baseline scene transform;
- raw-to-scene basis conversion;
- absolute per-frame transform;
- script-only output;
- optional execute mode;
- local marked integration tests.

### Phase 3 stabilization: display parity and Linux headless execution

**Implemented.**

- shared `visualize` display-policy helpers;
- matching threshold/top-fraction rows, legacy style, color scale, point order,
  projections, viewport, and colorbar;
- Linux `--offscreen --script` execution;
- renderer-completion artifact independent of process return code;
- parity and headless-failure regression tests.

### Phase 4: compositor and MP4

**Implemented.**

- structure/landscape frame validation;
- deterministic layout and annotations;
- Pillow composition;
- FFmpeg encoding;
- stage status, logs, and recovery behavior.

### Phase 5: examples and public documentation

- `examples/animation/ino80/README.md`;
- EA and RV waypoint examples;
- session setup guide without large MRC/`.cxs` files;
- script-only, execute, reverse, and ping-pong examples;
- CLI help and release notes;
- final validation report.

### Presentation follow-up

**Planned.** Tighten projection spacing and dynamic overlays, then add an
optional two-session structure view without changing trajectory, transform,
display-selection, or encoding semantics.

---

## 17. Test and acceptance plan

Normal tests cover:

- two-waypoint and multisegment SLERP;
- exact endpoints and nonduplicated joins;
- holds, reverse, and ping-pong;
- quaternion sign continuity and normalization;
- EA/RV round trips;
- deterministic frame counts;
- adjacent-frame SO(3) distance;
- canonical-to-raw back-mapping;
- identity and non-identity baselines;
- passive/active direction;
- pivot invariance;
- unchanged reference transform;
- projection marker coordinates;
- exact display-row/style parity with equivalent `visualize` inputs;
- Linux execute arguments include `--offscreen` and `--script`;
- return code zero with a missing/failed completion artifact still fails;
- script contents;
- manifest state transitions;
- frame counts and dimensions;
- failure messages.

The presentation refinement includes focused tests for compact panel
geometry, one coordinate-only annotation block, a shared colorbar, fixed
stacked destination rectangles, and `side_by_side` compatibility.

Focused tests cover the tighter gap, resolution-scaled marker and coordinate
text, unchanged display-row selection, optional secondary-session parsing,
matching initial model transforms, independent renderer completion,
synchronized frame counts, fixed dual destinations, and rejection of dual view
with `side_by_side`.

Focused tests cover crop validation, fixed crop/destination geometry, inward
alignment, unchanged source frames, single-view compatibility, and visible
`SLD` colorbar text while retaining the underlying `sld_display` field.

Optional `pytest -m chimerax` validation covers a small asymmetric density,
known 90-degree rotations, exact model resolution, scene composition order, and
rendered frame counts.

Do not use full PNG hashes as the primary assertion across platforms. Validate
dimensions, numerical transforms, frame indices, marker neighborhoods, and
manifest values.

The feature is accepted only when:

1. EA and RV paths produce one shared SO(3) trajectory truth.
2. Raw and canonical paths resolve to the correct raw physical RO.
3. Reference density stays fixed.
4. Moving density rotates in the validated direction about the declared pivot.
5. Every frame uses an absolute transform from the saved baseline.
6. Script-only requires no ChimeraX installation.
7. Execute mode stops on ChimeraX failure.
8. Composite and trajectory frame counts match exactly.
9. H.264 output is validated when encoding is requested.
10. Existing run, landscape, selection, and metadata artifacts remain unchanged.
11. The manifest makes every convention and frame assertion auditable.
12. User documentation distinguishes rigid-body rendering from reconstruction.

---

## 18. Implementation-blocking design gates

Before code work begins, resolve and test:

1. the existing centralized EA/RV/Rotation helpers to reuse;
2. whether every supported canonical frame is a proper, invertible
   change-of-basis for canonical-to-raw RO mapping;
3. the exact meaning and source of a physical raw-to-ChimeraX-scene transform;
4. the passive/active delta direction against current normalized RO semantics;
5. ChimeraX `Place` multiplication order;
6. the baseline RO represented by the prepared example session;
7. the resume/finalization workflow after script-only rendering;
8. closed-loop semantics, which remain outside the minimum first slice.

No implementation should bypass these gates with convention guessing or
visual-only validation.
