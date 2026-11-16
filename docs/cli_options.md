# All command options

This reference contains the public `--help` output from cryoROLE **2.0.0a1**, generated on 2026-11-16. Hidden compatibility aliases are intentionally omitted. Your installed `cryorole COMMAND --help` remains the source of truth if its version differs.

[CLI policies and behavior](cli_reference.md) · [Worked examples](command_workflow.md) · [README demo](../README.md#3-cli-walkthrough)

[align](#align) · [preflight](#preflight) · [run](#run) · [status](#status) · [next](#next) · [guide](#guide) · [explore](#explore) · [gui](#gui) · [visualize](#visualize) · [animate](#animate) · [canonical-views](#canonical-views) · [canonicalize](#canonicalize) · [select](#select) · [export](#export) · [manifest](#manifest)

## Command overview

```text
usage: cryorole [-h] [--version]
                {align,preflight,run,status,next,guide,explore,gui,visualize,animate,canonical-views,canonicalize,select,export,manifest}
                ...

cryoROLE: relative-orientation landscapes between two refined bodies.

positional arguments:
  {align,preflight,run,status,next,guide,explore,gui,visualize,animate,canonical-views,canonicalize,select,export,manifest}
    align               Align STAR metadata before cryorole run --row-aligned.
    preflight           Check inputs, matching, policies, and resources
                        without creating a run.
    run                 Generate RO/SLD landscape artifacts without selection
                        by default.
    status              Inspect actual run-bundle artifacts and integrity.
    next                Recommend conservative next commands from run
                        artifacts.
    guide               Plan or resume the cryoROLE workflow safely.
    explore             Explore draft radius selections in an offline local
                        browser.
    gui                 Open the desktop interface (install with: pip install
                        "cryorole[gui]").
    visualize           Render display-only views from an existing run bundle.
    animate             Export Phase 1-4 trajectory, rendered frames,
                        composites, and optional MP4.
    canonical-views     Export camera-only ChimeraX views along the canonical
                        axes.
    canonicalize        Canonicalize an existing landscape explicitly.
    select              Create a named scientific Selection from a saved
                        landscape.
    export              Export a selection from a run bundle; use --run-dir
                        RUN --selection-id ID.
    manifest            Write a provenance manifest from available workflow
                        sections.

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

Typical workflow:
  cryorole preflight --ref ref.cs --mov mov.cs
  cryorole run --ref ref.cs --mov mov.cs --output-dir my_run
  cryorole canonicalize --run-dir my_run
  cryorole visualize --run-dir my_run --space canonical
  cryorole select --run-dir my_run --selection-id region_01 --mode radius -c 0 0 0 -r 10
  cryorole export --run-dir my_run --selection-id region_01
Run `cryorole COMMAND --help` for details and examples.
```

## align

```text
usage: cryorole align [-h] [--ref REF] [--mov MOV] [--align-id ALIGN_ID]
                      [--output-dir OUTPUT_DIR] [--key KEY [KEY ...]]
                      [--key-pair REF_COL=MOV_COL] [--float-tol COL=TOL]
                      [--path-mode PATH_MODE]
                      [--duplicate-policy {exclude,first}] [--overwrite]
                      [--coordinate-match {unchanged,recentered-exact}]
                      [--recenter-shift X Y Z] [--ref-angpix REF_ANGPIX]
                      [--micrograph-angpix MICROGRAPH_ANGPIX]
                      [--via-extraction INPUT OUTPUT]
                      [--fix-subtract-coordinates SUBTRACT_JOB_OR_STAR]
                      [--subtract-input SUBTRACT_INPUT] [--center X Y Z]
                      [--model-angpix MODEL_ANGPIX] [--apply-to APPLY_TO]

Align STAR metadata before cryorole run --row-aligned. Every strategy is exact and auditable: outputs are verbatim row subsets plus match_table.csv and align_report.json, and the next `cryorole run --row-aligned` command is printed. Also corrects the stale coordinates of a recentred RELION subtraction (--fix-subtract-coordinates).

options:
  -h, --help            show this help message and exit
  --ref REF             Reference STAR metadata path.
  --mov MOV             Moving STAR metadata path.
  --align-id ALIGN_ID   Output id. Default: default.
  --output-dir OUTPUT_DIR
                        Output directory. Default: <directory of
                        --ref>/cryorole_alignments/<align-id>/.
  --key KEY [KEY ...]   Explicit STAR key columns for alignment (same names in
                        both files).
  --key-pair REF_COL=MOV_COL
                        Cross-column key, e.g.
                        _rlnImageName=_rlnImageOriginalName. May repeat
                        (composite key).
  --float-tol COL=TOL   Numeric key tolerance, e.g. _rlnCoordinateX=0.1. May
                        repeat.
  --path-mode PATH_MODE
                        Path normalization for path-like keys: exact,
                        basename, or suffix:N. Default: exact.
  --duplicate-policy {exclude,first}
                        Duplicate key handling. Default: exclude.
  --overwrite           Replace artifacts for this align id only; source STAR
                        files are unchanged.

RELION recentring (exact):
  --coordinate-match {unchanged,recentered-exact}
                        unchanged: micrograph + identical coordinates.
                        recentered-exact: --ref is a RELION re-extraction
                        input and --mov its output; rows are paired only if
                        predicted integer coordinates, residual origins,
                        angles and one-to-one assignment all verify.
  --recenter-shift X Y Z
                        RELION --recenter_x/y/z of the extraction (reference
                        pixels).
  --ref-angpix REF_ANGPIX
                        RELION --ref_angpix of the extraction, if it was
                        given.
  --micrograph-angpix MICROGRAPH_ANGPIX
                        Micrograph pixel size (Å). Default:
                        _rlnMicrographPixelSize, else
                        _rlnMicrographOriginalPixelSize.
  --via-extraction INPUT OUTPUT
                        Chain: --ref ↔ INPUT (particle names) → INPUT ↔ OUTPUT
                        (recentered-exact) → OUTPUT ↔ --mov (particle names).
                        Only rows verified on every link are paired.

RELION subtraction coordinate correction:
  --fix-subtract-coordinates SUBTRACT_JOB_OR_STAR
                        Write a copy of a recentred subtraction STAR with
                        corrected _rlnCoordinateX/Y (RELION leaves them
                        stale). Reads note.txt of the job; every value is
                        verified before writing.
  --subtract-input SUBTRACT_INPUT
                        Subtraction input STAR (the refinement's
                        run_data.star).
  --center X Y Z        Subtraction --center_x/y/z (reference-model pixels).
  --model-angpix MODEL_ANGPIX
                        Reference-model pixel size used by the subtraction.
  --apply-to APPLY_TO   Transfer the corrected coordinates to this STAR
                        (joined on _rlnImageName).

Examples:
  cryorole align --ref consensus.star --mov body.star --key-pair _rlnImageName=_rlnImageOriginalName
  cryorole align --ref Extract_in.star --mov Extract/job055/particles.star \
      --coordinate-match recentered-exact --recenter-shift -9 22 -102
  cryorole align --ref ref.star --mov mov.star --via-extraction Extract_in.star Extract/job055/particles.star \
      --recenter-shift -9 22 -102
  cryorole align --fix-subtract-coordinates Subtract/job041/

`cryorole preflight --ref REF --mov MOV` lists candidate keys when default matching fails.
```

## preflight

```text
usage: cryorole preflight [-h] --ref REF --mov MOV [--ref-domain REF_DOMAIN]
                          [--mov-domain MOV_DOMAIN] [--output-dir OUTPUT_DIR]
                          [--row-aligned] [--allow-low-overlap]
                          [--sld-metric {rotvec_euclidean,so3_geodesic}]
                          [--no-visualize] [--json [PATH]]

options:
  -h, --help            show this help message and exit
  --ref REF             Reference-domain pose metadata.
  --mov MOV             Moving-domain pose metadata.
  --ref-domain REF_DOMAIN
                        Reference domain label.
  --mov-domain MOV_DOMAIN
                        Moving domain label.
  --output-dir OUTPUT_DIR
                        Prospective run directory.
  --row-aligned         Assert row N matches row N.
  --allow-low-overlap   Explicitly allow overlap below the public 50%
                        threshold.
  --sld-metric {rotvec_euclidean,so3_geodesic}
                        Prospective SLD metric. Default: so3_geodesic.
  --no-visualize        Estimate a run without previews.
  --json [PATH]         Write JSON to PATH, or stdout when PATH is omitted.

Examples:
  cryorole preflight --ref ref.cs --mov mov.cs
  cryorole preflight --ref consensus.star --mov body.star --json report.json
```

## run

```text
usage: cryorole run [-h] --ref REF --mov MOV [--ref-domain REF_DOMAIN]
                    [--mov-domain MOV_DOMAIN] [--row-aligned]
                    [--allow-low-overlap] [--dry-run] [--json [PATH]]
                    [--sld-metric {rotvec_euclidean,so3_geodesic}]
                    [--output-dir OUTPUT_DIR] [--overwrite] [--no-visualize]
                    [--jobs JOBS] [--point-style {scatter,density}]

options:
  -h, --help            show this help message and exit
  --ref REF             Reference-domain pose metadata.
  --mov MOV             Moving-domain pose metadata.
  --ref-domain REF_DOMAIN
                        Reference domain name.
  --mov-domain MOV_DOMAIN
                        Moving domain name.
  --row-aligned         Assert that ref and mov row N refer to the same
                        particle.
  --allow-low-overlap   Explicitly allow key-match coverage below the public
                        50% safety threshold.
  --dry-run             Run the shared preflight only; do not create a run
                        bundle.
  --json [PATH]         With --dry-run, write the structured preflight report
                        (default: stdout).
  --sld-metric {rotvec_euclidean,so3_geodesic}
                        SLD kNN distance metric. Default: so3_geodesic;
                        rotvec_euclidean is a faster approximation for very
                        large datasets.
  --output-dir OUTPUT_DIR
                        Run bundle output directory. Default:
                        cryorole_outputs.
  --overwrite           Replace an existing --output-dir bundle (published
                        only after the new run completes).
  --no-visualize        Skip raw quick-look visualization.
  --jobs JOBS           Processes for the quick-look figures: 1-3 or auto.
                        Default: auto (up to 3, limited by CPU cores and free
                        memory; 1 below 100k particles).
  --point-style {scatter,density}
                        How 2D projection figures draw points: scatter
                        (default; one marker per point) or density (one image
                        cell per screen pixel holding the visible SLD value;
                        about 15x faster at 500k points, near-identical look,
                        slightly thinner in sparse regions).

Examples:
  cryorole run --ref ref.cs --mov mov.cs --output-dir my_run
  cryorole run --ref aligned_ref.star --mov aligned_mov.star --row-aligned --output-dir my_run
  cryorole run --ref ref.cs --mov mov.cs --dry-run   # same checks as preflight
```

## status

```text
usage: cryorole status [-h] [--run-dir RUN_DIR] [--json [PATH]]

Inspect a run bundle. Exits 0 for a completed or legacy bundle, 3 for a missing, failed or incomplete one.

options:
  -h, --help         show this help message and exit
  --run-dir RUN_DIR  Run bundle. Omitted: the current directory if it is a
                     bundle, else ./cryorole_outputs.
  --json [PATH]

Example:
  cryorole status --run-dir my_run
```

## next

```text
usage: cryorole next [-h] [--run-dir RUN_DIR] [--json [PATH]]

Recommend next commands. Exits 0 for a completed or legacy bundle, 3 for a missing, failed or incomplete one.

options:
  -h, --help         show this help message and exit
  --run-dir RUN_DIR  Run bundle. Omitted: the current directory if it is a
                     bundle, else ./cryorole_outputs.
  --json [PATH]

Example:
  cryorole next --run-dir my_run
```

## guide

```text
usage: cryorole guide [-h] [--ref REF] [--mov MOV] [--run-dir RUN_DIR]
                      [--output-dir OUTPUT_DIR] [--row-aligned]
                      [--allow-low-overlap] [--non-interactive]
                      [--execute-run] [--json [PATH]]

options:
  -h, --help            show this help message and exit
  --ref REF             Reference metadata for a new workflow.
  --mov MOV             Moving metadata for a new workflow.
  --run-dir RUN_DIR     Existing run bundle to resume.
  --output-dir OUTPUT_DIR
                        Prospective output for a new run.
  --row-aligned         Explicit row-aligned assertion.
  --allow-low-overlap   Explicit low-overlap override.
  --non-interactive     Print a plan and never wait for input.
  --execute-run         Explicitly execute a ready run; never creates
                        selections.
  --json [PATH]

Examples:
  cryorole guide --ref ref.cs --mov mov.cs --non-interactive
  cryorole guide --run-dir my_run
```

## explore

```text
usage: cryorole explore [-h] [--run-dir RUN_DIR] [--space {raw,canonical}]
                        [--canonical-id CANONICAL_ID]
                        [--selection-id SELECTION_ID]
                        [--threshold THRESHOLD | --top-fraction TOP_FRACTION]
                        [--max-display-points MAX_DISPLAY_POINTS]
                        [--colormap COLORMAP] [--port PORT] [--no-open]

options:
  -h, --help            show this help message and exit
  --run-dir RUN_DIR     Run bundle. Omitted: the current directory if it is a
                        bundle, else ./cryorole_outputs.
  --space {raw,canonical}
  --canonical-id CANONICAL_ID
                        Canonical frame (with --space canonical). Omitted: the
                        only frame in the bundle.
  --selection-id SELECTION_ID
                        Optional existing selection overlay.
  --threshold THRESHOLD
                        Initial display-only SLD threshold.
  --top-fraction TOP_FRACTION
                        Initial display-only top SLD fraction.
  --max-display-points MAX_DISPLAY_POINTS
  --colormap COLORMAP   Accessible default display colormap; legacy rainbow_r
                        remains available.
  --port PORT           Loopback port; 0 selects an available port.
  --no-open             Do not open the default browser automatically.

Examples:
  cryorole explore --run-dir my_run --space canonical
```

## gui

```text
usage: cryorole gui [-h] [--no-gl] [--run-dir RUN_DIR]

options:
  -h, --help         show this help message and exit
  --no-gl            Do not use OpenGL: no 3D view, and no OpenGL check at
                     start-up.
  --run-dir RUN_DIR  Run bundle to open in the View page. Omitted: open one
                     from File > Open run.

Examples:
  cryorole gui
  cryorole gui --run-dir my_run
  cryorole gui --no-gl   # 2D and 1D views only; skips the OpenGL check
  __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia cryorole gui   # monitor on an onboard chip
```

## visualize

```text
usage: cryorole visualize [-h] [--run-dir RUN_DIR] [--space {raw,canonical}]
                          [--canonical-id CANONICAL_ID]
                          [--selection-id SELECTION_ID]
                          [--use-selected-landscape] [--visual-id VISUAL_ID]
                          [--view VIEW] [--colormap COLORMAP]
                          [--representation {euler,rotvec,both}]
                          [--top-fraction TOP_FRACTION | --sld-threshold THRESHOLD | --all]
                          [--range AXIS:LOWER:UPPER] [--format FORMATS]
                          [--vmin VMIN] [--vmax VMAX]
                          [--point-size POINT_SIZE] [--opacity ALPHA]
                          [--axis-limit AXIS:LOWER:UPPER]
                          [--max-points MAX_POINTS] [--bins BINS]
                          [--hist-mode {count,percent}] [--kde]
                          [--kde-bandwidth KDE_BANDWIDTH]
                          [--3d-mode {interactive,static}] [--overwrite]

options:
  -h, --help            show this help message and exit
  --run-dir RUN_DIR     Run bundle. Omitted: the current directory if it is a
                        bundle, else ./cryorole_outputs.
  --space {raw,canonical}
                        Landscape space to visualize. Default: raw.
  --canonical-id CANONICAL_ID
                        Canonical frame (with --space canonical). Omitted: the
                        only frame in the bundle.
  --selection-id SELECTION_ID
                        Visualize only particles from
                        run_dir/selections/SELECTION_ID.
  --use-selected-landscape
                        Read selections/SELECTION_ID/selected_landscape/landsc
                        ape.npz directly.
  --visual-id VISUAL_ID
                        Visualization output id. Default: default.
  --view VIEW           Comma-separated views to generate: 2d,1d,3d. Default:
                        2d.
  --colormap COLORMAP   Matplotlib colormap for display colors. Default:
                        rainbow_r.
  --representation {euler,rotvec,both}
  --top-fraction TOP_FRACTION
                        Display the top fraction by SLD.
  --sld-threshold THRESHOLD
                        Display rows with sld_display >= T. Default: 1.
  --all                 Display all candidate particles instead of the default
                        SLD threshold.
  --range AXIS:LOWER:UPPER
                        Display-only range filter, e.g. --range alpha:-30:30
                        or --range x:-0.4:0.4.
  --format FORMATS      Static figure format; repeat or use commas. Default:
                        png.
  --vmin VMIN           Display-only color minimum.
  --vmax VMAX           Display-only color maximum.
  --point-size POINT_SIZE
                        Display-only scatter point size.
  --opacity ALPHA       Display-only point opacity: 0 transparent, 1 opaque.
                        Uses the style default when omitted.
  --axis-limit AXIS:LOWER:UPPER
                        Viewport limit only; may repeat and does not filter
                        rows.
  --max-points MAX_POINTS
                        Override the deterministic point cap for 2D/3D views.
  --bins BINS           1D histogram bins: auto or a positive integer.
                        Default: auto.
  --hist-mode {count,percent}
                        1D histogram y-axis mode. Default: percent.
  --kde                 Add coordinate KDE curves to 1D views.
  --kde-bandwidth KDE_BANDWIDTH
                        KDE bandwidth when --kde is used: scott, silverman, or
                        a positive float.
  --3d-mode {interactive,static}
                        3D output mode. Default: interactive offline HTML.
  --overwrite

Examples:
  cryorole visualize --run-dir my_run
  cryorole visualize --run-dir my_run --space canonical --view 2d,1d --top-fraction 0.4
  cryorole visualize --run-dir my_run --selection-id region_01 --visual-id region_01
```

## animate

```text
usage: cryorole animate [-h] --run-dir RUN_DIR
                        [--coordinate-set {raw,canonical}]
                        [--canonical-id CANONICAL_ID] --path-csv PATH_CSV
                        --path-space {ea,rv} --chimerax-session
                        CHIMERAX_SESSION
                        [--secondary-chimerax-session SECONDARY_CHIMERAX_SESSION]
                        [--tertiary-chimerax-session TERTIARY_CHIMERAX_SESSION]
                        --reference-model-id REFERENCE_MODEL_ID
                        --moving-model-id MOVING_MODEL_ID --pivot X Y Z
                        --baseline-ro BASELINE_RO --map-frame
                        {raw,canonical,explicit}
                        [--ro-direction {active,legacy-inverse}]
                        [--map-frame-transform MAP_FRAME_TRANSFORM]
                        [--euler-convention {auto,extrinsic_zyx,intrinsic_zyx}]
                        [--frames-per-segment FRAMES_PER_SEGMENT] [--fps FPS]
                        [--hold-frames HOLD_FRAMES] [--reverse] [--ping-pong]
                        [--threshold THRESHOLD | --top-fraction TOP_FRACTION]
                        [--colormap COLORMAP] [--vmin VMIN] [--vmax VMAX]
                        [--range AXIS:LOWER:UPPER] [--point-size POINT_SIZE]
                        [--landscape-alpha LANDSCAPE_ALPHA]
                        [--trail-frames TRAIL_FRAMES] [--show-current-values]
                        [--axis-limits AXIS_LIMITS]
                        [--render-mode {script-only,execute}]
                        [--chimerax-bin CHIMERAX_BIN]
                        [--canvas-width CANVAS_WIDTH]
                        [--canvas-height CANVAS_HEIGHT]
                        [--structure-width STRUCTURE_WIDTH]
                        [--structure-height STRUCTURE_HEIGHT]
                        [--structure-horizontal-crop FRACTION | --dual-structure-horizontal-crop FRACTION]
                        [--structure-vertical-crop FRACTION]
                        [--composite-width COMPOSITE_WIDTH]
                        [--composite-height COMPOSITE_HEIGHT]
                        [--layout {stacked,side_by_side}]
                        [--background-color BACKGROUND_COLOR] [--no-encode]
                        [--ffmpeg-bin FFMPEG_BIN] [--ffprobe-bin FFPROBE_BIN]
                        [--crf CRF] [--movie-name MOVIE_NAME] --output-dir
                        OUTPUT_DIR [--overwrite]

options:
  -h, --help            show this help message and exit
  --run-dir RUN_DIR     Existing cryoROLE run bundle.
  --coordinate-set {raw,canonical}
                        Waypoint and landscape coordinate set. Default: raw.
  --canonical-id CANONICAL_ID
  --path-csv PATH_CSV   EA or RV waypoint CSV.
  --path-space {ea,rv}
  --chimerax-session CHIMERAX_SESSION
                        Preconfigured .cxs session.
  --secondary-chimerax-session SECONDARY_CHIMERAX_SESSION
                        Optional second .cxs with the same model transforms
                        and a different view.
  --tertiary-chimerax-session TERTIARY_CHIMERAX_SESSION
                        Optional third .cxs with the same model transforms and
                        a different view; requires --secondary-chimerax-
                        session.
  --reference-model-id REFERENCE_MODEL_ID
                        Exact stationary ChimeraX model ID; repeat for a rigid
                        reference group.
  --moving-model-id MOVING_MODEL_ID
                        Exact ChimeraX model ID; repeat for a rigid moving
                        group.
  --pivot X Y Z         Rotation pivot in ChimeraX scene coordinates.
  --baseline-ro BASELINE_RO
                        Saved moving-model RO assertion: identity, first-
                        waypoint, ea:A,B,G, or rv:X,Y,Z.
  --map-frame {raw,canonical,explicit}
  --ro-direction {active,legacy-inverse}
                        How the moving model follows RO. active (default,
                        validated): turn by RO_target @ inverse(RO_baseline).
                        legacy-inverse: the inverse, as in movies made before
                        2026-10-12; only for reproducing them.
  --map-frame-transform MAP_FRAME_TRANSFORM
  --euler-convention {auto,extrinsic_zyx,intrinsic_zyx}
  --frames-per-segment FRAMES_PER_SEGMENT
  --fps FPS
  --hold-frames HOLD_FRAMES
  --reverse
  --ping-pong
  --threshold THRESHOLD
                        Display rows with sld_display >= threshold.
  --top-fraction TOP_FRACTION
  --colormap COLORMAP
  --vmin VMIN
  --vmax VMAX
  --range AXIS:LOWER:UPPER
  --point-size POINT_SIZE
  --landscape-alpha LANDSCAPE_ALPHA
  --trail-frames TRAIL_FRAMES
  --show-current-values
  --axis-limits AXIS_LIMITS
                        JSON file defining all three panel limits.
  --render-mode {script-only,execute}
  --chimerax-bin CHIMERAX_BIN
  --canvas-width CANVAS_WIDTH
  --canvas-height CANVAS_HEIGHT
  --structure-width STRUCTURE_WIDTH
  --structure-height STRUCTURE_HEIGHT
  --structure-horizontal-crop FRACTION
                        Crop this fraction from both horizontal sides of each
                        multi-view frame.
  --dual-structure-horizontal-crop FRACTION
                        Compatibility alias for --structure-horizontal-crop.
  --structure-vertical-crop FRACTION
                        Crop this fraction from both vertical sides of each
                        multi-view frame.
  --composite-width COMPOSITE_WIDTH
  --composite-height COMPOSITE_HEIGHT
  --layout {stacked,side_by_side}
                        Composite layout. Default: stacked.
  --background-color BACKGROUND_COLOR
  --no-encode           Write validated composite PNGs without encoding MP4.
  --ffmpeg-bin FFMPEG_BIN
  --ffprobe-bin FFPROBE_BIN
  --crf CRF
  --movie-name MOVIE_NAME
  --output-dir OUTPUT_DIR
  --overwrite
```

## canonical-views

```text
usage: cryorole canonical-views [-h] --run-dir RUN_DIR
                                [--canonical-id CANONICAL_ID]
                                --chimerax-session CHIMERAX_SESSION
                                --map-frame {raw,explicit}
                                [--map-frame-transform MAP_FRAME_TRANSFORM]
                                [--render-mode {script-only,execute}]
                                [--chimerax-bin CHIMERAX_BIN] [--width WIDTH]
                                [--height HEIGHT] [--save-sessions]
                                --output-dir OUTPUT_DIR [--overwrite]

options:
  -h, --help            show this help message and exit
  --run-dir RUN_DIR     Existing cryoROLE run bundle.
  --canonical-id CANONICAL_ID
  --chimerax-session CHIMERAX_SESSION
                        Preconfigured .cxs session.
  --map-frame {raw,explicit}
                        Scene basis of the input session.
  --map-frame-transform MAP_FRAME_TRANSFORM
  --render-mode {script-only,execute}
  --chimerax-bin CHIMERAX_BIN
  --width WIDTH
  --height HEIGHT
  --save-sessions       Also save one camera-only .cxs session for each
                        canonical axis.
  --output-dir OUTPUT_DIR
  --overwrite
```

## canonicalize

```text
usage: cryorole canonicalize [-h] [--run-dir RUN_DIR]
                             [--canonical-id CANONICAL_ID] [--overwrite]
                             [--write-csv] [--fit-top FIT_TOP_FRACTION]
                             [--positive-side {low,high}]
                             [--use-frame USE_FRAME] [--no-visualize]
                             [--projection-panels] [--jobs JOBS]
                             [--point-style {scatter,density}]

options:
  -h, --help            show this help message and exit
  --run-dir RUN_DIR     Run bundle. Omitted: the current directory if it is a
                        bundle, else ./cryorole_outputs.
  --canonical-id CANONICAL_ID
                        Canonical output id. Default: default.
  --overwrite           Replace canonical/<canonical-id>/ only; the raw
                        landscape and other frames are unchanged.
  --write-csv           Also write canonical_landscape.csv. Off by default: it
                        is large and canonical_landscape.npz holds the same
                        data.
  --fit-top FIT_TOP_FRACTION
                        Fraction of highest-sld_raw points used to fit
                        canonical axes. Default: 0.40.
  --positive-side {low,high}
                        Canonical axis direction: low or high density skew.
                        Default: low.
  --use-frame USE_FRAME
                        Existing canonical_frame.json to apply instead of
                        fitting.
  --no-visualize        Skip canonical quick-look visualization.
  --projection-panels   Also write each 2D projection panel as its own PNG (18
                        extra figures). The 3-view figures already show every
                        panel.
  --jobs JOBS           Processes for the preview figures: 1-3 or auto.
                        Default: auto (up to 3, limited by CPU cores and free
                        memory; 1 below 100k particles).
  --point-style {scatter,density}
                        How 2D projection figures draw points: scatter
                        (default; one marker per point) or density (one image
                        cell per screen pixel holding the visible SLD value;
                        about 15x faster at 500k points, near-identical look,
                        slightly thinner in sparse regions).

Examples:
  cryorole canonicalize --run-dir my_run
  cryorole canonicalize --run-dir my_run --canonical-id fit50 --fit-top 0.5
  cryorole canonicalize --run-dir other_run --use-frame my_run/canonical/default/canonical_frame.json
```

## select

```text
usage: cryorole select [-h] [--run-dir RUN_DIR] --selection-id NAME
                       [--space {raw,canonical}] [--canonical-id CANONICAL_ID]
                       [--mode {radius,threshold,range,random,metadata}]
                       [--center A B C]
                       [--center-representation {euler,rotvec}]
                       [--radius RADIUS] [--radius-rad RADIUS_RAD]
                       [--metric {so3,rotvec}] [--sld-min SLD_MIN]
                       [--sld-max SLD_MAX] [--range-bound AXIS:LOWER:UPPER]
                       [--fraction FRACTION] [--seed SEED]
                       [--metadata-domain {ref,mov}]
                       [--metadata-column METADATA_COLUMN]
                       [--metadata-value METADATA_VALUE] [--split-by-value]
                       [--write-selected-landscape] [--recompute-sld]
                       [--overwrite]

Save a scientific subset for inspection and export. Choose a selection name and a mode. Selection evaluates the full parent landscape, independently of display sampling or filters.

options:
  -h, --help            show this help message and exit

Common inputs and naming:
  --run-dir RUN_DIR     Run bundle. Omitted: the current directory if it is a
                        bundle, else ./cryorole_outputs.
  --selection-id NAME   Required user-chosen name, e.g. region_01. Saved under
                        RUN/selections/NAME/. No default name.
  --space {raw,canonical}
                        Parent coordinate space. Default: raw.
  --canonical-id CANONICAL_ID
                        Canonical frame; only with --space canonical. Omitted:
                        the only frame in the bundle.
  --mode {radius,threshold,range,random,metadata}
                        radius: nearby orientations; threshold: SLD bounds;
                        range: coordinate bounds; random: random subset;
                        metadata: source values. Default: radius.

radius mode:
  Select near a center using SO(3) geodesic distance by default. Requires --center and exactly one radius option. Euler conventions come from the parent landscape.

  --center A B C, -c A B C
                        Center in the chosen --space: Euler degrees or rotvec
                        radians.
  --center-representation {euler,rotvec}
                        Representation of --center. Default: euler.
  --radius RADIUS, -r RADIUS
                        Radius in degrees; exclusive with --radius-rad.
  --radius-rad RADIUS_RAD
                        Radius in radians; independent of center units.
  --metric {so3,rotvec}
                        so3: rotation geodesic; rotvec: Euclidean rotation-
                        vector distance. Default: so3. Neither uses Euler
                        Euclidean distance.

threshold mode:
  Select by sld_raw, including both boundaries. Requires at least one bound; display colors and filters do not control selection.

  --sld-min SLD_MIN     Inclusive lower SLD bound (>= 0).
  --sld-max SLD_MAX     Inclusive upper SLD bound (>= lower bound).

range mode:
  Select the intersection of coordinate-axis bounds. Use alpha/beta/gamma in degrees OR x/y/z in radians; do not mix representations. Euler lower > upper wraps across the periodic seam; rotvec bounds must be ordered.

  --range-bound AXIS:LOWER:UPPER
                        Required; repeat for multiple axes. Bounds are
                        inclusive. Blank sides are open, e.g. alpha::20. At
                        least one side must be constrained. Repeat of the same
                        axis uses the last bound.

random mode:
  Sample without replacement from all parent rows; ceil(fraction * row count) rows are selected. This creates a Selection, unlike visualization sampling.

  --fraction FRACTION   Required fraction, 0 < F <= 1.
  --seed SEED           Non-negative random seed. Omitted: a seed is
                        generated; either way it is recorded in
                        selection.json.

metadata mode:
  Use run-time CryoSPARC CS or RELION STAR metadata and recorded source-row provenance. CS supports scalar integers, booleans, and UTF-8 text; empty strings are excluded. Requires domain, column, and either values or split. No rematching or external annotation file.

  --metadata-domain {ref,mov}
                        Required source domain; never inferred.
  --metadata-column METADATA_COLUMN
                        Required source particle column, e.g.
                        alignments3D/class (CS) or rlnClassNumber (STAR).
  --metadata-value METADATA_VALUE
                        Values to include (union), e.g. 0,1. CS uses typed
                        exact matching; booleans accept true/false/1/0.
                        Exclusive with --split-by-value.
  --split-by-value      One standard Selection per matched value: NAME_VALUE,
                        with filename-safe components. CS limit: 100 groups;
                        collisions fail before writing. Exclusive with
                        --metadata-value.

Output controls:
  --write-selected-landscape
                        Also write selected_landscape/ for visualize
                        --selection-id NAME --use-selected-landscape.
  --recompute-sld       Requires --write-selected-landscape. Recompute subset
                        SLD and preserve parent SLD fields separately.
  --overwrite           Explicitly replace this selection ID only (generated
                        child IDs in split mode); parent landscapes and
                        unrelated selections remain unchanged.

Examples (replace RUN with your run directory):
  cryorole select --run-dir RUN --selection-id region_01 --mode radius --center 0 0 0 --radius 15
  cryorole select --run-dir RUN --selection-id high_sld --mode threshold --sld-min 2
  cryorole select --run-dir RUN --selection-id alpha_window --mode range --range-bound alpha:-20:20
  cryorole select --run-dir RUN --selection-id sample_10pct --mode random --fraction 0.1 --seed 7
  cryorole select --run-dir RUN --selection-id cs_classes --mode metadata --metadata-domain ref --metadata-column alignments3D/class --metadata-value 0,1
  cryorole select --run-dir RUN --selection-id ref_classes --mode metadata --metadata-domain ref --metadata-column rlnClassNumber --metadata-value 1,3
  cryorole select --run-dir RUN --selection-id ref_classes --mode metadata --metadata-domain ref --metadata-column rlnClassNumber --split-by-value
```

## export

```text
usage: cryorole export [-h] [--run-dir RUN_DIR] [--selection-id SELECTION_ID]
                       [--selection SELECTION] [--domain {ref,mov,both}]
                       [--format {auto,relion_star,cryosparc_cs,keys}]
                       [--output-dir OUTPUT_DIR] [--overwrite]
                       [--relocated-ref RELOCATED_REF]
                       [--relocated-mov RELOCATED_MOV]
                       [--allow-unverified-source]
                       [--cryosparc-group-ref CRYOSPARC_GROUP_REF]
                       [--cryosparc-group-mov CRYOSPARC_GROUP_MOV]
                       [--no-cryosparc-group]
                       {selection} ...

positional arguments:
  {selection}
    selection           Compatibility alias for cryorole export.

options:
  -h, --help            show this help message and exit
  --run-dir RUN_DIR     Run bundle. Omitted: the current directory if it is a
                        bundle, else ./cryorole_outputs.
  --selection-id SELECTION_ID
                        Selection under RUN/selections/. Omitted: the only
                        selection in the bundle.
  --selection SELECTION
                        Advanced: direct path to a selection.json artifact.
  --domain {ref,mov,both}
                        Source metadata domain to export: ref, mov, or both.
                        Default: both.
  --format {auto,relion_star,cryosparc_cs,keys}
                        Output format; auto resolves per source domain.
                        Default: auto.
  --output-dir OUTPUT_DIR
                        Advanced: override output directory. Default:
                        RUN/exports/<selection_id>/.
  --overwrite           Replace export artifacts in the output directory only;
                        source, run, and selection files are unchanged.
  --relocated-ref RELOCATED_REF
                        Explicit relocated ref source; content must match the
                        recorded SHA-256.
  --relocated-mov RELOCATED_MOV
                        Explicit relocated mov source; content must match the
                        recorded SHA-256.
  --allow-unverified-source
                        Advanced legacy override for bundles that predate
                        source hashes; recorded in export report.
  --cryosparc-group-ref CRYOSPARC_GROUP_REF
                        CryoSPARC result group (.csg) of the ref source job to
                        mirror, e.g. PROJECT/J75/J75_particles.csg. Omitted:
                        the only one next to the ref .cs, if any.
  --cryosparc-group-mov CRYOSPARC_GROUP_MOV
                        As --cryosparc-group-ref, for the mov source.
  --no-cryosparc-group  Write the .cs subsets alone, without mirroring a
                        CryoSPARC result group (.csg), e.g. when the source
                        job's passthrough files are not available.

Examples:
  cryorole export --run-dir my_run --selection-id region_01
  cryorole export --run-dir my_run --selection-id region_01 --domain mov --format relion_star
```

## manifest

```text
usage: cryorole manifest [-h] --output OUTPUT [--workflow-name WORKFLOW_NAME]
                         [--command-string COMMAND_STRING] [--overwrite]
                         [--compute-file-hashes]
                         [--hash-algorithm HASH_ALGORITHM]
                         [--include-selected-particle-keys]

options:
  -h, --help            show this help message and exit
  --output OUTPUT       Manifest JSON output path.
  --workflow-name WORKFLOW_NAME
  --command-string COMMAND_STRING
  --overwrite
  --compute-file-hashes
  --hash-algorithm HASH_ALGORITHM
  --include-selected-particle-keys
                        Include selected particle keys in the manifest.
```
