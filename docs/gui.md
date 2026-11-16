# Desktop GUI

The desktop workflow is **Inputs → Parameters → Run → View → Select → Export**.
You can also use **File → Open run** to review an existing CLI run bundle.
The GUI uses the same services, selection evaluator and standard artifacts as
the command line. For a first session, follow the [illustrated README](../README.md).

## First session with the demo

Follow the README's [three preparation commands](../README.md#3-cli-walkthrough)
before opening the demo. Its original STAR files need verified row-aligned
matching and its supplied display frame uses a non-default axis direction;
these two options are currently command-line only. `expected/` contains
comparison artifacts, not a complete run bundle.

Open `ino80hex_run` with **File → Open run...**, or start with
`cryorole gui --run-dir ino80hex_run`. In View, choose **Canonical: high**,
**Euler angles (deg)** and **All particles**. The README's
[GUI walkthrough](../README.md#4-gui-walkthrough)
selects a **Radius** of 10° around the typed centre α = 40, β = 0, γ = 0 and
saves it with **Confirm selection**: 5,009 demo particles, the same as the
command-line `select` step. The centre and radius are demo values, not a
scientific cutoff. Ordinary identity-matched inputs can start directly on Inputs.

## Install and start

Install the optional desktop dependencies in the Python environment you intend
to use:

```bash
python -m pip install ".[gui]"   # from the cryorole checkout
cryorole gui
cryorole gui --run-dir /path/to/run
```

Linux Qt 6.5+ needs `libxcb-cursor0` for the X11 platform plugin. Use a terminal
in the workstation's own graphical session. A missing display is different
from an unusable OpenGL driver: the startup probe reports the reason before
opening the main window. Do not invent a `DISPLAY` value.

The GUI probes the actual 3D widget in a child process. If OpenGL is unavailable,
2D and 1D still work and the **3D view...** button explains the limitation.
`cryorole gui --no-gl` deliberately uses this mode. Software OpenGL is reported
as such. For a Linux desktop driven by an onboard/BMC adapter with a separate
NVIDIA GPU, PRIME offload may be needed:

```bash
__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia cryorole gui
```

The intended platforms are Linux and macOS. Offscreen tests verify 2D/1D,
workflow and artifact behavior; they do not establish physical desktop,
driver, 3D-picking or high-particle-count performance acceptance.

## Start or resume a run

Inputs accepts a `.star` or `.cs` file, or a refinement job folder for each
domain. Reference and moving body have distinct meanings; do not swap them to
make the display look preferable. Ambiguous job folders list candidates for
you to choose. Both domains must come from the same software. Choose particle
metadata containing 3D orientations, not a map or optimiser/model STAR. See
[CryoSPARC file choice](cryosparc_workflow.md#which-output-should-i-use) and
[RELION file choice](relion_workflow.md#which-output-should-i-use).

The GUI uses default identity matching. It has no row-aligned, explicit mapping
or low-overlap override control. For cases requiring those CLI options, first
prepare and run the verified inputs with the documented CLI workflow, then open
the resulting bundle here. Equal counts or nearby micrograph coordinates alone
do not establish a pair of corresponding particles.

Preflight checks are automatic. Changing an input or checked option immediately
invalidates the old result. **Run** stays disabled until the current checks
finish and are usable. Review warnings before starting. Parameters includes
the output folder, SLD settings and optional figures/CSV.

Run uses a child process. Progress and errors appear in the page and Log panel.
**Cancel** requests a stop at a stage boundary; cancellation publishes no run.
If it cannot stop promptly, **Stop now** terminates the worker; the error may
identify a staging folder left for recovery. Closing during Run or
Canonicalize asks before stopping. Canonicalize writes a new frame rather than
replacing the raw landscape. **Canonicalize...** offers **Frame name**,
**Fit to top** (percentage) and preview figures. The GUI uses the default axis
sign policy; `--positive-side high` or `--use-frame` require the CLI. The fit
percentage controls frame fitting, not how many particles are retained.

## Inspect a landscape

Choose Raw or a canonical frame and Euler angles (degrees) or rotation vectors
(radians). The initial display filter is `sld_display >= 1`. Top SLD percent
and All particles are alternatives. Display filters, axis ranges, the colour
scale and the deterministic 3D sample change no scientific selection.

2D and 1D use every row passing the display filter. Only 3D uses its point cap.
The 3D view opens in its own window and follows the same settings. Standard,
Fit core, Fit all and Custom axis ranges apply consistently across views.

### Save figures

Click **Save figures...**, choose a new figure name and either 2D projections
or 2D plus 1D distributions. The action captures the current run, space,
representation, display filter, colour scale and **Axis range** settings.
Pending filter edits are applied before capturing. Changing run or controls
while saving does not redirect the export. The figures are written in a child
process, like a run, so the window stays responsive while they render.

Use Custom ranges to preserve a particular range in the figures. Mouse
pan/zoom, saved-selection overlays and draft highlights are not included.
These are standard landscape figures, not a screenshot or selection export.
The name is never overwritten. A rendering error publishes no partial figure
set; choose a new name or correct the error and retry. Let saving finish before
closing the window; the GUI keeps the window open while the writer is active.

The result is under `visualizations/raw/<name>/` or
`visualizations/canonical/<canonical-id>/<name>/`, with PNGs and
`visualization_report.json`. Its `gui_display_snapshot` records run/landscape
identity, the overlay/viewport exclusions and the equivalent CLI command.
Successful actions also appear in the Commands panel.

### Review an existing selection

**Show saved selection** loads a previously confirmed Selection from this run.
**Refresh list** includes selections created in another window. The overlay
shows its full selected count and how many rows pass the current display
filter. A zero visible count does not mean the saved selection is empty.

Overlaying is read-only: it creates no draft, enables no Confirm action and
writes no selection or metadata. It matches saved particle keys to the full
current landscape, including across raw/canonical views. A different run clears
the overlay. A late load after clearing or changing run is discarded.

An active new draft takes visual priority; the overlay label states when the
draft highlight is shown. Clearing the draft reveals the saved overlay again.
Missing keys, count mismatches and another run's identity are refused. Legacy
selections without a recorded parent run ID are not overlaid; open their
original compatible workflow or create a new verified Selection.

## Confirm and export particles

The Select page evaluates a new radius or other selection over the full parent
landscape, independent of display filters and samples. A click proposes a
centre at the nearest displayed particle; it is not a saved decision. Radius
selection is a neighbourhood in 3D rotation space, not a rectangle or lasso in
a 2D panel. Radius is in degrees even when centre coordinates use rotation-vector
radians. Wait for the exact count, choose a new
selection name and use **Confirm selection** to create the normal Selection.
For non-radius modes such as **Random fraction**, the action is labelled
**Create selection**. Fractions are entered as percentages in the GUI; the CLI
uses a fraction from 0 to 1. **0.10%** in the GUI is `--fraction 0.001`.

Switching Euler/RV coordinates converts the typed centre rather than
reinterpreting its units. The saved interaction records what was displayed.

Export consumes a confirmed selection and recorded ref/mov source-row IDs. It
does not reselect, rematch or rewrite poses. Choose domain and a new output
folder. Source hashes and row provenance are verified. Existing folders are
kept. Errors include a remedy and leave no partly published metadata export.

For CryoSPARC, a result group describes the source job's `.cs` and passthrough
files. The preview reports an unusable/missing group instead of offering it as
ready. Choose the correct `.csg` when needed. The completed export lists its
files and **Import Result Group** instructions; files must be copied inside
the CryoSPARC project directory by an account allowed to write there. The GUI
does not copy them into a project or submit an external job automatically.

## Acceptance checklist

Before release, check normal startup on Linux and macOS, no-GL fallback,
separate 3D window/picking, and representative particle-count performance on
real target desktops. Exercise current-input preflight, Run/Cancel/Stop now,
Canonicalize failure/close, run switches during background operations,
Selection confirmation and source-subset export. Check Save figures against
the displayed settings and a standard CLI render. Record screenshots, platform,
driver and limits. Synthetic/offscreen results do not complete this checklist.
