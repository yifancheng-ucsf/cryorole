#!/usr/bin/env python
"""cryoROLE ChimeraX landscape viewer.

Revision: 2026-07-28 public release filename.

Single-file ChimeraX command script for visualizing cryoROLE 2.0 landscapes.

How to load this script in ChimeraX
-----------------------------------
open /path/to/cryorole_chimerax_viewer.py

Main commands
-------------
cryorole open /path/to/display_table.csv
cryorole open /path/to/display_table.csv space canonical basis euler metric sld_display top_fraction 0.4
cryorole open /path/to/display_table.csv space canonical basis rv metric sld_display threshold 1.0
cryorole update #1 threshold 1.2
cryorole update #1 top_fraction 0.2
cryorole update #1 color_range 0.6,1.6
cryorole update #1 fade_below 2.0 fade_alpha 0.15
cryorole update #1 space raw basis euler metric sld_raw
cryorole info #1
cryorole map #1 full_bounds true grid_spacing 1.0 sdev 1.0

Design notes
------------
This keeps the cryoROLE 1.0 workflow:
    open python script in ChimeraX -> register commands -> load CSV -> display colored 3D cloud

But the parser is schema-aware for cryoROLE 2.0 display tables:
    canonical_ea_zyx_alpha_deg / beta / gamma
    raw_ea_zyx_alpha_deg / beta / gamma
    canonical_rv_x_rad / y / z
    raw_rv_x_rad / y / z
    compact display_table columns: euler_alpha / beta / gamma and rotvec_x / y / z
    sld_display / sld_raw

It also keeps legacy cryoROLE 1.0 fallbacks:
    Alpha / Beta / Gamma
    RotVec_X / RotVec_Y / RotVec_Z
    RND

This file intentionally does not perform cryoROLE analysis, canonicalization,
or particle export. It is only a lightweight ChimeraX viewer for 3D landscape
inspection, density coloring, and threshold/top-fraction filtering.
"""

import os
import csv
import re
import numpy as np


# ----------------------------
# Header / schema utilities
# ----------------------------

_GREEK_TRANSLATION = {
    "α": "alpha",
    "β": "beta",
    "γ": "gamma",
    "Α": "alpha",
    "Β": "beta",
    "Γ": "gamma",
}


def _normalize_name(name):
    """Normalize a table column or user-supplied column name for robust matching."""
    if name is None:
        return ""
    s = str(name).strip().strip('"').strip("'")
    for k, v in _GREEK_TRANSLATION.items():
        s = s.replace(k, v)
    s = s.casefold()
    # Make common spellings equivalent.
    s = s.replace("rotationvector", "rotvec")
    s = s.replace("rotation_vector", "rotvec")
    s = s.replace("rotation-vector", "rotvec")
    s = s.replace("eulerangle", "ea")
    s = s.replace("euler_angle", "ea")
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def _detect_delimiter(header_line):
    """Infer a simple delimiter from the header line."""
    if "\t" in header_line:
        return "\t"
    if "," in header_line:
        return ","
    if ";" in header_line:
        return ";"
    return None  # whitespace-delimited


def _read_header(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    with open(path, "r", newline="") as f:
        first = f.readline()
    if not first:
        raise ValueError(f"File is empty: {path}")

    first = first.lstrip("\ufeff").rstrip("\n\r")
    delimiter = _detect_delimiter(first)
    if delimiter is None:
        header = first.split()
    else:
        header = next(csv.reader([first], delimiter=delimiter))
    header = [h.strip() for h in header]

    norm_to_index = {}
    for i, h in enumerate(header):
        n = _normalize_name(h)
        if n and n not in norm_to_index:
            norm_to_index[n] = i

    return header, delimiter, norm_to_index


def _index_for(norm_to_index, aliases):
    """Return first matching column index from alias list."""
    for a in aliases:
        n = _normalize_name(a)
        if n in norm_to_index:
            return norm_to_index[n]
    return None


def _canonical_space(space):
    if space is None:
        return "auto"
    s = _normalize_name(space)
    if s in ("canonical", "canon", "display"):
        return "canonical"
    if s in ("raw", "original"):
        return "raw"
    if s in ("auto", "default"):
        return "auto"
    if s in ("legacy", "old"):
        return "legacy"
    raise ValueError("space must be 'canonical', 'raw', 'auto', or 'legacy'")


def _canonical_basis(basis):
    if basis is None:
        return "euler"
    b = _normalize_name(basis)
    if b in ("euler", "ea", "zyx", "euler_zyx"):
        return "euler"
    if b in ("rv", "rotvec", "rotation_vector", "rotationvector"):
        return "rv"
    raise ValueError("basis must be 'euler' or 'rv'/'rotvec'")


def _euler_aliases_for_space(space):
    return [
        [f"{space}_ea_zyx_alpha_deg", f"{space}_ea_zyx_beta_deg", f"{space}_ea_zyx_gamma_deg"],
        [f"{space}_euler_zyx_alpha_deg", f"{space}_euler_zyx_beta_deg", f"{space}_euler_zyx_gamma_deg"],
        [f"{space}_euler_alpha_deg", f"{space}_euler_beta_deg", f"{space}_euler_gamma_deg"],
        [f"{space}_ea_alpha_deg", f"{space}_ea_beta_deg", f"{space}_ea_gamma_deg"],
        [f"{space}_alpha_deg", f"{space}_beta_deg", f"{space}_gamma_deg"],
        [f"{space}_alpha", f"{space}_beta", f"{space}_gamma"],
    ]


def _rv_aliases_for_space(space):
    return [
        [f"{space}_rv_x_rad", f"{space}_rv_y_rad", f"{space}_rv_z_rad"],
        [f"{space}_rotvec_x_rad", f"{space}_rotvec_y_rad", f"{space}_rotvec_z_rad"],
        [f"{space}_rotation_vector_x_rad", f"{space}_rotation_vector_y_rad", f"{space}_rotation_vector_z_rad"],
        [f"{space}_rv_x", f"{space}_rv_y", f"{space}_rv_z"],
        [f"{space}_rotvec_x", f"{space}_rotvec_y", f"{space}_rotvec_z"],
    ]


def _compact_euler_aliases():
    """Aliases used by current cryoROLE 2.x display_table.csv files.

    These tables often store only one coordinate source at a time and use a
    separate 'coordinate_source' column with values such as 'canonical'.
    """
    return [
        ["euler_alpha", "euler_beta", "euler_gamma"],
        ["ea_alpha", "ea_beta", "ea_gamma"],
        ["alpha_deg", "beta_deg", "gamma_deg"],
        ["alpha", "beta", "gamma"],
    ]


def _compact_rv_aliases():
    """Aliases used by current cryoROLE 2.x display_table.csv files."""
    return [
        ["rotvec_x", "rotvec_y", "rotvec_z"],
        ["rv_x", "rv_y", "rv_z"],
        ["rotation_vector_x", "rotation_vector_y", "rotation_vector_z"],
        ["x", "y", "z"],
    ]


def _infer_coordinate_source(path, delimiter, norm_to_index, max_rows=1000):
    """Infer coordinate source from a compact display table, if present.

    Returns 'canonical', 'raw', another string value, or None.  This is only
    used for metadata/reporting and for choosing default metrics; coordinate
    columns are still resolved from the header.
    """
    idx = _index_for(norm_to_index, ["coordinate_source", "coord_source", "source", "space"])
    if idx is None:
        return None

    values = []
    try:
        with open(path, "r", newline="") as f:
            # skip header
            next(f, None)
            if delimiter is None:
                for _, line in zip(range(max_rows), f):
                    parts = line.strip().split()
                    if len(parts) > idx and parts[idx]:
                        values.append(_normalize_name(parts[idx]))
            else:
                reader = csv.reader(f, delimiter=delimiter)
                for _, row in zip(range(max_rows), reader):
                    if len(row) > idx and row[idx]:
                        values.append(_normalize_name(row[idx]))
    except Exception:
        return None

    values = [v for v in values if v]
    if not values:
        return None
    unique = sorted(set(values))
    if len(unique) == 1:
        v = unique[0]
        if v in ("canonical", "canon", "display"):
            return "canonical"
        if v in ("raw", "original"):
            return "raw"
        return v
    if "canonical" in unique and "raw" not in unique:
        return "canonical"
    if "raw" in unique and "canonical" not in unique:
        return "raw"
    return "mixed"


def _legacy_euler_aliases():
    return [
        ["Alpha", "Beta", "Gamma"],
        ["alpha", "beta", "gamma"],
        ["α", "β", "γ"],
        ["a", "b", "g"],
    ]


def _legacy_rv_aliases():
    return [
        ["RotVec_X", "RotVec_Y", "RotVec_Z"],
        ["RotVecX", "RotVecY", "RotVecZ"],
        ["RV_X", "RV_Y", "RV_Z"],
        ["rv_x", "rv_y", "rv_z"],
        ["x", "y", "z"],
    ]


def _match_triplet(norm_to_index, triplet_alias_groups):
    """Return (x_idx, y_idx, z_idx, matched_aliases) from possible triplet alias groups."""
    for group in triplet_alias_groups:
        idx = [_index_for(norm_to_index, [name]) for name in group]
        if all(i is not None for i in idx):
            return tuple(idx), tuple(group)
    return None, None


def _resolve_xyz_columns(norm_to_index, space="auto", basis="euler", coordinate_source=None):
    """Find coordinate columns for requested space/basis."""
    space = _canonical_space(space)
    basis = _canonical_basis(basis)

    spaces_to_try = [space] if space != "auto" else ["canonical", "raw"]

    triplet_groups = []
    if basis == "euler":
        for sp in spaces_to_try:
            if sp != "legacy":
                triplet_groups.extend(_euler_aliases_for_space(sp))
        triplet_groups.extend(_compact_euler_aliases())
        triplet_groups.extend(_legacy_euler_aliases())
    else:
        for sp in spaces_to_try:
            if sp != "legacy":
                triplet_groups.extend(_rv_aliases_for_space(sp))
        triplet_groups.extend(_compact_rv_aliases())
        triplet_groups.extend(_legacy_rv_aliases())

    triplet, matched_aliases = _match_triplet(norm_to_index, triplet_groups)
    if triplet is None:
        raise ValueError(
            f"Could not find {space}/{basis} coordinate columns. "
            "Expected cryoROLE 2.0 columns such as "
            "'canonical_ea_zyx_alpha_deg/beta/gamma', "
            "compact display_table columns 'euler_alpha/euler_beta/euler_gamma' "
            "or 'rotvec_x/rotvec_y/rotvec_z', or legacy columns "
            "'Alpha/Beta/Gamma' or 'RotVec_X/Y/Z'."
        )

    # Infer the resolved space from the matched names.  Current display_table.csv
    # files may use compact coordinate names plus coordinate_source=canonical/raw.
    matched_joined = "_".join(_normalize_name(a) for a in matched_aliases)
    resolved_space = "legacy"
    if "canonical" in matched_joined:
        resolved_space = "canonical"
    elif "raw" in matched_joined:
        resolved_space = "raw"
    elif coordinate_source in ("canonical", "raw"):
        resolved_space = coordinate_source
    elif space in ("canonical", "raw"):
        # User explicitly requested a space but the file uses compact columns.
        resolved_space = space
    elif any(_normalize_name(a) in ("euler_alpha", "rotvec_x", "rv_x") for a in matched_aliases):
        resolved_space = "display"

    return triplet, resolved_space, basis, matched_aliases


def _resolve_metric_column(norm_to_index, metric=None, space="auto"):
    """Find color/density metric column."""
    if metric is not None:
        idx = _index_for(norm_to_index, [metric])
        if idx is None:
            raise ValueError(f"Could not find metric column '{metric}'")
        return idx, metric

    # Default preference. sld_display is the cryoROLE 2.0 display metric.
    # RND is the cryoROLE 1.0 fallback.
    space = _canonical_space(space)
    candidates = []
    if space == "raw":
        candidates.extend(["sld_raw", "rnd_raw", "RND_raw"])
    candidates.extend([
        "sld_display",
        "rnd_display",
        "RND",
        "rnd",
        "sld_raw",
        "sld",
        "density",
        "color",
        "colour",
    ])

    idx = _index_for(norm_to_index, candidates)
    if idx is None:
        raise ValueError(
            "Could not find a density/color metric column. "
            "Expected one of: sld_display, sld_raw, RND, rnd."
        )

    # Report the first candidate that matched.
    for c in candidates:
        if _normalize_name(c) in norm_to_index and norm_to_index[_normalize_name(c)] == idx:
            return idx, c
    return idx, "metric"


def _read_numeric_columns(path, delimiter, usecols):
    """Read only selected numeric columns from a possibly wide CSV/TSV table."""
    data = np.genfromtxt(
        path,
        delimiter=delimiter,
        skip_header=1,
        usecols=list(usecols),
        dtype=np.float32,
        comments=None,
        invalid_raise=False,
    )
    if data.size == 0:
        raise ValueError("No numeric data rows were read.")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def _load_xyz_metric(path, space="auto", basis="euler", metric=None):
    header, delimiter, norm_to_index = _read_header(path)
    coordinate_source = _infer_coordinate_source(path, delimiter, norm_to_index)
    requested_space = _canonical_space(space)
    if (
        requested_space in ("canonical", "raw")
        and coordinate_source in ("canonical", "raw")
        and requested_space != coordinate_source
    ):
        raise ValueError(
            f"This compact display table reports coordinate_source={coordinate_source!r}, "
            f"but you requested space={requested_space!r}.  Load the matching "
            "raw/canonical display_table.csv, or omit the space keyword."
        )
    xyz_cols, resolved_space, resolved_basis, xyz_aliases = _resolve_xyz_columns(
        norm_to_index, space=space, basis=basis, coordinate_source=coordinate_source
    )
    metric_col, resolved_metric = _resolve_metric_column(norm_to_index, metric=metric, space=resolved_space)

    usecols = list(xyz_cols) + [metric_col]
    data = _read_numeric_columns(path, delimiter, usecols)
    xyz = data[:, 0:3].astype(np.float32, copy=False)
    color_values = data[:, 3].astype(np.float32, copy=False)

    ok = np.isfinite(xyz).all(axis=1) & np.isfinite(color_values)
    if not np.all(ok):
        xyz = xyz[ok]
        color_values = color_values[ok]

    if len(xyz) == 0:
        raise ValueError("No finite coordinate/color rows remain after reading the table.")

    meta = {
        "path": path,
        "header": header,
        "delimiter": delimiter,
        "xyz_columns": xyz_cols,
        "metric_column": metric_col,
        "space": resolved_space,
        "basis": resolved_basis,
        "metric": resolved_metric if metric is None else metric,
        "coordinate_source": coordinate_source,
        "xyz_aliases": xyz_aliases,
        "n_total_loaded": len(xyz),
    }
    return xyz, color_values, meta


# ----------------------------
# ChimeraX point cloud logic
# ----------------------------

def cryorole_open(
    session,
    open_path,
    space="auto",
    basis="euler",
    metric=None,
    palette=None,
    color_range=None,
    threshold=None,
    min_value=None,
    max_value=None,
    top_fraction=None,
    fade_below=None,
    fade_alpha=None,
):
    """Load a cryoROLE 2.0/legacy table as a colored 3D point cloud."""
    xyz_all, color_all, meta = _load_xyz_metric(open_path, space=space, basis=basis, metric=metric)

    from chimerax.core.models import Surface
    name = f"cryoROLE {meta['space']} {meta['basis']} ({len(xyz_all)} points)"
    s = Surface(name, session)
    s.SESSION_SAVE_DRAWING = True
    s.display_style = s.Dot
    s._cryorole_xyz_all = xyz_all
    s._cryorole_color_all = color_all
    s._cryorole_meta = meta
    s._cryorole_palette = palette
    s._cryorole_color_range = color_range
    s._cryorole_fade_below = fade_below
    s._cryorole_fade_alpha = fade_alpha
    session.models.add([s])

    _apply_filter_and_color(
        session,
        s,
        palette=palette,
        color_range=color_range,
        threshold=threshold,
        min_value=min_value,
        max_value=max_value,
        top_fraction=top_fraction,
        fade_below=fade_below,
        fade_alpha=fade_alpha,
    )
    return s


# Convenience alias.  ChimeraX command functions must explicitly have
# an initial argument named "session"; using *args/**kwargs here triggers
# ValueError("Missing initial 'session' argument") in recent ChimeraX builds.
def cryorole_load(
    session,
    open_path,
    space="auto",
    basis="euler",
    metric=None,
    palette=None,
    color_range=None,
    threshold=None,
    min_value=None,
    max_value=None,
    top_fraction=None,
    fade_below=None,
    fade_alpha=None,
):
    return cryorole_open(
        session,
        open_path,
        space=space,
        basis=basis,
        metric=metric,
        palette=palette,
        color_range=color_range,
        threshold=threshold,
        min_value=min_value,
        max_value=max_value,
        top_fraction=top_fraction,
        fade_below=fade_below,
        fade_alpha=fade_alpha,
    )


def cryorole_update(
    session,
    points_model,
    space=None,
    basis=None,
    metric=None,
    palette=None,
    color_range=None,
    threshold=None,
    min_value=None,
    max_value=None,
    top_fraction=None,
    fade_below=None,
    fade_alpha=None,
):
    """Update coordinates, filtering, and coloring of an existing cryoROLE point cloud."""
    if not hasattr(points_model, "_cryorole_meta"):
        raise ValueError("This model does not appear to be a cryoROLE point cloud.")

    meta = points_model._cryorole_meta

    # If coordinate space/basis/metric changes, re-read the source file.
    requested_space = space if space is not None else meta.get("space", "auto")
    requested_basis = basis if basis is not None else meta.get("basis", "euler")
    requested_metric = metric if metric is not None else meta.get("metric", None)

    need_reload = (
        space is not None
        or basis is not None
        or metric is not None
    )

    if need_reload:
        xyz_all, color_all, new_meta = _load_xyz_metric(
            meta["path"],
            space=requested_space,
            basis=requested_basis,
            metric=requested_metric,
        )
        points_model._cryorole_xyz_all = xyz_all
        points_model._cryorole_color_all = color_all
        points_model._cryorole_meta = new_meta
        meta = new_meta

    if palette is None:
        palette = getattr(points_model, "_cryorole_palette", None)
    if color_range is None:
        color_range = getattr(points_model, "_cryorole_color_range", None)

    points_model._cryorole_palette = palette
    points_model._cryorole_color_range = color_range

    _apply_filter_and_color(
        session,
        points_model,
        palette=palette,
        color_range=color_range,
        threshold=threshold,
        min_value=min_value,
        max_value=max_value,
        top_fraction=top_fraction,
        fade_below=fade_below,
        fade_alpha=fade_alpha,
    )


def _make_filter_mask(color_values, threshold=None, min_value=None, max_value=None, top_fraction=None):
    if threshold is not None:
        if min_value is not None and float(min_value) != float(threshold):
            raise ValueError("Use either threshold or min_value, not conflicting values for both.")
        min_value = threshold

    mask = np.ones(len(color_values), dtype=bool)

    if min_value is not None:
        mask &= (color_values >= float(min_value))
    if max_value is not None:
        mask &= (color_values <= float(max_value))

    if top_fraction is not None:
        tf = float(top_fraction)
        if tf <= 0 or tf > 1:
            raise ValueError("top_fraction must be in the range (0, 1].")
        cutoff = np.quantile(color_values, 1.0 - tf)
        mask &= (color_values >= cutoff)

    return mask


def _apply_filter_and_color(
    session,
    points_model,
    palette=None,
    color_range=None,
    threshold=None,
    min_value=None,
    max_value=None,
    top_fraction=None,
    fade_below=None,
    fade_alpha=None,
):
    xyz_all = points_model._cryorole_xyz_all
    color_all = points_model._cryorole_color_all
    meta = points_model._cryorole_meta

    mask = _make_filter_mask(
        color_all,
        threshold=threshold,
        min_value=min_value,
        max_value=max_value,
        top_fraction=top_fraction,
    )

    xyz = xyz_all[mask]
    color_values = color_all[mask]

    if len(xyz) == 0:
        raise ValueError("Filtering removed all points. Use a lower threshold or larger top_fraction.")

    colors = _point_colors(color_values, palette, color_range)

    if fade_below is not None:
        fb = float(fade_below)
        fa = 0.15 if fade_alpha is None else float(fade_alpha)
        fa = max(0.0, min(1.0, fa))
        low_mask = (color_values < fb)
        if np.any(low_mask):
            colors[low_mask, 3] = np.uint8(round(255.0 * fa))

    _update_point_cloud(points_model, xyz, colors)

    points_model._cryorole_last_mask = mask
    points_model._cryorole_last_filter = {
        "threshold": threshold,
        "min_value": min_value,
        "max_value": max_value,
        "top_fraction": top_fraction,
        "fade_below": fade_below,
        "fade_alpha": fade_alpha,
    }

    xyz_ranges = " ".join(
        f"({xyz[:, a].min():.3g}, {xyz[:, a].max():.3g})"
        for a in (0, 1, 2)
    )
    crange = f"({color_values.min():.5g}, {color_values.max():.5g})"
    total = len(color_all)
    shown = len(color_values)
    msg = (
        f"cryoROLE: shown {shown:,}/{total:,} points; "
        f"space={meta.get('space')} basis={meta.get('basis')} metric={meta.get('metric')}; "
        f"xyz ranges {xyz_ranges}; color range {crange}"
    )
    session.logger.info(msg)


def _update_point_cloud(points_model, xyz, colors):
    vertices = xyz.astype(np.float32, copy=True)
    from numpy import arange, int32
    dots = arange(len(vertices), dtype=int32).reshape((len(vertices), 1))
    points_model.set_geometry(vertices, None, dots)
    points_model.vertex_colors = colors


def _point_colors(color_values, palette, color_range):
    """Map scalar values to ChimeraX RGBA8 vertex colors."""
    from chimerax.surface.colorvol import _use_full_range, _colormap_with_range

    if _use_full_range(color_range, palette):
        color_range = (float(color_values.min()), float(color_values.max()))

    colormap = _colormap_with_range(palette, color_range, "rainbow")
    return colormap.interpolated_rgba8(color_values)


def cryorole_info(session, points_model):
    """Print metadata for an existing cryoROLE point cloud."""
    if not hasattr(points_model, "_cryorole_meta"):
        raise ValueError("This model does not appear to be a cryoROLE point cloud.")

    meta = points_model._cryorole_meta
    last_filter = getattr(points_model, "_cryorole_last_filter", {})
    n_total = len(points_model._cryorole_color_all)
    n_shown = len(points_model.vertices) if points_model.vertices is not None else 0

    lines = [
        "cryoROLE point cloud",
        f"  source: {meta.get('path')}",
        f"  total loaded rows: {n_total:,}",
        f"  currently shown: {n_shown:,}",
        f"  space: {meta.get('space')}",
        f"  basis: {meta.get('basis')}",
        f"  metric: {meta.get('metric')}",
        f"  xyz column indices: {meta.get('xyz_columns')}",
        f"  metric column index: {meta.get('metric_column')}",
        f"  matched xyz names: {meta.get('xyz_aliases')}",
        f"  fade_below: {getattr(points_model, '_cryorole_fade_below', None)}",
        f"  fade_alpha: {getattr(points_model, '_cryorole_fade_alpha', None)}",
        f"  last filter: {last_filter}",
    ]
    session.logger.info("\n".join(lines))


def cryorole_map(
    session,
    points_model,
    sdev=1.0,
    grid_spacing=1.0,
    cutoff_range=5,
    bounds=None,
    full_bounds=False,
):
    """Create a Gaussian density map from the currently visible cryoROLE points."""
    if not hasattr(points_model, "_cryorole_meta"):
        raise ValueError("This model does not appear to be a cryoROLE point cloud.")

    if full_bounds:
        basis = points_model._cryorole_meta.get("basis", "euler")
        if basis == "rv":
            bounds = (-3.14, 3.14, -3.14, 3.14, -3.14, 3.14)
        else:
            bounds = (-180, 180, -90, 90, -180, 180)

    xyz = points_model.vertices
    if xyz is None or len(xyz) == 0:
        raise ValueError("No visible points available for map creation.")

    from chimerax.map import molmap, volume_from_grid_data
    if bounds is None:
        grid = molmap.bounding_grid(xyz, grid_spacing, pad=0)
    else:
        if len(bounds) != 6:
            raise ValueError("bounds must contain 6 numbers: xmin xmax ymin ymax zmin zmax")
        grid = _bounding_grid(bounds, grid_spacing)

    from numpy import ones, float32
    weights = ones(len(xyz), float32)
    molmap.add_gaussians(grid, xyz, weights, sdev, cutoff_range)
    v = volume_from_grid_data(grid, session)
    v.name = f"cryoROLE map {len(xyz):,} points"
    return v


def _bounding_grid(bounds, spacing):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    nx = max(1, int(round((xmax - xmin) / spacing)))
    ny = max(1, int(round((ymax - ymin) / spacing)))
    nz = max(1, int(round((zmax - zmin) / spacing)))
    shape = (nz, ny, nx)
    origin = (xmin, ymin, zmin)
    from numpy import zeros, float32
    matrix = zeros(shape, float32)
    from chimerax.map_data import ArrayGridData
    return ArrayGridData(matrix, origin, (spacing, spacing, spacing))


# ----------------------------
# Command registration
# ----------------------------

def register_command(logger):
    from chimerax.core.commands import (
        register,
        CmdDesc,
        ModelArg,
        OpenFileNameArg,
        StringArg,
        ColormapArg,
        ColormapRangeArg,
        FloatArg,
        FloatsArg,
        BoolArg,
    )

    common_keywords = [
        ("space", StringArg),
        ("basis", StringArg),
        ("metric", StringArg),
        ("palette", ColormapArg),
        ("color_range", ColormapRangeArg),
        ("threshold", FloatArg),
        ("min_value", FloatArg),
        ("max_value", FloatArg),
        ("top_fraction", FloatArg),
        ("fade_below", FloatArg),
        ("fade_alpha", FloatArg),
    ]

    # IMPORTANT: ChimeraX CmdDesc instances are single-use.  In newer ChimeraX
    # versions (including 1.10 dev builds), reusing the same CmdDesc for an alias
    # raises: ValueError("Can not reuse CmdDesc instances").  Therefore, each
    # registered command gets its own fresh CmdDesc object, even if the schema is
    # identical.
    open_desc = CmdDesc(
        required=[("open_path", OpenFileNameArg)],
        keyword=list(common_keywords),
        synopsis="Load a cryoROLE 2.0 or legacy cryoROLE point-cloud table",
    )
    load_desc = CmdDesc(
        required=[("open_path", OpenFileNameArg)],
        keyword=list(common_keywords),
        synopsis="Alias of 'cryorole open'",
    )
    register("cryorole open", open_desc, cryorole_open, logger=logger)
    register("cryorole load", load_desc, cryorole_load, logger=logger)

    update_desc = CmdDesc(
        required=[("points_model", ModelArg)],
        keyword=common_keywords,
        synopsis="Update display coordinates, filtering, and coloring of a cryoROLE point cloud",
    )
    register("cryorole update", update_desc, cryorole_update, logger=logger)

    info_desc = CmdDesc(
        required=[("points_model", ModelArg)],
        synopsis="Show metadata for a cryoROLE point cloud",
    )
    register("cryorole info", info_desc, cryorole_info, logger=logger)

    map_desc = CmdDesc(
        required=[("points_model", ModelArg)],
        keyword=[
            ("sdev", FloatArg),
            ("grid_spacing", FloatArg),
            ("bounds", FloatsArg),
            ("full_bounds", BoolArg),
        ],
        synopsis="Create a Gaussian density map from currently visible cryoROLE points",
    )
    register("cryorole map", map_desc, cryorole_map, logger=logger)

    logger.info(
        "Registered cryoROLE 2.0 commands: "
        "cryorole open, cryorole load, cryorole update, cryorole info, cryorole map"
    )


register_command(session.logger)
