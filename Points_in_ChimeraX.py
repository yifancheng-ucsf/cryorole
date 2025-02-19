#!/usr/bin/env python
""" Points_in_chimerax.py
This script registers ChimeraX commands to load, update, and display point cloud data from a text file.
Usage examples: open Points_in_chimerax.py # execute in ChimeraX command line
points load euler test.csv # Loads file using Alpha/Beta/Gamma columns; # full_bounds = (-180, 180, -90, 90, -180, 180)
points load rotvec test.csv # Loads file using RotVec_X/RotVec_Y/RotVec_Z columns; # full_bounds = (-3.14, 3.14, -3.14, 3.14, -3.14, 3.14)
points update #1 min_value 0.3 # Updates filtering and coloring on existing point cloud model #1
points map #1 full_bounds true # Creates a map from model #1 using full_bounds based on its mode """

import os
import numpy as np

def open_points_file(session, path, mode):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File {path} not found.")

    # Read the file, assuming the first line contains headers (tab-delimited)
    with open(path, 'r') as f:
        header = f.readline().strip().split('\t')

    # Select columns based on mode
    if mode.lower() == 'rotvec':
        x_idx = header.index('RotVec_X')
        y_idx = header.index('RotVec_Y')
        z_idx = header.index('RotVec_Z')
    elif mode.lower() == 'euler':
        x_idx = header.index('Alpha')
        y_idx = header.index('Beta')
        z_idx = header.index('Gamma')
    else:
        raise ValueError("mode must be either 'euler' or 'rotvec'")

    rnd_index = header.index('RND')
    xyz_columns = (x_idx, y_idx, z_idx)
    color_column = rnd_index

    # Load the data (assuming tab-delimited values)
    points = np.loadtxt(path, delimiter='\t', skiprows=1, dtype=np.float32)

    from chimerax.core.models import Surface
    s = Surface(f'{len(points)} points', session)
    s.SESSION_SAVE_DRAWING = True
    s.points = points
    s.display_style = s.Dot
    session.models.add([s])

    # Save column indices and mode for later use
    s.xyz_columns = xyz_columns
    s.color_column = color_column
    s.mode = mode.lower()
    return s

def points_load(session, mode, open_path, palette=None, color_range=None, min_value=None, max_value=None):
    """ Load a new point cloud file with the specified mode ("euler" or "rotvec") and apply filtering and coloring. """
    points_model = open_points_file(session, open_path, mode)
    filter_and_color_points(points_model, palette=palette, color_range=color_range,
                            min_value=min_value, max_value=max_value)

def points_update(session, points_model, palette=None, color_range=None, min_value=None, max_value=None):
    """ Update filtering and coloring on an existing point cloud model. """
    filter_and_color_points(points_model, palette=palette, color_range=color_range,
                            min_value=min_value, max_value=max_value)

def filter_and_color_points(points_model, xyz_columns=None, color_column=None, palette=None, color_range=None,
                            min_value=None, max_value=None):
    xyz_columns = xyz_columns if xyz_columns is not None else points_model.xyz_columns
    color_column = color_column if color_column is not None else points_model.color_column
    points_model.xyz_columns = xyz_columns
    points_model.color_column = color_column

    xyz = points_model.points[:, xyz_columns]
    color_values = points_model.points[:, color_column]

    # Apply filtering based on min_value and/or max_value
    if min_value is not None or max_value is not None:
        mask = np.ones_like(color_values, dtype=bool)
        if min_value is not None:
            mask &= (color_values >= min_value)
        if max_value is not None:
            mask &= (color_values <= max_value)
        xyz = xyz[mask]
        color_values = color_values[mask]

    colors = point_colors(color_values, palette, color_range)
    update_point_cloud(points_model, xyz, colors)

    xyz_ranges = ' '.join(f'({"%.1f" % xyz[:, a].min()}, {"%.1f" % xyz[:, a].max()})'
                          for a in (0, 1, 2))
    crange = f'({color_values.min():.5g}, {color_values.max():.5g})'
    msg = f'{len(xyz)} points shown, xyz ranges {xyz_ranges}, color range {crange}'
    session.logger.info(msg)

def update_point_cloud(points_model, xyz, colors):
    vertices = xyz.copy()
    from numpy import arange, int32
    dots = arange(len(vertices), dtype=int32).reshape((len(vertices), 1))
    points_model.set_geometry(vertices, None, dots)
    points_model.vertex_colors = colors

def point_colors(color_values, palette, color_range):
    from chimerax.surface.colorvol import _use_full_range, _colormap_with_range
    if _use_full_range(color_range, palette):
        color_range = (color_values.min(), color_values.max())
    colormap = _colormap_with_range(palette, color_range, 'rainbow')
    colors = colormap.interpolated_rgba8(color_values)
    return colors

def points_map(session, points_model, sdev=1.0, grid_spacing=1.0, cutoff_range=5, bounds=None, full_bounds=False):
    # If full_bounds is True, set bounds based on the mode stored in the model.
    if full_bounds:
        if hasattr(points_model, 'mode'):
            if points_model.mode == 'rotvec':
                bounds = (-3.14, 3.14, -3.14, 3.14, -3.14, 3.14)
            else: bounds = (-180, 180, -90, 90, -180, 180)
        else: bounds = (-180, 180, -90, 90, -180, 180)


    xyz = points_model.vertices
    from chimerax.map import molmap, volume_from_grid_data
    if bounds is None:
        grid = molmap.bounding_grid(xyz, grid_spacing, pad=0)
    else:
        grid = bounding_grid(bounds, grid_spacing)

    from numpy import ones, float32
    weights = ones(len(xyz), float32)  # Weight each point equally.
    molmap.add_gaussians(grid, xyz, weights, sdev, cutoff_range)
    v = volume_from_grid_data(grid, session)
    v.name = f'map {len(xyz)} points'
    return v

def bounding_grid(bounds, spacing):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    nx = int((xmax - xmin) / spacing)
    ny = int((ymax - ymin) / spacing)
    nz = int((zmax - zmin) / spacing)
    shape = (nz, ny, nx)
    origin = (xmin, ymin, zmin)
    from numpy import zeros, float32
    matrix = zeros(shape, float32)
    from chimerax.map_data import ArrayGridData
    grid = ArrayGridData(matrix, origin, (spacing, spacing, spacing))
    return grid

def register_command(logger):
    from chimerax.core.commands import (register, CmdDesc, ModelArg, OpenFileNameArg, StringArg,
                                        ColormapArg, ColormapRangeArg, FloatArg, FloatsArg, BoolArg)

    # Register the "points load" command.
    load_desc = CmdDesc(required=[('mode', StringArg), ('open_path', OpenFileNameArg)],
                        keyword=[('palette', ColormapArg),
                                ('color_range', ColormapRangeArg),
                                ('min_value', FloatArg),
                                ('max_value', FloatArg)],
                        synopsis='Load a point cloud file with specified mode ("euler" or "rotvec")')
    register('points load', load_desc, points_load, logger=logger)

    # Register the "points update" command.
    update_desc = CmdDesc(required=[('points_model', ModelArg)],
                        keyword=[('palette', ColormapArg),
                                ('color_range', ColormapRangeArg),
                                ('min_value', FloatArg),
                                ('max_value', FloatArg)],
                        synopsis='Update filtering and coloring of an existing point cloud')
    register('points update', update_desc, points_update, logger=logger)

    # Register the "points map" command.
    map_desc = CmdDesc(required=[('points_model', ModelArg)],
                    keyword=[('sdev', FloatArg),
                                ('grid_spacing', FloatArg),
                                ('bounds', FloatsArg),
                                ('full_bounds', BoolArg)],
                    synopsis='Create a map by placing a Gaussian at each point of a point cloud')
    register('points map', map_desc, points_map, logger=logger)
register_command(session.logger)