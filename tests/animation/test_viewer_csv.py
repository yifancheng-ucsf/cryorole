"""Standalone viewer parser contracts, without a ChimeraX runtime."""
import ast
import csv
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cryorole.export.visualization import write_landscape_visualizations
from cryorole.models.landscape import Landscape


@pytest.fixture(scope="module")
def load_table():
    source = Path(__file__).resolve().parents[2] / "cryorole_chimerax_viewer.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    registration = [node for node in tree.body if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                    and isinstance(node.value.func, ast.Name) and node.value.func.id == "register_command"]
    assert len(registration) == 1
    tree.body = [node for node in tree.body if node not in registration]
    namespace = {"__name__": "viewer_parser_test"}
    exec(compile(tree, str(source), "exec"), namespace)
    return namespace["_load_xyz_metric"]


def test_actual_public_display_table_preserves_quoted_particle_names(tmp_path, load_table):
    keys = ["plain", "particle,comma", 'particle"quote', "particle\nnext line", "粒子 μ"]
    coordinates = np.random.default_rng(0).normal(0, .1, (len(keys), 3))
    density = np.arange(1, len(keys) + 1, dtype=float)
    landscape = Landscape(data=pd.DataFrame({
        "particle_key": keys, "coordinates_analysis": list(coordinates), "coordinates_display": list(coordinates),
        "sld_unfloored": density, "sld_raw": density, "sld_display": density,
        "sld_was_floored": False, "sld_local_k_mean": 1., "sld_effective_local_k_mean": 1., "sld_distance_floor": 1e-4,
    }))
    write_landscape_visualizations(landscape, tmp_path, coordinate_source="analysis", representation="rotvec",
                                   formats=("png",), artifact_layout="run_bundle", display_table_filename="display_table.csv")
    source = tmp_path / "display_table.csv"
    before = source.read_bytes()
    xyz, values, metadata = load_table(str(source), space="raw", basis="rv", metric="sld_display")
    np.testing.assert_allclose(xyz, coordinates, atol=1e-7)
    np.testing.assert_array_equal(values, density)
    assert metadata["n_total_loaded"] == len(keys)
    assert source.read_bytes() == before


@pytest.mark.parametrize("delimiter", [",", "\t", ";"])
def test_utf8_bom_quoted_delimited_single_row(tmp_path, load_table, delimiter):
    source = tmp_path / "table.txt"
    with source.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream, delimiter=delimiter)
        writer.writerow(["particle_key", "coordinate_source", "rotvec_x", "rotvec_y", "rotvec_z", "sld_display"])
        writer.writerow(['粒子,tab\tsemi;"quote"\nnext', "raw", .1, -.2, .3, 4.])
    xyz, values, metadata = load_table(str(source), space="raw", basis="rv", metric="sld_display")
    np.testing.assert_allclose(xyz, [[.1, -.2, .3]], atol=1e-7)
    np.testing.assert_array_equal(values, [4.])
    assert metadata["space"] == "raw"


def test_legacy_whitespace_single_row_still_loads(tmp_path, load_table):
    source = tmp_path / "legacy.txt"
    source.write_text("particle Alpha Beta Gamma RND\np1 10 20 30 4\n")
    xyz, values, metadata = load_table(str(source), space="legacy", basis="euler")
    np.testing.assert_array_equal(xyz, [[10, 20, 30]])
    np.testing.assert_array_equal(values, [4])
    assert metadata["basis"] == "euler"


def test_existing_nonfinite_display_filter_is_preserved(tmp_path, load_table):
    source = tmp_path / "finite.csv"
    source.write_text("particle_key,rotvec_x,rotvec_y,rotvec_z,sld_display\nok,1,2,3,4\nbad,nan,2,3,4\n")
    xyz, values, metadata = load_table(str(source), basis="rv", metric="sld_display")
    np.testing.assert_array_equal(xyz, [[1, 2, 3]])
    np.testing.assert_array_equal(values, [4])
    assert metadata["n_total_loaded"] == 1
