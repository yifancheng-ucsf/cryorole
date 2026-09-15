"""Executable regression coverage for the generated offline viewer."""

import json
import re
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

from cryorole.models.landscape import Landscape
from cryorole.visualize.service import _write_interactive_3d


def test_generated_viewer_script_parses_and_handles_large_display(tmp_path):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required for generated JavaScript execution")
    landscape = Landscape(data=pd.DataFrame({
        "particle_key": ['particle"\\line\n</script>', "p2"],
        "coordinates_analysis": [np.array([0.1, 0.2, 0.3]), np.array([-0.2, 0.3, 0.1])],
        "coordinates_display": [np.array([0.1, 0.2, 0.3]), np.array([-0.2, 0.3, 0.1])],
        "sld_display": [1.5, 2.5],
        "sld_unfloored": [1.5, 2.5], "sld_raw": [1.5, 2.5],
        "sld_was_floored": [False, False], "sld_local_k_mean": [1., 1.],
        "sld_effective_local_k_mean": [1., 1.], "sld_distance_floor": [0., 0.],
    }))
    _write_interactive_3d(
        landscape, output_dir=tmp_path, coordinate_source="analysis",
        representation="both", euler_sequence="zyx", colormap="rainbow_r",
        color_vmin=1.5, color_vmax=2.5, max_points=50_000, random_seed=0,
    )
    html = (tmp_path / "landscape_3d.html").read_text(encoding="utf-8")
    scripts = re.findall(r"<script[^>]*>([\s\S]*?)</script>", html)
    assert len(scripts) == 2
    payload = json.loads(scripts[0])
    assert payload["particle_keys"][0] == landscape.data.particle_key[0]
    assert payload["selection_enabled"] is False
    # A canvas stub verifies initialization and large-array drawing separately
    # from the real-browser sign-off. It is not a browser rendering substitute.
    harness = r'''
const vm = require('vm');
const payload = JSON.parse(process.argv[1]);
const script = process.argv[2];
let painted = 0;
const ctx = {clearRect(){}, beginPath(){}, arc(){}, fill(){painted++}};
const elements = {
  data: {textContent: JSON.stringify(payload)},
  plot: {getContext(){return ctx}, getBoundingClientRect(){return {left:0,top:80,width:800,height:500}}, setPointerCapture(){}},
  rep: {value:'rotvec', appendChild(){}}, tip:{style:{}}, reset:{}, summary:{}
};
const sandbox = {document:{getElementById:id=>elements[id],createElement:()=>({})},
  innerWidth:800,innerHeight:600,devicePixelRatio:1,addEventListener(){}};
vm.createContext(sandbox);
vm.runInContext(script,sandbox);
if(painted !== 2) throw Error('Initial drawing failed');
const before = vm.runInContext('JSON.stringify(screen)', sandbox);
elements.plot.onpointerdown({clientX:100,clientY:100,pointerId:1});
elements.plot.onpointermove({clientX:200,clientY:180});
elements.plot.onpointerup();
if(vm.runInContext('JSON.stringify(screen)', sandbox) === before) throw Error('Rotation failed');
elements.reset.onclick();
if(vm.runInContext('JSON.stringify(screen)', sandbox) !== before) throw Error('Reset failed');
elements.plot.onwheel({deltaY:100,preventDefault(){}});
if(vm.runInContext('JSON.stringify(screen)', sandbox) === before) throw Error('Zoom failed');
elements.reset.onclick();
const point = vm.runInContext('screen[0]', sandbox);
elements.plot.onpointermove({clientX:point[0],clientY:80+point[1]*500/elements.plot.height});
if(!elements.tip.textContent.includes(payload.particle_keys[point[3]])) throw Error('Hover key failed');
if(!elements.tip.textContent.includes('RV x (rad)')) throw Error('Hover coordinates failed');
elements.rep.value='euler'; elements.rep.onchange();
if(vm.runInContext('JSON.stringify(screen)', sandbox) === before) throw Error('Representation switch failed');
elements.rep.value='rotvec';
painted=0;
vm.runInContext("d.representations.rotvec.values = Array.from({length:150000}, () => [0.1,0.2,0.3]); d.colors = Array(150000).fill('#fff'); draw()",sandbox);
if(painted !== 150000) throw Error('Large display failed');
'''
    result = subprocess.run(
        [node, "-e", harness, json.dumps(payload), scripts[1]],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, result.stderr
