"""Run the real browser asset under a deterministic DOM when Node is available."""
from pathlib import Path
import shutil
import subprocess

import pytest


def test_browser_draft_contract():
    node=shutil.which('node')
    if node is None:
        pytest.skip('Node is unavailable; browser contract can also run with the audit V8 adapter')
    root=Path(__file__).resolve().parents[2]
    completed=subprocess.run([node,str(Path(__file__).with_name('browser_draft_contract.js')),
                              str(root/'cryorole/interactive/assets/app.js')],capture_output=True,text=True,timeout=30)
    assert completed.returncode==0,completed.stdout+completed.stderr
