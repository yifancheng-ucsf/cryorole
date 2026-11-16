"""Command lines for the GUI's child processes (worker, display probe).

``python -m MODULE`` puts the working directory first on ``sys.path``, so a
``json.py`` or a ``cryorole/`` folder where the GUI was started would be
imported instead of the installed modules. These commands run the module the
same way but drop that entry first (the ``-c`` form adds it as ``''``).
"""

from __future__ import annotations

import sys

_RUN_MODULE = (
    "import os, runpy, sys\n"
    "if sys.path and sys.path[0] in ('', os.getcwd()):\n"
    "    del sys.path[0]\n"
    "runpy.run_module({module!r}, run_name='__main__', alter_sys=True)\n"
)


def module_command(module: str, *args: str) -> list[str]:
    """``[python, ...]`` that runs ``module`` as ``__main__`` without the working directory on ``sys.path``."""

    return [sys.executable, "-c", _RUN_MODULE.format(module=module), *args]
