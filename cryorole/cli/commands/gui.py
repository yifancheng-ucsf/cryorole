"""CLI adapter for the desktop GUI."""

from __future__ import annotations


def gui_command(args) -> int:
    from cryorole.gui import launch

    return launch(no_gl=args.no_gl, run_dir=args.run_dir)
