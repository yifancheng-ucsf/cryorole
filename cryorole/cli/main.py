"""Root parser, dispatch, and public compatibility exports for cryoROLE."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from cryorole.cli.commands.artifacts import (
    align_command,
    export_metadata_command,
    export_selection_command,
    manifest_command,
)
from cryorole.cli.commands.downstream import (
    animate_command,
    canonical_views_command,
    canonicalize_command,
    select_command,
    visualize_command,
)
from cryorole.cli.commands.workflow import (
    explore_command,
    guide_command,
    next_command,
    preflight_command,
    run_command,
    status_command,
)
from cryorole.cli.parsers import build_command_parser


def build_parser() -> argparse.ArgumentParser:
    """Build the public parser from argparse-only registrations."""

    return build_command_parser(
        {
            "align": align_command,
            "run": run_command,
            "preflight": preflight_command,
            "status": status_command,
            "next": next_command,
            "guide": guide_command,
            "explore": explore_command,
            "visualize": visualize_command,
            "animate": animate_command,
            "canonical_views": canonical_views_command,
            "canonicalize": canonicalize_command,
            "select": select_command,
            "export": export_selection_command,
            "export_metadata": export_metadata_command,
            "manifest": manifest_command,
        }
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the cryoROLE CLI.

    Expected problems are printed as ``cryorole: error: <message>`` plus a
    ``what to do:`` line and exit with status 2. Ctrl-C exits with 130.
    Unexpected errors print a short report (exit 1); set ``CRYOROLE_DEBUG=1``
    to see the traceback instead.
    """

    import os

    from cryorole.errors import classify_exception
    from cryorole.logs import configure_cli_logging

    parser = build_parser()
    args = parser.parse_args(argv)
    configure_cli_logging(quiet=bool(getattr(args, "quiet", False)))
    try:
        return args.handler(args)
    except KeyboardInterrupt:
        parser.exit(130, "cryorole: interrupted; nothing further was written.\n")
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001 - the CLI boundary reports every failure cleanly
        if os.environ.get("CRYOROLE_DEBUG"):
            raise
        error = classify_exception(exc)
        text = f"cryorole: error: {error.message}\n"
        if error.remedy:
            text += f"  what to do: {error.remedy}\n"
        parser.exit(1 if error.code == "internal" else 2, text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
