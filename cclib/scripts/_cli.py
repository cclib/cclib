# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Helpers shared by cclib's command line scripts."""

import argparse
from pathlib import Path

import cclib


class VersionAction(argparse.Action):
    """Print the cclib version and install location, then exit.

    argparse's built-in ``version`` action re-wraps its text to the terminal
    width, which would split the install path across lines, so the message is
    written directly instead.
    """

    def __init__(
        self,
        option_strings,
        dest=argparse.SUPPRESS,
        default=argparse.SUPPRESS,
        help="show the cclib version and install location, then exit",
    ) -> None:
        super().__init__(
            option_strings=option_strings, dest=dest, default=default, nargs=0, help=help
        )

    def __call__(self, parser, namespace, values, option_string=None) -> None:
        # Deliberately stdout, not parser.exit(message=...), which writes to stderr.
        print(f"{parser.prog} (cclib {cclib.__version__})")
        print(f"installed at {Path(cclib.__file__).parent}")
        parser.exit()


def add_version_argument(parser: argparse.ArgumentParser) -> None:
    """Add a ``--version`` flag to a cclib script's argument parser."""
    parser.add_argument("--version", action=VersionAction)
