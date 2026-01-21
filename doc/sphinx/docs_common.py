# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import os
import sys


def check_cclib(cclib):
    """Make sure we are importing code from a subdirectory, which should exist
    and should have been updated just before running this script. Note that
    this script does not assume any version in the module and just takes
    what it finds... so an appropriate checkout should be done first."""
    if cclib.__file__[: len(os.getcwd())] != os.getcwd():
        print("Do not seem to be importing from current directory")
        sys.exit(1)
