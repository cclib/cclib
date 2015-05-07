# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""This script runs the regression framework in the cclib-data repostiory."""

from __future__ import print_function

import os
import sys


if __name__ == "__main__":

    # Assume the cclib-data repository is cloned in this directory.
    regression_dir = os.path.join("..", "data", "regression")
    sys.path.append(regression_dir)
    import regression

    opt_traceback = "--traceback" in sys.argv
    opt_status = "--status" in sys.argv

    # This can be used to limit the programs we want to run regressions for.
    which = [arg for arg in sys.argv[1:] if not arg in ["--status", "--traceback"]]

    regression.main(which, opt_traceback, opt_status, regression_dir)
