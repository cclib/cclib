# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

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
