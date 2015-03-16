# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Run parser tests for cclib."""

import sys

from data import *


if __name__ == "__main__":

    # These allow the parsers and modules tested to be filtered on the command line
    # with any number of arguments. No matching parsers/modules implies all of them.
    parsers = [p for p in parsers if p in sys.argv] or parsers
    modules = [m for m in test_modules if m in sys.argv] or test_modules

    # These options are used for Travis CI.
    status = "status" in sys.argv or "--status" in sys.argv
    terse = "terse" in sys.argv or "--terse" in sys.argv

    tests = testall(parsers, modules, status, terse)
    visualtests()