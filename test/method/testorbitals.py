# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2016 the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test the various population analyses (MPA, LPA, CSPA) in cclib"""

import sys
import os
import logging
import unittest

import numpy

from cclib.method import Orbitals
from cclib.parser import Gaussian
from cclib.parser import Psi4

sys.path.insert(1, "..")

from ..test_data import getdatafile


class RestrictedCalculationTest(unittest.TestCase):
    """Check retricted calculation."""
    def setUp(self):
        self.data, self.logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])

    def test_closed_shell(self):
        self.assertTrue(Orbitals(self.data).closed_shell())


class UnrestrictedCalculationTest(unittest.TestCase):
    """Check unrestricted calculation."""
    def setUp(self):
        self.data, self.logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_un_sp.log"])

    def test_closed_shell(self):
        self.assertFalse(Orbitals(self.data).closed_shell())


class RestrictedOpenShellCalculationTest(unittest.TestCase):
    """Check restricted open shell calcualtion."""
    def setUp(self):
        self.data, self.logfile = getdatafile(Psi4, "basicPsi4-1.0", ["dvb_sp_rohf.out"])

    def test_closed_shel(self):
        self.assertFalse(Orbitals(self.data).closed_shell())

# TODO: add a case (regression) with an unrestricted calculation for a closed shell system.
# For example, in regressions: Gaussian/Gaussian03/Mo4OSibdt2


tests = [RestrictedCalculationTest, UnrestrictedCalculationTest]


if __name__ == "__main__":
    for test in tests:
        thistest = unittest.makeSuite(test)
        unittest.TextTestRunner(verbosity=2).run(thistest)
