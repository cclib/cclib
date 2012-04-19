# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

__revision__ = "$Revision$"

import os
import logging
import unittest

import numpy

import bettertest
from testall import getfile
from cclib.method import MPA, LPA, CSPA
from cclib.parser import Gaussian


class GaussianMPATest(bettertest.TestCase):
    """Mulliken Population Analysis test"""
    def setUp(self):
        self.data, self.logfile = getfile(Gaussian, "basicGaussian03", "dvb_sp.out")
        self.analysis = MPA(self.data)
        self.analysis.logger.setLevel(0)
        self.analysis.calculate()
    def testsum(self):
        """Do the Mulliken charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos) - self.data.charge
        totalpopulation = sum(self.analysis.fragcharges)
        self.assertInside(totalpopulation, formalcharge, 0.001)

class GaussianLPATest(bettertest.TestCase):
    """Lowdin Population Analysis test"""
    def setUp(self):
        self.data, self.logfile = getfile(Gaussian, "basicGaussian03", "dvb_sp.out")
        self.analysis = LPA(self.data)
        self.analysis.logger.setLevel(0)
        self.analysis.calculate()
    def testsum(self):
        """Do the Lowdin charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos) - self.data.charge
        totalpopulation = sum(self.analysis.fragcharges)
        self.assertInside(totalpopulation, formalcharge, 0.001)

class GaussianCSPATest(bettertest.TestCase):
    """C-squared Population Analysis test"""
    def setUp(self):
        self.data, self.logfile = getfile(Gaussian, "basicGaussian03", "dvb_sp.out")
        self.analysis = CSPA(self.data)
        self.analysis.logger.setLevel(0)
        self.analysis.calculate()
    def testsum(self):
        """Do the CSPA charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos) - self.data.charge
        totalpopulation = sum(self.analysis.fragcharges)
        self.assertInside(totalpopulation, formalcharge, 0.001)


tests = [GaussianMPATest, GaussianLPATest, GaussianCSPATest]

       
if __name__ == "__main__":
    for test in tests:
        thistest = unittest.makeSuite(test)
        unittest.TextTestRunner(verbosity=2).run(thistest)
