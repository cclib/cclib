"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

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
        self.data = getfile(Gaussian, "basicGaussian03", "dvb_sp.out")
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
        self.data = getfile(Gaussian, "basicGaussian03", "dvb_sp.out")
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
        self.data = getfile(Gaussian, "basicGaussian03", "dvb_sp.out")
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
