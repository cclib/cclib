"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Rev$"

import os
import logging
import unittest

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

import bettertest
from testall import getfile
from cclib.method import MPA, CSPA
from cclib.parser import Gaussian


class GaussianMPATest(bettertest.TestCase):
    """Mulliken Population Analysis test"""
    def setUp(self):
        self.data = getfile(Gaussian, "basicGaussian03", "dvb_sp.out")
        self.analysis = MPA(self.data)
        self.analysis.logger.setLevel(0)
        self.analysis.calculate()
    def testsum(self):
        """Do the MPA charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos)
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
        formalcharge = sum(self.data.atomnos)
        totalpopulation = sum(self.analysis.fragcharges)
        self.assertInside(totalpopulation, formalcharge, 0.001)


tests = [GaussianMPATest, GaussianCSPATest]

       
if __name__ == "__main__":
    for test in tests:
        thistest = unittest.makeSuite(test)
        unittest.TextTestRunner(verbosity=2).run(thistest)
