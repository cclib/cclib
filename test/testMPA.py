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
from cclib.method import MPA
from cclib.parser import Gaussian


class GenericMPATest(bettertest.TestCase):
    """MPA test"""
    def testsum(self):
        """Do the Mulliken charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos)
        totalmpa = sum(self.MPA.fragcharges)
        self.assertInside(totalmpa, formalcharge, 0.001)
        
class GaussianMPATest(GenericMPATest):
    def setUp(self):
        self.data = getfile(Gaussian, "basicGaussian03", "dvb_sp.out")
        self.MPA = MPA(self.data)
        self.MPA.logger.setLevel(0)
        self.MPA.calculate()

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(GaussianMPATest)
        
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
