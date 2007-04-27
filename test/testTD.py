import os, unittest

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

from testall import getfile
from cclib.parser import ADF, GAMESS, Gaussian, Jaguar, GAMESSUK
import bettertest


class GenericTDTest(bettertest.TestCase):
    """Time-dependent HF or DFT calculations."""
    def testenergiesnumber(self):
        """Is the length of etenergies correct?"""
        self.assertEqual(len(self.data.etenergies), self.number)
    
    def testoscsnumber(self):
        """Is the length of eotscs correct?"""
        self.assertEqual(len(self.data.etoscs), self.number)

    def testrotatsnumber(self):
        """Is the length of etrotats correct?"""
        self.assertEqual(len(self.data.etrotats), self.number)

    def testsecsnumber(self):
        """Is the length of etsecs correct?"""
        self.assertEqual(len(self.data.etsecs), self.number)

    def testsymsnumber(self):
        """Is the length of etsyms correct?"""
        self.assertEqual(len(self.data.etsyms), self.number)

class GaussianTDDFTTest(GenericTDTest):
    def setUp(self):
        self.data = data[0]
        self.number = 24
    
    def testsyms(self):
        """Is etsyms populated by singlets and triplets 50/50?"""
        singlets = [sym for sym in self.data.etsyms if "Singlet" in sym]
        triplets = [sym for sym in self.data.etsyms if "Triplet" in sym]
        self.assertEqual(len(singlets), self.number/2)
        self.assertEqual(len(triplets), self.number/2)

names = [ "Gaussian" ]
tests = [ GaussianTDDFTTest ]
data = [ getfile(Gaussian, "basicGaussian03","CO_TD_delta.log") ]
              
if __name__=="__main__":
    total = errors = failures = 0

    for name,test in zip(names,tests):
        print "\n**** Testing %s ****" % name
        myunittest = unittest.makeSuite(test)
        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)

    print "\n\n********* SUMMARY OF TD TEST **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
