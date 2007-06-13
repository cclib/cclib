__revision__ = "$Revision$"

import os, unittest

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

import bettertest
from testall import gettestdata


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
        self.data = testdata[self.__class__.__name__]["data"]
        self.number = 24
    
    def testsyms(self):
        """Is etsyms populated by singlets and triplets 50/50?"""
        singlets = [sym for sym in self.data.etsyms if "Singlet" in sym]
        triplets = [sym for sym in self.data.etsyms if "Triplet" in sym]
        self.assertEqual(len(singlets), self.number/2)
        self.assertEqual(len(triplets), self.number/2)


# Load test data using information in file.
testdata = gettestdata(module="TD")
              
if __name__=="__main__":
    total = errors = failures = 0
    for test in testdata:
        module = testdata[test]["module"]
        print "\n**** test%s for %s ****" %(module, '/'.join(testdata[test]["location"]))
        test = eval(test)
        myunittest = unittest.makeSuite(test)
        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)

    print "\n\n********* SUMMARY OF TD TEST **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
