__revision__ = "$Revision$"

import os, unittest

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

import bettertest
from testall import gettestdata


class GenericCCTest(bettertest.TestCase):
    """Coupled-Cluster calculations."""

    def testsign(self):
        corrections = self.data.ccenergies - self.data.scfenergies
        self.failUnless(numpy.alltrue(corrections < 0.0))

class GenericCCDTest(GenericCCTest):
    """CCD calculations."""
    def testsign(self):
        """CCD: Are the Coupled-Cluster correction negative?"""
        super(GenericCCDTest, self).testsign()

class GenericCCSDTest(GenericCCTest):
    """CCSD calculations."""
    def testsign(self):
        """CCSD: Are the Coupled-Cluster correction negative?"""
        super(GenericCCSDTest, self).testsign()

class GenericCCSDTTest(GenericCCTest):
    """CCSD(T) calculations."""
    def testsign(self):
        """CCSD(T): Are the Coupled-Cluster correction negative?"""
        super(GenericCCSDTTest, self).testsign()

class GAMESSUSCCDTest(GenericCCDTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GAMESSUSCCSDTest(GenericCCSDTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GAMESSUSCCSDTTest(GenericCCSDTTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GaussianCCDTest(GenericCCDTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GaussianCCSDTest(GenericCCSDTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GaussianCCSDTTest(GenericCCSDTTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

              
# Load test data using information in file.
testdata = gettestdata(module="CC")

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

    print "\n\n********* SUMMARY OF CC TEST **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
