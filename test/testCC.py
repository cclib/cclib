import os, unittest
from Numeric import alltrue
from testall import getfile
from cclib.parser import GAMESS, Gaussian
import bettertest

class GenericCCTest(bettertest.TestCase):
    """Coupled-Cluster calculations."""
    def testsign(self):
        corrections = self.data.ccenergies - self.data.scfenergies
        self.failUnless(alltrue(corrections < 0.0))

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

class GaussianCCDTest(GenericCCDTest):
    def setUp(self):
        self.data = data[0]

class GaussianCCSDTest(GenericCCSDTest):
    def setUp(self):
        self.data = data[1]

class GaussianCCSDTTest(GenericCCSDTTest):
    def setUp(self):
        self.data = data[2]

class GAMESSCCDTest(GenericCCDTest):
    def setUp(self):
        self.data = data[3]

class GAMESSCCSDTest(GenericCCSDTest):
    def setUp(self):
        self.data = data[4]

class GAMESSCCSDTTest(GenericCCSDTTest):
    def setUp(self):
        self.data = data[5]

names = [ "Gaussian", "Gaussian", "Gaussian",
          "GAMESS", "GAMESS", "GAMESS" ]
tests = [ GaussianCCDTest, GaussianCCSDTest, GaussianCCSDTTest,
          GAMESSCCDTest, GAMESSCCSDTest, GAMESSCCSDTTest ]
data = [getfile(Gaussian, "basicGaussian03","water_ccd.log"),
        getfile(Gaussian, "basicGaussian03","water_ccsd.log"),
        getfile(Gaussian, "basicGaussian03","water_ccsd(t).log"),
        getfile(GAMESS, "basicGAMESS-US","water_ccd.out"),
        getfile(GAMESS, "basicGAMESS-US","water_ccsd.out"),
        getfile(GAMESS, "basicGAMESS-US","water_ccsd(t).out"),
        ]
              
if __name__=="__main__":
    total = errors = failures = 0

    for name,test in zip(names,tests):
        print "\n**** Testing %s ****" % name
        myunittest = unittest.makeSuite(test)
        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)

    print "\n\n********* SUMMARY OF MP TEST **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
