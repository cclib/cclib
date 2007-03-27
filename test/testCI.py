import os, unittest
from testall import getfile
from cclib.parser import Gaussian, GAMESS, Jaguar
import bettertest

class GenericCITest(bettertest.TestCase):
    """CI calculations."""
    nstates = 0
    def testetenergies(self):
        """ Are there nstate elements in etenergies?"""
        self.assertEqual(len(self.data.etenergies), self.nstates)

class GenericCISTest(GenericCITest):
    """CIS calculations."""

class GaussianCISTest(GenericCISTest):
    def setUp(self):
        self.data = data[0]
        self.nstates = 10
    def testetsecs(self):
        """ Are there nstate elements in etsecs?"""
        self.assertEqual(len(self.data.etsecs), self.nstates)

class GAMESSCISTest(GenericCISTest):
    def setUp(self):
        self.data = data[1]
        self.nstates = 10
    def testetsecs(self):
        """ Are there nstate elements in etsecs?"""
        self.assertEqual(len(self.data.etsecs), self.nstates)

class Jaguar65CISTest(GenericCISTest):
    def setUp(self):
        self.data = data[2]
        self.nstates = 5

names = [ "Gaussian", "GAMESS", "Jaguar" ]
tests = [ GaussianCISTest, GAMESSCISTest, Jaguar65CISTest ]
data = [getfile(Gaussian, "basicGaussian03", "water_cis.log"),
        getfile(GAMESS, "basicGAMESS-US", "water_cis.out"),
        getfile(Jaguar, "basicJaguar6.5", "water_cis.out")
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
