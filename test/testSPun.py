import os
import unittest
import bettertest
import Numeric

from testall import getfile
from cclib.parser import ADF, GAMESS, Gaussian, Jaguar, GAMESSUK

class GenericSPunTest(bettertest.TestCase):
    """Restricted single point calculations with MO coeffs and overlap info."""
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        self.assertEquals(self.data.aooverlaps.shape,(self.data.nbasis,self.data.nbasis))

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        self.assertEquals(self.data.mocoeffs.shape,(2,self.data.nmo,self.data.nbasis))

    def testhomos(self):
        """Are the homos correct?"""
        self.assertArrayEquals(self.data.homos,Numeric.array([34,33],"i"),"%s != array([34,33],'i')" % Numeric.array_repr(self.data.homos))

class GaussianSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = data[0]

class GamessUSSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = data[1]

class PCGamessSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = data[2]

class ADFSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = data[3]

    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        #ADF uses fooverlaps
        self.assertEquals(self.data.fooverlaps.shape,(self.data.nbasis,self.data.nbasis))

class Jaguar42SPunTest(GenericSPunTest):
    def setUp(self):
        self.data = data[4]

class Jaguar65SPunTest(GenericSPunTest):
    def setUp(self):
        self.data = data[5]
        
class GamessUKSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = data[6]

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x (max(homos)+5) x nbasis?"""
        self.assertEquals(self.data.mocoeffs.shape,(2,max(self.data.homos)+6,self.data.nbasis))


names = [ "Gaussian", "PCGamess", "GAMESS", "ADF", "Jaguar 4.2",
          "Jaguar 6.5", "GAMESS UK"]
tests = [ GaussianSPunTest, PCGamessSPunTest,
          GamessUSSPunTest, ADFSPunTest,
          Jaguar42SPunTest, Jaguar65SPunTest,
          GamessUKSPunTest ]
data = [ getfile(Gaussian,"basicGaussian03","dvb_un_sp.out"),
         getfile(GAMESS,"basicGAMESS-US","dvb_un_sp.out"),
         getfile(GAMESS,"basicPCGAMESS","dvb_un_sp.out"),
         getfile(ADF,"basicADF2004.01","dvb_un_sp.adfout"),
         getfile(Jaguar, "basicJaguar4.2", "dvb_un_sp.out"),
         getfile(Jaguar, "basicJaguar6.5", "dvb_un_sp.out"),
         getfile(GAMESSUK, "basicGAMESS-UK", "dvb_un_sp_b.out")]
              
if __name__=="__main__":
    total = errors = failures = 0

    for name,test in zip(names,tests):
        print "\n**** Testing %s SPun ****" % name
        myunittest = unittest.makeSuite(test)
        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)

    print "\n\n********* SUMMARY OF SPun **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
