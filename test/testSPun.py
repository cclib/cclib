import os, unittest
from cclib.parser import GAMESS,G03,ADF,Jaguar
from Numeric import array
from testall import getfile

class GenericSPunTest(unittest.TestCase):
    """Restricted single point calculations with MO coeffs and overlap info."""
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        self.assertEquals(self.data.aooverlaps.shape,(self.data.nbasis,self.data.nbasis))

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nindep x nbasis?"""
        self.assertEquals(self.data.mocoeffs.shape,(2,self.data.nindep,self.data.nbasis))

    def testhomos(self):
        """What are the homos?"""
        self.assertEquals(type(self.data.homos),type(array([])))
        self.assertEquals(self.data.homos,array([34,33],"i"))

class GaussianSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = getfile(G03,"basicGaussian03","dvb_un_sp.out")

class GamessUSSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = getfile(GAMESS,"basicGAMESS-US","dvb_un_sp.out")

class PCGamessSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = getfile(GAMESS,"basicPCGAMESS","dvb_un_sp.out")

class ADFSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = getfile(ADF,"basicADF2004.01","dvb_un_sp.adfout")

names = [ "Gaussian", "PCGamess", "GAMESS", "ADF" ]
tests = [ GaussianSPunTest, PCGamessSPunTest,
          GamessUSSPunTest, ADFSPunTest ]
              
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
