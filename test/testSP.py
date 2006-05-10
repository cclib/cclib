import os, unittest
from cclib.parser import GAMESS,G03,ADF,Jaguar
from Numeric import array
from testall import getfile

class GenericSPTest(unittest.TestCase):
    """Restricted single point calculations with MO coeffs and overlap info."""
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        self.assertEquals(self.data.aooverlaps.shape,(self.data.nbasis,self.data.nbasis))

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nindep x nbasis?"""
        self.assertEquals(self.data.mocoeffs.shape,(1,self.data.nindep,self.data.nbasis))

class GaussianSPTest(GenericSPTest):
    def setUp(self):
        self.data = getfile(G03,"basicGaussian03","dvb_sp.out")

class GamessUSSPTest(GenericSPTest):
    def setUp(self):
        self.data = getfile(GAMESS,"basicGAMESS-US","dvb_sp.out")

class PCGamessSPTest(GenericSPTest):
    def setUp(self):
        self.data = getfile(GAMESS,"basicPCGAMESS","dvb_sp.out")

class ADFSPTest(GenericSPTest):
    def setUp(self):
        self.data = getfile(ADF,"basicADF2004.01","dvb_sp_b.adfout")
    
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        #ADF uses fooverlaps
        self.assertEquals(self.data.fooverlaps.shape,(self.data.nbasis,self.data.nbasis))

names = [ "Gaussian", "PCGamess", "GAMESS", "ADF", "Jaguar" ]
tests = [ GaussianSPTest, PCGamessSPTest,
          GamessUSSPTest, ADFSPTest ]
              
if __name__=="__main__":
    total = errors = failures = 0

    for name,test in zip(names,tests):
        print "\n**** Testing %s SP ****" % name
        myunittest = unittest.makeSuite(test)
        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)

    print "\n\n********* SUMMARY OF SP **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
