import os, unittest
from Numeric import array
from testall import getfile
from cclib.parser import ADF, GAMESS, Gaussian, Jaguar, GAMESSUK

class GenericSPTest(unittest.TestCase):
    """Restricted single point calculations with MO coeffs and overlap info."""
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        self.assertEquals(self.data.aooverlaps.shape,(self.data.nbasis,self.data.nbasis))

    def testatomcoords(self):
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        self.assertEquals(self.data.atomcoords.shape,(1,self.data.natom,3))
    
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        self.assertEquals(self.data.mocoeffs.shape,(1,self.data.nmo,self.data.nbasis))

class GaussianSPTest(GenericSPTest):
    def setUp(self):
        self.data = data[0]

class GamessUSSPTest(GenericSPTest):
    def setUp(self):
        self.data = data[1]

class PCGamessSPTest(GenericSPTest):
    def setUp(self):
        self.data = data[2]

class ADFSPTest(GenericSPTest):
    def setUp(self):
        self.data = data[3]
    
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        #ADF uses fooverlaps
        self.assertEquals(self.data.fooverlaps.shape,(self.data.nbasis,self.data.nbasis))

class JaguarSPTest(GenericSPTest):
    def setUp(self):
        self.data = data[4]

class GamessUKSPTest(GenericSPTest):
    def setUp(self):
        self.data = data[5]

names = [ "Gaussian", "PCGamess", "GAMESS", "ADF", "Jaguar", "GAMESS UK"]
tests = [ GaussianSPTest, PCGamessSPTest,
          GamessUSSPTest, ADFSPTest,
          JaguarSPTest, GamessUKSPTest]
data = [getfile(Gaussian, "basicGaussian03","dvb_sp.out"),
        getfile(GAMESS, "basicGAMESS-US","dvb_sp.out"),
        getfile(GAMESS, "basicPCGAMESS","dvb_sp.out"),
        getfile(ADF, "basicADF2004.01","dvb_sp_b.adfout"),
        getfile(Jaguar, "basicJaguar", "eg02", "dvb_sp.out"),
        getfile(GAMESSUK, "basicGAMESS-UK", "dvb_sp.out")]
              
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
