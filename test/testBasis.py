import os, unittest
from Numeric import array
from testall import getfile
from cclib.parser import ADF, GAMESS, Gaussian, GAMESSUK

class GenericBasisTest(unittest.TestCase):
    """Some type of calculation so long as it has basis set information."""
    def testgbasis(self):
        """Is gbasis the right length?"""
        self.assertEquals(self.data.natom, len(self.data.gbasis))
    
    def testnames(self):
        """Test the names of the basis sets."""
        for atom in self.data.gbasis:
            for fns in atom:
                self.assert_(fns[0] in ['S', 'P'],
                             "%s not one of S or P" % fns[0])

    def testsizeofbasis(self):
        """Test the basis set size."""
        total = 0
        multiple = {'S':1, 'P':3}
        for atom in self.data.gbasis:
            for fns in atom:
                total += multiple[fns[0]] # Add 3 for P
        self.assertEquals(self.data.nbasis, total)
    
    def testcoeffs(self):
        """Test the coeffs of the basis sets."""
        for atom in self.data.gbasis:
            if len(atom)==1: # i.e. a 'H'
                coeffs = atom[0][1]
                self.assertAlmostEqual(coeffs[0][0], 3.42525, 5)
                self.assertAlmostEqual(coeffs[0][1], 0.15433, 5)
            else: # i.e. a 'C'
                self.assertEquals(len(atom), 3)
                s_coeffs = atom[1][1]
                p_coeffs = atom[2][1]
                self.assertAlmostEqual(s_coeffs[0][0], 2.9412, 4)
                self.assertAlmostEqual(p_coeffs[0][0], 2.9412, 4)
                self.assertAlmostEqual(s_coeffs[0][1], -0.1000, 4)
                self.assertAlmostEqual(p_coeffs[0][1], 0.1559, 4)

class GaussianBasisTest(GenericBasisTest):
    def setUp(self):
        self.data = data[0]

class GamessUSBasisTest(GenericBasisTest):
    def setUp(self):
        self.data = data[1]

class PCGamessBasisTest(GenericBasisTest):
    def setUp(self):
        self.data = data[2]

class GamessUKBasisTest(GenericBasisTest):
    def setUp(self):
        self.data = data[3]

names = [ "Gaussian", "PCGamess", "GAMESS", "GAMESS UK"]
tests = [ GaussianBasisTest, PCGamessBasisTest,
          GamessUSBasisTest, 
          GamessUKBasisTest]
data = [getfile(Gaussian, "basicGaussian03","dvb_sp_basis.log"),
        getfile(GAMESS, "basicGAMESS-US","dvb_sp.out"),
        getfile(GAMESS, "basicPCGAMESS","dvb_sp.out"),
        getfile(GAMESSUK, "basicGAMESS-UK", "dvb_sp.out")]
              
if __name__=="__main__":
    total = errors = failures = 0

    for name,test in zip(names,tests):
        print "\n**** Testing %s for Gbasis ****" % name
        myunittest = unittest.makeSuite(test)
        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)

    print "\n\n********* SUMMARY OF GBASIS **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
