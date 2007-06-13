__revision__ = "$Revision$"

import os, unittest

import bettertest
from testall import gettestdata


class GenericBasisTest(bettertest.TestCase):
    """Any type of calculation that has basis set information."""

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

class GamessUKBasisTest(GenericBasisTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GamessUSBasisTest(GenericBasisTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GaussianBasisTest(GenericBasisTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class PCGamessBasisTest(GenericBasisTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]


# Load test data using information in file.
testdata = gettestdata(module="Basis")
              
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

    print "\n\n********* SUMMARY OF BASIS TEST **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
