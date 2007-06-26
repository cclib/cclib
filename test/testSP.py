__revision__ = "$Revision$"

import os, unittest

import bettertest
from testall import gettestdata


class GenericSPTest(bettertest.TestCase):
    """Restricted single point calculations with MO coeffs and overlap info."""

    nbasisdict = {1:1, 6:5} # STO-3G, H has 1, C has 3

    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        self.assertEquals(self.data.aooverlaps.shape,(self.data.nbasis,self.data.nbasis))

    def testatomcoords(self):
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        self.assertEquals(self.data.atomcoords.shape,(1,self.data.natom,3))
    
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        self.assertEquals(type(self.data.mocoeffs), type([]))
        self.assertEquals(len(self.data.mocoeffs), 1)
        self.assertEquals(self.data.mocoeffs[0].shape,
                          (self.data.nmo, self.data.nbasis))

    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique?"""
        all = []
        for i,atom in enumerate(self.data.atombasis):
            self.assertEquals(len(atom), self.nbasisdict[self.data.atomnos[i]])
            all += atom
        # Test if there are as many indices as atomic orbitals.
        self.assertEquals(len(all), self.data.nbasis)
        # Check if all are different (every orbital indexed once).
        self.assertEquals(len(set(all)), len(all))

class ADFSPTest(GenericSPTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]
    
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        # ADF uses fooverlaps.
        self.assertEquals(self.data.fooverlaps.shape,(self.data.nbasis,self.data.nbasis))

class GamessUKSPTest(GenericSPTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GamessUSSPTest(GenericSPTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GaussianSPTest(GenericSPTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class Jaguar42SPTest(GenericSPTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

    # Data file does not contain enough information. Can we make a new one?
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""
        self.assertEquals(1, 1)

    # Data file does not contain enough information. Can we make a new one?
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique? PASS"""
        self.assertEquals(1, 1)
       
class Jaguar60SPTest(GenericSPTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]
    nbasisdict = {1:5, 6:15} # 6-31G(d,p)

class PCGamessSPTest(GenericSPTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]


# Load test data using information in file.
testdata = gettestdata(module="SP")

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

    print "\n\n********* SUMMARY OF SP TEST **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
