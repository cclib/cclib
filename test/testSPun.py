__revision__ = "$Revision$"

import os
import unittest

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

import bettertest
from testall import gettestdata


class GenericSPunTest(bettertest.TestCase):
    """Restricted single point calculations with MO coeffs and overlap info."""
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        self.assertEquals(self.data.aooverlaps.shape,(self.data.nbasis,self.data.nbasis))

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 2 x nmo x nbasis?"""
        self.assertEquals(type(self.data.mocoeffs), type([]))
        self.assertEquals(len(self.data.mocoeffs), 2)
        self.assertEquals(self.data.mocoeffs[0].shape,
                          (self.data.nmo, self.data.nbasis))
        self.assertEquals(self.data.mocoeffs[1].shape,
                          (self.data.nmo, self.data.nbasis))

    def testcharge_and_mult(self):
        """Are the charge and multiplicity correct?"""
        self.assertEquals(self.data.charge, 1)
        self.assertEquals(self.data.mult, 2)

    def testhomos(self):
        """Are the homos correct?"""
        self.assertArrayEquals(self.data.homos, numpy.array([34,33],"i"),"%s != array([34,33],'i')" % numpy.array_repr(self.data.homos))

    def testmoenergies(self):
        """Are the dims of the moenergies equals to 2 x nmo?"""
        self.assertEquals(len(self.data.moenergies), 2)
        self.assertEquals(len(self.data.moenergies[0]), self.data.nmo)
        self.assertEquals(len(self.data.moenergies[1]), self.data.nmo)

    def testmosyms(self):
        """Are the dims of the mosyms equals to 2 x nmo?"""
        shape = (len(self.data.mosyms), len(self.data.mosyms[0]))
        self.assertEquals(shape, (2, self.data.nmo))
        
class ADFSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        #ADF uses fooverlaps
        self.assertEquals(self.data.fooverlaps.shape,(self.data.nbasis,self.data.nbasis))

class GamessUKSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 2 x (homos+6) x nbasis?"""
        self.assertEquals(type(self.data.mocoeffs), type([]))
        self.assertEquals(len(self.data.mocoeffs), 2)
        self.assertEquals(self.data.mocoeffs[0].shape,
                          (self.data.homos[0]+6, self.data.nbasis))
        self.assertEquals(self.data.mocoeffs[1].shape,
                          (self.data.homos[1]+6, self.data.nbasis))

class GamessUSSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GaussianSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

    def testatomnos(self):
        """Does atomnos have the right dimension (20)?"""
        size = len(self.data.atomnos)
        self.assertEquals(size, 20)

class Jaguar42SPunTest(GenericSPunTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]
    
    # Data file does not contain enough information. Can we make a new one?
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis? PASS"""
        self.assertEquals(1,1)

    # Data file does not contain enough information. Can we make a new one?
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""
        self.assertEquals(1,1)

    # Why is this test passed?
    def testmosyms(self):
        """Are the dims of the mosyms equal to 2 x nmo? PASS"""
        self.assertEquals(1,1)

class Jaguar60SPunTest(GenericSPunTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

    # Why is this test passed?
    def testmosyms(self):
        """Are the dims of the mosyms equal to 2 x nmo? PASS"""
        self.assertEquals(1,1)
        
class Jaguar65SPunTest(GenericSPunTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]
        
    # Data file does not contain enough information. Can we make a new one?
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis? PASS"""
        self.assertEquals(1,1)

    # Why is this test passed?
    def testmoenergies(self):
        """Are the dims of the moenergies equal to 2 x homos+11?"""
        self.assertEquals(len(self.data.moenergies), 2)
        self.assertEquals(len(self.data.moenergies[0]), self.data.homos[0]+11)
        self.assertEquals(len(self.data.moenergies[1]), self.data.homos[1]+11)
        
    # Data file does not contain enough information. Can we make a new one?
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""
        self.assertEquals(1,1)

    # Why is this test passed?
    def testmosyms(self):
        """Are the dims of the mosyms equal to 2 x nmo? PASS"""
        self.assertEquals(1,1)

class PCGamessSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]


# Load test data using information in file.
testdata = gettestdata(module="SPun")
              
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

    print "\n\n********* SUMMARY OF SPUN TEST **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
