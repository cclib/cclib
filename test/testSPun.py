__revision__ = "$Revision$"

import os
import unittest
import bettertest

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

from testall import getfile
from cclib.parser import ADF, GAMESS, Gaussian, Jaguar, GAMESSUK


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
        
class GaussianSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = data[0]

    def testatomnos(self):
        """Does atomnos have the right dimension (20)?"""
        size = len(self.data.atomnos)
        self.assertEquals(size, 20)

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
    
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis? PASS"""
        self.assertEquals(1,1)

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""
        self.assertEquals(1,1)

    def testmosyms(self):
        """Are the dims of the mosyms equal to 2 x nmo? PASS"""
        self.assertEquals(1,1)

class Jaguar60SPunTest(GenericSPunTest):
    def setUp(self):
        self.data = data[5]

    def testmosyms(self):
        """Are the dims of the mosyms equal to 2 x nmo? PASS"""
        self.assertEquals(1,1)
        
class Jaguar65SPunTest(GenericSPunTest):
    def setUp(self):
        self.data = data[6]
        
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis? PASS"""
        self.assertEquals(1,1)

    def testmoenergies(self):
        """Are the dims of the moenergies equal to 2 x homos+11?"""
        self.assertEquals(len(self.data.moenergies), 2)
        self.assertEquals(len(self.data.moenergies[0]), self.data.homos[0]+11)
        self.assertEquals(len(self.data.moenergies[1]), self.data.homos[1]+11)
        
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""
        self.assertEquals(1,1)

    def testmosyms(self):
        """Are the dims of the mosyms equal to 2 x nmo? PASS"""
        self.assertEquals(1,1)

class GamessUKSPunTest(GenericSPunTest):
    def setUp(self):
        self.data = data[7]

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 2 x (homos+6) x nbasis?"""
        self.assertEquals(type(self.data.mocoeffs), type([]))
        self.assertEquals(len(self.data.mocoeffs), 2)
        self.assertEquals(self.data.mocoeffs[0].shape,
                          (self.data.homos[0]+6, self.data.nbasis))
        self.assertEquals(self.data.mocoeffs[1].shape,
                          (self.data.homos[1]+6, self.data.nbasis))

names = [ "Gaussian", "PCGamess", "GAMESS", "ADF",
          "Jaguar 4.2", "Jaguar 6.0", "Jaguar 6.5",
          "GAMESS UK"]
tests = [ GaussianSPunTest, PCGamessSPunTest,
          GamessUSSPunTest, ADFSPunTest,
          Jaguar42SPunTest, Jaguar60SPunTest, Jaguar65SPunTest,
          GamessUKSPunTest ]
data = [ getfile(Gaussian,"basicGaussian03","dvb_un_sp_b.log"),
         getfile(GAMESS,"basicGAMESS-US","dvb_un_sp.out"),
         getfile(GAMESS,"basicPCGAMESS","dvb_un_sp.out"),
         getfile(ADF,"basicADF2004.01","dvb_un_sp.adfout"),
         getfile(Jaguar, "basicJaguar4.2", "dvb_un_sp.out"),
         getfile(Jaguar, "basicJaguar6.0", "dvb_un_sp.out"),
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

    print "\n\n********* SUMMARY OF SPUN TEST **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
