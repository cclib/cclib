import os, unittest
from Numeric import alltrue, array, resize
from testall import getfile
from cclib.parser import ADF, GAMESS, Gaussian, Jaguar, GAMESSUK
import bettertest

class GenericMPTest(bettertest.TestCase):
    """MPx calculations."""

    def testsizeandshape(self):
        """Are the dimensions of mpenergies correct?"""
        self.assertEqual(self.data.mpenergies.shape, (len(self.data.scfenergies),
                                                      self.level-1))

    def testchange(self):
        """Are Moller-Plesset corrections negative?"""
        if self.level == 2:
          corrections = self.data.mpenergies[:,0] - self.data.scfenergies
        else:
          corrections = self.data.mpenergies[:,self.level-2] - self.data.mpenergies[:,self.level-3]
        self.failUnless(alltrue(corrections < 0.0))

class GAMESSMP2Test(GenericMPTest):
    def setUp(self):
        self.data = data[0]
        self.level = 2

class GaussianMP2Test(GenericMPTest):
    def setUp(self):
        self.data = data[1]
        self.level = 2

class GaussianMP3Test(GenericMPTest):
    def setUp(self):
        self.data = data[2]
        self.level = 3

class GaussianMP4Test(GenericMPTest):
    def setUp(self):
        self.data = data[3]
        self.level = 4

class GaussianMP5Test(GenericMPTest):
    def setUp(self):
        self.data = data[4]
        self.level = 5

names = [ "GAMESS MP2",
          "Gaussian MP2", "Gaussian MP3", "Gaussian MP4", "Gaussian MP5" ]
tests = [ GAMESSMP2Test, 
          GaussianMP2Test, GaussianMP3Test, GaussianMP4Test, GaussianMP5Test ]
data = [getfile(GAMESS, "basicGAMESS-US", "water_mp2.out"),
        getfile(Gaussian, "basicGaussian03","water_mp2.log"),
        getfile(Gaussian, "basicGaussian03","water_mp3.log"),
        getfile(Gaussian, "basicGaussian03","water_mp4.log"),
        getfile(Gaussian, "basicGaussian03","water_mp5.log"),
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
