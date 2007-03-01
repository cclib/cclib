import os, unittest
from Numeric import alltrue, array, resize
from testall import getfile
from cclib.parser import ADF, GAMESS, Gaussian, Jaguar, GAMESSUK
import bettertest

class GenericMP2Test(bettertest.TestCase):
    """MP2 calculations."""
    level = 2
    def testsizeandshape(self):
        """(MP2) Are the dimensions of mpenergies correct?"""
        self.assertEqual(self.data.mpenergies.shape, (len(self.data.scfenergies),
                                                      self.level-1))

    def testchange(self):
        """(MP2) Are Moller-Plesset corrections negative?"""
        if self.level == 2:
            corrections = self.data.mpenergies[:,0] - self.data.scfenergies
        else:
            corrections = self.data.mpenergies[:,self.level-2] - self.data.mpenergies[:,self.level-3]
        self.failUnless(alltrue(corrections < 0.0))

class GenericMP3Test(GenericMP2Test):
    """MP3 calculations."""
    level = 3
    def testsizeandshape(self):
        """(MP3) Are the dimensions of mpenergies correct?"""
        super(GenericMP3Test,self).testsizeandshape()
    def testchange(self):
        """(MP3) Are Moller-Plesset corrections negative?"""
        super(GenericMP3Test,self).testchange()

class GenericMP4Test(GenericMP2Test):
    """MP4 calculations."""
    level = 4
    def testsizeandshape(self):
        """(MP4) Are the dimensions of mpenergies correct?"""
        super(GenericMP4Test,self).testsizeandshape()
    def testchange(self):
        """(MP4) Are Moller-Plesset corrections negative?"""
        super(GenericMP4Test,self).testchange()

class GenericMP5Test(GenericMP2Test):
    """MP5 calculations."""
    level = 5
    def testsizeandshape(self):
        """(MP5) Are the dimensions of mpenergies correct?"""
        super(GenericMP5Test,self).testsizeandshape()
    def testchange(self):
        """(MP5) Are Moller-Plesset corrections negative?"""
        super(GenericMP5Test,self).testchange()

class GAMESSMP2Test(GenericMP2Test):
    def setUp(self):
        self.data = data[0]
        
class GaussianMP2Test(GenericMP2Test):
    def setUp(self):
        self.data = data[1]
        
class GaussianMP3Test(GenericMP3Test):
    def setUp(self):
        self.data = data[2]

class GaussianMP4Test(GenericMP4Test):
    def setUp(self):
        self.data = data[3]

class GaussianMP5Test(GenericMP5Test):
    def setUp(self):
        self.data = data[4]

class GAMESSUKMP2Test(GenericMP2Test):
    def setUp(self):
        self.data = data[5]

class GAMESSUKMP3Test(GenericMP3Test):
    def setUp(self):
        self.data = data[6]

class PCGAMESSMP2Test(GenericMP2Test):
    def setUp(self):
        self.data = data[7]

class PCGAMESSMP3Test(GenericMP3Test):
    def setUp(self):
        self.data = data[8]

class PCGAMESSMP4Test(GenericMP4Test):
    def setUp(self):
        self.data = data[9]

class PCGAMESSMP4SDTQTest(GenericMP4Test):
    def setUp(self):
        self.data = data[10]

names = [ "GAMESS",
          "Gaussian", "Gaussian", "Gaussian", "Gaussian",
          "GAMESS UK", "GAMESS UK",
          "PC GAMESS", "PC GAMESS", "PC GAMESS", "PC GAMESS" ]
tests = [ GAMESSMP2Test, 
          GaussianMP2Test, GaussianMP3Test, GaussianMP4Test, GaussianMP5Test,
          GAMESSUKMP2Test, GAMESSUKMP3Test, PCGAMESSMP2Test, PCGAMESSMP3Test,
          PCGAMESSMP4Test, PCGAMESSMP4SDTQTest]
data = [getfile(GAMESS, "basicGAMESS-US", "water_mp2.out"),
        getfile(Gaussian, "basicGaussian03","water_mp2.log"),
        getfile(Gaussian, "basicGaussian03","water_mp3.log"),
        getfile(Gaussian, "basicGaussian03","water_mp4.log"),
        getfile(Gaussian, "basicGaussian03","water_mp5.log"),
        getfile(GAMESSUK, "basicGAMESS-UK","water_mp2.out"),
        getfile(GAMESSUK, "basicGAMESS-UK","water_mp3.out"),
        getfile(GAMESS, "basicPCGAMESS", "water_mp2.out"),
        getfile(GAMESS, "basicPCGAMESS", "water_mp3.out"),
        getfile(GAMESS, "basicPCGAMESS", "water_mp4.out"),
        getfile(GAMESS, "basicPCGAMESS", "water_mp4_sdtq.out"),
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
