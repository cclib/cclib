__revision__ = "$Revision$"

import os, unittest

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

import bettertest
from testall import gettestdata


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
        self.failUnless(numpy.alltrue(corrections < 0.0))

class GenericMP3Test(GenericMP2Test):
    """MP3 calculations."""
    level = 3
    
    def testsizeandshape(self):
        """(MP3) Are the dimensions of mpenergies correct?"""
        super(GenericMP3Test,self).testsizeandshape()
    def testchange(self):
        """(MP3) Are Moller-Plesset corrections negative?"""
        super(GenericMP3Test,self).testchange()

class GenericMP4SDQTest(GenericMP2Test):
    """MP4(SDQ) calculations."""
    level = 4
    
    def testsizeandshape(self):
        """(MP4-SDQ) Are the dimensions of mpenergies correct?"""
        super(GenericMP4SDQTest,self).testsizeandshape()
    def testchange(self):
        """(MP4-SDQ) Are Moller-Plesset corrections negative?"""
        super(GenericMP4SDQTest,self).testchange()

class GenericMP4SDTQTest(GenericMP2Test):
    """MP4(SDTQ) calculations."""
    level = 4
    
    def testsizeandshape(self):
        """(MP4-SDTQ) Are the dimensions of mpenergies correct?"""
        super(GenericMP4SDTQTest,self).testsizeandshape()
    def testchange(self):
        """(MP4-SDTQ) Are Moller-Plesset corrections negative?"""
        super(GenericMP4SDTQTest,self).testchange()

class GenericMP5Test(GenericMP2Test):
    """MP5 calculations."""
    level = 5
    
    def testsizeandshape(self):
        """(MP5) Are the dimensions of mpenergies correct?"""
        super(GenericMP5Test,self).testsizeandshape()
    def testchange(self):
        """(MP5) Are Moller-Plesset corrections negative?"""
        super(GenericMP5Test,self).testchange()

class GAMESSUKMP2Test(GenericMP2Test):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GAMESSUKMP3Test(GenericMP3Test):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GAMESSUSMP2Test(GenericMP2Test):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]
        
class GaussianMP2Test(GenericMP2Test):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]
        
class GaussianMP3Test(GenericMP3Test):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GaussianMP4SDTQTest(GenericMP4SDTQTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GaussianMP5Test(GenericMP5Test):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class Jaguar65LMP2Test(GenericMP2Test):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class PCGAMESSMP2Test(GenericMP2Test):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class PCGAMESSMP3Test(GenericMP3Test):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class PCGAMESSMP4SDQTest(GenericMP4SDQTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class PCGAMESSMP4SDTQTest(GenericMP4SDTQTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]


# Load test data using information in file.
testdata = gettestdata(module="MP")
              
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

    print "\n\n********* SUMMARY OF MP TEST **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
