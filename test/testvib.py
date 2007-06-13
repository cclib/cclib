__revision__ = "$Revision$"

import os, unittest

import bettertest
from testall import gettestdata


class GenericIRTest(bettertest.TestCase):
    """Vibrational frequency calculations."""

    def testvibdisps(self):
        """Are the dimensions of vibdisps consistent with 3N-6 x N x 3"""
        numvib = 3*len(self.data.atomnos) - 6
        self.assertEqual(self.data.vibdisps.shape,
                        (numvib, len(self.data.atomnos), 3))

    def testlengths(self):
        """Are the lengths of vibfreqs and vibirs correct?"""
        numvib = 3*len(self.data.atomnos) - 6
        self.assertEqual(len(self.data.vibfreqs), numvib)
        self.assertEqual(len(self.data.vibirs), numvib)

    def testfreqval(self):
        """Is the highest freq value 3630 +/- 200 cm-1?"""
        self.assertInside(self.data.vibfreqs[-1], 3630, 200)

    def testirintens(self):
        """Is the maximum IR intensity 100 +/- 10 km mol-1?"""
        self.assertInside(max(self.data.vibirs), 100, 10)

##    def testmaxvibdisps(self):
##        """What is the maximum value of displacement for a H vs a C?"""
##        Cvibdisps = compress(self.data.atomnos==6, self.data.vibdisps, 1)
##        Hvibdisps = compress(self.data.atomnos==1, self.data.vibdisps, 1)
##        self.assertEqual(max(abs(Cvibdisps).flat), 1.0)
        

class ADFIRTest(GenericIRTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]
    
class GamessUKIRTest(GenericIRTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GamessUSIRTest(GenericIRTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GaussianIRTest(GenericIRTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

    def testvibsyms(self):
        """Is the length of vibsyms correct?"""
        numvib = 3*len(self.data.atomnos) - 6        
        self.assertEqual(len(self.data.vibsyms), numvib)
       
class Jaguar42IRTest(GenericIRTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

    def testvibsyms(self):
            """Is the length of vibsyms correct?"""
            numvib = 3*len(self.data.atomnos) - 6        
            self.assertEqual(len(self.data.vibsyms), numvib)

class Jaguar65IRTest(GenericIRTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

    def testvibsyms(self):
            """Is the length of vibsyms correct?"""
            numvib = 3*len(self.data.atomnos) - 6        
            self.assertEqual(len(self.data.vibsyms), numvib)

class PCGamessIRTest(GenericIRTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

    def testirintens(self):
        """Is the maximum IR intensity 135 +/- 5 km mol-1?"""
        self.assertInside(max(self.data.vibirs), 135, 5)     

class GenericRamanTest(bettertest.TestCase):
    """Raman calculations."""

    def testlengths(self):
        """Is the length of vibramans correct?"""
        numvib = 3*len(self.data.atomnos) - 6
        self.assertEqual(len(self.data.vibramans), numvib)

    def testramanintens(self):
        """Is the maximum Raman intensity 575 +/- 5 A**4/amu?"""
        self.assertInside(max(self.data.vibramans), 575, 5)

class GamessUKRamanTest(GenericRamanTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

class GaussianRamanTest(GenericRamanTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]

    def testramanintens(self):
        """Is the maximum Raman intensity 1066 +/- 5 A**4/amu?"""
        self.assertInside(max(self.data.vibramans), 1066, 5)

class PCGamessRamanTest(GenericRamanTest):
    def setUp(self):
        self.data = testdata[self.__class__.__name__]["data"]


# Load test data using information in file.
testdata = gettestdata(module="vib")

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

    print "\n\n********* SUMMARY OF VIB TEST **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
