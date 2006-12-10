import os, unittest
from Numeric import array, compress
from testall import getfile
from cclib.parser import ADF, GAMESS, Gaussian, Jaguar, GAMESSUK
import bettertest

class GenericVibTest(bettertest.TestCase):
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
        

class GaussianVibTest(GenericVibTest):
    def setUp(self):
        self.data = data[0]

    def testvibsyms(self):
        """Is the length of vibsyms correct?"""
        numvib = 3*len(self.data.atomnos) - 6        
        self.assertEqual(len(self.data.vibsyms), numvib)
       
class GamessUSVibTest(GenericVibTest):
    def setUp(self):
        self.data = data[1]

class PCGamessVibTest(GenericVibTest):
    def setUp(self):
        self.data = data[2]

    def testirintens(self):
        """Is the maximum IR intensity 135 +/- 5 km mol-1?"""
        self.assertInside(max(self.data.vibirs), 135, 5)     

class ADFVibTest(GenericVibTest):
    def setUp(self):
        self.data = data[3]
    
class Jaguar42VibTest(GenericVibTest):
    def setUp(self):
        self.data = data[4]

    def testvibsyms(self):
            """Is the length of vibsyms correct?"""
            numvib = 3*len(self.data.atomnos) - 6        
            self.assertEqual(len(self.data.vibsyms), numvib)

class Jaguar65VibTest(GenericVibTest):
    def setUp(self):
        self.data = data[4]

    def testvibsyms(self):
            """Is the length of vibsyms correct?"""
            numvib = 3*len(self.data.atomnos) - 6        
            self.assertEqual(len(self.data.vibsyms), numvib)

class GamessUKVibTest(GenericVibTest):
    def setUp(self):
        self.data = data[6]

class GenericRamanTest(bettertest.TestCase):
    """Raman calculations."""

    def testlengths(self):
        """Is the length of vibramans correct?"""
        numvib = 3*len(self.data.atomnos) - 6
        self.assertEqual(len(self.data.vibramans), numvib)

    def testramanintens(self):
        """Is the maximum Raman intensity 575 +/- 5 A**4/amu?"""
        self.assertInside(max(self.data.vibramans), 575, 5)

class GaussianRamanTest(GenericRamanTest):
    def setUp(self):
        self.data = data[7]

    def testramanintens(self):
        """Is the maximum Raman intensity 1066 +/- 5 A**4/amu?"""
        self.assertInside(max(self.data.vibramans), 1066, 5)

class GamessUKRamanTest(GenericRamanTest):
    def setUp(self):
        self.data = data[8]

class PCGamessRamanTest(GenericRamanTest):
    def setUp(self):
        self.data = data[9]

names = [ "Gaussian", "PCGamess", "GAMESS", "ADF", "Jaguar 4.2",
          "Jaguar 6.5", "GAMESS UK",
          "Gaussian", "GAMESS UK", "PCGamess"]
tests = [ GaussianVibTest, PCGamessVibTest,
          GamessUSVibTest, ADFVibTest,
          Jaguar42VibTest, Jaguar65VibTest,
          GamessUKVibTest,
          GaussianRamanTest, GamessUKRamanTest,
          PCGamessRamanTest ]
data = [getfile(Gaussian, "basicGaussian03","dvb_ir.out"),
        getfile(GAMESS, "basicGAMESS-US","dvb_ir.out"),
        getfile(GAMESS, "basicPCGAMESS","dvb_ir.out"),
        getfile(ADF, "basicADF2004.01","dvb_ir.adfout"),
        getfile(Jaguar, "basicJaguar4.2", "dvb_ir.out"),
        getfile(Jaguar, "basicJaguar6.5", "dvb_ir.out"),
        getfile(GAMESSUK, "basicGAMESS-UK", "dvb_ir.out"),
        getfile(Gaussian, "basicGaussian03", "dvb_raman.out"),
        getfile(GAMESSUK, "basicGAMESS-UK","dvb_raman.out"),
        getfile(GAMESS, "basicPCGAMESS","dvb_raman.out"),
        ]
              
if __name__=="__main__":
    total = errors = failures = 0

    for name,test in zip(names,tests):
        print "\n**** Testing %s vibrations ****" % name
        myunittest = unittest.makeSuite(test)
        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)

    print "\n\n********* SUMMARY OF SP **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
