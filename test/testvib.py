import os, unittest
from Numeric import array
from testall import getfile
from cclib.parser import ADF, GAMESS, Gaussian, Jaguar, GAMESSUK

class GenericVibTest(unittest.TestCase):
    """Vibrational frequency calculations."""

    def testvibcarts(self):
        """Are the dimensions of vibcarts consistent with 3N-6 x N x 3"""
        numvib = 3*len(self.data.atomnos) - 6
        self.assertEqual(self.data.vibcarts.shape,
                        (numvib, len(self.data.atomnos), 3))

class GaussianVibTest(GenericVibTest):
    def setUp(self):
        self.data = data[0]

class GamessUSVibTest(GenericVibTest):
    def setUp(self):
        self.data = data[1]

class PCGamessVibTest(GenericVibTest):
    def setUp(self):
        self.data = data[2]

class ADFVibTest(GenericVibTest):
    def setUp(self):
        self.data = data[3]
    
class Jaguar42VibTest(GenericVibTest):
    def setUp(self):
        self.data = data[4]

class Jaguar65VibTest(GenericVibTest):
    def setUp(self):
        self.data = data[4]

class GamessUKVibTest(GenericVibTest):
    def setUp(self):
        self.data = data[6]

names = [ "Gaussian", "PCGamess", "GAMESS", "ADF", "Jaguar 4.2",
          "Jaguar 6.5", "GAMESS UK"]
tests = [ GaussianVibTest, PCGamessVibTest,
          GamessUSVibTest, ADFVibTest,
          Jaguar42VibTest, Jaguar65VibTest,
          GamessUKVibTest]
data = [getfile(Gaussian, "basicGaussian03","dvb_ir.out"),
        getfile(GAMESS, "basicGAMESS-US","dvb_ir.out"),
        getfile(GAMESS, "basicPCGAMESS","dvb_ir.out"),
        getfile(ADF, "basicADF2004.01","dvb_ir.adfout"),
        getfile(Jaguar, "basicJaguar4.2", "dvb_ir.out"),
        getfile(Jaguar, "basicJaguar6.5", "dvb_ir.out"),
        getfile(GAMESSUK, "basicGAMESS-UK", "dvb_ir.out")]
              
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
