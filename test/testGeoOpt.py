import os, unittest
from Numeric import array
from testall import getfile
from cclib.parser import ADF, GAMESS, Gaussian, Jaguar

class GenericGeoOptTest(unittest.TestCase):
    def testhomos(self):
        """Is the index of the homo equal to 34?"""
        self.assertEquals(self.data.homos,array([34]))

    def testnatom(self):
        """Is the number of atoms equal to 20?"""
        self.assertEquals(self.data.natom,20)

    def testnbasis(self):
        """Is the number of basis set functions equal to 60?"""
        self.assertEquals(self.data.nbasis,60)

    def testscfenergy(self):
        """Is the SCF energy within 3eV(?) of -382.3?"""
        self.assert_(self.data.scfenergies[-1]+382.3<3)

    def testnormalisesym(self):
        """Did this subclasses overwrite normalisesym?"""
        self.assertNotEquals(self.data.normalisesym("A"),"ERROR: This should be overwritten by this subclass")

    def testlengthmoenergies(self):
        """Is the number of evalues equal to 60?"""
        self.assertEquals(60,len(self.data.moenergies[0]))

    def testsymlabels(self):
        """Are all the symmetry labels either Ag/u or Bg/u?"""
        sumwronglabels = sum([x not in ['Ag','Bu','Au','Bg'] for x in self.data.mosyms[0]])
        self.assertEquals(sumwronglabels,0)

    def testscfvaluetype(self):
        """Do the scf values have the right type?"""
        self.assert_(type(self.data.scfvalues[0])==type(array([])) and type(self.data.scfvalues)==type([]))

    def testscfvaluedim(self):
        """Do the scf values have the right dimensions?"""
        self.assert_(len(self.data.scfvalues)==len(self.data.geovalues) and len(self.data.scfvalues[0])==len(self.data.scftargets))

class GaussianGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = getfile(Gaussian,"basicGaussian03","dvb_gopt.out")

class GamessUSGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = getfile(GAMESS,"basicGAMESS-US","dvb_gopt_a.out")

class PCGamessGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = getfile(GAMESS,"basicPCGAMESS","dvb_gopt_a.out")

class ADFGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = getfile(ADF,"basicADF2004.01","dvb_gopt_b.adfout")

    def testscfvaluedim(self):
        """Do the scf values have the right dimensions? 
           ADF calculations one more SCF cycle after the geometry is converged"""
        self.assert_(len(self.data.scfvalues)==len(self.data.geovalues)+1 and len(self.data.scfvalues[0])==len(self.data.scftargets))

class JaguarGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = getfile(Jaguar,"basicJaguar","eg01","dvb_gopt.out")

names = [ "Gaussian", "PCGamess", "GAMESS", "ADF", "Jaguar" ]
tests = [ GaussianGeoOptTest, PCGamessGeoOptTest,
          GamessUSGeoOptTest, ADFGeoOptTest,
          JaguarGeoOptTest ]

if __name__=="__main__":
    total = errors = failures = 0
    for name,test in zip(names,tests):
        print "\n**** Testing %s Geo Opt ****" % name
        myunittest = unittest.makeSuite(test)
        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)
    
    print "\n\n********* SUMMARY OF Geo Opt **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)
