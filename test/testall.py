import os, unittest
from cclib.parser import GAMESS,G03,ADF,Jaguar
from Numeric import array

class GenericSPTest(unittest.TestCase):
    """Restricted single point calculations with MO coeffs and overlap info."""
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        self.assertEquals(self.data.aooverlaps.shape,(self.data.nbasis,self.data.nbasis))

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nindep x nbasis?"""
        self.assertEquals(self.data.mocoeffs.shape,(1,self.data.nindep,self.data.nbasis))

class GaussianSPTest(GenericSPTest):
    def setUp(self):
        self.data = getfile(G03,"basicGaussian03","dvb_sp.out")

class GamessUSSPTest(GenericSPTest):
    def setUp(self):
        self.data = getfile(GAMESS,"basicGAMESS-US","dvb_sp.out")

class PCGamessSPTest(GenericSPTest):
    def setUp(self):
        self.data = getfile(GAMESS,"basicPCGAMESS","dvb_sp.out")

class ADFSPTest(GenericSPTest):
    def setUp(self):
        self.data = getfile(ADF,"basicADF2004.01","dvb_sp_b.adfout")

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
        self.data = getfile(G03,"basicGaussian03","dvb_gopt.out")

class GamessUSGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = getfile(GAMESS,"basicGAMESS-US","dvb_gopt_a.out")

class PCGamessGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = getfile(GAMESS,"basicPCGAMESS","dvb_gopt_a.out")

class ADFGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = getfile(ADF,"basicADF2004.01","dvb_gopt.adfout")

    def testscfvaluedim(self):
        """Do the scf values have the right dimensions? 
           ADF calculations one more SCF cycle after the geometry is converged"""
        self.assert_(len(self.data.scfvalues)==len(self.data.geovalues)+1 and len(self.data.scfvalues[0])==len(self.data.scftargets))

class JaguarGeoOptTest(GenericGeoOptTest):
    def setUp(self):
        self.data = getfile(Jaguar,"basicJaguar","eg01","dvb_gopt.out")

def getfile(parser,*location):
    """Returns a parsed logfile."""
    if parser.__name__ in ['GAMESS','ADF','Jaguar']:
        fullpath = ("..","data",parser.__name__) + location
    elif parser.__name__=="G03":
        fullpath = ("..","data","Gaussian") + location
    logfile = parser(os.path.join(*fullpath))
    logfile.logger.setLevel(0)
    logfile.parse()
    return logfile

def visualtests():
    """These are not formal tests -- but they should be eyeballed."""
    logfiles = [ getfile(G03,"basicGaussian03","dvb_gopt.out"),
                 getfile(GAMESS,"basicPCGAMESS","dvb_gopt_a.out"),
                 getfile(GAMESS,"basicGAMESS-US","dvb_gopt_a.out"),
                 getfile(ADF,"basicADF2004.01","dvb_gopt.adfout"),
                 getfile(Jaguar,"basicJaguar","eg01","dvb_gopt.out")]

    print "\n\nMO energies of optimised dvb"
    print "    ","".join(["%8s" % x for x in ['Gaussian','PCGAMESS','GAMESS-US','ADF','Jaguar']])
    print "HOMO", "  ".join(["%+2.4f" % x.moenergies[0,x.homos[0]] for x in logfiles])
    print "LUMO", "  ".join(["%+2.4f" % x.moenergies[0,x.homos[0]+1] for x in logfiles])
    print "H-L  ", "  ".join(["%2.4f" % (x.moenergies[0,x.homos[0]+1]-x.moenergies[0,x.homos[0]],) for x in logfiles])

    
    
if __name__=="__main__":
    names = [ "Gaussian", "PCGamess", "GAMESS", "ADF", "Jaguar" ]
    tests = [ GaussianGeoOptTest, PCGamessGeoOptTest,
              GamessUSGeoOptTest, ADFGeoOptTest,
              JaguarGeoOptTest ]
    total = errors = failures = 0
    for name,test in zip(names,tests):
        print "\n**** Testing %s Geo Opt ****" % name
        myunittest = unittest.makeSuite(test)
        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)

    tests = [ GaussianSPTest, PCGamessSPTest,
              GamessUSSPTest, ADFSPTest ]
    for name,test in zip(names,tests):
        print "\n**** Testing %s SP ****" % name
        myunittest = unittest.makeSuite(test)
        a = unittest.TextTestRunner(verbosity=2).run(myunittest)
        total += a.testsRun
        errors += len(a.errors)
        failures += len(a.failures)

    print "\n\n********* SUMMARY OF EVERYTHING **************"
    print "TOTAL: %d\tPASSED: %d\tFAILED: %d\tERRORS: %d" % (total,total-(errors+failures),failures,errors)

    print "\n\n*** Visual tests ***"
    visualtests()
