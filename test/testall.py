import os, unittest
from cclib.parser import GAMESS,G03,ADF
from Numeric import array

class GenericGeoOptTest(unittest.TestCase):
    def testhomos(self):
        """Is the index of the homo equal to 34?"""
        self.assertEquals(self.data.homos,array([34]))

    def testnatom(self):
        """Is the number of atoms equal to 20?"""
        self.assertEquals(self.data.natom,20)

    def testnbasis(self):
        """Is the number of basis set function equal to 60?"""
        self.assertEquals(self.data.nbasis,60)

    def testscfenergy(self):
        """Is the SCF energy within 3eV(?) of -382.3?"""
        self.assert_(self.data.scfenergies[-1]+382.3<3)

    def testnormalisesym(self):
        """Did this subclasses overwrite normalisesym?"""
        self.assertNotEquals(self.data.normalisesym("A"),"ERROR: This should be overwritten by this subclass")

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

def getfile(parser,*location):
    """Returns a parsed logfile."""
    if parser.__name__ in ['GAMESS','ADF']:
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
                 getfile(ADF,"basicADF2004.01","dvb_gopt.adfout") ]

    print "\n\nMO energies of optimised dvb"
    print "    ","".join(["%8s" % x for x in ['Gaussian','PCGAMESS','GAMESS-US','ADF']])
    print "HOMO", "  ".join(["%+2.4f" % x.moenergies[0,x.homos[0]] for x in logfiles])
    print "LUMO", "  ".join(["%+2.4f" % x.moenergies[0,x.homos[0]+1] for x in logfiles])
    print "H-L  ", "  ".join(["%2.4f" % (x.moenergies[0,x.homos[0]+1]-x.moenergies[0,x.homos[0]],) for x in logfiles])

    for x in logfiles:
        print x.mosyms
    
    
if __name__=="__main__":
    gaussiantests = unittest.makeSuite(GaussianGeoOptTest)
    pcgamesstests = unittest.makeSuite(PCGamessGeoOptTest)
    gamessustests = unittest.makeSuite(GamessUSGeoOptTest)
    adftests = unittest.makeSuite(ADFGeoOptTest)
    print "\n*** Testing Gaussian dvb_gopt.out ***"
    unittest.TextTestRunner(verbosity=2).run(gaussiantests)
    print "\n\n*** Testing PCGAMESS dvb_gopt_a.out ***"    
    unittest.TextTestRunner(verbosity=2).run(pcgamesstests)
    print "\n\n*** Testing GAMESS-US dvb_gopt_a.out ***"    
    unittest.TextTestRunner(verbosity=2).run(gamessustests)
    print "\n\n*** Testing ADF dvb_gopt.adfout ***"
    unittest.TextTestRunner(verbosity=2).run(adftests)
    print "\n\n*** Visual tests ***"
    visualtests()
