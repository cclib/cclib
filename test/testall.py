import os, unittest
from cclib.parser import GAMESS,G03,ADF
from Numeric import array

class GenericTest(unittest.TestCase):
    def testhomos(self):
        self.assertEquals(self.data.homos,array([34]))

    def testnatom(self):
        self.assertEquals(self.data.natom,20)

    def testnbasis(self):
        self.assertEquals(self.data.nbasis,60)

    def testscfenergy(self):
        self.assert_(self.data.scfenergies[-1]+382.3<3)

    def testhomoenergy(self):
        self.assert_(self.data.moenergies[0,self.data.homos[0]]+4.165<0.5,"HOMO energy is %f (for G03 it's -4.165)" % self.data.moenergies[0,self.data.homos[0]])

class GaussianTest(GenericTest):
    def setUp(self):
        self.data = getfile(G03,"basicGaussian03","dvb_gopt.out")

class GamessTest(GenericTest):
    def setUp(self):
        self.data = getfile(GAMESS,"basicPCGAMESS","dvb_gopt_a.out")

class ADFTest(GenericTest):
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
    pass
    
    
if __name__=="__main__":
    gaussiantests = unittest.makeSuite(GaussianTest)
    gamesstests = unittest.makeSuite(GamessTest)
    adftests = unittest.makeSuite(ADFTest)
    alltests = unittest.TestSuite((gaussiantests,gamesstests,adftests))
    unittest.TextTestRunner(verbosity=2).run(alltests)
