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
        self.data = G03(os.path.join("..","data","Gaussian","basicGaussian03","dvb_gopt.out"))
        self.data.logger.setLevel(0)
        self.data.parse()

class GamessTest(GenericTest):
    def setUp(self):
        self.data = GAMESS(os.path.join("..","data","GAMESS","basicPCGAMESS","dvb_gopt_a.out"))
        self.data.logger.setLevel(0)
        self.data.parse()

class ADFTest(GenericTest):
    def setUp(self):
        self.data = ADF(os.path.join("..","data","ADF","basicADF2004.01","dvb_gopt.adfout"))
        self.data.logger.setLevel(0)
        self.data.parse()
    
if __name__=="__main__":
    gaussiantests = unittest.makeSuite(GaussianTest)
    gamesstests = unittest.makeSuite(GamessTest)
    adftests = unittest.makeSuite(ADFTest)
    alltests = unittest.TestSuite((gaussiantests,gamesstests,adftests))
    unittest.TextTestRunner(verbosity=2).run(alltests)
