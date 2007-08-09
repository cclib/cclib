__revision__ = "$Revision$"

import bettertest


class GenericSPTest(bettertest.TestCase):
    """Restricted single point unittest."""

    nbasisdict = {1:1, 6:5} # STO-3G, H has 1, C has 3.

    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        self.assertEquals(self.data.aooverlaps.shape,(self.data.nbasis,self.data.nbasis))

    def testatomcoords(self):
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        self.assertEquals(self.data.atomcoords.shape,(1,self.data.natom,3))
    
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        self.assertEquals(type(self.data.mocoeffs), type([]))
        self.assertEquals(len(self.data.mocoeffs), 1)
        self.assertEquals(self.data.mocoeffs[0].shape,
                          (self.data.nmo, self.data.nbasis))

    def testaooverlaps(self):
        """Are the first row and colm of the overlap matrix identical?"""
        self.assertEquals(sum(self.data.aooverlaps[0,:] -
                              self.data.aooverlaps[:,0]),
                          0)

    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique?"""
        all = []
        for i,atom in enumerate(self.data.atombasis):
            self.assertEquals(len(atom), self.nbasisdict[self.data.atomnos[i]])
            all += atom
        # Test if there are as many indices as atomic orbitals.
        self.assertEquals(len(all), self.data.nbasis)
        # Check if all are different (every orbital indexed once).
        self.assertEquals(len(set(all)), len(all))

class ADFSPTest(GenericSPTest):
    """ADF restricted single point unittest."""

    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        # ADF uses fooverlaps.
        self.assertEquals(self.data.fooverlaps.shape,(self.data.nbasis,self.data.nbasis))

class GamessUKSPTest(GenericSPTest):
    """GAMESS-UK restricted single point unittest."""

class GamessUSSPTest(GenericSPTest):
    """GAMESS-US restricted single point unittest."""

class GaussianSPTest(GenericSPTest):
    """Gaussian restricted single point unittest."""

class Jaguar42SPTest(GenericSPTest):
    """Jaguar4.2 restricted single point unittest."""

    # Data file does not contain enough information. Can we make a new one?
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""
        self.assertEquals(1, 1)

    # Data file does not contain enough information. Can we make a new one?
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique? PASS"""
        self.assertEquals(1, 1)
       
class Jaguar60SPTest(GenericSPTest):
    """Jaguar6.0 restricted single point unittest."""

    nbasisdict = {1:5, 6:15} # 6-31G(d,p)

class MolproSPTest(GenericSPTest):
    """Molpro restricted single point unittest."""

class PCGamessSPTest(GenericSPTest):
    """PC-GAMESS restricted single point unittest."""


if __name__=="__main__":

    from testall import testmodule
    testmodule("SP")
