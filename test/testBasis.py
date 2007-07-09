__revision__ = "$Revision$"

import bettertest


class GenericBasisTest(bettertest.TestCase):
    """Basis set unittest."""

    def testgbasis(self):
        """Is gbasis the right length?"""
        self.assertEquals(self.data.natom, len(self.data.gbasis))
    
    def testnames(self):
        """Test the names of the basis sets."""
        for atom in self.data.gbasis:
            for fns in atom:
                self.assert_(fns[0] in ['S', 'P'],
                             "%s not one of S or P" % fns[0])

    def testsizeofbasis(self):
        """Test the basis set size."""
        total = 0
        multiple = {'S':1, 'P':3}
        for atom in self.data.gbasis:
            for fns in atom:
                total += multiple[fns[0]] # Add 3 for P
        self.assertEquals(self.data.nbasis, total)
    
    def testcoeffs(self):
        """Test the coeffs of the basis sets."""
        for atom in self.data.gbasis:
            if len(atom)==1: # i.e. a 'H'
                coeffs = atom[0][1]
                self.assertAlmostEqual(coeffs[0][0], 3.42525, 5)
                self.assertAlmostEqual(coeffs[0][1], 0.15433, 5)
            else: # i.e. a 'C'
                self.assertEquals(len(atom), 3)
                s_coeffs = atom[1][1]
                p_coeffs = atom[2][1]
                self.assertAlmostEqual(s_coeffs[0][0], 2.9412, 4)
                self.assertAlmostEqual(p_coeffs[0][0], 2.9412, 4)
                self.assertAlmostEqual(s_coeffs[0][1], -0.1000, 4)
                self.assertAlmostEqual(p_coeffs[0][1], 0.1559, 4)

class GamessUKBasisTest(GenericBasisTest):
    """GAMESS-UK basis set unittest."""

class GamessUSBasisTest(GenericBasisTest):
    """GAMESS-US basis set unittest."""

class GaussianBasisTest(GenericBasisTest):
    """Gaussian basis set unittest."""

class MolproBasisTest(GenericBasisTest):
    """Molpro basis set unittest."""

class PCGamessBasisTest(GenericBasisTest):
    """PC-GAMESS basis set unittest."""

              
if __name__=="__main__":

    from testall import testmodule
    testmodule("Basis")
