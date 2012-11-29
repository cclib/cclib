# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

__revision__ = "$Revision$"

import bettertest


class GenericBasisTest(bettertest.TestCase):
    """Basis set unittest."""

    names = ['S', 'P', 'D', 'F', 'G']
    multiple = {'S':1, 'P':3, 'D':6, 'F':10, 'G':15}
    multiple_spher = {'S':1, 'P':3, 'D':5, 'F':7, 'G':9}
    spherical = False

    def testgbasis(self):
        """Is gbasis the right length?"""
        self.assertEquals(self.data.natom, len(self.data.gbasis))
    
    def testnames(self):
        """Are the name of basis set functions acceptable?"""
        for atom in self.data.gbasis:
            for fns in atom:
                self.assert_(fns[0] in self.names,
                             "%s not one of S or P" % fns[0])

    def testsizeofbasis(self):
        """Is the basis set the correct size?"""
        total = 0
        multiple = self.multiple
        if self.spherical:
            multiple = self.multiple_spher
        for atom in self.data.gbasis:
            for fns in atom:
                 # Add 3 for P, 5 or 6 for D, and so forth.
                total += multiple[fns[0]]
        self.assertEquals(self.data.nbasis, total)
    
    def testcoeffs(self):
        """Are the basis set coefficients correct?"""
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


class GenericBigBasisTest(GenericBasisTest):
    """Big basis set unittest."""
    
    # Write up a new test, and/or revise the one inherited.
    def testcoeffs(self):
        """Are the basis set coefficients correct? PASS"""
        self.assertEqual(1, 1)


class GamessUKBasisTest(GenericBasisTest):
    """GAMESS-UK basis set unittest."""


class GamessUSBasisTest(GenericBasisTest):
    """GAMESS-US basis set unittest."""


class GamessUSBigBasisTest(GenericBigBasisTest):
    """GAMESS-US big basis set unittest."""

    old_tests = ["GAMESS/GAMESS-US/MoOCl4-sp_2005.06.27.r3.out.bz2"]

class GaussianBasisTest(GenericBasisTest):
    """Gaussian basis set unittest."""


class GaussianBigBasisTest(GenericBigBasisTest):
    """Gaussian big basis set unittest."""
    
    spherical = True


class MolproBasisTest(GenericBasisTest):
    """Molpro basis set unittest."""


class MolproBigBasisTest(GenericBigBasisTest):
    """Molpro big basis set unittest."""

    spherical = True


class PCGamessBasisTest(GenericBasisTest):
    """PC-GAMESS basis set unittest."""


class PCGamessBigBasisTest(GenericBigBasisTest):
    """PC-GAMESS big basis set unittest."""

              
if __name__=="__main__":

    from testall import testall
    testall(modules=["Basis"])
