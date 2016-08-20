# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006,2007,2012,2014,2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test logfiles related to basis sets"""

import os
import unittest


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericBasisTest(unittest.TestCase):
    """Generic basis set unittest"""

    # The number of contraction per atom, by atom number.
    contractions = { 1: 1, 6: 3 }

    # Number of components in each contraction by subshell type,
    # so that we can infer nbasis from gbasis. Note how we assume
    # the basis set is not is spherical representation.
    names = ['S', 'P', 'D', 'F', 'G']
    multiple = {'S': 1, 'P': 3, 'D': 6, 'F': 10, 'G': 15}
    multiple_spher = {'S': 1, 'P': 3, 'D': 5, 'F': 7, 'G': 9}
    spherical = False

    # These are the expected exponents and coefficients for the first
    # Gaussians in particular shells for hydrogen and carbon atoms.
    gbasis_H_1s_func0 = [3.42525, 0.15433]
    gbasis_C_2s_func0 = [2.9412, -0.1000]
    gbasis_C_2p_func0 = [2.9412, 0.1559]

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
        multiple = self.multiple_spher if self.spherical else self.multiple
        for atom in self.data.gbasis:
            for (ftype, contraction) in atom:
                total += multiple[ftype]

        self.assertEquals(self.data.nbasis, total)

    def testcontractions(self):
        """Are the number of contractions on all atoms correct?"""
        for iatom, atom in enumerate(self.data.gbasis):
            atomno = self.data.atomnos[iatom]
            self.assertEquals(len(atom), self.contractions[atomno])

    def testprimitives(self):
        """Are all primitives 2-tuples?"""
        for atom in self.data.gbasis:
            for ftype, contraction in atom:
                for primitive in contraction:
                    self.assertEquals(len(primitive), 2)

    def testcoeffs(self):
        """Are the atomic basis set exponents and coefficients correct?"""

        for iatom,atom in enumerate(self.data.gbasis):
            if self.data.atomnos[iatom] == 1:
                coeffs = atom[0][1]
                self.assertAlmostEqual(coeffs[0][0], self.gbasis_H_1s_func0[0], 4)
                self.assertAlmostEqual(coeffs[0][1], self.gbasis_H_1s_func0[1], 4)
            else:
                s_coeffs = atom[1][1]
                p_coeffs = atom[2][1]
                self.assertAlmostEqual(s_coeffs[0][0], self.gbasis_C_2s_func0[0], 4)
                self.assertAlmostEqual(p_coeffs[0][0], self.gbasis_C_2p_func0[0], 4)
                self.assertAlmostEqual(s_coeffs[0][1], self.gbasis_C_2s_func0[1], 4)
                self.assertAlmostEqual(p_coeffs[0][1], self.gbasis_C_2p_func0[1], 4)


class JaguarBasisTest(GenericBasisTest):
    """Customized basis set unittest"""

    # For some reason, Jaguar seems to use slightly different coefficients for
    # contractions in the STO-3G basis set. Or perhaps we don't understand something.
    gbasis_H_1s_func0 = [3.42525, 0.24050]
    gbasis_C_2s_func0 = [2.941249, -0.29565]
    gbasis_C_2p_func0 = [2.941249, 0.22135]


class GenericBigBasisTest(GenericBasisTest):
    """Generic big basis set unittest"""

    contractions = { 6: 20 }

    @unittest.skip('Write up a new test, and/or revise the one inherited.')
    def testcoeffs(self):
        """Are the basis set coefficients correct?"""
        self.assertEqual(1, 1)

    @unittest.skip('# of contractions is 20 for VQZ, but 29 for CVQZ; unify files first.')
    def testcontractions(self):
        """"""
        self.assertEqual(1, 1)


class DALTONBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""
    spherical = True


class GaussianBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""
    spherical = True


class JaguarBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""

    spherical = True

    # Jaguar only goes up to F functions.
    names = ['S', 'P', 'D', 'F']


class MolproBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""
    spherical = True


class PsiBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""
    spherical = True


class QChemBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""
    spherical = True


if __name__=="__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['Basis'])
    suite.testall()
