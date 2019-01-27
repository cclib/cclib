# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test logfiles related to basis sets"""

import os
import unittest

from skip import skipForParser

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

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testgbasis(self):
        """Is gbasis the right length?"""
        self.assertEqual(self.data.natom, len(self.data.gbasis))

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testnames(self):
        """Are the name of basis set functions acceptable?"""
        for atom in self.data.gbasis:
            for fns in atom:
                self.assertTrue(fns[0] in self.names,
                             "%s not one of S or P" % fns[0])

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testsizeofbasis(self):
        """Is the basis set the correct size?"""

        total = 0
        multiple = self.multiple_spher if self.spherical else self.multiple
        for atom in self.data.gbasis:
            for (ftype, contraction) in atom:
                total += multiple[ftype]

        self.assertEqual(self.data.nbasis, total)

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testcontractions(self):
        """Are the number of contractions on all atoms correct?"""
        for iatom, atom in enumerate(self.data.gbasis):
            atomno = self.data.atomnos[iatom]
            self.assertEqual(len(atom), self.contractions[atomno])

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testprimitives(self):
        """Are all primitives 2-tuples?"""
        for atom in self.data.gbasis:
            for ftype, contraction in atom:
                for primitive in contraction:
                    self.assertEqual(len(primitive), 2)

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
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


class MolcasBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""
    spherical = True


class MolproBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""
    spherical = True


class Psi4BigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""
    spherical = True


class QChemBigBasisTest(GenericBigBasisTest):
    """Customized big basis set unittest"""
    spherical = True


if __name__=="__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['Basis'])
    suite.testall()
