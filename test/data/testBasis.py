# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test logfiles related to basis sets"""

import pytest


class GenericBasisTest:
    """Generic basis set unittest"""

    # The number of contraction per atom, by atom number.
    contractions = {1: 1, 6: 3}

    # Number of components in each contraction by subshell type,
    # so that we can infer nbasis from gbasis. Note how we assume
    # the basis set is not is spherical representation.
    names = ["S", "P", "D", "F", "G"]
    multiple = {"S": 1, "P": 3, "D": 6, "F": 10, "G": 15}
    multiple_spher = {"S": 1, "P": 3, "D": 5, "F": 7, "G": 9}
    spherical = False

    # These are the expected exponents and coefficients for the first
    # Gaussians in particular shells for hydrogen and carbon atoms.
    gbasis_H_1s_func0 = [3.42525, 0.15433]
    gbasis_C_2s_func0 = [2.9412, -0.1000]
    gbasis_C_2p_func0 = [2.9412, 0.1559]

    def testgbasis(self, data) -> None:
        """Is gbasis the right length?"""
        assert data.natom == len(data.gbasis)

    def testnames(self, data) -> None:
        """Are the name of basis set functions acceptable?"""
        for atom in data.gbasis:
            for fns in atom:
                assert fns[0] in self.names, f"{fns[0]} not one of S or P"

    def testsizeofbasis(self, data) -> None:
        """Is the basis set the correct size?"""

        total = 0
        multiple = self.multiple_spher if self.spherical else self.multiple
        for atom in data.gbasis:
            for ftype, contraction in atom:
                total += multiple[ftype]

        assert data.nbasis == total

    def testcontractions(self, data) -> None:
        """Are the number of contractions on all atoms correct?"""
        for iatom, atom in enumerate(data.gbasis):
            atomno = data.atomnos[iatom]
            assert len(atom) == self.contractions[atomno]

    def testprimitives(self, data) -> None:
        """Are all primitives 2-tuples?"""
        for atom in data.gbasis:
            for ftype, contraction in atom:
                for primitive in contraction:
                    assert len(primitive) == 2

    def testcoeffs(self, data) -> None:
        """Are the atomic basis set exponents and coefficients correct?"""

        for iatom, atom in enumerate(data.gbasis):
            if data.atomnos[iatom] == 1:
                coeffs = atom[0][1]
                assert round(abs(coeffs[0][0] - self.gbasis_H_1s_func0[0]), 4) == 0
                assert round(abs(coeffs[0][1] - self.gbasis_H_1s_func0[1]), 4) == 0
            else:
                s_coeffs = atom[1][1]
                p_coeffs = atom[2][1]
                assert round(abs(s_coeffs[0][0] - self.gbasis_C_2s_func0[0]), 4) == 0
                assert round(abs(p_coeffs[0][0] - self.gbasis_C_2p_func0[0]), 4) == 0
                assert round(abs(s_coeffs[0][1] - self.gbasis_C_2s_func0[1]), 4) == 0
                assert round(abs(p_coeffs[0][1] - self.gbasis_C_2p_func0[1]), 4) == 0

    def testatomcoords(self, data) -> None:
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        expected_shape = (1, data.natom, 3)
        assert data.atomcoords.shape == expected_shape


class JaguarBasisTest(GenericBasisTest):
    """Customized basis set unittest"""

    # For some reason, Jaguar seems to use slightly different coefficients for
    # contractions in the STO-3G basis set. Or perhaps we don't understand something.
    gbasis_H_1s_func0 = [3.42525, 0.24050]
    gbasis_C_2s_func0 = [2.941249, -0.29565]
    gbasis_C_2p_func0 = [2.941249, 0.22135]


class GenericBigBasisTest(GenericBasisTest):
    """Generic big basis set unittest"""

    contractions = {6: 20}

    @pytest.mark.skip("Write up a new test, and/or revise the one inherited.")
    def testcoeffs(self, data):
        """Are the basis set coefficients correct?"""
        assert True

    @pytest.mark.skip("# of contractions is 20 for VQZ, but 29 for CVQZ; unify files first.")
    def testcontractions(self, data):
        """"""
        assert True


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
    names = ["S", "P", "D", "F"]


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
