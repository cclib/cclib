# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test coupled cluster logfiles"""

import numpy
import pytest


class GenericCCTest:
    """Generic coupled cluster unittest"""

    rel_thresh = 1.0e-4

    def testsizeandshape(self, data) -> None:
        """Are the dimensions of ccenergies correct?"""
        assert data.ccenergies.shape == (len(data.scfenergies),)

    def testsign(self, data) -> None:
        """Are the coupled cluster corrections negative?"""
        corrections = data.ccenergies - data.scfenergies
        assert numpy.all(corrections < 0.0)


class GenericCC2Test(GenericCCTest):
    # Turbomole 7.4
    corr_energy = -0.0422913114

    def testenergycc2(self, data) -> None:
        """Is the CC2 correlation energy within the target?"""
        e_scf = data.scfenergies[0]
        e_cc = data.ccenergies[0]
        e_corr = e_cc - e_scf
        assert pytest.approx(e_corr, rel=self.rel_thresh) == self.corr_energy


class GenericCCDTest(GenericCCTest):
    # Q-Chem 5.4
    corr_energy = -0.05304913

    def testenergyccd(self, data) -> None:
        """Is the CCD correlation energy within the target?"""
        e_scf = data.scfenergies[0]
        e_cc = data.ccenergies[0]
        e_corr = e_cc - e_scf
        assert pytest.approx(e_corr, rel=self.rel_thresh) == self.corr_energy


class GenericCCSDTest(GenericCCTest):
    # Q-Chem 5.4
    corr_energy = -0.05335475

    def testenergyccsd(self, data) -> None:
        """Is the CCSD correlation energy within the target?"""
        e_scf = data.scfenergies[0]
        e_cc = data.ccenergies[0]
        e_corr = e_cc - e_scf
        assert pytest.approx(e_corr, rel=self.rel_thresh) == self.corr_energy


class GenericCCSDPTTest(GenericCCTest):
    # Q-Chem 5.4
    corr_energy = -0.05335475 + -0.00007679

    def testenergyccsdpt(self, data) -> None:
        """Is the CCSD(T) correlation energy within the target?"""
        e_scf = data.scfenergies[0]
        e_cc = data.ccenergies[0]
        e_corr = e_cc - e_scf
        assert pytest.approx(e_corr, rel=self.rel_thresh) == self.corr_energy


class DALTONCCSDPTTest(GenericCCSDPTTest):
    # DALTON 2015
    corr_energy = -73.4491137256 - -73.4345879538


class GamessCCDTest(GenericCCDTest):
    # GAMESS-US 2018
    corr_energy = -0.2070122255


class GamessCCSDTest(GenericCCSDTest):
    # GAMESS-US 2018
    corr_energy = -0.2078353755


class GamessCCSDPTTest(GenericCCSDPTTest):
    # GAMESS-US 2018
    corr_energy = -0.2107439309


class MolcasCCSDTest(GenericCCSDTest):
    # OpenMolcas 18.0
    corr_energy = -0.0507029554992


class MolproCCDTest(GenericCCDTest):
    # Molpro 2012
    corr_energy = -0.212003853457


class MolproCCSDTest(GenericCCSDTest):
    # Molpro 2012
    corr_energy = -0.212805593281


class MolproCCSDPTTest(GenericCCSDPTTest):
    # Molpro 2012
    corr_energy = -0.215948162128


class NWChemCCSDPTTest(GenericCCSDPTTest):
    # NWChem 7.0
    corr_energy = -0.053506090792006


class OrcaCCSDTest(GenericCCSDTest):
    # ORCA 5.0
    corr_energy = -0.049913571


class OrcaCCSDPTTest(GenericCCSDPTTest):
    # ORCA 5.0
    corr_energy = -0.049982067


class Psi4CCSDTest(GenericCCSDTest):
    # Psi4 1.3.1
    corr_energy = -0.053429410909519


class Psi4CCSDPTTest(GenericCCSDPTTest):
    # Psi4 1.3.1
    corr_energy = -0.053429410909519 + -0.000076703052954


class TurbomoleCCSDTest(GenericCCSDTest):
    # TODO The two Turbomole CCSD(T) energies aren't similar...
    rel_thresh = 5.0e-3

    # Turbomole 7.4
    corr_energy = -0.0510005307
