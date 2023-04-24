# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test coupled cluster logfiles"""

import os
import unittest

import numpy
import pytest

__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericCCTest(unittest.TestCase):
    """Generic coupled cluster unittest"""

    rel_thresh = 1.0e-4

    def testsizeandshape(self):
        """Are the dimensions of ccenergies correct?"""
        assert self.data.ccenergies.shape == \
                         (len(self.data.scfenergies),)

    def testsign(self):
        """Are the coupled cluster corrections negative?"""
        corrections = self.data.ccenergies - self.data.scfenergies
        assert numpy.alltrue(corrections < 0.0)


class GenericCC2Test(GenericCCTest):

    # Turbomole 7.4
    corr_energy = -1.1508051574141973

    def testenergycc2(self):
        """Is the CC2 correlation energy within the target?"""
        e_scf = self.data.scfenergies[0]
        e_cc = self.data.ccenergies[0]
        e_corr = e_cc - e_scf
        assert pytest.approx(e_corr, rel=self.rel_thresh) == self.corr_energy


class GenericCCDTest(GenericCCTest):

    # Q-Chem 5.4
    corr_energy = -1.4435403900740766

    def testenergyccd(self):
        """Is the CCD correlation energy within the target?"""
        e_scf = self.data.scfenergies[0]
        e_cc = self.data.ccenergies[0]
        e_corr = e_cc - e_scf
        assert pytest.approx(e_corr, rel=self.rel_thresh) == self.corr_energy


class GenericCCSDTest(GenericCCTest):

    # Q-Chem 5.4
    corr_energy = -1.4518567335733223

    def testenergyccsd(self):
        """Is the CCSD correlation energy within the target?"""
        e_scf = self.data.scfenergies[0]
        e_cc = self.data.ccenergies[0]
        e_corr = e_cc - e_scf
        assert pytest.approx(e_corr, rel=self.rel_thresh) == self.corr_energy


class GenericCCSDPTTest(GenericCCTest):

    # Q-Chem 5.4
    corr_energy = -1.4539460237174353

    def testenergyccsdpt(self):
        """Is the CCSD(T) correlation energy within the target?"""
        e_scf = self.data.scfenergies[0]
        e_cc = self.data.ccenergies[0]
        e_corr = e_cc - e_scf
        assert pytest.approx(e_corr, rel=self.rel_thresh) == self.corr_energy


class DALTONCCSDPTTest(GenericCCSDPTTest):

    # DALTON 2015
    corr_energy = -0.39526637011522325


class GamessCCDTest(GenericCCDTest):

    # GAMESS-US 2018
    corr_energy = -5.633089378137811


class GamessCCSDTest(GenericCCSDTest):

    # GAMESS-US 2018
    corr_energy = -5.655488429741581


class GamessCCSDPTTest(GenericCCSDPTTest):

    # GAMESS-US 2018
    corr_energy = -5.734634247949089


class MolcasCCSDTest(GenericCCSDTest):

    # OpenMolcas 18.0
    corr_energy = -1.3796976452617855


class MolproCCDTest(GenericCCDTest):

    # Molpro 2012
    corr_energy = -5.76891848850255


class MolproCCSDTest(GenericCCSDTest):

    # Molpro 2012
    corr_energy = -5.790734939563208


class MolproCCSDPTTest(GenericCCSDPTTest):

    # Molpro 2012
    corr_energy = -5.876248590504929


class NWChemCCSDPTTest(GenericCCSDPTTest):

    # NWChem 7.0
    corr_energy = -1.4559748390679488


class OrcaCCSDTest(GenericCCSDTest):

    # ORCA 5.0
    corr_energy = -1.3582212909295777


class OrcaCCSDPTTest(GenericCCSDPTTest):

    # ORCA 5.0
    corr_energy = -1.360085161959887


class Psi4CCSDTest(GenericCCSDTest):

    # Psi4 1.3.1
    corr_energy = -1.4538882732524598


class Psi4CCSDPTTest(GenericCCSDPTTest):

    # Psi4 1.3.1
    corr_energy = -1.4559754695608262


class TurbomoleCCSDTest(GenericCCSDTest):

    # TODO The two Turbomole CCSD(T) energies aren't similar...
    rel_thresh = 5.0e-3

    # Turbomole 7.4
    corr_energy = -1.3877950772714485


if __name__ == "__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['CC'])
    suite.testall()
