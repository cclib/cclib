# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test configuration interaction (CI) logfiles in cclib"""

import os
import unittest

import numpy

from skip import skipForParser

__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericCISTest(unittest.TestCase):
    """Generic CIS(RHF)/STO-3G water unittest"""

    nstates = 5

    # First four singlet/triplet state excitation energies [wavenumber].
    # Based on output in GAMESS test.
    etenergies0 = numpy.array([98614.56, 114906.59, 127948.12, 146480.64])
    etenergies1 = numpy.array([82085.34,  98999.11, 104077.89, 113978.37])

    # First four singlet/triplet state excitation orbitals and coefficients.
    # Tuples: (from MO, to MO, coefficient) - don't need spin indices.
    # Based on output in GAMESS test (using the "dets" algorithm).
    # Note that only coefficients larger than 0.1 are included here, as
    #   the Gaussian test does not contain smaller ones.
    # The simple example of water should yield the same first 4 states in all programs.
    # Note that the GAMESS test "water_cis_dets" calculates also triplet states,
    #  but the resulting transition dipoles and oscillator strengths are not correct
    #  and this file is not tested here (although it probably should be).
    etsecs0 = [ [(4, 5, -1.0)],
                [(4, 6, -1.0)],
                [(3, 5,  0.96689)],
                [(2, 5,  0.4407), (3, 6, 0.89770)] ]
    etsecs1 = [ [(4, 5,  1.0)],
                [(2, 6, -0.231), (3, 5, -0.9676)],
                [(4, 6,  1.0)],
                [(2, 5, -0.536), (3, 6, -0.843)] ]

    etsecs_precision = 0.0005

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testetenergiesvalues(self):
        """ Are etenergies within 50 wavenumbers of the correct values?"""
        indices0 = [i for i in range(self.nstates) if self.data.etsyms[i][0] == "S"]
        indices1 = [i for i in range(self.nstates) if self.data.etsyms[i][0] == "T"]
        singlets = [self.data.etenergies[i] for i in indices0]
        triplets = [self.data.etenergies[i] for i in indices1]
        # All programs do singlets.
        singletdiff = singlets[:4] - self.etenergies0
        self.assertTrue(numpy.alltrue(singletdiff < 50))
        # Not all programs do triplets (i.e. Jaguar).
        if len(triplets) >= 4:
            tripletdiff = triplets[:4] - self.etenergies1
            self.assertTrue(numpy.alltrue(tripletdiff < 50))

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testsecs(self):
        """Is the sum of etsecs close to 1?"""
        etsec = self.data.etsecs[2] # Pick one with several contributors
        sumofsec = sum([z*z for (x, y, z) in etsec])
        self.assertAlmostEqual(sumofsec, 1.0, delta=0.02)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testetsecsvalues(self):
        """ Are etsecs correct and coefficients close to the correct values?"""
        indices0 = [i for i in range(self.nstates) if self.data.etsyms[i][0] == "S"]
        indices1 = [i for i in range(self.nstates) if self.data.etsyms[i][0] == "T"]
        singlets = [self.data.etsecs[i] for i in indices0]
        triplets = [self.data.etsecs[i] for i in indices1]
        # All programs do singlets.
        for i in range(4):
            for exc in self.etsecs0[i]:
                found = False
                for s in singlets[i]:
                    if s[0][0] == exc[0] and s[1][0] == exc[1]:
                        found = True
                        self.assertAlmostEqual(abs(s[2]), abs(exc[2]), delta=self.etsecs_precision)
                if not found:
                    self.fail("Excitation %i->%s not found (singlet state %i)" %(exc[0], exc[1], i))
        # Not all programs do triplets (i.e. Jaguar).
        if len(triplets) >= 4:
            for i in range(4):
                for exc in self.etsecs1[i]:
                    found = False
                    for s in triplets[i]:
                        if s[0][0] == exc[0] and s[1][0] == exc[1]:
                            found = True
                            self.assertAlmostEqual(abs(s[2]), abs(exc[2]), delta=self.etsecs_precision)
                    if not found:
                        self.fail("Excitation %i->%s not found (triplet state %i)" %(exc[0], exc[1], i))


class GAMESSCISTest(GenericCISTest):
    """Customized CIS(RHF)/STO-3G water unittest"""

    def testnocoeffs(self):
        """Are natural orbital coefficients the right size?"""
        self.assertEqual(self.data.nocoeffs.shape, (self.data.nmo, self.data.nbasis))

    def testnooccnos(self):
        """Are natural orbital occupation numbers the right size?"""
        self.assertEqual(self.data.nooccnos.shape, (self.data.nmo, ))


class GaussianCISTest(GenericCISTest):
    """Customized CIS(RHF)/STO-3G water unittest"""
    nstates = 10

    def testnocoeffs(self):
        """Are natural orbital coefficients the right size?"""
        self.assertEqual(self.data.nocoeffs.shape, (self.data.nmo, self.data.nbasis))

    def testnooccnos(self):
        """Are natural orbital occupation numbers the right size?"""
        self.assertEqual(self.data.nooccnos.shape, (self.data.nmo, ))



class Jaguar83CISTest(GenericCISTest):
    """Customized CIS(RHF)/STO-3G water unittest"""

    # The Jaguar8.3 job was created using 6-31G instead of STO-3G.
    etsecs0 = [ [(4, 5,  0.99186)],
                [(4, 6, -0.98594)],
                [(3, 5, -0.98321)],
                [(2, 5,  0.19240), (3, 6, -0.97090)] ]
    etsecs1 = [ [(4, 5,  1.0)],
                [(2, 6, -0.231), (3, 5, -0.9676)],
                [(4, 6,  1.0)],
                [(2, 5, -0.536), (3, 6, -0.843)] ]


class QChemCISTest(GenericCISTest):
    """Customized CIS(RHF)/STO-3G water unittest"""
    nstates = 10


if __name__=="__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['CI'])
    suite.testall()
