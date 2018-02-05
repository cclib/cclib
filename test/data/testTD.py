# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test single point time-dependent logfiles in cclib"""

import os
import unittest

import numpy

from skip import skipForParser


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericTDTest(unittest.TestCase):
    """Generic time-dependent HF/DFT unittest"""

    number = 5
    expected_l_max = 41000

    @skipForParser('DALTON', 'etoscs are not parsed')
    def testenergies(self):
        """Is the l_max reasonable?"""

        self.assertEqual(len(self.data.etenergies), self.number)

        # Note that if all oscillator strengths are zero (like for triplets)
        # then this will simply pick out the first energy.
        idx_lambdamax = numpy.argmax(self.data.etoscs)
        self.assertAlmostEqual(self.data.etenergies[idx_lambdamax], self.expected_l_max, delta=5000)

    @skipForParser('DALTON', 'Oscillator strengths will have to be calculated, not just parsed.')
    def testoscs(self):
        """Is the maximum of etoscs in the right range?"""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertAlmostEqual(max(self.data.etoscs), 0.67, delta=0.1)

    @skipForParser('DALTON', '???')
    def testsecs(self):
        """Is the sum of etsecs close to 1?"""
        self.assertEqual(len(self.data.etsecs), self.number)
        lowestEtrans = self.data.etsecs[1]
        sumofsec = sum([z*z for (x, y, z) in lowestEtrans])
        self.assertAlmostEqual(sumofsec, 1.0, delta=0.16)

    @skipForParser('DALTON', '???')
    def testsecs_transition(self):
        """Is the lowest E transition from the HOMO or to the LUMO?"""
        idx_minenergy = numpy.argmin(self.data.etoscs)
        sec = self.data.etsecs[idx_minenergy]
        t = [(c*c, s, e) for (s, e, c) in sec]
        t.sort()
        t.reverse()
        self.assert_(t[0][1][0] == self.data.homos[0] or
                     t[0][2][0] == self.data.homos[0]+1, t[0])

    def testsymsnumber(self):
        """Is the length of etsyms correct?"""
        self.assertEqual(len(self.data.etsyms), self.number)


class ADFTDDFTTest(GenericTDTest):
    """Customized time-dependent DFT unittest"""
    number = 5

    def testsecs(self):
        """Is the sum of etsecs close to 1?"""
        self.assertEqual(len(self.data.etsecs), self.number)
        lowestEtrans = self.data.etsecs[1]

        #ADF squares the etsecs
        sumofsec = sum([z for (x, y, z) in lowestEtrans])
        self.assertAlmostEqual(sumofsec, 1.0, delta=0.16)


class DALTONTDTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    number = 20


class GaussianTDDFTTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    expected_l_max = 48000

    def testrotatsnumber(self):
        """Is the length of etrotats correct?"""
        self.assertEqual(len(self.data.etrotats), self.number)


class GAMESSUSTDDFTTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""
    number = 10


class JaguarTDDFTTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    expected_l_max = 48000

    def testoscs(self):
        """Is the maximum of etoscs in the right range?"""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertAlmostEqual(max(self.data.etoscs), 1.0, delta=0.2)


class OrcaTDDFTTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    number = 10
    expected_l_max = 48000

    def testoscs(self):
        """Is the maximum of etoscs in the right range?"""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertAlmostEqual(max(self.data.etoscs), .09, delta=0.01)

class OrcaROCISTest(GenericTDTest):
    """Customized test for ROCIS"""
    number = 57
    expected_l_max = 2316970.8
    n_spectra = 9

    def testoscs(self):
        """Is the maximum of etoscs in the right range?"""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertAlmostEqual(max(self.data.etoscs), 0.015, delta=0.1)

    def testTransprop(self):
        """Check the number of spectra parsed"""
        self.assertEqual(len(self.data.transprop), self.n_spectra)
        tddft_length = 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS'
        self.assertIn(tddft_length, self.data.transprop)

    def testsymsnumber(self):
        """ORCA ROCIS has no symmetry"""
        pass

    def testsecs(self):
        """ROCIS does not form singly excited configurations (secs)"""
        pass

    def testsecs_transition(self):
        """ROCIS does not form singly excited configurations (secs)"""
        pass


class QChemTDDFTTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    number = 10
    expected_l_max = 48000

    def testoscs(self):
        """Is the maximum of etoscs in the right range?"""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertAlmostEqual(max(self.data.etoscs), .9, delta=0.1)


class GenericTDDFTtrpTest(GenericTDTest):
    """Generic time-dependent HF/DFT (triplet) unittest"""

    number = 5
    expected_l_max = 24500

    def testoscs(self):
        """Triplet excitations should be disallowed."""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertAlmostEqual(max(self.data.etoscs), 0.0, delta=0.01)


if __name__=="__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['TD'])
    suite.testall()
