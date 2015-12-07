# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2007,2008,2012,2014,2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

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
        idx_lambdamax = [i for i, x in enumerate(self.data.etoscs)
                         if x == max(self.data.etoscs)][0]
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
        idx_minenergy = [i for i, x in enumerate(self.data.etenergies)
                         if x == min(self.data.etenergies)][0]
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
        self.assertAlmostEqual(max(self.data.etoscs), 1.0, delta=0.1)


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
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['TD'])
    suite.testall()
