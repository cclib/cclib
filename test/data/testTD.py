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
from skip import skipForLogfile


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericTDTest(unittest.TestCase):
    """Generic time-dependent HF/DFT unittest"""

    number = 5
    expected_l_max = 41000

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForLogfile('Turbomole/basicTurbomole7.4/CO_cc2_TD_trip', 'Oscillator strengths are not available for Turbomole triplets using ricc2 but are required for testenergies()')
    def testenergies(self):
        """Is the l_max reasonable?"""

        self.assertEqual(len(self.data.etenergies), self.number)

        # Note that if all oscillator strengths are zero (like for triplets)
        # then this will simply pick out the first energy.
        idx_lambdamax = numpy.argmax(self.data.etoscs)
        self.assertAlmostEqual(self.data.etenergies[idx_lambdamax], self.expected_l_max, delta=5000)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForLogfile("Turbomole/basicTurbomole7.4/CO_cc2_TD_trip", "Oscillator strengths are not available for triplets with Turbomole's ricc2")
    def testoscs(self):
        """Is the maximum of etoscs in the right range?"""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertAlmostEqual(max(self.data.etoscs), 0.67, delta=0.1)

    @skipForParser('FChk','The parser is still being developed so we skip this test')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    def testsecs(self):
        """Is the sum of etsecs close to 1?"""
        self.assertEqual(len(self.data.etsecs), self.number)
        lowestEtrans = self.data.etsecs[numpy.argmin(self.data.etenergies)]
        sumofsec = sum([z*z for (x, y, z) in lowestEtrans])
        self.assertAlmostEqual(sumofsec, 1.0, delta=0.16)

    @skipForParser('FChk', 'This is true for calculations without symmetry, but not with?')
    @skipForParser('DALTON', 'This is true for calculations without symmetry, but not with?')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    def testsecs_transition(self):
        """Is the lowest E transition from the HOMO or to the LUMO?"""
        lowestEtrans = self.data.etsecs[numpy.argmin(self.data.etenergies)]
        t = list(reversed(sorted([(c*c, s, e) for (s, e, c) in lowestEtrans])))
        self.assertTrue(t[0][1][0] == self.data.homos[0] or
                        t[0][2][0] == self.data.homos[0] + 1, t[0])

    @skipForParser('Molcas','The parser is still being developed so we skip this test')    
    def testsymsnumber(self):
        """Is the length of etsyms correct?"""
        self.assertEqual(len(self.data.etsyms), self.number)

    @skipForParser('ADF', 'etrotats are not yet implemented')
    @skipForParser('DALTON', 'etrotats are not yet implemented')
    @skipForParser('FChk', 'etrotats are not yet implemented')
    @skipForParser('GAMESS', 'etrotats are not yet implemented')
    @skipForParser('GAMESSUK', 'etrotats are not yet implemented')
    @skipForParser('Jaguar', 'etrotats are not yet implemented')
    @skipForParser('NWChem', 'etrotats are not yet implemented')
    @skipForParser('QChem', 'Q-Chem cannot calculate rotatory strengths')
    @skipForLogfile("Turbomole/basicTurbomole7.4/CO_cc2_TD", "Rotatory strengths are not currently available for ricc2")
    @skipForLogfile("Turbomole/basicTurbomole7.4/CO_adc2_TD", "Rotatory strengths are not currently available for ricc2")
    def testrotatsnumber(self):
        """Is the length of etrotats correct?"""
        self.assertEqual(len(self.data.etrotats), self.number)

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

    def testetdipsshape(self):
        """Is the shape of etdips correct?"""
        self.assertEqual(numpy.shape(self.data.etdips), (self.number, 3))

    def testetveldipsshape(self):
        """Is the shape of etveldips correct?"""
        self.assertEqual(numpy.shape(
            self.data.etveldips), (self.number, 3))

    def testetmagdipsshape(self):
        """Is the shape of etmagdips correct?"""
        self.assertEqual(numpy.shape(self.data.etmagdips), (self.number, 3))

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
        self.assertAlmostEqual(max(self.data.etoscs), 1.17, delta=0.01)


class QChemTDDFTTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    number = 10
    expected_l_max = 48000

    def testoscs(self):
        """Is the maximum of etoscs in the right range?"""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertAlmostEqual(max(self.data.etoscs), 0.9, delta=0.1)


class GenericTDDFTtrpTest(GenericTDTest):
    """Generic time-dependent HF/DFT (triplet) unittest"""

    number = 5
    expected_l_max = 24500

    def testoscs(self):
        """Triplet excitations should be disallowed."""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertAlmostEqual(max(self.data.etoscs), 0.0, delta=0.01)


class OrcaROCISTest(GenericTDTest):
    """Customized test for ROCIS"""
    number = 57
    expected_l_max = 2316970.8
    n_spectra = 2

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

    def testrotatsnumber(self):
        """ROCIS does not calculate rotatory strengths"""
        pass


class OrcaROCIS40Test(OrcaROCISTest):
    """Customized test for ROCIS"""
    # In ORCA 4.0, an additional spectrum ("COMBINED ELECTRIC DIPOLE +
    # MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (Origin Independent,
    # Length Representation)") was present that is not in ORCA 4.1.
    # Changed via 1085. Only the Electric Dipole are parsed
    n_spectra = 2


class TurbomoleTDTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""
    
    number = 10
    expected_l_max = 91432
    
    def testoscs(self):
        """Is the maximum of etoscs in the right range?"""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertAlmostEqual(max(self.data.etoscs), 0.19, delta=0.1)

class TurbomoleTDADC2Test(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""
    
    number = 10
    expected_l_max = 136329
    
    def testoscs(self):
        """Is the maximum of etoscs in the right range?"""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertAlmostEqual(max(self.data.etoscs), 0.80, delta=0.1)

class TurbomoleTDTripTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""
    
    number = 10
    expected_l_max = 51530

    def testoscs(self):
        """Is the maximum of etoscs in the right range?"""
        self.assertEqual(len(self.data.etoscs), self.number)
        self.assertAlmostEqual(max(self.data.etoscs), 0.84, delta=0.1)
        
class TurbomoleTDCC2TripTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""
    # This test is for triplets with ricc2, which does not support oscillator strengths.
    
    number = 10

    def testenergies(self):
        """Is the l_max reasonable?"""
        self.assertEqual(len(self.data.etenergies), self.number)
        

if __name__ =="__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['TD'])
    suite.testall()
    
