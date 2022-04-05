# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test logfiles with vibration output in cclib"""

import os
import unittest

from skip import skipForParser, skipForLogfile


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericIRTest(unittest.TestCase):
    """Generic vibrational frequency unittest"""

    # Unit tests should normally give this value for the largest IR intensity.
    max_IR_intensity = 100

    # Unit tests may give these values for the largest force constant and reduced mass, respectively.
    max_force_constant = 10.0
    max_reduced_mass = 6.9

    # reference zero-point correction from Gaussian 16 dvb_ir.out
    zpve = 0.1771

    def setUp(self):
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        self.numvib = 3*len(self.data.atomnos) - 6

    def testbasics(self):
        """Are basic attributes correct?"""
        self.assertEqual(self.data.natom, 20)

    @skipForParser('NWChem', 'Not implemented for this parser')
    def testvibdisps(self):
        """Are the dimensions of vibdisps consistent with numvib x N x 3"""
        self.assertEqual(len(self.data.vibfreqs), self.numvib)
        self.assertEqual(self.data.vibdisps.shape,
                         (self.numvib, len(self.data.atomnos), 3))

    def testlengths(self):
        """Are the lengths of vibfreqs and vibirs (and if present, vibsyms, vibfconnsts and vibrmasses) correct?"""
        self.assertEqual(len(self.data.vibfreqs), self.numvib)
        if hasattr(self.data, 'vibirs'):
            self.assertEqual(len(self.data.vibirs), self.numvib)
        if hasattr(self.data, 'vibsyms'):
            self.assertEqual(len(self.data.vibsyms), self.numvib)
        if hasattr(self.data, 'vibfconsts'):
            self.assertEqual(len(self.data.vibfconsts), self.numvib)
        if hasattr(self.data, 'vibrmasses'):
            self.assertEqual(len(self.data.vibrmasses), self.numvib)

    def testfreqval(self):
        """Is the highest freq value 3630 +/- 200 wavenumber?"""
        self.assertAlmostEqual(max(self.data.vibfreqs), 3630, delta=200)

    @skipForParser('Psi4', 'Psi cannot print IR intensities')
    def testirintens(self):
        """Is the maximum IR intensity 100 +/- 10 km/mol?"""
        self.assertAlmostEqual(max(self.data.vibirs), self.max_IR_intensity, delta=10)

    @skipForParser('ADF', 'ADF cannot print force constants')
    @skipForParser('DALTON', 'DALTON cannot print force constants')
    @skipForParser('GAMESS', 'GAMESS-US cannot print force constants')
    @skipForParser('GAMESSUK', 'GAMESS-UK cannot print force constants')
    @skipForParser('Molcas', 'Molcas cannot print force constants')
    @skipForParser('Molpro', 'Molpro cannot print force constants')
    @skipForParser('NWChem', 'Not implemented for this parser')
    @skipForParser('ORCA', 'ORCA cannot print force constants')
    @skipForParser('Turbomole', 'Turbomole cannot print force constants')
    @skipForLogfile('Jaguar/Jaguar4.2', 'Data file does not contain force constants')
    @skipForLogfile('Psi4/Psi4-1.0', 'Data file contains vibrational info with cartesian coordinates')
    def testvibfconsts(self):
        """Is the maximum force constant 10. +/- 0.1 mDyn/angstrom?"""
        self.assertAlmostEqual(max(self.data.vibfconsts), self.max_force_constant, delta=0.1)

    @skipForParser('ADF', 'ADF cannot print reduced masses')
    @skipForParser('DALTON', 'DALTON cannot print reduced masses')
    @skipForParser('GAMESSUK', 'GAMESSUK cannot print reduced masses')
    @skipForParser('Molpro', 'Molpro cannot print reduced masses')
    @skipForParser('NWChem', 'Not implemented for this parser')
    @skipForParser('ORCA', 'ORCA cannot print reduced masses')
    @skipForLogfile('GAMESS/PCGAMESS', 'Data file does not contain reduced masses')
    @skipForLogfile('Psi4/Psi4-1.0', 'Data file does not contain reduced masses')
    def testvibrmasses(self):
        """Is the maximum reduced mass 6.9 +/- 0.1 daltons?"""
        self.assertAlmostEqual(max(self.data.vibrmasses), self.max_reduced_mass, delta=0.1)

    @skipForParser('Psi3', 'not implemented yet')
    def testzeropointcorrection(self):
        """Is the zero-point correction correct?"""
        self.assertAlmostEqual(self.data.zpve, self.zpve, delta=1.0e-3)


class ADFIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    # ???
    def testzeropointcorrection(self):
        """Is the zero-point correction correct?"""
        self.assertAlmostEqual(self.data.zpve, self.zpve, delta=1.0e-2)


class FireflyIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    max_IR_intensity = 135
    # ???
    zpve = 0.1935


class GaussianIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    def testvibsyms(self):
        """Is the length of vibsyms correct?"""
        self.assertEqual(len(self.data.vibsyms), self.numvib)

    def testzeropointcorrection(self):
        # reference zero-point correction from dvb_ir.out
        zpve = 0.1771
        """Is the zero-point correction correct?"""
        self.assertAlmostEqual(self.data.zpve, zpve, delta=0.001)

    entropy_places = 6
    enthalpy_places = 3
    freeenergy_places = 3

    def testtemperature(self):
        """Is the temperature 298.15 K?"""
        self.assertAlmostEqual(298.15, self.data.temperature)

    def testpressure(self):
        """Is the pressure 1 atm?"""
        self.assertAlmostEqual(1, self.data.pressure)

    def testentropy(self):
         """Is the entropy reasonable"""
         self.assertAlmostEqual(0.0001462623335480945, self.data.entropy, self.entropy_places)

    def testenthalpy(self):
         """Is the enthalpy reasonable"""
         self.assertAlmostEqual(-382.12130688525264, self.data.enthalpy, self.enthalpy_places)

    def testfreeenergy(self):
         """Is the freeenergy reasonable"""
         self.assertAlmostEqual(-382.164915, self.data.freeenergy, self.freeenergy_places)

    def testfreeenergyconsistency(self):
        """Does G = H - TS hold"""
        self.assertAlmostEqual(self.data.enthalpy - self.data.temperature * self.data.entropy, self.data.freeenergy, self.freeenergy_places)

class JaguarIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    # Jagur outputs vibrational info with cartesian coordinates
    max_force_constant = 3.7
    max_reduced_mass = 2.3

    def testvibsyms(self):
        """Is the length of vibsyms correct?"""
        self.assertEqual(len(self.data.vibsyms), self.numvib)


class MolcasIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    max_IR_intensity = 65
    zpve = 0.1783

    entropy_places = 6
    enthalpy_places = 3
    freeenergy_places = 3

    def testtemperature(self):
        """Is the temperature 298.15 K?"""
        self.assertAlmostEqual(298.15, self.data.temperature)

    def testpressure(self):
        """Is the pressure 1 atm?"""
        self.assertAlmostEqual(1, self.data.pressure)

    def testentropy(self):
         """Is the entropy reasonable"""
         self.assertAlmostEqual(0.00013403320476271246, self.data.entropy, self.entropy_places)

    def testenthalpy(self):
         """Is the enthalpy reasonable"""
         self.assertAlmostEqual(-382.11385, self.data.enthalpy, self.enthalpy_places)

    def testfreeenergy(self):
         """Is the freeenergy reasonable"""
         self.assertAlmostEqual(-382.153812, self.data.freeenergy, self.freeenergy_places)

    def testfreeenergyconsistency(self):
        """Does G = H - TS hold"""
        self.assertAlmostEqual(self.data.enthalpy - self.data.temperature * self.data.entropy, self.data.freeenergy, self.freeenergy_places)

class NWChemIRTest(GenericIRTest):
    """Generic imaginary vibrational frequency unittest"""

    def setUp(self):
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        self.numvib = 3*len(self.data.atomnos)


class OrcaIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    # ORCA has a bug in the intensities for version < 4.0
    max_IR_intensity = 215
    zpve = 0.1921

    enthalpy_places = 3
    entropy_places = 6
    freeenergy_places = 3

    def testtemperature(self):
        """Is the temperature 298.15 K?"""
        self.assertAlmostEqual(298.15, self.data.temperature)

    def testpressure(self):
        """Is the pressure 1 atm?"""
        self.assertAlmostEqual(1, self.data.pressure)

    def testenthalpy(self):
         """Is the enthalpy reasonable"""
         self.assertAlmostEqual(-381.85224835, self.data.enthalpy, self.enthalpy_places)

    def testentropy(self):
         """Is the entropy reasonable"""
         self.assertAlmostEqual(0.00012080325339594164, self.data.entropy, self.entropy_places)

    def testfreeenergy(self):
         """Is the freeenergy reasonable"""
         self.assertAlmostEqual(-381.88826585, self.data.freeenergy, self.freeenergy_places)

    def testfreeenergyconsistency(self):
        """Does G = H - TS hold"""
        self.assertAlmostEqual(self.data.enthalpy - self.data.temperature * self.data.entropy, self.data.freeenergy, self.freeenergy_places)


class QChemIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    enthalpy_places = 3
    entropy_places = 6
    freeenergy_places = 3

    def testtemperature(self):
        """Is the temperature 298.15 K?"""
        self.assertEqual(298.15, self.data.temperature)

    def testpressure(self):
        """Is the pressure 1 atm?"""
        self.assertAlmostEqual(1, self.data.pressure)

    def testenthalpy(self):
         """Is the enthalpy reasonable"""
         self.assertAlmostEqual(0.1871270552135131, self.data.enthalpy, self.enthalpy_places)

    def testentropy(self):
         """Is the entropy reasonable"""
         self.assertAlmostEqual(0.00014667348271900577, self.data.entropy, self.entropy_places)

    def testfreeenergy(self):
         """Is the freeenergy reasonable"""
         self.assertAlmostEqual(0.14339635634084155, self.data.freeenergy, self.freeenergy_places)

    # Molecular mass of DVB in mD.
    molecularmass = 130078.25

    def testatommasses(self):
        """Do the atom masses sum up to the molecular mass (130078.25+-0.1mD)?"""
        mm = 1000 * sum(self.data.atommasses)
        self.assertAlmostEqual(
            mm, 130078.25, delta=0.1, msg=f"Molecule mass: {mm:f} not 130078 +- 0.1mD"
        )

    def testhessian(self):
        """Do the frequencies from the Hessian match the printed frequencies?"""

    def testfreeenergyconsistency(self):
        """Does G = H - TS hold"""
        self.assertAlmostEqual(self.data.enthalpy - self.data.temperature * self.data.entropy, self.data.freeenergy, self.freeenergy_places)


class GamessIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""
    # Molecular mass of DVB in mD.
    molecularmass = 130078.25
    enthalpy_places = 3
    entropy_places = 6
    freeenergy_places = 3

    def testatommasses(self):
        """Do the atom masses sum up to the molecular mass (130078.25+-0.1mD)?"""
        mm = 1000 * sum(self.data.atommasses)
        self.assertAlmostEqual(
            mm, 130078.25, delta=0.1, msg=f"Molecule mass: {mm:f} not 130078 +- 0.1mD"
        )

    def testtemperature(self):
        """Is the temperature 298.15 K?"""
        self.assertAlmostEqual(298.15, self.data.temperature)

    def testpressure(self):
        """Is the pressure 1 atm?"""
        self.assertAlmostEqual(1, self.data.pressure)

    def testenthalpy(self):
         """Is the enthalpy reasonable"""
         self.assertAlmostEqual(-381.86372805188300, self.data.enthalpy, self.enthalpy_places)

    def testentropy(self):
         """Is the entropy reasonable"""
         self.assertAlmostEqual(0.00014875961938, self.data.entropy, self.entropy_places)

    def testfreeenergy(self):
         """Is the freeenergy reasonable"""
         self.assertAlmostEqual(-381.90808120060200, self.data.freeenergy, self.freeenergy_places)

    def testfreeenergyconsistency(self):
        """Does G = H - TS hold"""
        self.assertAlmostEqual(self.data.enthalpy - self.data.temperature * self.data.entropy, self.data.freeenergy, self.freeenergy_places)


class Psi4IRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    # RHF is used for Psi4 IR test data instead of B3LYP
    max_force_constant = 9.37
    zpve = 0.1917


class TurbomoleIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    # ???
    zpve = 0.1725


class GenericIRimgTest(unittest.TestCase):
    """Generic imaginary vibrational frequency unittest"""

    def setUp(self):
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        self.numvib = 3*len(self.data.atomnos) - 6

    def testvibdisps(self):
        """Are the dimensions of vibdisps consistent with numvib x N x 3"""
        self.assertEqual(self.data.vibdisps.shape,
                         (self.numvib, len(self.data.atomnos), 3))

    def testlengths(self):
        """Are the lengths of vibfreqs and vibirs correct?"""
        self.assertEqual(len(self.data.vibfreqs), self.numvib)
        self.assertEqual(len(self.data.vibirs), self.numvib)

    def testfreqval(self):
        """Is the lowest freq value negative?"""
        self.assertTrue(self.data.vibfreqs[0] < 0)

##    def testmaxvibdisps(self):
##        """What is the maximum value of displacement for a H vs a C?"""
##        Cvibdisps = compress(self.data.atomnos==6, self.data.vibdisps, 1)
##        Hvibdisps = compress(self.data.atomnos==1, self.data.vibdisps, 1)
##        self.assertEqual(max(abs(Cvibdisps).flat), 1.0)


class GenericRamanTest(unittest.TestCase):
    """Generic Raman unittest"""

    # This value is in amu.
    max_raman_intensity = 575

    def setUp(self):
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        self.numvib = 3*len(self.data.atomnos) - 6

    def testlengths(self):
        """Is the length of vibramans correct?"""
        self.assertEqual(len(self.data.vibramans), self.numvib)

    # The tolerance for this number has been increased, since ORCA
    # failed to make it inside +/-5, but it would be nice in the future
    # to determine is it's not too much work whether this is due to
    # algorithmic differences, or to differences in the input basis set
    # or coordinates. The first would be OK, but in the second case the
    # unit test jobs should be made more comparable. With cclib, we first
    # of all want to succeed in parsing, but would also like to remain
    # as comparable between programs as possible (for these tests).
    # Note also that this value is adjusted for Gaussian and DALTON - why?
    def testramanintens(self):
        """Is the maximum Raman intensity correct?"""
        self.assertAlmostEqual(max(self.data.vibramans), self.max_raman_intensity, delta=8)

        # We used to test this, but it seems to vary wildly between
        # programs... perhaps we could use it if we knew why...
        #self.assertInside(self.data.vibramans[1], 2.6872, 0.0001)

    def testvibdisps(self):
        """Is the length and value of vibdisps correct?"""
        assert hasattr(self.data, "vibdisps")
        assert len(self.data.vibdisps) == 54


class DALTONRamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    max_raman_intensity = 745


class GaussianRamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    max_raman_intensity = 1066


class OrcaRamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    max_raman_intensity = 1045


class QChemRamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    max_raman_intensity = 588

if __name__=="__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['vib'])
    suite.testall()
