# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test logfiles with vibration output in cclib"""

import pytest
from skip import skipForLogfile, skipForParser


class GenericIRTest:
    """Generic vibrational frequency unittest"""

    highest_freq = 3630
    highest_freq_thresh = 200

    # Unit tests should normally give this value for the largest IR intensity.
    max_IR_intensity = 100

    # Unit tests may give these values for the largest force constant and reduced mass, respectively.
    max_force_constant = 10.0
    force_constant_thresh = 0.1
    max_reduced_mass = 6.9

    # reference zero-point correction from Gaussian 16 dvb_ir.out
    zpve = 0.1771
    entropy = 0.0001462623335480945
    enthalpy = -382.12130688525264
    freeenergy = -382.164915

    zpve_thresh = 1.0e-3
    # TODO refactor from places to thresh
    enthalpy_places = 3
    entropy_places = 6
    freeenergy_places = 3

    # Molecular mass of DVB in mD.
    molecularmass = 130078.25
    molecularmass_thresh = 0.25

    @pytest.fixture
    def numvib(self, data) -> int:
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        return 3 * len(data.atomnos) - 6

    def testbasics(self, data) -> None:
        """Are basic attributes correct?"""
        assert data.natom == 20

    @skipForLogfile("FChk/basicGaussian09", "not printed in older versions than 16")
    @skipForLogfile("FChk/basicQChem5.4", "not printed")
    @skipForParser("xTB", "Custom treatment")
    def testvibdisps(self, data, numvib) -> None:
        """Are the dimensions of vibdisps consistent with numvib x N x 3"""
        assert len(data.vibfreqs) == numvib
        assert data.vibdisps.shape == (numvib, len(data.atomnos), 3)

    @skipForLogfile("FChk/basicGaussian09", "not printed in older versions than 16")
    @skipForLogfile("FChk/basicQChem5.4", "not printed")
    def testlengths(self, data, numvib) -> None:
        """Are the lengths of vibfreqs and vibirs (and if present, vibsyms, vibfconnsts and vibrmasses) correct?"""
        assert len(data.vibfreqs) == numvib
        if hasattr(data, "vibirs"):
            assert len(data.vibirs) == numvib
        if hasattr(data, "vibsyms"):
            assert len(data.vibsyms) == numvib
        if hasattr(data, "vibfconsts"):
            assert len(data.vibfconsts) == numvib
        if hasattr(data, "vibrmasses"):
            assert len(data.vibrmasses) == numvib

    @skipForLogfile("FChk/basicGaussian09", "not printed in older versions than 16")
    @skipForLogfile("FChk/basicQChem5.4", "not printed")
    def testfreqval(self, data) -> None:
        """Does the highest frequency value match?"""
        assert abs(max(data.vibfreqs) - self.highest_freq) < self.highest_freq_thresh

    @skipForLogfile("FChk/basicGaussian09", "not printed in older versions than 16")
    @skipForLogfile("FChk/basicQChem5.4", "not printed")
    @skipForLogfile(
        "Psi4/basicPsi4-1.2.1/dvb_ir_rhf.out", "not implemented in versions older than 1.7"
    )
    @skipForLogfile(
        "Psi4/basicPsi4-1.3.1/dvb_ir_rhf.out", "not implemented in versions older than 1.7"
    )
    def testirintens(self, data) -> None:
        """Is the maximum IR intensity 100 +/- 10 km/mol?"""
        assert abs(max(data.vibirs) - self.max_IR_intensity) < 10

    @skipForParser("ADF", "ADF cannot print force constants")
    @skipForParser("DALTON", "DALTON cannot print force constants")
    @skipForParser("GAMESS", "GAMESS-US cannot print force constants")
    @skipForParser("GAMESSUK", "GAMESS-UK cannot print force constants")
    @skipForParser("Molcas", "Molcas cannot print force constants")
    @skipForParser("Molpro", "Molpro cannot print force constants")
    @skipForParser("NWChem", "Not implemented for this parser")
    @skipForParser("ORCA", "ORCA cannot print force constants")
    @skipForParser("Turbomole", "Turbomole cannot print force constants")
    @skipForParser("xTB", "xTB does not print force constants")
    @skipForLogfile("FChk/basicGaussian09", "not printed in older versions than 16")
    @skipForLogfile("FChk/basicQChem5.4", "not printed")
    def testvibfconsts(self, data) -> None:
        """Is the maximum force constant 10. +/- 0.1 mDyn/angstrom?"""
        assert abs(max(data.vibfconsts) - self.max_force_constant) < self.force_constant_thresh

    @skipForParser("ADF", "ADF cannot print reduced masses")
    @skipForParser("DALTON", "DALTON cannot print reduced masses")
    @skipForParser("GAMESSUK", "GAMESSUK cannot print reduced masses")
    @skipForParser("Molpro", "Molpro cannot print reduced masses")
    @skipForParser("NWChem", "Not implemented for this parser")
    @skipForParser("ORCA", "ORCA cannot print reduced masses")
    @skipForLogfile("FChk/basicGaussian09", "not printed in older versions than 16")
    @skipForLogfile("FChk/basicQChem5.4", "not printed")
    @skipForLogfile("GAMESS/PCGAMESS", "Data file does not contain reduced masses")
    def testvibrmasses(self, data) -> None:
        """Is the maximum reduced mass 6.9 +/- 0.1 daltons?"""
        assert abs(max(data.vibrmasses) - self.max_reduced_mass) < 0.1

    @skipForParser("FChk", "not printed")
    @skipForParser("Psi3", "not implemented yet")
    def testzeropointcorrection(self, data) -> None:
        """Is the zero-point correction correct?"""
        assert abs(data.zpve - self.zpve) < self.zpve_thresh

    @skipForParser("ADF", "not implemented yet")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Gaussian", "not implemented yet")
    @skipForParser("Jaguar", "not implemented yet")
    @skipForParser("Molcas", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("ORCA", "not implemented yet")
    @skipForParser("Psi4", "not implemented yet")
    @skipForLogfile(
        "QChem/basicQChem5.4/dvb_ir.out", "needs to be rerun with print level turned up"
    )
    @skipForParser("Turbomole", "not implemented yet")
    @skipForParser("xTB", "not implemented yet")
    def testhessian(self, data) -> None:
        """Are the dimensions of the molecular Hessian correct?"""
        assert data.hessian.shape == (3 * data.natom, 3 * data.natom)

    def testhessian_frequencies(self, data) -> None:
        """Do the frequencies from the Hessian match the printed frequencies?"""

    @skipForParser("DALTON", "not implemented yet")
    @skipForParser("FChk", "not printed")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    def testtemperature(self, data) -> None:
        """Is the temperature 298.15 K?"""
        assert round(abs(298.15 - data.temperature), 7) == 0

    @skipForParser("DALTON", "not implemented yet")
    @skipForParser("FChk", "not printed")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("Psi4", "not implemented yet")
    @skipForParser("xTB", "not printed")
    @skipForParser("Turbomole", "not implemented yet")
    def testpressure(self, data) -> None:
        """Is the pressure 1 atm?"""
        assert round(abs(1 - data.pressure), 7) == 0

    @skipForParser("DALTON", "not implemented yet")
    @skipForParser("FChk", "not printed")
    @skipForParser("Jaguar", "not implemented yet")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    def testentropy(self, data) -> None:
        """Is the entropy reasonable"""
        assert round(abs(self.entropy - data.entropy), self.entropy_places) == 0

    @skipForParser("ADF", "not implemented yet")
    @skipForParser("DALTON", "not implemented yet")
    @skipForParser("FChk", "not printed")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    def testenthalpy(self, data) -> None:
        """Is the enthalpy reasonable"""
        assert round(abs(self.enthalpy - data.enthalpy), self.enthalpy_places) == 0

    @skipForParser("ADF", "not implemented yet")
    @skipForParser("DALTON", "not implemented yet")
    @skipForParser("FChk", "not printed")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    def testfreeenergy(self, data) -> None:
        """Is the freeenergy reasonable"""
        assert round(abs(self.freeenergy - data.freeenergy), self.freeenergy_places) == 0

    @skipForParser("ADF", "not implemented yet")
    @skipForParser("DALTON", "not implemented yet")
    @skipForParser("FChk", "not printed")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    def testfreeenergyconsistency(self, data) -> None:
        """Does G = H - TS hold"""
        assert (
            round(
                abs(data.enthalpy - data.temperature * data.entropy - data.freeenergy),
                self.freeenergy_places,
            )
            == 0
        )

    @skipForParser("FChk", "not printed")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Jaguar", "not implemented yet")
    @skipForParser("Molcas", "not implemented yet")
    @skipForParser("Molpro", "not implemented yet")
    @skipForParser("GAMESSUK", "not implemented yet")
    @skipForParser("Turbomole", "not implemented yet")
    @skipForParser("xTB", "not implemented yet")
    def testatommasses(self, data) -> None:
        """Do the atom masses sum up to the molecular mass?"""
        mm = 1000 * sum(data.atommasses)
        assert abs(mm - self.molecularmass) < self.molecularmass_thresh, (
            f"Molecule mass: {mm:f} not {self.molecularmass:f} +- {self.molecularmass_thresh:f} mD"
        )


class ADFIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    zpve = 0.1759
    entropy = 0.00013953096105839138

    zpve_thresh = 1.1e-3
    entropy_places = 4


class FireflyIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    max_IR_intensity = 135

    zpve = 0.1935
    enthalpy = -379.5751787863937
    freeenergy = -379.61838132136285

    entropy_places = 5


class GaussianIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    @skipForParser("FChk", "not printed")
    def testvibsyms(self, data, numvib) -> None:
        """Is the length of vibsyms correct?"""
        assert len(data.vibsyms) == numvib


class JaguarIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    # Jagur outputs vibrational info with cartesian coordinates
    max_force_constant = 3.7
    max_reduced_mass = 2.3

    freeenergy_places = 2

    def testvibsyms(self, data, numvib) -> None:
        """Is the length of vibsyms correct?"""
        assert len(data.vibsyms) == numvib


class MolcasIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    max_IR_intensity = 65

    zpve = 0.1783
    entropy = 0.00013403320476271246
    enthalpy = -382.11385
    freeenergy = -382.153812


class NWChemIRTest(GenericIRTest):
    """Generic imaginary vibrational frequency unittest"""

    @pytest.fixture
    def numvib(self, data) -> int:
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        return 3 * len(data.atomnos)


class GamessIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    entropy = 0.00014875961938
    enthalpy = -381.86372805188300
    freeenergy = -381.90808120060200


class OrcaIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    # ORCA 5.0
    entropy = 0.00014384698977024988
    enthalpy = -381.86823907
    freeenergy = -381.91112705

    enthalpy_places = 2
    entropy_places = 5
    freeenergy_places = 2

    molecularmass = 130190


class Psi4HFIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    max_IR_intensity = 146
    max_force_constant = 9.37

    zpve = 0.1917
    entropy = 0.00013229523
    enthalpy = -379.57027841
    freeenergy = -379.60972224


class Psi4KSIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    enthalpy_places = 2
    freeenergy_places = 2


class TurbomoleIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    # ???
    zpve = 0.1725


class XTBIRTest(GenericIRTest):
    """Customized vibrational frequency unittest"""

    highest_freq = 3131.43
    highest_freq_thresh = 10

    zpve = 0.161236858990
    entropy = 0.0001493609321651285
    enthalpy = -26.266467652754
    freeenergy = -26.310999637373
    max_reduced_mass = 11.43

    @pytest.fixture
    def numvib(self, data) -> int:
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        return 3 * len(data.atomnos)


class GenericIRimgTest:
    """Generic imaginary vibrational frequency unittest"""

    @pytest.fixture
    def numvib(self, data) -> int:
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        return 3 * len(data.atomnos) - 6

    def testvibdisps(self, data, numvib) -> None:
        """Are the dimensions of vibdisps consistent with numvib x N x 3"""
        assert data.vibdisps.shape == (numvib, len(data.atomnos), 3)

    def testlengths(self, data, numvib) -> None:
        """Are the lengths of vibfreqs and vibirs correct?"""
        assert len(data.vibfreqs) == numvib
        assert len(data.vibirs) == numvib

    def testfreqval(self, data) -> None:
        """Is the lowest freq value negative?"""
        assert data.vibfreqs[0] < 0


##    def testmaxvibdisps(self, data) -> None:
##        """What is the maximum value of displacement for a H vs a C?"""
##        Cvibdisps = compress(data.atomnos==6, data.vibdisps, 1)
##        Hvibdisps = compress(data.atomnos==1, data.vibdisps, 1)
##        self.assertEqual(max(abs(Cvibdisps).flat), 1.0)


class GenericRamanTest:
    """Generic Raman unittest"""

    # This value is in amu.
    max_raman_intensity = 575

    @pytest.fixture
    def numvib(self, data) -> int:
        """Initialize the number of vibrational frequencies on a per molecule basis"""
        return 3 * len(data.atomnos) - 6

    def testlengths(self, data, numvib) -> None:
        """Is the length of vibramans correct?"""
        assert len(data.vibramans) == numvib

    # The tolerance for this number has been increased, since ORCA
    # failed to make it inside +/-5, but it would be nice in the future
    # to determine is it's not too much work whether this is due to
    # algorithmic differences, or to differences in the input basis set
    # or coordinates. The first would be OK, but in the second case the
    # unit test jobs should be made more comparable. With cclib, we first
    # of all want to succeed in parsing, but would also like to remain
    # as comparable between programs as possible (for these tests).
    # Note also that this value is adjusted for Gaussian and DALTON - why?
    def testramanintens(self, data) -> None:
        """Is the maximum Raman intensity correct?"""
        assert abs(max(data.vibramans) - self.max_raman_intensity) < 8

        # We used to test this, but it seems to vary wildly between
        # programs... perhaps we could use it if we knew why...
        # self.assertInside(data.vibramans[1], 2.6872, 0.0001)

    def testvibdisps(self, data) -> None:
        """Is the length and value of vibdisps correct?"""
        assert hasattr(data, "vibdisps")
        assert len(data.vibdisps) == 54


class DALTONRamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    max_raman_intensity = 745


class GaussianRamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    max_raman_intensity = 1066


class OrcaRamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    max_raman_intensity = 1045


class Orca6RamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    # This value has changed again in Orca 6 for some reason...
    max_raman_intensity = 1037


class QChemRamanTest(GenericRamanTest):
    """Customized Raman unittest"""

    max_raman_intensity = 588
