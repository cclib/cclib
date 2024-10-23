# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test single point time-dependent logfiles in cclib"""

from cclib.parser import utils

import numpy
from skip import skipForLogfile, skipForParser


class GenericTDTest:
    """Generic time-dependent HF/DFT unittest"""

    number = 5
    # ???
    expected_l_max = 0.18680974536292055
    expected_f_max = 0.67
    expected_f_max_thresh = 0.1
    symmetries = ["Singlet-Bu", "Singlet-Bu", "Singlet-Ag", "Singlet-Bu", "Singlet-Ag"]
    sumofsec = 1.0
    method = "TD-DFT"

    @skipForParser("ADF", "excited_states_method not yet implemented")
    @skipForParser("DALTON", "excited_states_method not yet implemented")
    @skipForParser("FChk", "excited_states_method not yet implemented")
    @skipForParser("GAMESS", "excited_states_method not yet implemented")
    @skipForParser("GAMESSUK", "excited_states_method not yet implemented")
    @skipForParser("Jaguar", "excited_states_method not yet implemented")
    @skipForParser("NWChem", "excited_states_method not yet implemented")
    @skipForParser("QChem", "excited_states_method not yet implemented")
    def testmetadata(self, data) -> None:
        """Did we parse an excited states method?"""
        assert data.metadata["excited_states_method"] == self.method

    @skipForLogfile(
        "Turbomole/basicTurbomole7.4/CO_cc2_TD_trip",
        "Oscillator strengths are not available for Turbomole triplets using ricc2 but are required for testenergies()",
    )
    def testenergies(self, data) -> None:
        """Is the l_max reasonable?"""

        assert len(data.etenergies) == self.number

        # Note that if all oscillator strengths are zero (like for triplets)
        # then this will simply pick out the first energy.
        idx_lambdamax = numpy.argmax(data.etoscs)
        assert abs(
            data.etenergies[idx_lambdamax]
            - utils.convertor(self.expected_l_max, "hartree", "wavenumber")
        ) < utils.convertor(0.022781676263770798, "hartree", "wavenumber")

    @skipForLogfile(
        "Turbomole/basicTurbomole7.4/CO_cc2_TD_trip",
        "Oscillator strengths are not available for triplets with Turbomole's ricc2",
    )
    def testoscs(self, data) -> None:
        """Is the maximum of etoscs in the right range?"""
        assert len(data.etoscs) == self.number
        assert abs(max(data.etoscs) - self.expected_f_max) < self.expected_f_max_thresh

    @skipForParser("FChk", "The parser is still being developed so we skip this test")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForLogfile(
        "Gaussian/basicGaussian16/dvb_eomccsd.log", "Transitions are not yet parsed for EOM-CCSD"
    )
    def testsecs(self, data) -> None:
        """Is the sum of etsecs close to 1?"""
        assert len(data.etsecs) == self.number
        lowestEtrans = data.etsecs[numpy.argmin(data.etenergies)]
        sumofsec = sum([z * z for (x, y, z) in lowestEtrans])
        assert abs(sumofsec - self.sumofsec) < 0.16

    @skipForParser("FChk", "This is true for calculations without symmetry, but not with?")
    @skipForParser("DALTON", "This is true for calculations without symmetry, but not with?")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForLogfile(
        "Gaussian/basicGaussian16/dvb_eomccsd.log", "Transitions are not yet parsed for EOM-CCSD"
    )
    def testsecs_transition(self, data) -> None:
        """Is the lowest E transition from the HOMO or to the LUMO?"""
        lowestEtrans = data.etsecs[numpy.argmin(data.etenergies)]
        t = list(reversed(sorted([(c * c, s, e) for (s, e, c) in lowestEtrans])))
        assert t[0][1][0] == data.homos[0] or t[0][2][0] == data.homos[0] + 1

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForLogfile("FChk/basicQChem5.4", "etsyms are not yet implemented")
    @skipForLogfile("ORCA/basicORCA5.0/dvb_adc2.log", "etsyms are not available for this method")
    @skipForLogfile(
        "ORCA/basicORCA5.0/dvb_eom_ccsd.log", "etsyms are not available for this method"
    )
    @skipForLogfile(
        "ORCA/basicORCA5.0/dvb_pno_eom_ccsd.log", "etsyms are not available for this method"
    )
    def testsymsnumber(self, data) -> None:
        """Is the length of etsyms correct?"""
        assert len(data.etsyms) == self.number

    @skipForParser("ADF", "etrotats are not yet implemented")
    @skipForParser("DALTON", "etsyms are not yet implemented")
    @skipForParser("FChk", "etsyms are not yet implemented")
    @skipForParser("GAMESS", "etsyms are not yet implemented")
    @skipForParser("GAMESSUK", "etsyms are not yet implemented")
    @skipForParser("Jaguar", "etsyms are not yet implemented")
    @skipForParser("NWChem", "etsyms are not yet implemented")
    @skipForParser("QChem", "etsyms are not yet implemented")
    @skipForLogfile("ORCA/basicORCA4.2", "etsyms are only available in ORCA >= 5.0")
    @skipForLogfile("ORCA/basicORCA4.1", "etsyms are only available in ORCA >= 5.0")
    @skipForLogfile("Gaussian/basicGaussian09", "symmetry is missing for this log file")
    @skipForLogfile("FChk/basicQChem5.4", "etsyms are not yet implemented")
    def testsyms(self, data) -> None:
        """Are the values of etsyms correct?"""
        assert data.etsyms == self.symmetries

    @skipForParser("ADF", "etrotats are not yet implemented")
    @skipForParser("DALTON", "etrotats are not yet implemented")
    @skipForParser("FChk", "etrotats are not yet implemented")
    @skipForParser("GAMESS", "etrotats are not yet implemented")
    @skipForParser("GAMESSUK", "etrotats are not yet implemented")
    @skipForParser("Jaguar", "etrotats are not yet implemented")
    @skipForParser("NWChem", "etrotats are not yet implemented")
    @skipForParser("QChem", "Q-Chem cannot calculate rotatory strengths")
    @skipForLogfile("FChk/basicQChem5.4", "Q-Chem cannot calculate rotatory strengths")
    @skipForLogfile(
        "Turbomole/basicTurbomole7.4/CO_cc2_TD",
        "Rotatory strengths are not currently available for ricc2",
    )
    @skipForLogfile(
        "Turbomole/basicTurbomole7.4/CO_adc2_TD",
        "Rotatory strengths are not currently available for ricc2",
    )
    @skipForLogfile(
        "ORCA/basicORCA6.0/dvb_eom_ccsd.log", "etrotats don't seem to be available for EOM-CCSD"
    )
    def testrotatsnumber(self, data) -> None:
        """Is the length of etrotats correct?"""
        assert len(data.etrotats) == self.number

    @skipForParser("ADF", "optstate is not yet implemented")
    @skipForParser("DALTON", "optstate are not yet implemented")
    @skipForParser("FChk", "optstate are not yet implemented")
    @skipForParser("GAMESS", "optstate are not yet implemented")
    @skipForParser("GAMESSUK", "optstate are not yet implemented")
    @skipForParser("Jaguar", "optstate are not yet implemented")
    @skipForParser("NWChem", "optstate are not yet implemented")
    @skipForParser("ORCA", "optstate are not yet implemented")
    @skipForParser("QChem", "optstate are not yet implemented")
    @skipForParser("Turbomole", "optstate are not yet implemented")
    def testoptstate(self, data) -> None:
        # All our examples have a default state-of-interest of 1 (index 0).
        assert data.metadata["opt_state"] == 0


class ADFTDDFTTest(GenericTDTest):
    """Customized time-dependent DFT unittest"""

    number = 5

    def testsecs(self, data) -> None:
        """Is the sum of etsecs close to 1?"""
        assert len(data.etsecs) == self.number
        lowestEtrans = data.etsecs[1]

        # ADF squares the etsecs
        sumofsec = sum([z for (x, y, z) in lowestEtrans])
        assert abs(sumofsec - 1.0) < 0.16


class DALTONTDTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    number = 20


class GaussianTDDFTTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    expected_l_max = 0.21870409213219966

    @skipForLogfile(
        "FChk/basicGaussian09", "etrotats are not available in fchk, only the main logfile"
    )
    @skipForLogfile(
        "FChk/basicGaussian16", "etrotats are not available in fchk, only the main logfile"
    )
    def testrotatsnumber(self, data) -> None:
        """Is the length of etrotats correct?"""
        assert len(data.etrotats) == self.number

    def testetdipsshape(self, data) -> None:
        """Is the shape of etdips correct?"""
        assert numpy.shape(data.etdips) == (self.number, 3)

    def testetveldipsshape(self, data) -> None:
        """Is the shape of etveldips correct?"""
        assert numpy.shape(data.etveldips) == (self.number, 3)

    def testetmagdipsshape(self, data) -> None:
        """Is the shape of etmagdips correct?"""
        assert numpy.shape(data.etmagdips) == (self.number, 3)


class GAMESSUSTDDFTTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    number = 10


class JaguarTDDFTTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    expected_l_max = 0.21870409213219966
    expected_f_max = 1.2


class OrcaTDDFTTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    number = 10
    expected_l_max = 0.21870409213219966
    symmetries = [
        "Triplet-Bu",
        "Triplet-Ag",
        "Triplet-Bu",
        "Triplet-Bu",
        "Triplet-Bu",
        "Singlet-Bu",
        "Singlet-Bu",
        "Singlet-Ag",
        "Singlet-Bu",
        "Singlet-Ag",
    ]
    method = "TDA"

    def testoscs(self, data) -> None:
        """Is the maximum of etoscs in the right range?"""
        assert len(data.etoscs) == self.number
        assert abs(max(data.etoscs) - 1.17) < 0.01


class QChemTDDFTTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    number = 10
    expected_l_max = 0.21870409213219966
    expected_f_max = 0.9


class GenericTDDFTtrpTest(GenericTDTest):
    """Generic time-dependent HF/DFT (triplet) unittest"""

    number = 5
    expected_l_max = 0.11163021369247692

    def testoscs(self, data) -> None:
        """Triplet excitations should be disallowed."""
        assert len(data.etoscs) == self.number
        assert abs(max(data.etoscs)) < 0.01


class OrcaROCISTest(GenericTDTest):
    """Customized test for ROCIS"""

    number = 4
    expected_l_max = 10.55689555
    expected_f_max = 0.015
    # per 1085, no VELOCITY DIPOLE MOMENTS are parsed
    n_spectra = 7

    # Do we want to parse ROCIS as its own method?
    method = "CIS"

    def testTransprop(self, data) -> None:
        """Check the number of spectra parsed"""
        assert len(data.transprop) == self.n_spectra
        tddft_length = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
        assert tddft_length in data.transprop

    def testsymsnumber(self, data) -> None:
        """ORCA ROCIS has no symmetry"""
        pass

    def testsecs(self, data) -> None:
        """ROCIS does not form singly excited configurations (secs)"""
        pass

    def testsecs_transition(self, data) -> None:
        """ROCIS does not form singly excited configurations (secs)"""
        pass

    def testrotatsnumber(self, data) -> None:
        """ROCIS does not calculate rotatory strengths"""
        pass

    def testsyms(self, data) -> None:
        """ROCIS does not show symmetries"""
        pass


class Orca6ROCISTest(OrcaROCISTest):
    # The energies are different in 5 Vs 6 for some reason I can't fathom,
    # perhaps a setting has changed from some old default value?
    expected_l_max = 10.718302
    n_spectra = 4


class TurbomoleTDTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    number = 10
    expected_l_max = 0.4165966608093103
    expected_f_max = 0.19
    symmetries = ["Singlet-A"] * 10

    @skipForLogfile("Turbomole/basicTurbomole7.4/CO_cc2_TD", "There are no dipole moments in ricc2")
    def testetmagdipsshape(self, data) -> None:
        """Is the shape of etmagdips correct?"""
        assert numpy.shape(data.etmagdips) == (self.number, 3)


class TurbomoleTDADC2Test(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    number = 10
    expected_l_max = 0.6211579
    expected_f_max = 0.8
    symmetries = ["Singlet-A"] * 10
    method = "ADC(2)"


class TurbomoleTDCC2Test(TurbomoleTDTest):
    """Customized time-dependent HF/DFT unittest"""

    method = "CC2"


class TurbomoleTDTripTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    number = 10
    expected_l_max = 0.2347872184128124
    expected_f_max = 0.84
    symmetries = ["Triplet-A"] * 10
    method = "RPA"


class TurbomoleTDCC2TripTest(GenericTDTest):
    """Customized time-dependent HF/DFT unittest"""

    # This test is for triplets with ricc2, which does not support oscillator strengths.

    number = 10
    symmetries = ["Triplet-A"] * 10
    method = "CC2"

    def testenergies(self, data) -> None:
        """Is the l_max reasonable?"""
        assert len(data.etenergies) == self.number


class OrcaETPostHFTest(GenericTDTest):
    """Tests for post-HF excited states with ORCA."""

    number = 2
    symmetries = ["Singlet", "Singlet"]
    expected_f_max = 1.0
    expected_l_max = 0.27338011516524957
    # Not sure why this value != 1 for these methods?
    # Perhaps remaining contributions were too small to print?
    sumofsec = 0.43
    method = "EOM-CCSD"


class OrcaADC2Test(OrcaETPostHFTest):
    method = "ADC(2)"


class OrcaSTEOMCCSDTest(OrcaETPostHFTest):
    """Test for STEOM-CCSD with Orca."""

    sumofsec = 1.0
    method = "STEOM-CCSD"


class GaussianEOMCCSDTest(GenericTDTest):
    """Test for EOM-CCSD with Gaussian."""

    number = 10
    expected_l_max = 0.2802797426880702
    expected_f_max = 0.9802
    symmetries = [
        "Triplet-Bu",
        "Triplet-Ag",
        "Triplet-Bu",
        "Singlet-Bu",
        "Triplet-Bu",
        "Triplet-Bu",
        "Singlet-Bu",
        "Triplet-Ag",
        "Triplet-Bu",
        "Triplet-Ag",
    ]
    method = "EOM-CCSD"
