# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test configuration interaction (CI) logfiles in cclib"""

from typing import TYPE_CHECKING

import numpy

if TYPE_CHECKING:
    from cclib.parser.data import ccData


class GenericCISTest:
    """Generic CIS(RHF)/STO-3G water unittest"""

    nstates = 5

    # First four singlet/triplet state excitation energies [wavenumber].
    # Based on output in GAMESS test.
    etenergies0 = numpy.array([98614.56, 114906.59, 127948.12, 146480.64])
    etenergies1 = numpy.array([82085.34, 98999.11, 104077.89, 113978.37])

    # First four singlet/triplet state excitation orbitals and coefficients.
    # Tuples: (from MO, to MO, coefficient) - don't need spin indices.
    # Based on output in GAMESS test (using the "dets" algorithm).
    # Note that only coefficients larger than 0.1 are included here, as
    #   the Gaussian test does not contain smaller ones.
    # The simple example of water should yield the same first 4 states in all programs.
    # Note that the GAMESS test "water_cis_dets" calculates also triplet states,
    #  but the resulting transition dipoles and oscillator strengths are not correct
    #  and this file is not tested here (although it probably should be).
    etsecs0 = [[(4, 5, -1.0)], [(4, 6, -1.0)], [(3, 5, 0.96689)], [(2, 5, 0.4407), (3, 6, 0.89770)]]
    etsecs1 = [
        [(4, 5, 1.0)],
        [(2, 6, -0.231), (3, 5, -0.9676)],
        [(4, 6, 1.0)],
        [(2, 5, -0.536), (3, 6, -0.843)],
    ]

    etsecs_precision = 0.0005

    def testetenergiesvalues(self, data: "ccData") -> None:
        """Are etenergies within 50 wavenumbers of the correct values?"""
        indices0 = [i for i in range(self.nstates) if data.etsyms[i][0] == "S"]
        indices1 = [i for i in range(self.nstates) if data.etsyms[i][0] == "T"]
        singlets = [data.etenergies[i] for i in indices0]
        triplets = [data.etenergies[i] for i in indices1]
        # All programs do singlets.
        singletdiff = singlets[:4] - self.etenergies0
        assert numpy.all(singletdiff < 50)
        # Not all programs do triplets (i.e. Jaguar).
        if len(triplets) >= 4:
            tripletdiff = triplets[:4] - self.etenergies1
            assert numpy.all(tripletdiff < 50)

    def testsecs(self, data: "ccData") -> None:
        """Is the sum of etsecs close to 1?"""
        etsec = data.etsecs[2]  # Pick one with several contributors
        sumofsec = sum([z * z for (x, y, z) in etsec])
        assert abs(sumofsec - 1.0) < 0.02

    def testetsecsvalues(self, data: "ccData") -> None:
        """Are etsecs correct and coefficients close to the correct values?"""
        indices0 = [i for i in range(self.nstates) if data.etsyms[i][0] == "S"]
        indices1 = [i for i in range(self.nstates) if data.etsyms[i][0] == "T"]
        singlets = [data.etsecs[i] for i in indices0]
        triplets = [data.etsecs[i] for i in indices1]
        # All programs do singlets.
        for i in range(4):
            for exc in self.etsecs0[i]:
                found = False
                for s in singlets[i]:
                    if s[0][0] == exc[0] and s[1][0] == exc[1]:
                        found = True
                        assert abs(abs(s[2]) - abs(exc[2])) < self.etsecs_precision
                assert found, (
                    f"Excitation {int(exc[0])}->{exc[1]} not found (singlet state {int(i)})"
                )
        # Not all programs do triplets (i.e. Jaguar).
        if len(triplets) >= 4:
            for i in range(4):
                for exc in self.etsecs1[i]:
                    found = False
                    for s in triplets[i]:
                        if s[0][0] == exc[0] and s[1][0] == exc[1]:
                            found = True
                            assert abs(abs(s[2]) - abs(exc[2])) < self.etsecs_precision
                    assert found, (
                        f"Excitation {int(exc[0])}->{exc[1]} not found (triplet state {int(i)})"
                    )


class GAMESSCISTest(GenericCISTest):
    """Customized CIS(RHF)/STO-3G water unittest"""

    def testnocoeffs(self, data: "ccData") -> None:
        """Are natural orbital coefficients the right size?"""
        assert data.nocoeffs.shape == (data.nmo, data.nbasis)

    def testnooccnos(self, data: "ccData") -> None:
        """Are natural orbital occupation numbers the right size?"""
        assert data.nooccnos.shape == (data.nmo,)


class GaussianCISTest(GenericCISTest):
    """Customized CIS(RHF)/STO-3G water unittest"""

    nstates = 10

    def testnocoeffs(self, data: "ccData") -> None:
        """Are natural orbital coefficients the right size?"""
        assert data.nocoeffs.shape == (data.nmo, data.nbasis)

    def testnooccnos(self, data: "ccData") -> None:
        """Are natural orbital occupation numbers the right size?"""
        assert data.nooccnos.shape == (data.nmo,)


class Jaguar83CISTest(GenericCISTest):
    """Customized CIS(RHF)/STO-3G water unittest"""

    # The Jaguar8.3 job was created using 6-31G instead of STO-3G.
    etsecs0 = [
        [(4, 5, 0.99186)],
        [(4, 6, -0.98594)],
        [(3, 5, -0.98321)],
        [(2, 5, 0.19240), (3, 6, -0.97090)],
    ]
    etsecs1 = [
        [(4, 5, 1.0)],
        [(2, 6, -0.231), (3, 5, -0.9676)],
        [(4, 6, 1.0)],
        [(2, 5, -0.536), (3, 6, -0.843)],
    ]


class QChemCISTest(GenericCISTest):
    """Customized CIS(RHF)/STO-3G water unittest"""

    nstates = 10
