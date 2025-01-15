# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the Nuclear method in cclib"""

import logging
import re
import sys
from typing import Sequence

from cclib.method import Nuclear
from cclib.parser import DALTON, Gaussian, Molcas, QChem, ccData, utils

import numpy as np

sys.path.insert(1, "..")

from ..test_data import getdatafile


class NuclearTest:
    def test_stoichiometry(self) -> None:
        """Testing stoichoimetry generation."""
        data = ccData()

        def check(atomnos: Sequence[int], formula: str, charge: int = 0) -> None:
            data.natom = len(atomnos)
            data.atomnos = np.array(atomnos)
            data.atomcoords = np.zeros((data.natom, 3))
            data.charge = charge
            assert Nuclear(data).stoichiometry() == formula

        # Basics and permutations.
        check([], "")
        check([6, 1, 6, 1, 1, 1], "C2H4")
        check([1, 1, 1, 6, 1, 6], "C2H4")

        # Charges.
        check([8], "O", charge=0)
        check([8], "O(+1)", charge=1)
        check([8], "O(-1)", charge=-1)
        check([8], "O(+2)", charge=2)
        check([8], "O(+9)", charge=9)

        # Element counts.
        check([6, 1], "CH")
        check([6] * 60, "C60")

        # Test the Hill system.
        check([8, 1, 1], "H2O")
        check([6, 8, 8, 1, 1], "CH2O2")
        check([16, 16, 8, 8], "O2S2")

    def test_repulsion_energy(self) -> None:
        """Testing nuclear repulsion energy for one logfile where it is printed."""

        data, logfile = getdatafile(QChem, "basicQChem5.4", ["water_mp4sdq.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        with open(logfile.filename) as f:
            output = f.read()
        line = re.search("Nuclear Repulsion Energy = .* hartrees", output).group()
        nre = float(line.split()[4])
        nre = utils.convertor(nre, "hartree", "eV")
        assert round(abs(nuclear.repulsion_energy() - nre), 5) == 0

    def test_principal_moments_of_inertia(self) -> None:
        """Testing principal moments of inertia and the principal axes for one
        logfile where it is printed.
        """

        data, logfile = getdatafile(DALTON, "basicDALTON-2015", ["dvb_sp_hf.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        ref_pmoi = []
        ref_axes = []
        with open(logfile.filename) as f:
            for line in f:
                if line.strip() == "Principal moments of inertia (u*A**2) and principal axes":
                    next(f)
                    next(f)
                    for _ in range(3):
                        tokens = [float(x) for x in next(f).split()[1:]]
                        ref_pmoi.append(tokens[0])
                        ref_axes.append(tokens[1:])
        pmoi, axes = nuclear.principal_moments_of_inertia("amu_angstrom_2")
        np.testing.assert_allclose(pmoi, ref_pmoi, rtol=0, atol=1.0e-4)
        # The phases of the eigenvectors may be different, but they
        # are still orthonormal within each set.
        np.testing.assert_allclose(np.abs(axes), np.abs(ref_axes), rtol=0, atol=1.0e-4)

    def test_principal_moments_of_inertia_2(self) -> None:
        """Testing principal moments of inertia and the principal axes for one
        logfile where it is printed.

        This test was added as a follow-up to PR #790.
        """
        data, _ = getdatafile(Gaussian, "basicGaussian16", ["dvb_ir.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        ref_pmoi = [390.07633792, 2635.01850639, 3025.09484431]
        pmoi, _ = nuclear.principal_moments_of_inertia("amu_bohr_2")
        np.testing.assert_allclose(pmoi, ref_pmoi, rtol=0, atol=1.0e-4)

    def test_rotational_constants(self) -> None:
        """Testing rotational constants for logfiles where they are printed."""

        data, logfile = getdatafile(DALTON, "basicDALTON-2015", ["dvb_sp_hf.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        ref_mhz = [0.0 for _ in range(3)]
        ref_ghz = [0.0 for _ in range(3)]
        ref_invcm = [0.0 for _ in range(3)]

        with open(logfile.filename) as f:
            for line in f:
                if line.strip() == "Rotational constants":
                    while line.split() != ["A", "B", "C"]:
                        line = next(f)
                    line = next(f)
                    ref_mhz = [float(x) for x in next(f).split()[:-1]]
                    ref_invcm = [float(x) for x in next(f).split()[:-1]]
                    break
        rotconsts_ghz = nuclear.rotational_constants("ghz")
        rotconsts_invcm = nuclear.rotational_constants("invcm")
        np.testing.assert_allclose(rotconsts_ghz * 1.0e3, ref_mhz, rtol=0, atol=1.0e-4)
        np.testing.assert_allclose(rotconsts_invcm, ref_invcm, rtol=0, atol=1.0e-4)

        data, logfile = getdatafile(Gaussian, "basicGaussian16", ["dvb_sp.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        with open(logfile.filename) as f:
            for line in f:
                if "Rotational constants (GHZ):" in line:
                    ref_ghz = [float(x) for x in line.split()[3:]]
                    break
        rotconsts_ghz = nuclear.rotational_constants("ghz")
        np.testing.assert_allclose(rotconsts_ghz, ref_ghz, rtol=0, atol=1.0e-5)

        data, logfile = getdatafile(Molcas, "basicOpenMolcas18.0", ["dvb_ir.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        with open(logfile.filename) as f:
            for line in f:
                if line.startswith(" Rotational Constants (cm-1)"):
                    # sort because they don't seem to come out in the same
                    # order, and a moment of inertia tensor isn't printed for
                    # comparison
                    ref_invcm = sorted([float(x) for x in line.split()[-3:]])
                    line = next(f)
                    ref_ghz = sorted([float(x) for x in line.split()[-3:]])
        rotconsts_ghz = sorted(nuclear.rotational_constants("ghz"))
        rotconsts_invcm = sorted(nuclear.rotational_constants("invcm"))
        np.testing.assert_allclose(rotconsts_ghz, ref_ghz, rtol=0, atol=1.0e-4)
        np.testing.assert_allclose(rotconsts_invcm, ref_invcm, rtol=0, atol=1.0e-4)
