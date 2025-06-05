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
from cclib.parser import DALTON, GAMESS, Gaussian, Molcas, QChem, ccData, utils

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

    def test_center_of_mass(self) -> None:
        data, logfile = getdatafile(DALTON, "basicDALTON-2015", ["dvb_sp_hf.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        np.testing.assert_allclose(
            nuclear.center_of_mass(), np.array([0.0, 0.0, 0.0]), rtol=0.0, atol=1.0e-13
        )

    def test_moment_of_inertia_tensor_molcas(self) -> None:
        """Testing the full moment of inertia tensor for a Molcas output."""

        # https://gitlab.com/Molcas/OpenMolcas/-/blob/805c93ebfa2cd5fa185ca03cf965580ebdfeb352/src/system_util/constants.F90
        # CONV_AMU_TO_AU_ -- Convert 1AMU to au.
        # CONST_AMU_IN_SI_ -- Atomic mass unit (1/12*m[C12]) in SI units.
        # CONST_ELECTRON_MASS_IN_SI_ -- Mass of the electron in SI units.

        # For this particular file.
        CODATA_SET = 2014

        if CODATA_SET == 2022:
            CONST_AMU_IN_SI_ = 1.66053906892e-27
            CONST_ELECTRON_MASS_IN_SI_ = 9.1093837139e-31
        elif CODATA_SET == 2018:
            CONST_AMU_IN_SI_ = 1.66053906660e-27
            CONST_ELECTRON_MASS_IN_SI_ = 9.1093837015e-31
        elif CODATA_SET == 2014:
            CONST_AMU_IN_SI_ = 1.660539040e-27
            CONST_ELECTRON_MASS_IN_SI_ = 9.10938356e-31
        elif CODATA_SET == 2010:
            CONST_AMU_IN_SI_ = 1.660538921e-27
            CONST_ELECTRON_MASS_IN_SI_ = 9.10938291e-31
        CONV_AMU_TO_AU_ = CONST_AMU_IN_SI_ / CONST_ELECTRON_MASS_IN_SI_
        uToau = CONV_AMU_TO_AU_

        # The lines before
        # https://gitlab.com/Molcas/OpenMolcas/-/blob/805c93ebfa2cd5fa185ca03cf965580ebdfeb352/src/gateway_util/rigrot.F90#L127
        # and the "Cartesian Coordinates / Bohr, Angstrom" part of the output
        # file indicate that `Coor` is in units of bohr and `rM` is in units
        # of au.  The section "Coordinates and Masses of Atoms, in au and A"
        # has masses printed in amu.

        data, logfile = getdatafile(Molcas, "basicOpenMolcas18.0", ["dvb_sp.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        ref_moi_tensor = np.zeros(shape=(3, 3))

        with open(logfile.filename) as inputfile:
            for line in inputfile:
                if line.strip() == "The Moment of Inertia Tensor / au":
                    line = next(inputfile)
                    for i in range(3):
                        line = next(inputfile)
                        ref_moi_tensor[i, : i + 1] = [utils.float(x) for x in line.split()[1:]]

        ref_moi_tensor = utils.symmetrize(ref_moi_tensor, use_triangle="lower")
        # TODO why the sign difference for some elements?
        np.testing.assert_allclose(
            abs(nuclear.moment_of_inertia_tensor(units="amu_bohr_2")),
            ref_moi_tensor / uToau,
            rtol=1.2e-4,
        )

    def test_principal_moments_of_inertia_dalton(self) -> None:
        """Testing principal moments of inertia and the principal axes for a
        DALTON output.
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
        np.testing.assert_allclose(pmoi, ref_pmoi)
        # The phases of the eigenvectors may be different, but they
        # are still orthonormal within each set.
        np.testing.assert_allclose(np.abs(axes), np.abs(ref_axes), rtol=5.4e-7)

    def test_principal_moments_of_inertia_gamessus(self) -> None:
        """Testing principal moments of inertia and the principal axes for a
        GAMESS-US output.
        """

        data, _ = getdatafile(GAMESS, "basicGAMESS-US2018", ["dvb_ir.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        # line 174
        ref_amu_angstrom_2 = [109.440, 736.960, 846.400]
        # line 5085
        ref_amu_bohr_2 = [390.81900, 2631.73138, 3022.55038]

        pmoi_amu_angstrom_2, _ = nuclear.principal_moments_of_inertia("amu_angstrom_2")
        pmoi_amu_bohr_2, _ = nuclear.principal_moments_of_inertia("amu_bohr_2")

        np.testing.assert_allclose(pmoi_amu_angstrom_2, ref_amu_angstrom_2, rtol=4.3e-6)
        np.testing.assert_allclose(pmoi_amu_bohr_2, ref_amu_bohr_2, rtol=1.3e-6)

    def test_principal_moments_of_inertia_gaussian(self) -> None:
        """Testing principal moments of inertia and the principal axes for a
        Gaussian output.
        """

        data, _ = getdatafile(Gaussian, "basicGaussian16", ["dvb_ir.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        ref_pmoi = [390.07633792, 2635.01850639, 3025.09484431]
        pmoi, _ = nuclear.principal_moments_of_inertia("amu_bohr_2")
        np.testing.assert_allclose(pmoi, ref_pmoi)

    def test_principal_moments_of_inertia_molcas(self) -> None:
        """Testing principal moments of inertia and the principal axes for a
        Molcas output.
        """

        data, _ = getdatafile(Molcas, "basicOpenMolcas18.0", ["dvb_sp.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        # CODATA 2014 from OpenMolcas
        CONST_AMU_IN_SI_ = 1.660539040e-27
        CONST_ELECTRON_MASS_IN_SI_ = 9.10938356e-31
        CONV_AMU_TO_AU_ = CONST_AMU_IN_SI_ / CONST_ELECTRON_MASS_IN_SI_

        # line 488
        ref_pmoi = sorted(np.asarray([0.5368e07, 0.4655e07, 0.7127e06]) / CONV_AMU_TO_AU_)
        pmoi, _ = nuclear.principal_moments_of_inertia("amu_bohr_2")
        np.testing.assert_allclose(pmoi, ref_pmoi, rtol=5.1e-5)

    def test_principal_moments_of_inertia_qchem(self) -> None:
        """Testing principal moments of inertia and the principal axes for a
        Q-Chem output.
        """

        data, _ = getdatafile(QChem, "basicQChem5.4", ["dvb_ir.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        # line 936
        ref_pmoi = [390.81956, 2631.73358, 3022.55314]
        pmoi, _ = nuclear.principal_moments_of_inertia("amu_bohr_2")
        np.testing.assert_allclose(pmoi, ref_pmoi, rtol=1.3e-6)

    def test_rotational_constants_dalton(self) -> None:
        """Testing rotational constants for a DALTON output."""

        data, logfile = getdatafile(DALTON, "basicDALTON-2015", ["dvb_sp_hf.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        # line 642
        ref_mhz = [4617.8434, 685.7618, 597.0921]
        # line 643
        ref_invcm = [0.154035, 0.022875, 0.019917]
        rotconsts_ghz = nuclear.rotational_constants("ghz")
        rotconsts_invcm = nuclear.rotational_constants("invcm")
        np.testing.assert_allclose(rotconsts_ghz * 1.0e3, ref_mhz, rtol=2.0e-5)
        np.testing.assert_allclose(rotconsts_invcm, ref_invcm, rtol=2.2e-3)

    def test_rotational_constants_gamess(self) -> None:
        """Testing rotational constants for a GAMESS-US output."""

        data, _ = getdatafile(GAMESS, "basicGAMESS-US2018", ["dvb_ir.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        # line 5088
        ref_ghz = [4.61361, 0.68513, 0.59655]
        rotconsts_ghz = nuclear.rotational_constants("ghz")
        # should be 4 but lower likely due to masses
        np.testing.assert_allclose(rotconsts_ghz, ref_ghz, rtol=9.3e-4)

    def test_rotational_constants_gaussian(self) -> None:
        """Testing rotational constants for a Gaussian output."""

        data, logfile = getdatafile(Gaussian, "basicGaussian16", ["dvb_sp.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        # line 362
        ref_ghz = [4.6266363, 0.6849065, 0.5965900]
        rotconsts_ghz = nuclear.rotational_constants("ghz")
        np.testing.assert_allclose(rotconsts_ghz, ref_ghz)

    def test_rotational_constants_molcas(self) -> None:
        """Testing rotational constants for a Molcas output."""

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
        np.testing.assert_allclose(rotconsts_ghz, ref_ghz, rtol=4.1e-5)
        np.testing.assert_allclose(rotconsts_invcm, ref_invcm, rtol=2.2e-3)
