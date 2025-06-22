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
import pytest

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

        with open(logfile.filename) as inputfile:
            output = inputfile.read()
        line = re.search("Nuclear Repulsion Energy = .* hartrees", output).group()
        nre = float(line.split()[4])
        nre = utils.convertor(nre, "hartree", "eV")
        assert nuclear.repulsion_energy(atomcoords_index=0) == pytest.approx(nre)

    def test_center_of_mass_dalton(self) -> None:
        """Testing the center of mass calculation for logfile where it is
        printed.

        For the divinylbenzene molecule used in most of our test data, the
        center of mass should be at the origin (zero).
        """

        data, _ = getdatafile(DALTON, "basicDALTON-2015", ["dvb_sp_hf.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        np.testing.assert_allclose(
            nuclear.center_of_mass(), np.array([0.0, 0.0, 0.0]), rtol=0.0, atol=1.0e-13
        )

    def test_moment_of_inertia_tensor_molcas(self) -> None:
        """Testing the full moment of inertia tensor for a Molcas output."""

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

        amu_to_au = _openmolcas_mass_amu_to_au(2014)
        np.testing.assert_allclose(
            nuclear.moment_of_inertia_tensor(units="amu_bohr_2"),
            utils.symmetrize(ref_moi_tensor, use_triangle="lower") / amu_to_au,
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
        with open(logfile.filename) as inputfile:
            for line in inputfile:
                if line.strip() == "Principal moments of inertia (u*A**2) and principal axes":
                    next(inputfile)
                    next(inputfile)
                    for _ in range(3):
                        tokens = [float(x) for x in next(inputfile).split()[1:]]
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
        np.testing.assert_allclose(
            nuclear.principal_moments_of_inertia("amu_angstrom_2")[0],
            [109.440, 736.960, 846.400],
            rtol=4.3e-6,
        )
        # line 5085
        np.testing.assert_allclose(
            nuclear.principal_moments_of_inertia("amu_bohr_2")[0],
            [390.81900, 2631.73138, 3022.55038],
            rtol=1.3e-6,
        )

    def test_principal_moments_of_inertia_gaussian(self) -> None:
        """Testing principal moments of inertia and the principal axes for a
        Gaussian output.
        """

        data, _ = getdatafile(Gaussian, "basicGaussian16", ["dvb_ir.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        # line 2030
        pmoi, moi_axes = nuclear.principal_moments_of_inertia("amu_bohr_2")
        np.testing.assert_allclose(pmoi, [390.07633792, 2635.01850639, 3025.09484431])
        np.testing.assert_allclose(
            np.abs(moi_axes),
            np.abs(
                np.array(
                    [
                        [0.02413, 0.99971, 0.00000],
                        [0.99971, -0.02413, 0.00000],
                        [-0.00000, -0.00000, 1.00000],
                    ]
                )
            ),
            rtol=1.8e-4,
        )

    def test_principal_moments_of_inertia_molcas(self) -> None:
        """Testing principal moments of inertia and the principal axes for a
        Molcas output.
        """

        data, _ = getdatafile(Molcas, "basicOpenMolcas18.0", ["dvb_sp.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        amu_to_au = _openmolcas_mass_amu_to_au(2014)

        # line 488
        np.testing.assert_allclose(
            nuclear.principal_moments_of_inertia("amu_bohr_2")[0],
            np.asarray([0.7127e06, 0.4655e07, 0.5368e07]) / amu_to_au,
            rtol=5.1e-5,
        )

    def test_principal_moments_of_inertia_qchem(self) -> None:
        """Testing principal moments of inertia and the principal axes for a
        Q-Chem output.
        """

        data, _ = getdatafile(QChem, "basicQChem5.4", ["dvb_ir.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        # line 936
        np.testing.assert_allclose(
            nuclear.principal_moments_of_inertia("amu_bohr_2")[0],
            [390.81956, 2631.73358, 3022.55314],
            rtol=1.3e-6,
        )

    def test_rotational_constants_dalton(self) -> None:
        """Testing rotational constants for a DALTON output."""

        data, _ = getdatafile(DALTON, "basicDALTON-2015", ["dvb_sp_hf.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        # line 642
        np.testing.assert_allclose(
            nuclear.rotational_constants("ghz") * 1.0e3,
            [4617.8434, 685.7618, 597.0921],
            rtol=2.0e-5,
        )
        # line 643
        np.testing.assert_allclose(
            nuclear.rotational_constants("invcm"), [0.154035, 0.022875, 0.019917], rtol=2.2e-3
        )

    def test_rotational_constants_gamessus(self) -> None:
        """Testing rotational constants for a GAMESS-US output."""

        data, _ = getdatafile(GAMESS, "basicGAMESS-US2018", ["dvb_ir.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        # line 5088
        # should be 4 but lower likely due to masses
        np.testing.assert_allclose(
            nuclear.rotational_constants("ghz"), [4.61361, 0.68513, 0.59655], rtol=9.3e-4
        )

    def test_rotational_constants_gaussian(self) -> None:
        """Testing rotational constants for a Gaussian output."""

        data, _ = getdatafile(Gaussian, "basicGaussian16", ["dvb_sp.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        # line 362
        np.testing.assert_allclose(
            nuclear.rotational_constants("ghz"), [4.6266363, 0.6849065, 0.5965900]
        )

    def test_rotational_constants_molcas(self) -> None:
        """Testing rotational constants for a Molcas output."""

        data, _ = getdatafile(Molcas, "basicOpenMolcas18.0", ["dvb_ir.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        # line 12670
        np.testing.assert_allclose(
            nuclear.rotational_constants("ghz"), [4.6160, 0.7067, 0.6129], atol=2.9e-5
        )
        # line 12669
        np.testing.assert_allclose(
            nuclear.rotational_constants("invcm"), [0.1540, 0.0236, 0.0204], rtol=2.2e-3
        )


def _openmolcas_mass_amu_to_au(codata_set: int) -> float:
    # https://gitlab.com/Molcas/OpenMolcas/-/blob/805c93ebfa2cd5fa185ca03cf965580ebdfeb352/src/system_util/constants.F90
    # CONV_AMU_TO_AU_ -- Convert 1AMU to au.
    # CONST_AMU_IN_SI_ -- Atomic mass unit (1/12*m[C12]) in SI units.
    # CONST_ELECTRON_MASS_IN_SI_ -- Mass of the electron in SI units.
    if codata_set == 2022:
        const_amu_in_si_ = 1.66053906892e-27
        const_electron_mass_in_si_ = 9.1093837139e-31
    elif codata_set == 2018:
        const_amu_in_si_ = 1.66053906660e-27
        const_electron_mass_in_si_ = 9.1093837015e-31
    elif codata_set == 2014:
        const_amu_in_si_ = 1.660539040e-27
        const_electron_mass_in_si_ = 9.10938356e-31
    elif codata_set == 2010:
        const_amu_in_si_ = 1.660538921e-27
        const_electron_mass_in_si_ = 9.10938291e-31
    else:
        raise ValueError(f"Unknown value for argument codata_set: {codata_set}")
    return const_amu_in_si_ / const_electron_mass_in_si_
