# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the DDEC6 in cclib"""

import os
import sys
from typing import Optional

from cclib.io import ccread
from cclib.method import DDEC6, volume
from cclib.method.calculationmethod import MissingAttributeError
from cclib.parser import Psi4

import numpy
import pytest
from numpy.testing import assert_allclose

from ..test_data import getdatafile


class DDEC6Test:
    """DDEC6 method tests."""

    def setup_method(self) -> None:
        self.parse()

    def parse(self, molecule_name: Optional[str] = None) -> None:
        if molecule_name is None:
            self.data, self.logfile = getdatafile(Psi4, "basicPsi4-1.2.1", ["water_mp2.out"])
        else:
            self.data = ccread(
                os.path.join(os.path.dirname(os.path.realpath(__file__)), f"{molecule_name}.out")
            )

    def testmissingrequiredattributes(self) -> None:
        """Is an error raised when required attributes are missing?"""
        for missing_attribute in DDEC6.required_attrs:
            self.parse()
            vol = volume.Volume((-4, -4, -4), (4, 4, 4), (0.2, 0.2, 0.2))
            delattr(self.data, missing_attribute)
            with pytest.raises(MissingAttributeError):
                trial = DDEC6(self.data, vol, os.path.dirname(os.path.realpath(__file__)))  # noqa: F841

    def test_proatom_read(self) -> None:
        """Are proatom densities imported correctly?"""

        self.parse()
        vol = volume.Volume((-4, -4, -4), (4, 4, 4), (0.2, 0.2, 0.2))

        self.analysis = DDEC6(self.data, vol, os.path.dirname(os.path.realpath(__file__)))

        refH_den = [
            2.66407645e-01,
            2.66407645e-01,
            2.66407643e-01,
            2.66407612e-01,
            2.66407322e-01,
        ]  # Hydrogen first five densities
        refH_r = [  # noqa: F841
            1.17745807e-07,
            4.05209491e-06,
            3.21078677e-05,
            1.39448474e-04,
            4.35643929e-04,
        ]  # Hydrogen first five radii
        refO_den = [
            2.98258510e02,
            2.98258510e02,
            2.98258509e02,
            2.98258487e02,
            2.98258290e02,
        ]  # Oxygen first five densities
        refO_r = [  # noqa: F841
            5.70916728e-09,
            1.97130512e-07,
            1.56506399e-06,
            6.80667366e-06,
            2.12872046e-05,
        ]  # Oxygen first five radii

        assert_allclose(self.analysis.proatom_density[0][0:5], refO_den, rtol=1e-3)
        assert_allclose(self.analysis.proatom_density[1][0:5], refH_den, rtol=1e-3)
        assert_allclose(self.analysis.proatom_density[2][0:5], refH_den, rtol=1e-3)

    def test_water_charges(self) -> None:
        """Are charges and quantities in each step of DDEC6 algorithm calculated correctly
        for water?

        Here, values are compared against `chargemol` calculations.
        Due to the differences in basis set used for calculation and slightly different integration
        grid, some discrepancy is inevitable in the comparison.
        """

        self.parse()
        # use precalculated fine cube file
        imported_vol = volume.read_from_cube(
            os.path.join(os.path.dirname(os.path.realpath(__file__)), "water_fine.cube")
        )

        analysis = DDEC6(self.data, imported_vol, os.path.dirname(os.path.realpath(__file__)))
        analysis.calculate()

        radial_indices = []
        for atomi in range(len(self.data.atomnos)):
            lst = []
            for radius in [0.05, 0.10, 0.15, 0.20, 0.25]:
                # find closest radius index
                lst.append(numpy.abs(analysis.radial_grid_r[atomi] - radius).argmin())
            radial_indices.append(lst)

        # values from `chargemol` calculation
        # which is based on proatomic densities calculated with different basis set.
        # discrepancy comes from the fact that `chargemol` grid & `horton` grid don't exactly match
        # (rtol is adjusted to account for this inevitable discrepancy)
        # STEP 1
        # Check assigned charges.
        assert_allclose(analysis.reference_charges[0], [-0.513006, 0.256231, 0.256775], rtol=0.10)
        # STEP 2
        # Check assigned charges.
        assert_allclose(analysis.reference_charges[1], [-0.831591, 0.415430, 0.416161], rtol=0.20)
        # STEP 3
        # Check integrated charge density (rho^cond(r)) on grid with integrated values (=nelec).
        assert abs(analysis.charge_density.integrate() - analysis.rho_cond.integrate()) < 1
        for atomi in range(len(analysis.data.atomnos)):
            assert (
                abs(
                    analysis._integrate_from_radial([analysis._cond_density[atomi]], [atomi])
                    + analysis.reference_charges[-1][atomi]
                    - analysis.data.atomnos[atomi]
                )
                < 0.5
            )
        # Also compare with data from `chargemol`
        # discrepancy comes from the fact that `chargemol` grid and `horton` grid do not exactly match
        assert_allclose(
            analysis.tau[0][radial_indices[0]],
            [0.999846160, 0.999739647, 0.999114037, 0.997077942, 0.994510889],
            rtol=0.10,
        )
        assert_allclose(
            analysis.tau[1][radial_indices[1]],
            [0.864765882, 0.848824620, 0.805562019, 0.760402501, 0.736949861],
            rtol=0.10,
        )
        assert_allclose(
            analysis.tau[2][radial_indices[2]],
            [0.845934391, 0.839099407, 0.803699493, 0.778428137, 0.698628724],
            rtol=0.10,
        )
        # STEP 4-7
        # Check values assigned to u_A
        assert_allclose(
            analysis.u_A,
            [
                [0.572349429, 0.296923935, 0.296520531],
                [0.563154399, 0.291919678, 0.291376710],
                [0.563475132, 0.292007655, 0.291508794],
                [0.565816045, 0.293131322, 0.292902112],
            ],
            atol=0.05,
        )
        # Check assigned charges
        assert_allclose(analysis.fragcharges, [-0.757097, 0.378410, 0.378687], atol=0.2)

    @pytest.mark.skipif(
        sys.version_info > (3, 8),
        reason="This test doesn't converge with newer psi4 versions availiable with python >3.8",
    )
    def test_chgsum_h2(self) -> None:
        """Are DDEC6 charges for hydrogen atoms in nonpolar H2 small as expected?

        Using much denser grid (spacing of 0.1 rather than 0.2 which is the cube file included
        in the test) gives [0.00046066, 0.00046066].
        """

        self.parse("h2")
        vol = volume.Volume((-2, -2, -2), (2, 2, 2), (0.2, 0.2, 0.2))
        analysis = DDEC6(self.data, vol, os.path.dirname(os.path.realpath(__file__)))
        analysis.calculate()

        assert abs(analysis.fragcharges[0] - analysis.fragcharges[1]) < 1e-12

    def test_chgsum_co(self) -> None:
        """Are DDEC6 charges for carbon monoxide reported as expected?

        Deviation from a total of zero (-0.00682) occurs because the integrated value of total
        density (14.006876594937234) is slightly larger than # of electrons.

        Using a finer grid reduces this discrepancy.
        """

        self.parse("co")
        imported_vol = volume.read_from_cube(
            os.path.join(os.path.dirname(os.path.realpath(__file__)), "co.cube")
        )
        analysis = DDEC6(self.data, imported_vol, os.path.dirname(os.path.realpath(__file__)))
        analysis.calculate()

        assert abs(numpy.sum(analysis.fragcharges) - 0) < 1e-2
        assert_allclose(analysis.fragcharges, [0.13221636, -0.13903595], atol=1e-3)

    def test_chg_nh3(self) -> None:
        """Are DDEC6 charges for ammonia reported as expected?

        Deviation from a total of zero (0.026545) occurs because the integrated value of total
        density (9.973453129261163) is slightly smaller than number of electrons.

        Using a finer grid reduces this discrepancy.
        """

        self.parse("nh3")
        imported_vol = volume.read_from_cube(
            os.path.join(os.path.dirname(os.path.realpath(__file__)), "nh3.cube")
        )
        analysis = DDEC6(self.data, imported_vol, os.path.dirname(os.path.realpath(__file__)))
        analysis.calculate()

        assert_allclose(
            analysis.fragcharges, [-0.7824003, 0.26854388, 0.26959206, 0.27081123], atol=1e-3
        )
