# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the Hirshfeld Method in cclib"""

import os

from cclib.io import ccread
from cclib.method import Hirshfeld, volume
from cclib.method.calculationmethod import MissingAttributeError
from cclib.parser import Psi4

import numpy
import pytest
from numpy.testing import assert_allclose

from ..test_data import getdatafile


class HirshfeldTest:
    """Hirshfeld method tests."""

    def setup_method(self) -> None:
        self.parse()

    def parse(self) -> None:
        self.data, self.logfile = getdatafile(Psi4, "basicPsi4-1.2.1", ["water_mp2.out"])
        self.volume = volume.Volume((-4, -4, -4), (4, 4, 4), (0.2, 0.2, 0.2))

    def testmissingrequiredattributes(self):
        """Is an error raised when required attributes are missing?"""
        for missing_attribute in Hirshfeld.required_attrs:
            self.parse()
            delattr(self.data, missing_attribute)
            with pytest.raises(MissingAttributeError):
                trialBader = Hirshfeld(  # noqa: F841
                    self.data, self.volume, os.path.dirname(os.path.realpath(__file__))
                )

    def test_proatom_read(self):
        """Are proatom densities imported correctly?"""

        self.parse()
        self.analysis = Hirshfeld(
            self.data, self.volume, os.path.dirname(os.path.realpath(__file__))
        )

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

    def test_water_charges(self):
        """Are Hirshfeld charges calculated correctly for water?

        Note. Table 1 in doi:10.1007/BF01113058 reports Hirshfeld charge for Hydrogen atom as
              0.11 when STO-3G basis set was used and
              0.18 when 6-311G** basis set was used.
              Here, Psi4 calculation was done using STO-3G.
        """

        self.parse()
        # use precalculated fine cube file
        imported_vol = volume.read_from_cube(
            os.path.join(os.path.dirname(os.path.realpath(__file__)), "water_fine.cube")
        )

        analysis = Hirshfeld(self.data, imported_vol, os.path.dirname(os.path.realpath(__file__)))
        analysis.calculate()

        # Check assigned charges
        assert_allclose(analysis.fragcharges, [-0.29084274, 0.14357639, 0.14357639], atol=0.1)

    def test_chgsum_h2(self):
        """Are Hirshfeld charges for hydrogen atoms in nonpolar H2 small as expected?"""

        h2path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "h2.out")
        data = ccread(h2path)
        vol = volume.Volume((-3, -3, -3), (3, 3, 3), (0.1, 0.1, 0.1))
        analysis = Hirshfeld(data, vol, os.path.dirname(os.path.realpath(__file__)))
        analysis.calculate()

        assert abs(numpy.sum(analysis.fragcharges) - 0) < 1e-2
        assert abs(analysis.fragcharges[0] - analysis.fragcharges[1]) < 1e-6

    def test_chgsum_co(self):
        """Are Hirshfeld charges for carbon monoxide reported as expected?

        Note. Table 1 in doi:10.1007/BF01113058 reports Hirshfeld charge for Carbon atom as
              0.06 when STO-3G basis set was used and
              0.14 when 6-311G** basis set was used.
              Here, Psi4 calculation was done using STO-3G.
        """

        copath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "co.out")
        data = ccread(copath)
        vol = volume.read_from_cube(
            os.path.join(os.path.dirname(os.path.realpath(__file__)), "co.cube")
        )
        analysis = Hirshfeld(data, vol, os.path.dirname(os.path.realpath(__file__)))
        analysis.calculate()

        assert abs(numpy.sum(analysis.fragcharges)) < 1e-2
        assert_allclose(analysis.fragcharges, [0.10590126, -0.11277786], atol=1e-3)
