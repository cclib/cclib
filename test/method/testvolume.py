# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the Volume and related methods in cclib"""

import sys
from pathlib import Path

from cclib.method import volume
from cclib.parser import Gaussian, Psi4

import numpy

sys.path.insert(1, "..")

from numpy.testing import assert_allclose

from ..test_data import getdatafile


class VolumeTest:
    def test_scinotation(self):
        """Does the scientific notation writer work as expected?"""

        assert volume.scinotation(1.0 / 654) == " 1.52905E-03"
        assert volume.scinotation(-1.0 / 654) == "-1.52905E-03"

    def test_wavefunction(self):
        """Does the volume occupied by the HOMO integrate to the correct
        values?
        """

        data_basis, _ = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])
        data_sp, _ = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])

        vol = volume.Volume((-3.0, -6.0, -2.0), (3.0, 6.0, 2.0), (0.25, 0.25, 0.25))

        wavefn = volume.wavefunction(data_sp, vol, data_sp.mocoeffs[0][data_sp.homos[0]])
        integral = wavefn.integrate()
        integral_square = wavefn.integrate_square()

        assert abs(integral) < 1e-6  # not necessarily true for all wavefns
        assert abs(integral_square - 1.00) < 1e-2  # true for all wavefns
        print(integral, integral_square)

    def test_density(self):
        """Does the volume occupied by the combined electron density integrate
        to the correct value?
        """

        data_basis, _ = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])
        data_sp, _ = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])

        vol = volume.Volume((-3.0, -6.0, -2.0), (3.0, 6.0, 2.0), (0.25, 0.25, 0.25))

        frontierorbs = [data_sp.mocoeffs[0][(data_sp.homos[0] - 3) : (data_sp.homos[0] + 1)]]
        density = volume.electrondensity(data_sp, vol, frontierorbs)
        integral = density.integrate()

        assert abs(integral - 8.00) < 1e-2
        print("Combined Density of 4 Frontier orbitals=", integral)

    def test_cube(self):
        """Does the cube file written out for electron density on a Cartesian grid match
        expected values?
        """

        data, logfile = getdatafile(Psi4, "basicPsi4-1.2.1", ["water_mp2.out"])

        # Reference values were calculated using cubegen method in Psi4.
        # First six rows are information about the coordinates of the grid and comments.
        tmp = []

        with open(Path(__file__).resolve().parent / "water_mp2.cube") as f:
            lines = f.readlines()
            for line in lines[6 : len(lines)]:
                tmp.extend(line.split())
        tmp = numpy.asanyarray(tmp, dtype=float)
        refcube = numpy.resize(tmp, (13, 13, 13))

        # Values for the grid below are constructed to match Psi4 cube file.
        vol = volume.Volume(
            (-1.587532, -1.587532, -1.356299),
            (1.58754, 1.58754, 1.81877),
            (0.26458860545, 0.26458860545, 0.26458860545),
        )
        density = volume.electrondensity(data, vol, [data.mocoeffs[0][: data.homos[0]]])

        assert_allclose(density.data, refcube, atol=0.5, rtol=0.1)

    def test_roundtrip_cube(self):
        """Write a cube file and then read it back. Check if the volume object contains
        identical information on each grid point"""

        data, logfile = getdatafile(Psi4, "basicPsi4-1.2.1", ["water_mp2.out"])
        vol = volume.Volume((-1, -1, -1), (1, 1, 1), (0.4, 0.4, 0.4))
        density = volume.electrondensity(data, vol, [data.mocoeffs[0][: data.homos[0]]])

        density.writeascube("coarsewater.cube")
        density_recovered = volume.read_from_cube("coarsewater.cube")

        assert_allclose(density.data, density_recovered.data, rtol=0.05)

    def test_zip_cube(self):
        """Check we can read from a zipped file."""
        data = volume.read_from_cube(Path(__file__).resolve().parent / "co.cube.zip")
        assert len(data.data) > 0
