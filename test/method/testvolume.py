# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the Volume and related methods in cclib"""

from __future__ import print_function

import os, sys
import unittest
import numpy

from cclib.method import volume
from cclib.parser import Gaussian, Psi4

sys.path.insert(1, "..")

from ..test_data import getdatafile
from numpy.testing import assert_allclose


class VolumeTest(unittest.TestCase):
    def test_scinotation(self):
        """Does the scientific notation writer work as expected?"""

        self.assertEqual(volume.scinotation(1.0 / 654), " 1.52905E-03")
        self.assertEqual(volume.scinotation(-1.0 / 654), "-1.52905E-03")

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

        self.assertAlmostEqual(integral, 0, delta=1e-6)  # not necessarily true for all wavefns
        self.assertAlmostEqual(integral_square, 1.00, delta=1e-2)  # true for all wavefns
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

        self.assertTrue(abs(integral - 8.00) < 1e-2)
        print("Combined Density of 4 Frontier orbitals=", integral)

    def test_cube(self):
        """Does the cube file written out for electron density on a Cartesian grid match
        expected values?
        """

        data, logfile = getdatafile(Psi4, "basicPsi4-1.2.1", ["water_mp2.out"])

        # Reference values were calculated using cubegen method in Psi4.
        # First six rows are information about the coordinates of the grid and comments.
        tmp = []
        
        with open(os.path.dirname(os.path.realpath(__file__)) + "/water_mp2.cube") as f:
            lines = f.readlines()
            for line in lines[6 : len(lines)]:
                tmp.extend(line.split())
        tmp = numpy.asanyarray(tmp, dtype = float)
        refcube = numpy.resize(tmp, (13, 13, 13))

        # Values for the grid below are constructed to match Psi4 cube file.
        vol = volume.Volume(
            (-1.587532, -1.587532, -1.356299), 
            (1.58754, 1.58754, 1.81877), 
            (0.26458860545, 0.26458860545, 0.26458860545)
        )
        density = volume.electrondensity(data, vol, [data.mocoeffs[0][: data.homos[0]]])

        assert_allclose(density.data, refcube, atol=.5, rtol=.1)
