# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the Volume and related methods in cclib"""

from __future__ import print_function

import sys
import unittest

from cclib.method import volume
from cclib.parser import Gaussian

sys.path.insert(1, "..")

from ..test_data import getdatafile


class VolumeTest(unittest.TestCase):

    def test_scinotation(self):
        """Does the scientific notation writer work as expected?"""

        self.assertEqual(volume.scinotation(1./654), ' 1.52905E-03')
        self.assertEqual(volume.scinotation(-1./654), '-1.52905E-03')

    def test_wavefunction(self):
        """Does the volume occupied by the HOMO integrate to the correct
        values?
        """

        data_basis, _ = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])
        data_sp, _ = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])

        vol = volume.Volume((-3.0, -6.0, -2.0), (3.0, 6.0, 2.0), (0.25, 0.25, 0.25))

        wavefn = volume.wavefunction(data_sp.atomcoords[0],
                                     data_sp.mocoeffs[0][data_sp.homos[0]],
                                     data_basis.gbasis, vol)
        integral = wavefn.integrate()
        integral_square = wavefn.integrate_square()

        self.assertTrue(abs(integral) < 1e-6) # not necessarily true for all wavefns
        self.assertTrue(abs(integral_square - 1.00) < 1e-3) #   true for all wavefns
        print(integral, integral_square)

    def test_density(self):
        """Does the volume occupied by the combined electron density integrate
        to the correct value?
        """

        data_basis, _ = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])
        data_sp, _ = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])

        vol = volume.Volume((-3.0, -6.0, -2.0), (3.0, 6.0, 2.0), (0.25, 0.25, 0.25))

        frontierorbs = [data_sp.mocoeffs[0][(data_sp.homos[0] - 3):(data_sp.homos[0] + 1)]]
        density = volume.electrondensity(data_sp.atomcoords[0],
                                         frontierorbs, data_basis.gbasis, vol)
        integral = density.integrate()

        self.assertTrue(abs(integral - 8.00) < 1e-2)
        print("Combined Density of 4 Frontier orbitals=", integral)
