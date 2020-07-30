# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the DDEC6 in cclib"""

from __future__ import print_function

import sys
import os
import logging
import unittest

import numpy

from cclib.method import DDEC6, volume
from cclib.parser import Psi4
from cclib.method.calculationmethod import MissingAttributeError

from numpy.testing import assert_allclose

from ..test_data import getdatafile


class DDEC6Test(unittest.TestCase):
    """DDEC6 method tests."""

    def setUp(self):
        super(DDEC6Test, self).setUp()
        self.parse()

    def parse(self):
        self.data, self.logfile = getdatafile(Psi4, "basicPsi4-1.2.1", ["water_mp2.out"])
        self.volume = volume.Volume((-4, -4, -4), (4, 4, 4), (0.2, 0.2, 0.2))

    def testmissingrequiredattributes(self):
        """Is an error raised when required attributes are missing?"""
        for missing_attribute in DDEC6.required_attrs:
            self.parse()
            delattr(self.data, missing_attribute)
            with self.assertRaises(MissingAttributeError):
                trialBader = DDEC6(self.data, self.volume)

    def test_proatom_read(self):
        """Are proatom densities imported correctly?"""

        self.parse()
        self.analysis = DDEC6(self.data, self.volume, os.path.dirname(os.path.realpath(__file__)))

        refH_den = [
            2.66407645e-01,
            2.66407645e-01,
            2.66407643e-01,
            2.66407612e-01,
            2.66407322e-01,
        ]  # Hydrogen first five densities
        refH_r = [
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
        refO_r = [
            5.70916728e-09,
            1.97130512e-07,
            1.56506399e-06,
            6.80667366e-06,
            2.12872046e-05,
        ]  # Oxygen first five radii

        assert_allclose(self.analysis.proatom_density[0][0:5], refO_den, rtol=1e-3)
        assert_allclose(self.analysis.proatom_density[1][0:5], refH_den, rtol=1e-3)
        assert_allclose(self.analysis.proatom_density[2][0:5], refH_den, rtol=1e-3)

    def test_step1_charges(self):
        """Are step 1 charges calculated correctly?"""

        self.parse()
        # use precalculated fine cube file
        imported_vol = volume.read_from_cube(os.path.join(os.path.dirname(os.path.realpath(__file__)), "water_fine.cube"))
        
        analysis = DDEC6(self.data, imported_vol, os.path.dirname(os.path.realpath(__file__)))
        analysis.calculate()
        
        # values from `chargemol` calculation
        assert_allclose(analysis.refcharges, [-0.513006, 0.256231, 0.256775], rtol=.1)
