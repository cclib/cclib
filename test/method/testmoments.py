# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the Moments method in cclib"""

from __future__ import print_function

import sys
import unittest

import numpy
from numpy.testing import assert_almost_equal

from cclib.method import Moments
from cclib.parser import Gaussian, NWChem

sys.path.insert(1, "..")

from ..test_data import getdatafile


class MomentsTest(unittest.TestCase):
    def test_results(self):
        data, _ = getdatafile(Gaussian, "basicGaussian16", ["water_mp2.log"])
        x = Moments(data).calculate()
        assert_almost_equal(x[0], [0, 0, 0], 5)
        assert_almost_equal(x[1], [0, 0, -0.91543], 5)
        assert_almost_equal(x[2][0] + x[2][3] + x[2][5], 0)

    def test_origin_displacement(self):
        data, _ = getdatafile(Gaussian, "basicGaussian16", ["water_mp2.log"])
        x = Moments(data).calculate()
        y = Moments(data).calculate(origin=(2e7,18,28))
        assert_almost_equal(x[1], y[1])

    def test_origin_at_center_of_nuclear_charge(self):
        data, _ = getdatafile(Gaussian, "basicGaussian16", ["water_mp2.log"])
        x = Moments(data).calculate(origin="nuccharge")
        assert_almost_equal(x[0], [0, 0, 0], 6)

    def test_origin_at_center_of_mass(self):
        data, _ = getdatafile(Gaussian, "basicGaussian16", ["water_mp2.log"])
        x = Moments(data).calculate(origin="mass")
        assert_almost_equal(x[0], [0, 0, 0.0524806])

    def test_user_provided_origin(self):
        data, _ = getdatafile(Gaussian, "basicGaussian16", ["water_mp2.log"])
        x = Moments(data).calculate(origin=(1,1,1))
        assert_almost_equal(x[0], [1, 1, 1])
                            
    def test_user_provided_masses(self):
        data, _ = getdatafile(Gaussian, "basicGaussian16", ["water_mp2.log"])
        x = Moments(data).calculate(masses=[1,1,1], origin="mass")
        assert_almost_equal(x[0], [0, 0, -0.2780383])

    def test_inplace_modifying(self):
        data, _ = getdatafile(NWChem, "basicNWChem6.5", ["water_mp2.out"])
        assert not hasattr(data, "moments")        
        Moments(data).calculate(overwrite=True)
        assert hasattr(data, "moments")

    def test_inplace_modifying2(self):
        data, _ = getdatafile(Gaussian, "basicGaussian16", ["water_mp2.log"])
        Moments(data).calculate(overwrite=True)
        assert_almost_equal(data.moments[1], [0, 0, -0.91543], 5)

        
if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(unittest.makeSuite(MomentsTest))
