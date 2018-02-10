# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the MBO method in cclib"""

from __future__ import print_function

import sys
import os
import logging
import unittest

import numpy

from cclib.method import MBO
from cclib.parser import Gaussian

sys.path.insert(1, "..")

from ..test_data import getdatafile


class MBOTest(unittest.TestCase):

    def test_mbo_sp(self):
        """Testing Mayer bond orders for restricted single point."""

        data, logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])
        mbo = MBO(data)
        mbo.logger.setLevel(logging.ERROR)
        mbo.calculate()

        e_mbo = numpy.loadtxt(os.path.dirname(os.path.realpath(__file__)) + "/dvb_sp.mbo")
        self.assertTrue(numpy.all(mbo.fragresults[0] >= e_mbo - 0.25))
        self.assertTrue(numpy.all(mbo.fragresults[0] <= e_mbo + 0.25))

    def test_mbo_un_sp(self):
        """Testing Mayer bond orders for unrestricted single point."""

        data, logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_un_sp.log"])
        mbo = MBO(data)
        mbo.logger.setLevel(logging.ERROR)
        mbo.calculate()

        e_mbo = numpy.loadtxt(os.path.dirname(os.path.realpath(__file__)) + "/dvb_un_sp.mbo")
        bond_orders = mbo.fragresults[0] + mbo.fragresults[1]
        self.assertTrue(numpy.all(bond_orders >= e_mbo - 0.30))
        self.assertTrue(numpy.all(bond_orders <= e_mbo + 0.30))

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(unittest.makeSuite(MBOTest))
