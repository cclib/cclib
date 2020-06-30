# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the various population analyses (MPA, LPA, CSPA, Bickelhaupt) in cclib"""

from __future__ import print_function

import sys
import os
import logging
import unittest

import numpy

from cclib.method import Bader, Volume
from cclib.parser import Psi4
from cclib.method.calculationmethod import MissingAttributeError

from numpy.testing import assert_allclose

from ..test_data import getdatafile


class BaderTest(unittest.TestCase):
    """Bader's QTAIM method tests."""

    def setUp(self):
        super(BaderTest, self).setUp()
        self.parse()

    def parse(self):
        self.data, self.logfile = getdatafile(Psi4, "basicPsi4-1.2.1", ["water_mp2.out"])
        self.volume = Volume((-4, -4, -4), (4, 4, 4), (0.2, 0.2, 0.2))

    def testmissingrequiredattributes(self):
        """Is an error raised when required attributes are missing?"""
        for missing_attribute in Bader.required_attrs:
            self.parse()
            delattr(self.data, missing_attribute)
            with self.assertRaises(MissingAttributeError):
                trialBader = Bader(self.data, self.volume)

    def test_val(self):
        """Do the calculated values match with known values?
        """

        self.data, logfile = getdatafile(Psi4, "basicPsi4-1.2.1", ["water_mp2.out"])
        self.volume = Volume((-4, -4, -4), (4, 4, 4), (0.2, 0.2, 0.2))
        self.analysis = Bader(self.data, self.volume)
        self.analysis.calculate()

        refData = [6.9378, 0.3639, 0.3827]  # values from `bader` package

        assert_allclose(self.analysis.fragcharges, refData, atol=0.15)
