# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the various population analyses (MPA, LPA, CSPA) in cclib"""

from __future__ import print_function

import sys
import os
import logging
import unittest

import numpy

from cclib.method import MPA, LPA, CSPA
from cclib.parser import Gaussian

sys.path.insert(1, "..")

from ..test_data import getdatafile


class GaussianMPATest(unittest.TestCase):
    """Mulliken Population Analysis test"""

    def setUp(self):
        self.data, self.logfile = getdatafile(Gaussian, "basicGaussian03", ["dvb_un_sp.out"])
        self.analysis = MPA(self.data)
        self.analysis.logger.setLevel(0)
        self.analysis.calculate()

    def testsumcharges(self):
        """Do the Mulliken charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos) - self.data.charge
        totalpopulation = sum(self.analysis.fragcharges)
        self.assertAlmostEqual(totalpopulation, formalcharge, delta=1.0e-3)

    def testsumspins(self):
        """Do the Mulliken spins sum up to the total formal spin?"""
        formalspin = self.data.homos[0] - self.data.homos[1]
        totalspin = sum(self.analysis.fragspins)
        self.assertAlmostEqual(totalspin, formalspin, delta=1.0e-3)


class GaussianLPATest(unittest.TestCase):
    """Lowdin Population Analysis test"""

    def setUp(self):
        self.data, self.logfile = getdatafile(Gaussian, "basicGaussian03", ["dvb_un_sp.out"])
        self.analysis = LPA(self.data)
        self.analysis.logger.setLevel(0)
        self.analysis.calculate()

    def testsumcharges(self):
        """Do the Lowdin charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos) - self.data.charge
        totalpopulation = sum(self.analysis.fragcharges)
        self.assertAlmostEqual(totalpopulation, formalcharge, delta=0.001)

    def testsumspins(self):
        """Do the Lowdin spins sum up to the total formal spin?"""
        formalspin = self.data.homos[0] - self.data.homos[1]
        totalspin = sum(self.analysis.fragspins)
        self.assertAlmostEqual(totalspin, formalspin, delta=1.0e-3)


class GaussianCSPATest(unittest.TestCase):
    """C-squared Population Analysis test"""

    def setUp(self):
        self.data, self.logfile = getdatafile(Gaussian, "basicGaussian03", ["dvb_un_sp.out"])
        self.analysis = CSPA(self.data)
        self.analysis.logger.setLevel(0)
        self.analysis.calculate()

    def testsumcharges(self):
        """Do the CSPA charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos) - self.data.charge
        totalpopulation = sum(self.analysis.fragcharges)
        self.assertAlmostEqual(totalpopulation, formalcharge, delta=1.0e-3)

    def testsumspins(self):
        """Do the CSPA spins sum up to the total formal spin?"""
        formalspin = self.data.homos[0] - self.data.homos[1]
        totalspin = sum(self.analysis.fragspins)
        self.assertAlmostEqual(totalspin, formalspin, delta=1.0e-3)


tests = [GaussianMPATest, GaussianLPATest, GaussianCSPATest]


if __name__ == "__main__":
    for test in tests:
        thistest = unittest.makeSuite(test)
        unittest.TextTestRunner(verbosity=2).run(thistest)
