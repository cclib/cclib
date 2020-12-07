#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the various population analyses (MPA, LPA, CSPA, Bickelhaupt) in cclib"""

import sys
import os
import logging
import unittest

import numpy

from cclib.method import CSPA, LPA, MPA, OPA, Bickelhaupt
from cclib.method.calculationmethod import MissingAttributeError
from cclib.parser import Gaussian

sys.path.insert(1, "..")

from ..test_data import getdatafile


class PopulationTest(unittest.TestCase):
    """Generic population method tests."""
    
    methods = (CSPA, LPA, MPA, OPA, Bickelhaupt)

    def parse(self):
        self.data, self.logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_un_sp.log"])

    def calculate(self, method_class):
        if not hasattr(self, 'data'):
            self.parse()
        self.analysis = method_class(self.data)
        self.analysis.logger.setLevel(0)
        self.analysis.calculate()

    def testmissingrequiredattributes(self):
        """Is an error raised when required attributes are missing?"""
        for missing_attribute in MPA.required_attrs:
            self.parse()
            delattr(self.data, missing_attribute)
            for method_class in self.methods:
                with self.assertRaises(MissingAttributeError):
                    self.calculate(method_class)

    def testmissingoverlaps(self):
        """Is an error raised when no overlaps are available?"""
        self.parse()
        for overlap_attribute in MPA.overlap_attributes:
            if hasattr(self.data, overlap_attribute):
                delattr(self.data, overlap_attribute)
        for method_class in self.methods:
            if method_class.overlap_attributes:
                with self.assertRaises(MissingAttributeError):
                    self.calculate(method_class)


class GaussianMPATest(unittest.TestCase):
    """Mulliken Population Analysis test"""

    def setUp(self):
        self.data, self.logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_un_sp.log"])
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
        self.data, self.logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_un_sp.log"])
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
        self.data, self.logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_un_sp.log"])
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

class GaussianBickelhauptTest(unittest.TestCase):
    """Bickelhaupt Population Analysis test"""
    
    def setUp(self):
        super().setUp()
        self.data, self.logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_un_sp.log"])
        self.analysis = Bickelhaupt(self.data)
        self.analysis.logger.setLevel(0)
        self.analysis.calculate()
    
    def testsumcharges(self):
        """Do the Bickelhaupt charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos) - self.data.charge
        totalpopulation = sum(self.analysis.fragcharges)
        self.assertAlmostEqual(totalpopulation, formalcharge, delta=1.0e-3)
        
    def testsumspins(self):
        """Do the Bickelhaupt spins sum up to the total formal spin?"""
        formalspin = self.data.homos[0] - self.data.homos[1]
        totalspin = sum(self.analysis.fragspins)
        self.assertAlmostEqual(totalspin, formalspin, delta=1.0e-3)
    
    def test_dvb_sp(self):
        """Testing Bickelhaupt charges (restricted) against outputs from Multiwfn."""
        data, logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])
        bpa = Bickelhaupt(data)
        bpa.logger.setLevel(logging.ERROR)
        bpa.calculate()
        
        e_bpa = numpy.loadtxt(os.path.dirname(os.path.realpath(__file__)) + "/dvb_sp.bpa")
        self.assertTrue(numpy.all(bpa.fragcharges >= e_bpa - 0.05))
        self.assertTrue(numpy.all(bpa.fragcharges <= e_bpa + 0.05))
        
    def test_dvb_un_sp(self):
        """Testing Bickelhaupt charges (unrestricted) against outputs from Multiwfn."""
        data, logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_un_sp.log"])
        bpa = Bickelhaupt(data)
        bpa.logger.setLevel(logging.ERROR)
        bpa.calculate()
        
        e_bpaalpha = numpy.loadtxt(os.path.dirname(os.path.realpath(__file__)) + "/dvb_un_sp.bpa")
        e_bpaspin = numpy.loadtxt(os.path.dirname(os.path.realpath(__file__)) + "/dvb_un_sp.bpaspin")
        
        self.assertTrue(numpy.all(bpa.fragcharges >= e_bpaalpha - 0.05))
        self.assertTrue(numpy.all(bpa.fragcharges <= e_bpaalpha + 0.05))
        self.assertTrue(numpy.all(bpa.fragspins >= e_bpaspin - 0.05))
        self.assertTrue(numpy.all(bpa.fragspins <= e_bpaspin + 0.05))

tests = [GaussianMPATest, GaussianLPATest, GaussianCSPATest, GaussianBickelhauptTest]


if __name__ == "__main__":
    for test in tests:
        thistest = unittest.makeSuite(test)
        unittest.TextTestRunner(verbosity=2).run(thistest)
