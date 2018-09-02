# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the Electrons method in cclib"""

from __future__ import print_function

import sys
import os
import re
import logging
import unittest

import numpy

from cclib.method import Electrons
from cclib.parser import Gaussian
from cclib.parser import QChem

sys.path.insert(1, "..")

from ..test_data import getdatafile


class ElectronsTest(unittest.TestCase):

    def test_count(self):
        """Test electron count for logfiles without pseudopotentials."""
        data, logfile = getdatafile(QChem, "basicQChem4.2", ["water_mp4sdq.out"])
        self.assertEqual(Electrons(data).count(), 10)
        self.assertEqual(Electrons(data).count(core=True), 10)
        data, logfile = getdatafile(Gaussian, "basicGaussian09", ["water_cis.log"])
        self.assertEqual(Electrons(data).count(), 10)
        self.assertEqual(Electrons(data).count(core=True), 10)

    def test_count_pseudopotential(self):
        """Test electron count for logfiles with pseudopotentials."""
        data, logfile = getdatafile(Gaussian, "basicGaussian09", ["Mo4OCl4-sp.log"])
        self.assertEqual(Electrons(data).count(), 120)
        self.assertEqual(Electrons(data).count(core=True), 188)

    def test_alpha_beta(self):
        """Test number of alpha and beta electrons"""
        # Systems with alpha == beta
        data, logfile = getdatafile(QChem, "basicQChem4.2", ["water_mp4sdq.out"])
        self.assertEqual(Electrons(data).alpha(), 5)
        self.assertEqual(Electrons(data).beta(), 5)
        data, logfile = getdatafile(Gaussian, "basicGaussian09", ["water_cis.log"])
        self.assertEqual(Electrons(data).alpha(), 5)
        self.assertEqual(Electrons(data).beta(), 5)
        # Systems with alpha != beta
        data, logfile = getdatafile(QChem, "basicQChem4.2", ["dvb_sp_un.out"])
        self.assertEqual(Electrons(data).alpha(), 35)
        self.assertEqual(Electrons(data).beta(), 34)
        data, logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_un_sp.log"])
        self.assertEqual(Electrons(data).alpha(), 35)
        self.assertEqual(Electrons(data).beta(), 34)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(unittest.makeSuite(ElectronsTest))
