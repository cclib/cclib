# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2016 the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test the Electrons method in cclib"""

from __future__ import print_function

import sys
import os
import re
import logging
import unittest

import numpy

sys.path.append("..")
from test_data import getdatafile

from cclib.method import Electrons
from cclib.parser import Gaussian
from cclib.parser import QChem


class ElectronsTest(unittest.TestCase):

    def test_count(self):
        """Test electron count for logfiles without pseudopotentials."""
        data, logfile = getdatafile(QChem, "basicQChem4.2", "water_mp4sdq.out")
        self.assertEqual(Electrons(data).count(), 10)
        self.assertEqual(Electrons(data).count(core=True), 10)
        data, logfile = getdatafile(Gaussian, "basicGaussian09", "water_cis.log")
        self.assertEqual(Electrons(data).count(), 10)
        self.assertEqual(Electrons(data).count(core=True), 10)

    def test_count_pseudopotential(self):
        """Test electron count for logfiles with pseudopotentials."""
        data, logfile = getdatafile(Gaussian, "basicGaussian09", "Mo4OCl4-sp.log")
        self.assertEqual(Electrons(data).count(), 120)
        self.assertEqual(Electrons(data).count(core=True), 188)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(unittest.makeSuite(ElectronsTest))
