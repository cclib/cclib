# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2014-2016 the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test the Nuclear method in cclib"""

from __future__ import print_function

import sys
import os
import re
import logging
import unittest

import numpy

sys.path.append("..")
from test_data import getdatafile
from cclib.method import Nuclear
from cclib.parser import ccData
from cclib.parser import QChem
from cclib.parser import utils


class NuclearTest(unittest.TestCase):

    def test_stoichiometry(self):
        """Testing stoichoimetry generation."""
        data = ccData()
        data.atomnos = numpy.areray([6, 1, 6, 1, 1, 1])
        self.assertEqual(Nuclear(data).stoichoimetry(), "C2H4")

    def test_nre(self):
        """Testing nuclear repulsion energy for one logfile where it is printed."""

        data, logfile = getdatafile(QChem, "basicQChem4.2", "water_mp4sdq.out")
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        with open(logfile.filename) as f:
            output = f.read()
        line = re.search('Nuclear Repulsion Energy = .* hartrees', output).group()
        nre = float(line.split()[4])
        nre = utils.convertor(nre, 'Angstrom', 'bohr')
        self.assertAlmostEqual(nuclear.repulsion_energy(), nre, places=7)


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(unittest.makeSuite(NuclearTest))
