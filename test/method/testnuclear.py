# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the Nuclear method in cclib"""

import sys
import os
import re
import logging
import unittest

import numpy

from cclib.method import Nuclear
from cclib.parser import ccData
from cclib.parser import QChem
from cclib.parser import utils

sys.path.insert(1, "..")

from ..test_data import getdatafile


class NuclearTest(unittest.TestCase):

    def test_stoichiometry(self):
        """Testing stoichoimetry generation."""
        data = ccData()

        def check(atomnos, formula, charge=0):
            data.natom = len(atomnos)
            data.atomnos = numpy.array(atomnos)
            data.atomcoords = numpy.zeros((data.natom, 3))
            data.charge = charge
            self.assertEqual(Nuclear(data).stoichiometry(), formula)

        # Basics and permutations.
        check([], "")
        check([6, 1, 6, 1, 1, 1], "C2H4")
        check([1, 1, 1, 6, 1, 6], "C2H4")
       
        # Charges.
        check([8], "O", charge=0)
        check([8], "O(+1)", charge=1)
        check([8], "O(-1)", charge=-1)
        check([8], "O(+2)", charge=2)
        check([8], "O(+9)", charge=9)

        # Element counts.
        check([6, 1], "CH")
        check([6] * 60, "C60")

        # Test the Hill system.
        check([8, 1, 1], "H2O")
        check([6, 8, 8, 1, 1], "CH2O2")
        check([16, 16, 8, 8], "O2S2")

    def test_nre(self):
        """Testing nuclear repulsion energy for one logfile where it is printed."""

        data, logfile = getdatafile(QChem, "basicQChem4.2", ["water_mp4sdq.out"])
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
