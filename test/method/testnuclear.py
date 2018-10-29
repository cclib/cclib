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

import numpy as np

from cclib.method import Nuclear
from cclib.parser import ccData
from cclib.parser import DALTON
from cclib.parser import Gaussian
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
            data.atomnos = np.array(atomnos)
            data.atomcoords = np.zeros((data.natom, 3))
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

    def test_repulsion_energy(self):
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

    @unittest.skipIf(sys.version_info < (3, 0), "The periodictable package doesn't work in Python2.")
    def test_principal_moments_of_inertia(self):
        """Testing principal moments of inertia and the principal axes for one
        logfile where it is printed.
        """

        data, logfile = getdatafile(DALTON, "basicDALTON-2015", ["dvb_sp_hf.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        ref_pmoi = []
        ref_axes = []
        with open(logfile.filename) as f:
            for line in f:
                if line.strip() == "Principal moments of inertia (u*A**2) and principal axes":
                    next(f)
                    next(f)
                    for _ in range(3):
                        tokens = [float(x) for x in next(f).split()[1:]]
                        ref_pmoi.append(tokens[0])
                        ref_axes.append(tokens[1:])
        pmoi, axes = nuclear.principal_moments_of_inertia()
        np.testing.assert_allclose(pmoi, ref_pmoi, rtol=0, atol=1.0e-4)
        # The phases of the eigenvectors may be different, but they
        # are still orthonormal within each set.
        np.testing.assert_allclose(np.abs(axes), np.abs(ref_axes), rtol=0, atol=1.0e-4)

    @unittest.skipIf(sys.version_info < (3, 0), "The periodictable package doesn't work in Python2.")
    def test_rotational_constants(self):
        """Testing rotational constants for two logfiles where they are
        printed.
        """

        data, logfile = getdatafile(DALTON, "basicDALTON-2015", ["dvb_sp_hf.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        with open(logfile.filename) as f:
            for line in f:
                if line.strip() == "Rotational constants":
                    while line.split() != ['A', 'B', 'C']:
                        line = next(f)
                    line = next(f)
                    ref_mhz = [float(x) for x in next(f).split()[:-1]]
                    ref_invcm = [float(x) for x in next(f).split()[:-1]]
        rotconsts_ghz = nuclear.rotational_constants('ghz')
        rotconsts_invcm = nuclear.rotational_constants('invcm')
        np.testing.assert_allclose(rotconsts_ghz * 1.0e3, ref_mhz, rtol=0, atol=1.0e-4)
        np.testing.assert_allclose(rotconsts_invcm, ref_invcm, rtol=0, atol=1.0e-4)

        data, logfile = getdatafile(Gaussian, "basicGaussian16", ["dvb_sp.out"])
        nuclear = Nuclear(data)
        nuclear.logger.setLevel(logging.ERROR)

        with open(logfile.filename) as f:
            for line in f:
                if "Rotational constants (GHZ):" in line:
                    ref_ghz = [float(x) for x in line.split()[3:]]
        rotconsts_ghz = nuclear.rotational_constants('ghz')
        np.testing.assert_allclose(rotconsts_ghz, ref_ghz, rtol=0, atol=1.0e-5)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(unittest.makeSuite(NuclearTest))
