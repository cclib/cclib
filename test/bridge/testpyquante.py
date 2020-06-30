# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import unittest

import numpy

from cclib.bridge import cclib2pyquante
from ..test_data import getdatafile
from cclib.parser.utils import find_package

from numpy.testing import assert_array_almost_equal


class PyquanteTest(unittest.TestCase):
    """Tests for the cclib2pyquante bridge in cclib."""

    def setUp(self):
        super(PyquanteTest, self).setUp()
        self._found_pyquante = find_package("PyQuante")
        self.data, self.logfile = getdatafile("Gaussian", "basicGaussian16", ["water_ccsd.log"])

    def test_makepyquante(self):
        # Test older PyQuante bridge
        from PyQuante.hartree_fock import hf
        from PyQuante.Molecule import Molecule

        reference = Molecule(
            "h2o",
            [(8, (0, 0, 0.119159)), (1, (0, 0.790649, -0.476637)), (1, (0, -0.790649, -0.476637)),],
            units="Angstroms",
        )
        refen, reforbe, reforbs = hf(reference)

        pyqmol = cclib2pyquante.makepyquante(self.data)
        en, orbe, orbs = hf(pyqmol)

        self.assertAlmostEqual(en, refen, delta=1.0e-6)


class pyquante2Test(unittest.TestCase):
    """Tests for the cclib2pyquante bridge in cclib."""

    def setUp(self):
        super(pyquante2Test, self).setUp()
        self._found_pyquante2 = find_package("pyquante2")
        self.data, self.logfile = getdatafile("Gaussian", "basicGaussian16", ["water_ccsd.log"])

    def test_makepyquante(self):
        # Test pyquante2 bridge
        from pyquante2 import molecule, rhf, h2o, basisset

        bfs = basisset(h2o)
        # Copied from water_ccsd.log
        refmol = molecule(
            [(8, 0.0, 0.0, 0.119159), (1, 0, 0.790649, -0.476637), (1, 0, -0.790649, -0.476637)],
            units="Angstroms",
        )
        refsolver = rhf(refmol, bfs)
        refsolver.converge()

        pyqmol = cclib2pyquante.makepyquante(self.data)
        pyqsolver = rhf(pyqmol, bfs)
        pyqsolver.converge()

        assert_array_almost_equal(refsolver.energies[-1], pyqsolver.energies[-1], decimal=6)


if __name__ == "__main__":
    unittest.main()
