# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import unittest

from cclib.bridge import cclib2pyquante
from cclib.parser.utils import find_package

import numpy
from numpy.testing import assert_array_almost_equal

from ..test_data import getdatafile


class pyquante2Test(unittest.TestCase):
    """Tests for the cclib2pyquante bridge in cclib."""

    def setUp(self):
        super(pyquante2Test, self).setUp()
        if not find_package("pyquante2"):
            raise ImportError("Must install pyquante2 to run this test")

        self.data, self.logfile = getdatafile("Gaussian", "basicGaussian16", ["water_ccsd.log"])

    def test_makepyquante(self):
        # Test pyquante2 bridge
        from pyquante2 import basisset, h2o, molecule, rhf

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
