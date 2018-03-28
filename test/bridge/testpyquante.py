# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import unittest

import numpy

from cclib.bridge import cclib2pyquante


class PyquanteTest(unittest.TestCase):
    """Tests for the cclib2pyquante bridge in cclib."""

    def test_makepyquante(self):
        from PyQuante.hartree_fock import hf
        atomnos = numpy.array([1, 8, 1],"i")
        a = numpy.array([[-1, 1, 0], [0, 0, 0], [1, 1, 0]], "f")
        pyqmol = cclib2pyquante.makepyquante(a, atomnos)
        en, orbe, orbs = hf(pyqmol)
        ref = -75.824754
        assert abs(en - ref) < 1.0e-6


if __name__ == "__main__":
    unittest.main()
