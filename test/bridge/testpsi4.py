# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import unittest

import numpy as np

from cclib.bridge import cclib2psi4


class Psi4Test(unittest.TestCase):
    """Tests for the cclib2psi4 bridge in cclib."""

    def test_makepsi4(self):
        import psi4
        from psi4 import energy

        psi4.core.set_output_file("psi4_output.dat", False)

        atomnos = np.array([1, 8, 1], "i")
        a = np.array([[-1, 1, 0], [0, 0, 0], [1, 1, 0]], "f")
        psi4mol = cclib2psi4.makepsi4(a, atomnos)
        en = energy("scf/6-31G**", molecule=psi4mol)
        ref = -75.82474605514503
        assert abs(en - ref) < 1.0e-6


if __name__ == "__main__":
    unittest.main()
