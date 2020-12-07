#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import unittest

import numpy as np

from cclib.bridge import cclib2pyscf
from cclib.parser.utils import find_package


class PyscfTest(unittest.TestCase):
    """Tests for the cclib2pyscf bridge in cclib."""

    def setUp(self):
        super(PyscfTest, self).setUp()        
        if not find_package('pyscf'):
            raise ImportError('Must install pyscf to run this test')

    def test_makepyscf(self):
        import pyscf
        from pyscf import scf

        atomnos = np.array([1, 8, 1], "i")
        atomcoords = np.array([[-1, 1, 0], [0, 0, 0], [1, 1, 0]], "f")
        pyscfmol = cclib2pyscf.makepyscf(atomcoords, atomnos)
        pyscfmol.basis = "6-31G**"
        pyscfmol.cart = True
        pyscfmol.verbose = 0
        pyscfmol.build()

        mhf = pyscfmol.HF(conv_tol=1e-6)
        en = mhf.kernel()
        ref = -75.824754602
        assert abs(en - ref) < 1.0e-6


if __name__ == "__main__":
    unittest.main()
