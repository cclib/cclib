# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import unittest
import numpy as np
from test.test_data import getdatafile
from cclib.bridge import cclib2pyscf
from cclib.parser.utils import find_package, convertor


class PyscfTest(unittest.TestCase):
    """Tests for the cclib2pyscf bridge in cclib."""

    def setUp(self):
        super(PyscfTest, self).setUp()
        if not find_package("pyscf"):
            raise ImportError("Must install pyscf to run this test")
        self.data, self.logfile = getdatafile(
            "GAMESS", "basicGAMESS-US2018", ["water_mp2.out"]
        )

    def test_makepyscf(self):
        import pyscf
        from pyscf import scf

        refen = convertor(self.data.scfenergies[-1],"eV","hartree")  # value in eVs
        pyscfmol = cclib2pyscf.makepyscf(self.data)
        pyscfmol.cart = True
        pyscfmol.verbose = 0
        pyscfmol.build()

        mhf = pyscfmol.HF(conv_tol=1e-9)
        en = mhf.kernel()
        assert abs(en - refen) < 1.0e-5
        # check that default basis is returned if basis is not present.
        del self.data.gbasis
        pyscfmol2 = cclib2pyscf.makepyscf(self.data)
        assert pyscfmol2.basis == "sto-3g"


if __name__ == "__main__":
    unittest.main()
