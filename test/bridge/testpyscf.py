# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import unittest
from test.test_data import getdatafile

from cclib.bridge import cclib2pyscf
from cclib.parser.utils import convertor, find_package

import numpy as np


class PyscfTest(unittest.TestCase):
    def setUp(self) -> None:
        super(PyscfTest, self).setUp()
        if not find_package("pyscf"):
            raise ImportError("Must install pyscf to run this test")
        self.data, self.logfile = getdatafile("GAMESS", "basicGAMESS-US2018", ["dvb_sp.out"])
        self.udata, self.ulogfile = getdatafile("GAMESS", "basicGAMESS-US2018", ["dvb_un_sp.out"])

    def test_makepyscf(self) -> None:
        import pyscf
        from pyscf import dft

        refen = convertor(self.data.scfenergies[-1], "eV", "hartree")  # value in eVs
        pyscfmol = cclib2pyscf.makepyscf(self.data)
        pyscfmol.cart = True
        pyscfmol.verbose = 0
        pyscfmol.build()

        mhf = dft.RKS(pyscfmol)
        mhf.xc = "b3lyp"
        en = mhf.kernel()
        assert abs(en - refen) < 1.0e-4
        # check that default basis is returned if basis is not present.
        del self.data.gbasis
        pyscfmol2 = cclib2pyscf.makepyscf(self.data)
        assert pyscfmol2.basis == "sto-3g"

    def test_makepyscf_mos(self) -> None:
        pyscfmol = cclib2pyscf.makepyscf(self.data)
        mo_coeff, mo_occ, mo_syms, mo_energies = cclib2pyscf.makepyscf_mos(self.data, pyscfmol)
        assert np.allclose(mo_energies, convertor(np.array(self.data.moenergies), "eV", "hartree"))
        # check first MO coefficient
        assert np.allclose(mo_coeff[0][0], self.data.mocoeffs[0][0][0])
        # check a random middle MO coefficient
        assert np.allclose(mo_coeff[0][10], self.data.mocoeffs[0][10][0])
        # test unrestricted code.
        pyscfmol = cclib2pyscf.makepyscf(self.udata)
        mo_coeff, mo_occ, mo_syms, mo_energies = cclib2pyscf.makepyscf_mos(self.udata, pyscfmol)
        assert np.allclose(mo_energies, convertor(np.array(self.udata.moenergies), "eV", "hartree"))
        # check first MO coefficient
        assert np.allclose(mo_coeff[0][0][0], self.udata.mocoeffs[0][0][0])
        # check a random middle MO coefficient
        assert np.allclose(mo_coeff[0][0][10], self.udata.mocoeffs[0][10][0])


if __name__ == "__main__":
    unittest.main()
    PyscfTest.test_makepyscf_from_mos()
