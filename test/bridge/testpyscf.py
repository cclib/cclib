# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

from test.test_data import getdatafile

from cclib.bridge import cclib2pyscf
from cclib.parser.utils import find_package

import numpy as np


class PyscfTest:
    @classmethod
    def setup_class(cls) -> None:
        if not find_package("pyscf"):
            raise ImportError("Must install pyscf to run this test")
        cls.data, cls.logfile = getdatafile("Gaussian", "basicGaussian16", ["dvb_sp.out"])
        cls.udata, cls.ulogfile = getdatafile("Gaussian", "basicGaussian16", ["dvb_un_sp.log"])

    def test_makepyscf(self) -> None:
        from pyscf import dft

        refen = self.data.scfenergies[-1]
        pyscfmol = cclib2pyscf.makepyscf(self.data)
        pyscfmol.cart = True
        pyscfmol.verbose = 0
        pyscfmol.build()

        mhf = dft.RKS(pyscfmol)
        mhf.xc = "b3lyp"
        en = mhf.kernel()
        assert abs(en - refen) < 5.0e-5
        # check that default basis is returned if basis is not present.
        del self.data.gbasis
        pyscfmol2 = cclib2pyscf.makepyscf(self.data)
        assert pyscfmol2.basis == "sto-3g"

    def test_makepyscf_mos(self) -> None:
        pyscfmol = cclib2pyscf.makepyscf(self.data)
        mo_coeff, mo_occ, mo_syms, mo_energies = cclib2pyscf.makepyscf_mos(self.data, pyscfmol)
        assert np.allclose(mo_energies, self.data.moenergies)
        # check first MO coefficient
        assert np.allclose(mo_coeff[0][0], self.data.mocoeffs[0][0][0])
        # check a random middle MO coefficient
        assert np.allclose(mo_coeff[0][10], self.data.mocoeffs[0][10][0])
        # test unrestricted code.
        pyscfmol = cclib2pyscf.makepyscf(self.udata)
        mo_coeff, mo_occ, mo_syms, mo_energies = cclib2pyscf.makepyscf_mos(self.udata, pyscfmol)
        assert np.allclose(mo_energies, self.udata.moenergies)
        # check first MO coefficient
        assert np.allclose(mo_coeff[0][0][0], self.udata.mocoeffs[0][0][0])
        # check a random middle MO coefficient
        assert np.allclose(mo_coeff[0][0][10], self.udata.mocoeffs[0][10][0])
