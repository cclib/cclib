# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the Electrons method in cclib"""

import sys
import unittest

from cclib.method import Electrons
from cclib.parser import Gaussian
from cclib.parser import QChem

sys.path.insert(1, "..")

from ..test_data import getdatafile


class ElectronsTest(unittest.TestCase):
    def test_count(self):
        """Test electron count for logfiles without pseudopotentials."""
        data, logfile = getdatafile(QChem, "basicQChem5.4", ["water_mp4sdq.out"])
        assert Electrons(data).count() == 10
        assert Electrons(data).count(core=True) == 10
        data, logfile = getdatafile(Gaussian, "basicGaussian09", ["water_cis.log"])
        assert Electrons(data).count() == 10
        assert Electrons(data).count(core=True) == 10

    def test_count_pseudopotential(self):
        """Test electron count for logfiles with pseudopotentials."""
        data, logfile = getdatafile(Gaussian, "basicGaussian09", ["Mo4OCl4-sp.log"])
        assert Electrons(data).count() == 120
        assert Electrons(data).count(core=True) == 188

    def test_alpha_beta(self):
        """Test number of alpha and beta electrons"""
        # Systems with alpha == beta
        data, logfile = getdatafile(QChem, "basicQChem5.4", ["water_mp4sdq.out"])
        assert Electrons(data).alpha() == 5
        assert Electrons(data).beta() == 5
        data, logfile = getdatafile(Gaussian, "basicGaussian09", ["water_cis.log"])
        assert Electrons(data).alpha() == 5
        assert Electrons(data).beta() == 5
        # Systems with alpha != beta
        data, logfile = getdatafile(QChem, "basicQChem5.4", ["dvb_sp_un.out"])
        assert Electrons(data).alpha() == 35
        assert Electrons(data).beta() == 34
        data, logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_un_sp.log"])
        assert Electrons(data).alpha() == 35
        assert Electrons(data).beta() == 34


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(unittest.makeSuite(ElectronsTest))
