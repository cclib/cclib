# -*- coding: utf-8 -*-
#
# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the various population analyses (MPA, LPA, CSPA) in cclib"""

import sys

from cclib.method import Orbitals
from cclib.parser import Gaussian, Psi4

sys.path.insert(1, "..")

from ..test_data import getdatafile


class RestrictedCalculationTest:
    """Check retricted calculation."""

    def setup_method(self) -> None:
        self.data, self.logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])

    def test_closed_shell(self) -> None:
        assert Orbitals(self.data).closed_shell()


class UnrestrictedCalculationTest:
    """Check unrestricted calculation."""

    def setup_method(self) -> None:
        self.data, self.logfile = getdatafile(Gaussian, "basicGaussian09", ["dvb_un_sp.log"])

    def test_closed_shell(self) -> None:
        assert not Orbitals(self.data).closed_shell()


class RestrictedOpenShellCalculationTest:
    """Check restricted open shell calcualtion."""

    def setup_method(self) -> None:
        self.data, self.logfile = getdatafile(Psi4, "basicPsi4-1.3.1", ["dvb_sp_rohf.out"])

    def test_closed_shel(self) -> None:
        assert not Orbitals(self.data).closed_shell()


# TODO: add a case (regression) with an unrestricted calculation for a closed shell system.
# For example, in regressions: Gaussian/Gaussian03/Mo4OSibdt2


tests = [RestrictedCalculationTest, UnrestrictedCalculationTest]
