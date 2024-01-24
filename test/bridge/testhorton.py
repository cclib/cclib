# -*- coding: utf-8 -*-
#
# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import os

from cclib.bridge import cclib2horton
from cclib.parser.utils import find_package

import numpy
from numpy.testing import assert_array_almost_equal

from ..test_data import getdatafile


class HortonTest:
    """Tests for the Horton bridge in cclib"""

    # Both horton and cclib can read in fchk files. The test routine utilizes this fact
    # and compares the attributes that were directly loaded from each package
    # and the attributes that were converted.

    def setup_method(self) -> None:
        self.data, self.logfile = getdatafile("FChk", "basicGaussian16", ["dvb_un_sp.fchk"])
        datadir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "data"))
        inputfile = os.path.join(datadir, "FChk", "basicGaussian16", "dvb_un_sp.fchk")

        if not find_package("iodata"):
            raise ImportError("Must install iodata to run this test")

        from iodata.api import load_one

        self._hortonver = 3

        self.iodat = load_one(filename=inputfile)

    def test_makehorton(self) -> None:
        """Check that the bridge from cclib to horton works correctly"""
        # First use `makehorton` function to generate IOData object converted from cclib ccData
        hortonequiv = cclib2horton.makehorton(self.data)

        # Identify attributes that should be verified
        check = ["spinpol"]  # float or int
        checkArr = ["atcoords", "atnums", "atcorenums"]  # one dimensional arrays

        # IOData in horton 3 uses attr and initializes all possible attributes by default value
        # of None. Thus, hasattr cannot be used here on its own.
        for attr in check:
            if (
                hasattr(self.iodat, attr)
                and hasattr(hortonequiv, attr)
                and getattr(self.iodat, attr) is not None
                and getattr(hortonequiv, attr) is not None
            ):
                assert abs(getattr(self.iodat, attr) - getattr(hortonequiv, attr)) < 1.0e-3

        for attr in checkArr:
            if (
                hasattr(self.iodat, attr)
                and hasattr(hortonequiv, attr)
                and isinstance(getattr(self.iodat, attr), numpy.ndarray)
                and isinstance(getattr(hortonequiv, attr), numpy.ndarray)
            ):
                assert_array_almost_equal(
                    getattr(self.iodat, attr), getattr(hortonequiv, attr), decimal=3
                )

    def test_makecclib(self) -> None:
        """Check that the bridge from horton to cclib works correctly"""
        # First use `makecclib` function to generate ccData object converted from horton IOData
        cclibequiv = cclib2horton.makecclib(self.iodat)

        # Identify attributes that should be verified
        check = ["mult"]  # float or int
        checkArr = ["atomcoords", "atomnos", "coreelectrons", "mocoeffs"]  # one dimensional arrays
        checkArrArr = ["polarizability"]  # two dimensional arrays
        checkChg = ["mulliken", "natural"]  # atomcharges attribute is dictionary with these keys

        for attr in check:
            if hasattr(self.data, attr) and hasattr(cclibequiv, attr):
                assert abs(getattr(self.data, attr) - getattr(cclibequiv, attr)) < 1.0e-3

        for attr in checkArr:
            if hasattr(self.data, attr) and hasattr(cclibequiv, attr):
                assert_array_almost_equal(
                    getattr(self.data, attr), getattr(cclibequiv, attr), decimal=3
                )

        for attr in checkArrArr:
            if hasattr(self.data, attr) and hasattr(cclibequiv, attr):
                assert_array_almost_equal(
                    getattr(self.data, attr)[0], getattr(cclibequiv, attr)[0], decimal=3
                )

        if hasattr(self.data, "atomcharges") and hasattr(cclibequiv, "atomcharges"):
            for chg in checkChg:
                if chg in self.data.atomcharges and chg in cclibequiv.atomcharges:
                    assert_array_almost_equal(
                        self.data.atomcharges[chg][0], cclibequiv.atomcharges[chg][0], decimal=3
                    )
