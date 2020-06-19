# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import unittest
import os, sys
import numpy

from cclib.bridge import cclib2horton
from ..test_data import getdatafile
from cclib.parser.utils import find_package

from numpy.testing import assert_array_almost_equal


class HortonTest(unittest.TestCase):
    """ Tests for the horton bridge in cclib """

    # Both horton and cclib can read in fchk files. The test routine utilizes this fact
    # and compares the attributes that were directly loaded from each package
    # and the attributes that were converted.

    def setUp(self):
        super(HortonTest, self).setUp()

        self.data, self.logfile = getdatafile(
            "Gaussian", "basicGaussian16", ["dvb_un_sp.fchk"]
        )
        datadir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..", "data")
        )
        inputfile = os.path.join(
            datadir, "Gaussian", "basicGaussian16", "dvb_un_sp.fchk"
        )

        self._old_horton = False
        self._found_horton = find_package("horton")
        self._found_iodata = find_package("iodata")

        if self._found_horton:
            try:
                from horton import __version__
            except:
                self._old_horton = (
                    True  # Old horton versions do not have this __version__ attribute
                )
            else:
                if __version__[0] == "2":
                    from horton.io.iodata import IOData
                    from horton import log

                    log.set_level(
                        0
                    )  # This suppresses horton outputs so that they do not flood the test log.
                    self._hortonver = 2
                    self.iodat = IOData.from_file(inputfile)
        if self._found_iodata:
            from iodata import IOData
            from iodata.orbitals import MolecularOrbitals
            from iodata.api import load_one

            self._hortonver = 3

            self.iodat = load_one(filename=inputfile)

    def test_makehorton(self):
        """ Check that the bridge from cclib to horton works correctly """
        # First use `makehorton` function to generate IOData object converted from cclib ccData
        hortonequiv = cclib2horton.makehorton(self.data)

        if self._hortonver == 2:
            # Identify attributes that should be verified
            check = ["pseudo_numbers", "ms2"]  # float or int
            checkArr = [
                "coordinates",
                "numbers",
                "orb_alpha",
                "orb_beta",
                "mulliken_charges",
                "npa_charges",
            ]  # one dimensional arrays
            checkArrArr = ["polar"]  # two dimensional arrays

            for attr in check:
                if hasattr(self.iodat, attr) and hasattr(hortonequiv, attr):
                    self.assertAlmostEqual(
                        getattr(self.iodat, attr),
                        getattr(hortonequiv, attr),
                        delta=1.0e-3,
                    )

            for attr in checkArr:
                if hasattr(self.iodat, attr) and hasattr(hortonequiv, attr):
                    assert_array_almost_equal(
                        getattr(self.iodat, attr), getattr(hortonequiv, attr), decimal=3
                    )

            for attr in checkArrArr:
                if hasattr(self.iodat, attr) and hasattr(hortonequiv, attr):
                    assert_array_almost_equal(
                        getattr(self.iodat, attr)[0],
                        getattr(hortonequiv, attr)[0],
                        decimal=3,
                    )

        elif self._hortonver == 3:
            # Identify attributes that should be verified
            check = ["spinpol"]  # float or int
            checkArr = ["atcoords", "atnums", "atcorenums"]  # one dimensional arrays

            # IOData in horton 3 uses attr and initializes all possible attributes by default value of None. Thus, hasattr cannot be used here by its own.
            for attr in check:
                if (
                    hasattr(self.iodat, attr)
                    and hasattr(hortonequiv, attr)
                    and getattr(self.iodat, attr) != None
                    and getattr(hortonequiv, attr) != None
                ):
                    self.assertAlmostEqual(
                        getattr(self.iodat, attr),
                        getattr(hortonequiv, attr),
                        delta=1.0e-3,
                    )

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

    def test_makecclib(self):
        """ Check that the bridge from horton to cclib works correctly """
        # First use `makecclib` function to generate ccData object converted from horton IOData
        cclibequiv = cclib2horton.makecclib(self.iodat)

        # Identify attributes that should be verified
        check = ["mult", "coreelectrons"]  # float or int
        checkArr = ["atomcoords", "atomnos", "mocoeffs"]  # one dimensional arrays
        checkArrArr = ["polarizability"]  # two dimensional arrays
        checkChg = ["mulliken", "natural"]  # partial charges

        for attr in check:
            if hasattr(self.data, attr) and hasattr(cclibequiv, attr):
                self.assertAlmostEqual(
                    getattr(self.data, attr), getattr(cclibequiv, attr), delta=1.0e-3
                )

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

        if hasattr(self.data, "atomcharges") and hasattr(cclibequiv, "atcomcharges"):
            for chg in checkChg:
                if chg in self.data.atomcharges and chg in cclibequiv.atomcharges:
                    assert_array_almost_equal(
                        self.data.atomcharges[chg][0],
                        cclibequiv.atomcharges[chg][0],
                        decimal=3,
                    )


if __name__ == "__main__":
    unittest.main()
