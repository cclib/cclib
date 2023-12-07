# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test scan logfiles in cclib"""

import os
import unittest

import cclib

import numpy
from skip import skipForParser

__filedir__ = os.path.realpath(os.path.dirname(__file__))


OPT_DONE = cclib.parser.data.ccData.OPT_DONE
OPT_NEW = cclib.parser.data.ccData.OPT_NEW


def is_optnew(optstatus_value) -> bool:
    return optstatus_value & OPT_NEW == OPT_NEW


def is_optdone(optstatus_value) -> bool:
    return optstatus_value & OPT_DONE == OPT_DONE


class GenericUnrelaxedScanTest(unittest.TestCase):
    """Generic unrelaxed potential energy surface scan unittest."""

    # extra indices
    extra = 0

    @skipForParser("Jaguar", "Not implemented")
    def testscannames(self):
        assert isinstance(self.data.scannames, list)

    @skipForParser("ORCA", "Not implemented")
    @skipForParser("Jaguar", "Not implemented")
    def testscanenergies(self):
        assert isinstance(self.data.scanenergies, list)

        # This checks the order of magnitude, and unit conversion if nothing else.
        numpy.testing.assert_array_less(numpy.array(self.data.scanenergies), -10000)

    @skipForParser("ORCA", "Not implemented")
    @skipForParser("Jaguar", "Not implemented")
    def testscanparm(self):
        assert isinstance(self.data.scanparm, list)

        # Each parameters should have as many values as there are scan
        # energies, or optimized point on the PES.
        for parm in self.data.scanparm:
            assert len(parm) == len(self.data.scanenergies)


class GenericRelaxedScanTest(GenericUnrelaxedScanTest):
    """Generic relaxed potential energy surface scan unittest."""

    # extra indices
    extra = 0

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testnumindices(self):
        """Do the number of indices match number of scan points."""
        assert len(self.data.optdone) == 12 + self.extra

    @skipForParser("Jaguar", "Does not work as expected")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("ORCA", "Does not work as expected")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testindices(self):
        """Do the indices match the results from geovalues."""
        indexes = self.data.optdone
        geovalues_from_index = self.data.geovalues[indexes]
        temp = numpy.all(self.data.geovalues <= self.data.geotargets, axis=1)
        geovalues = self.data.geovalues[temp]
        numpy.testing.assert_array_equal(geovalues, geovalues_from_index)

    @skipForParser("Jaguar", "Not implemented")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("ORCA", "Not implemented")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testoptstatus(self):
        """Does optstatus contain expected values?"""
        # The input coordinates were at a stationary point.
        assert is_optdone(self.data.optstatus[0])

        assert len(self.data.converged_geometries) == len(self.data.optdone)
        for idone in self.data.optdone:
            assert is_optdone(self.data.optstatus[idone])
            if idone != len(self.data.optstatus) - 1:
                assert is_optnew(self.data.optstatus[idone + 1])

    @skipForParser("Jaguar", "Not implemented")
    def testscannames(self):
        assert isinstance(self.data.scannames, list)

    @skipForParser("Jaguar", "Not implemented")
    @skipForParser("ORCA", "Not implemented")
    def testscanenergies(self):
        assert isinstance(self.data.scanenergies, list)

        # This checks the order of magnitude, and unit conversion if nothing else.
        numpy.testing.assert_array_less(numpy.array(self.data.scanenergies), -10000)

    @skipForParser("Jaguar", "Not implemented")
    @skipForParser("ORCA", "Not implemented")
    def testscanparm(self):
        assert isinstance(self.data.scanparm, list)

        # Each parameters should have as many values as there are scan
        # energies, or optimized point on the PES.
        for parm in self.data.scanparm:
            assert len(parm) == len(self.data.scanenergies)


class GaussianUnrelaxedScanTest(GenericUnrelaxedScanTest):
    """Customized unrelaxed potential energy surface scan unittest"""

    extra = 1


class GaussianRelaxedScanTest(GenericRelaxedScanTest):
    """Customized relaxed potential energy surface scan unittest"""

    extra = 1


class JaguarRelaxedScanTest(GenericRelaxedScanTest):
    """Customized relaxed potential energy surface scan unittest"""

    extra = 1


class OrcaRelaxedScanTest(GenericRelaxedScanTest):
    """Customized relaxed potential energy surface scan unittest"""

    extra = 1


if __name__ == "__main__":
    import sys

    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite

    suite = DataSuite(["Scan"])
    suite.testall()
