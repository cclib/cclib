# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test scan logfiles in cclib"""

import cclib

import numpy
import pytest
from skip import skipForParser

OPT_DONE = cclib.parser.data.ccData.OPT_DONE
OPT_NEW = cclib.parser.data.ccData.OPT_NEW


def is_optnew(optstatus_value) -> bool:
    return optstatus_value & OPT_NEW == OPT_NEW


def is_optdone(optstatus_value) -> bool:
    return optstatus_value & OPT_DONE == OPT_DONE


class GenericUnrelaxedScanTest:
    """Generic unrelaxed potential energy surface scan unittest."""

    @pytest.fixture
    def extra(self) -> int:
        """extra indices"""
        return 0

    @skipForParser("Jaguar", "Not implemented")
    def testscannames(self, data) -> None:
        assert isinstance(data.scannames, list)

    @skipForParser("ORCA", "Not implemented")
    @skipForParser("Jaguar", "Not implemented")
    def testscanenergies(self, data) -> None:
        assert isinstance(data.scanenergies, list)

        # This checks the order of magnitude, and unit conversion if nothing else.
        numpy.testing.assert_array_less(numpy.array(data.scanenergies), -378)

    @skipForParser("ORCA", "Not implemented")
    @skipForParser("Jaguar", "Not implemented")
    def testscanparm(self, data) -> None:
        assert isinstance(data.scanparm, list)

        # Each parameters should have as many values as there are scan
        # energies, or optimized point on the PES.
        for parm in data.scanparm:
            assert len(parm) == len(data.scanenergies)


class GenericRelaxedScanTest(GenericUnrelaxedScanTest):
    """Generic relaxed potential energy surface scan unittest."""

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testnumindices(self, data, extra) -> None:
        """Do the number of indices match number of scan points."""
        assert len(data.optdone) == 12 + extra

    @skipForParser("Jaguar", "Does not work as expected")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("ORCA", "Does not work as expected")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testindices(self, data) -> None:
        """Do the indices match the results from geovalues."""
        indexes = data.optdone
        geovalues_from_index = data.geovalues[indexes]
        temp = numpy.all(data.geovalues <= data.geotargets, axis=1)
        geovalues = data.geovalues[temp]
        numpy.testing.assert_array_equal(geovalues, geovalues_from_index)

    @skipForParser("Jaguar", "Not implemented")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("ORCA", "Not implemented")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testoptstatus(self, data) -> None:
        """Does optstatus contain expected values?"""
        # The input coordinates were at a stationary point.
        assert is_optdone(data.optstatus[0])

        assert len(data.converged_geometries) == len(data.optdone)
        for idone in data.optdone:
            assert is_optdone(data.optstatus[idone])
            if idone != len(data.optstatus) - 1:
                assert is_optnew(data.optstatus[idone + 1])

    @skipForParser("Jaguar", "Not implemented")
    def testscannames(self, data) -> None:
        assert isinstance(data.scannames, list)

    @skipForParser("Jaguar", "Not implemented")
    @skipForParser("ORCA", "Not implemented")
    def testscanparm(self, data) -> None:
        assert isinstance(data.scanparm, list)

        # Each parameters should have as many values as there are scan
        # energies, or optimized point on the PES.
        for parm in data.scanparm:
            assert len(parm) == len(data.scanenergies)


class GaussianUnrelaxedScanTest(GenericUnrelaxedScanTest):
    """Customized unrelaxed potential energy surface scan unittest"""

    @pytest.fixture
    def extra(self) -> int:
        """extra indices"""
        return 1


class GaussianRelaxedScanTest(GenericRelaxedScanTest):
    """Customized relaxed potential energy surface scan unittest"""

    @pytest.fixture
    def extra(self) -> int:
        """extra indices"""
        return 1


class JaguarRelaxedScanTest(GenericRelaxedScanTest):
    """Customized relaxed potential energy surface scan unittest"""

    @pytest.fixture
    def extra(self) -> int:
        """extra indices"""
        return 1
