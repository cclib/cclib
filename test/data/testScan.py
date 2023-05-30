# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test scan logfiles in cclib"""

import cclib

import numpy
from skip import skipForParser

OPT_DONE = cclib.parser.data.ccData.OPT_DONE
OPT_NEW = cclib.parser.data.ccData.OPT_NEW


def is_optnew(optstatus_value) -> bool:
    return optstatus_value & OPT_NEW == OPT_NEW


def is_optdone(optstatus_value) -> bool:
    return optstatus_value & OPT_DONE == OPT_DONE


class GenericRelaxedScanTest_optdone_bool:
    """Generic relaxed potential energy surface scan unittest."""

    datatype = cclib.parser.data.ccData_optdone_bool

    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    def testoptdone(self, data) -> None:
        """Is the optimization finished?"""
        assert isinstance(data.optdone, bool)
        assert data.optdone

    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    def testindices(self, data) -> None:
        """Do the indices match the results from geovalues."""
        assert data.optdone and numpy.all(data.geovalues[-1] <= data.geotargets)

    @skipForParser("Jaguar", "Not implemented")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("ORCA", "Not implemented")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testoptstatus(self, data) -> None:
        """Does optstatus contain expected values?"""

        # The input and final coordinates were at a stationary points.
        assert is_optnew(data.optstatus[0])
        assert is_optdone(data.optstatus[0])
        assert is_optdone(data.optstatus[-1])


class GenericUnrelaxedScanTest:
    """Generic unrelaxed potential energy surface scan unittest."""

    # extra indices
    extra = 0

    @skipForParser("Jaguar", "Not implemented")
    def testscannames(self, data) -> None:
        assert isinstance(data.scannames, list)

    @skipForParser("ORCA", "Not implemented")
    @skipForParser("Jaguar", "Not implemented")
    def testscanenergies(self, data) -> None:
        assert isinstance(data.scanenergies, list)

        # This checks the order of magnitude, and unit conversion if nothing else.
        numpy.testing.assert_array_less(numpy.array(data.scanenergies), -10000)

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

    # extra indices
    extra = 0

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testnumindices(self, data) -> None:
        """Do the number of indices match number of scan points."""
        assert len(data.optdone) == 12 + self.extra

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
    def testscanenergies(self, data) -> None:
        assert isinstance(data.scanenergies, list)

        # This checks the order of magnitude, and unit conversion if nothing else.
        numpy.testing.assert_array_less(numpy.array(data.scanenergies), -10000)

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
