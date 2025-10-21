# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test scan logfiles in cclib"""

from typing import TYPE_CHECKING

import cclib

import numpy
import pytest
from common import is_optdone, is_optnew
from skip import skipForParser

if TYPE_CHECKING:
    from cclib.parser.data import ccData


class GenericUnrelaxedScanTest:
    """Generic unrelaxed potential energy surface scan unittest."""

    @pytest.fixture
    def extra(self) -> int:
        """extra indices"""
        return 0

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    def testscannames(self, data: "ccData") -> None:
        assert isinstance(data.scannames, list)

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    def testscanenergies(self, data: "ccData") -> None:
        assert isinstance(data.scanenergies, list)

        # This checks the order of magnitude, and unit conversion if nothing else.
        numpy.testing.assert_array_less(
            numpy.array(data.scanenergies), cclib.parser.utils.convertor(-378, "hartree", "eV")
        )

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    def testscanparm(self, data: "ccData") -> None:
        assert isinstance(data.scanparm, list)

        assert len(data.scanparm) == len(data.scannames)

        # Each parameters should have as many values as there are scan
        # energies, or optimized point on the PES.
        for parm in data.scanparm:
            assert len(parm) == len(data.scanenergies)

    @skipForParser("Jaguar", "Not implemented")
    def testscancoords(self, data) -> None:
        """Are the final coordinates for each scan point consistent?"""

        assert isinstance(data.scancoords, numpy.ndarray)
        if hasattr(data, "scanenergies"):
            assert len(data.scancoords) == len(data.scanenergies)
        # In an unrelaxed scan, the only coordinates present in the file
        # should be the ones associated with the coordinates at each scan
        # point.
        assert data.scancoords.shape == data.atomcoords.shape
        numpy.testing.assert_array_equal(data.scancoords, data.atomcoords)


class GenericRelaxedScanTest(GenericUnrelaxedScanTest):
    """Generic relaxed potential energy surface scan unittest."""

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testnumindices(self, data: "ccData", extra: int) -> None:
        """Do the number of indices match number of scan points?"""
        assert len(data.optdone) == 12 + extra

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testindices(self, data: "ccData") -> None:
        """Do the optdone indices match the results from geovalues?"""
        mask_converged_geovalues = numpy.all(data.geovalues <= data.geotargets, axis=1)
        indices_converged_geovalues = [i for i, v in enumerate(mask_converged_geovalues) if v]
        # Depending on the program, it's possible that convergence may be
        # triggered even if all convergence criteria are not met, but if they
        # are met, the program should consider it converged and we should have
        # set optdone.
        assert set(indices_converged_geovalues) - set(data.optdone) == set()

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testoptstatus(self, data: "ccData") -> None:
        """Does optstatus contain expected values?"""
        # The input coordinates were at a stationary point.
        #
        # This depends on the input, and not all currently obey this
        # requirement.
        #
        # assert is_optdone(data.optstatus[0])

        assert len(data.converged_geometries) == len(data.optdone)
        for idone in data.optdone:
            assert is_optdone(data.optstatus[idone])
            if idone != len(data.optstatus) - 1:
                assert is_optnew(data.optstatus[idone + 1])

    @skipForParser("Jaguar", "Not implemented")
    def testscancoords(self, data: "ccData") -> None:
        """Are the final coordinates for each scan point consistent?"""

        assert isinstance(data.scancoords, numpy.ndarray)
        if hasattr(data, "scanenergies"):
            assert len(data.scancoords) == len(data.scanenergies)
        # A relaxed scan is where a geometry optimization is performed for
        # each set of parameters along the scan rather than taking the
        # geometry as-is.  That means each point on the scan is considered
        # done when its geometry optimization has converged.
        assert data.scancoords.shape == (len(data.optdone), data.natom, 3)
        for iscan, idone in enumerate(data.optdone):
            numpy.testing.assert_array_equal(data.scancoords[iscan], data.atomcoords[idone])


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
