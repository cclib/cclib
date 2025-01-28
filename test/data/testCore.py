# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test logfiles with core electron data in cclib"""

from cclib.parser.utils import PeriodicTable

import numpy
from skip import skipForParser


class GenericCoreTest:
    """Generic core electrons unittest"""

    coredict = {"Mo": 28, "O": 0, "Cl": 10}
    charge = -2

    @skipForParser("CFOUR", "This parser is still being developed")
    @skipForParser("FChk", "Core electrons do not seem to be available")
    def testcorrect(self, data) -> None:
        """Is coreelectrons equal to what it should be?"""
        pt = PeriodicTable()
        ans = []
        for x in data.atomnos:
            ans.append(self.coredict[pt.element[x]])
        ans = numpy.array(ans, "i")
        numpy.testing.assert_array_equal(data.coreelectrons, ans)

    def testcharge(self, data) -> None:
        """Is the total charge correct?"""
        assert data.charge == self.charge


class ADFCoreTest(GenericCoreTest):
    """Customized core electrons unittest"""

    # For some reason ADF does not have a core in this test for chlorine atoms.
    # This might be fixable in the input.
    coredict = {"Mo": 28, "O": 0, "Cl": 0}

    # These calculations were run on the neutral complex.
    charge = 0


class JaguarCoreTest(GenericCoreTest):
    """Customized core electrons unittest"""

    # This test was done using LanL2DZ instead of the smaller variant.
    coredict = {"Mo": 36, "O": 0, "Cl": 10}
