# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2014,2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test scan logfiles in cclib"""

import os
import unittest

import numpy

from skip import skipForParser


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericScanTest(unittest.TestCase):
    """Generic relaxed potential energy surfance scan unittest"""

    # extra indices
    extra = 0

    def testnumindices(self):
        """Do the number of indices match number of scan points."""

        if self.data._attrtypes["optdone"] is bool:
            self.assertEquals(self.data.optdone, True)
        else:
            self.assertEquals(len(self.data.optdone), 12 + self.extra)

    def testindices(self):
        """Do the indices match the results from geovalues."""

        if self.data._attrtypes["optdone"] is bool:
            assert self.data.optdone and numpy.all(self.data.geovalues[-1] <= self.data.geotargets)
        else:
            indexes = self.data.optdone
            geovalues_from_index = self.data.geovalues[indexes]
            temp = numpy.all(self.data.geovalues <= self.data.geotargets, axis=1)
            geovalues = self.data.geovalues[temp]
            numpy.testing.assert_array_equal(geovalues, geovalues_from_index)


class GaussianScanTest(GenericScanTest):
    """Customized relaxed potential energy surface scan unittest"""
    extra = 1


class JaguarScanTest(GenericScanTest):
    """Customized relaxed potential energy surface scan unittest"""
    extra = 1


class OrcaScanTest(GenericScanTest):
    """Customized relaxed potential energy surface scan unittest"""
    extra = 1


if __name__=="__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['Scan'])
    suite.testall()
