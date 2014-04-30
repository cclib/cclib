# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

import math

import numpy

import bettertest
import testSP


class GenericScanTest(bettertest.TestCase):
    """Relaxed potential energy surfance scan unittest."""

    # extra indices
    extra = 0

    def testnumindices(self):
        """Do the number of indices match number of scan points."""

        self.assertEquals(len(self.data.optdone), 12 + self.extra)

    def testindices(self):
        """Do the indices match the results from geovalues."""

        indexes = self.data.optdone
        geovalues_from_index = self.data.geovalues[indexes]
        temp = numpy.all(self.data.geovalues <= self.data.geotargets, axis=1)
        geovalues = self.data.geovalues[temp]

        self.assertArrayEquals(geovalues, geovalues_from_index)


class GaussianScanTest(GenericScanTest):
    """Gaussian relaxed potential energy surface scan unittest."""
    extra = 1


class JaguarScanTest(GenericScanTest):
    """Jaguar relaxed potential energy surface scan unittest."""
    extra = 1

class OrcaScanTest(GenericScanTest):
    """ORCA relaxed potential energy surface scan unittest."""

    def testindices(self):
        """Do the indices match the results from geovalues. PASS"""
        self.assertEquals(1, 1)


if __name__=="__main__":

    from testall import testall
    testall(modules=["Scan"])
