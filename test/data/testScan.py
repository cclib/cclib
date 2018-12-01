# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test scan logfiles in cclib"""

import os
import unittest

import numpy
import cclib

from skip import skipForParser


__filedir__ = os.path.realpath(os.path.dirname(__file__))


OPT_DONE = cclib.parser.data.ccData.OPT_DONE
OPT_NEW = cclib.parser.data.ccData.OPT_NEW


class GenericScanTestBase(unittest.TestCase):
    """Base relaxed potential energy surface scan unittest."""

    def assertOptNew(self, optstatus_value):
        return optstatus_value & OPT_NEW == OPT_NEW

    def assertOptDone(self, optstatus_value):
        return optstatus_value & OPT_DONE == OPT_DONE


class GenericScanTest_optdone_bool(GenericScanTestBase):
    """Generic relaxed potential energy surface scan unittest."""

    datatype = cclib.parser.data.ccData_optdone_bool

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    def testoptdone(self):
        """Is the optimization finished?"""
        self.assertIsInstance(self.data.optdone, bool)
        self.assertEqual(self.data.optdone, True)

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    def testindices(self):
        """Do the indices match the results from geovalues."""
        assert self.data.optdone and numpy.all(self.data.geovalues[-1] <= self.data.geotargets)

    @skipForParser("Jaguar", "Not implemented")
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser("ORCA", "Not implemented")
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testoptstatus(self):
        """Does optstatus contain expected values?"""

        # The input and final coordinates were at a stationary points.
        self.assertOptNew(self.data.optstatus[0])
        self.assertOptDone(self.data.optstatus[0])
        self.assertOptDone(self.data.optstatus[-1])


class GenericScanTest(GenericScanTestBase):
    """Generic relaxed potential energy surface scan unittest."""

    # extra indices
    extra = 0
    
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testnumindices(self):
        """Do the number of indices match number of scan points."""
        self.assertEqual(len(self.data.optdone), 12 + self.extra)

    @skipForParser("Jaguar", "Does not work as expected")    
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser("ORCA", "Does not work as expected")
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testindices(self):
        """Do the indices match the results from geovalues."""
        indexes = self.data.optdone
        geovalues_from_index = self.data.geovalues[indexes]
        temp = numpy.all(self.data.geovalues <= self.data.geotargets, axis=1)
        geovalues = self.data.geovalues[temp]
        numpy.testing.assert_array_equal(geovalues, geovalues_from_index)

    @skipForParser("Jaguar", "Not implemented")
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser("ORCA", "Not implemented")
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testoptstatus(self):
        """Does optstatus contain expected values?"""
        OPT_NEW = self.data.OPT_NEW
        OPT_DONE = self.data.OPT_DONE
        # The input coordinates were at a stationary point.
        self.assertOptDone(self.data.optstatus[0])

        self.assertEqual(len(self.data.converged_geometries), len(self.data.optdone))
        for idone in self.data.optdone:
            self.assertOptDone(self.data.optstatus[idone])
            if idone != len(self.data.optstatus) - 1:
                self.assertOptNew(self.data.optstatus[idone+1])


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
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['Scan'])
    suite.testall()
