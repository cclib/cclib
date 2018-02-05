# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test logfiles with core electron data in cclib"""

import os
import unittest

import numpy

from cclib.parser.utils import PeriodicTable


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericCoreTest(unittest.TestCase):
    """Generic core electrons unittest"""

    coredict = {'Mo': 28, 'O':0, 'Cl':10}

    def testcorrect(self):
        """Is coreelectrons equal to what it should be?"""
        pt = PeriodicTable()
        ans = []
        for x in self.data.atomnos:
            ans.append(self.coredict[pt.element[x]])
        ans = numpy.array(ans, "i")
        numpy.testing.assert_array_equal(self.data.coreelectrons, ans)


class ADFCoreTest(GenericCoreTest):
    """Customized core electrons unittest"""

    # For some reason ADF does not have a core in this test for chlorine atoms.
    # This might be fixable in the input.
    coredict = {'Mo': 28, 'O':0, 'Cl':0}


class JaguarCoreTest(GenericCoreTest):
    """Customized core electrons unittest"""

    # This test was done using LanL2DZ instead of the smaller variant.
    coredict = {'Mo': 36, 'O':0, 'Cl':10}

           
if __name__=="__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['Suite'])
    suite.testall()
