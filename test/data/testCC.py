# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test coupled cluster logfiles"""

import os
import unittest

import numpy

from skip import skipForParser

__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericCCTest(unittest.TestCase):
    """Generic coupled cluster unittest"""

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testsign(self):
        """Are the coupled cluster corrections negative?"""
        corrections = self.data.ccenergies - self.data.scfenergies
        self.assertTrue(numpy.alltrue(corrections < 0.0))


if __name__ == "__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['CC'])
    suite.testall()
