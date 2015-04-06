# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2007,2012,2014,2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test coupled cluster logfiles"""

import os
import unittest

import numpy


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericCCTest(unittest.TestCase):
    """Generic coupled cluster unittest"""

    def testsign(self):
        """Are the coupled cluster corrections negative?"""
        corrections = self.data.ccenergies - self.data.scfenergies
        self.failUnless(numpy.alltrue(corrections < 0.0))


if __name__ == "__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['CC'])
    suite.testall()
