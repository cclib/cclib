# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Run parser unit tests for cclib."""

from __future__ import print_function

import unittest

import numpy

import cclib


class LogfileTest(unittest.TestCase):
    """Unit tests for Logfile class."""

    logfile = cclib.parser.logfileparser.Logfile('')

    def test_float(self):
        """Are floats converted from strings correctly?"""
        self.assertEqual(self.logfile.float("0.0"), 0.0)
        self.assertEqual(self.logfile.float("1.0"), 1.0)
        self.assertEqual(self.logfile.float("-1.0"), -1.0)
        self.assertEqual(self.logfile.float("1.2345E+02"), 123.45)
        self.assertEqual(self.logfile.float("1.2345D+02"), 123.45)
        self.assertTrue(numpy.isnan(self.logfile.float("*")))
        self.assertTrue(numpy.isnan(self.logfile.float("*****")))

    def test_normalisesym(self):
        """Does this method return ERROR in base class?"""
        self.assertTrue("ERROR" in self.logfile.normalisesym(""))


if __name__ == "__main__":

    unittest.main()
