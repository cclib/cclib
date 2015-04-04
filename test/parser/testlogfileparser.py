# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Unit tests for parser logfileparser module."""

import unittest

import numpy

import cclib


class LogfileTest(unittest.TestCase):
    """Unit tests for the Logfile class."""

    def test_float_basic(self):
        """Are floats converted from strings correctly?"""
        float = cclib.parser.logfileparser.Logfile('').float
        self.assertEqual(float("0.0"), 0.0)
        self.assertEqual(float("1.0"), 1.0)
        self.assertEqual(float("-1.0"), -1.0)

    def test_float_numeric_format(self):
        """Does numeric formatting get converted correctly?"""
        float = cclib.parser.logfileparser.Logfile('').float
        self.assertEqual(float("1.2345E+02"), 123.45)
        self.assertEqual(float("1.2345D+02"), 123.45)

    def test_float_stars(self):
        """Does the function return nan for stars?"""
        float = cclib.parser.logfileparser.Logfile('').float
        self.assertTrue(numpy.isnan(float("*")))
        self.assertTrue(numpy.isnan(float("*****")))

    def test_normalisesym_base_class_error(self):
        """Does this method return ERROR in base class?"""
        normalisesym = cclib.parser.logfileparser.Logfile('').normalisesym
        self.assertTrue("ERROR" in normalisesym(""))


if __name__ == "__main__":
    unittest.main()
