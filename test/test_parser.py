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
    """Unit tests for the logfile module."""

    def test_float(self):
        """Are floats converted from strings correctly?"""
        float = cclib.parser.logfileparser.Logfile('').float
        self.assertEqual(float("0.0"), 0.0)
        self.assertEqual(float("1.0"), 1.0)
        self.assertEqual(float("-1.0"), -1.0)
        self.assertEqual(float("1.2345E+02"), 123.45)
        self.assertEqual(float("1.2345D+02"), 123.45)
        self.assertTrue(numpy.isnan(float("*")))
        self.assertTrue(numpy.isnan(float("*****")))

    def test_normalisesym(self):
        """Does this method return ERROR in base class?"""
        normalisesym = cclib.parser.logfileparser.Logfile('').normalisesym
        self.assertTrue("ERROR" in normalisesym(""))


class UtilsTest(unittest.TestCase):
    """Unit tests for the utils module."""

    def test_convertor_basic(self):
        """Are some basic conversions correct?"""
        convertor = cclib.parser.utils.convertor
        self.assertAlmostEqual(convertor(1.89, "bohr", "Angstrom"), 1.0, places=3)
        self.assertAlmostEqual(convertor(0.529, "Angstrom", "bohr"), 1.0, places=3)
        self.assertAlmostEqual(convertor(627.5, "kcal", "hartree"), 1.0, places=3)

    def test_convertor_pairs(self):
        """Do flipped conversions correspond to each other?"""

        convertor = cclib.parser.utils.convertor

        pairs_proportional = (
            ("Angstrom", "bohr"),
            ("cm-1", "eV"), ("cm-1", "kcal"), ("eV", "kJmol-1"),
            ("coulomb", "e"),
        )
        pairs_inverse = (("nm", "cm-1"),)

        for unit1, unit2 in pairs_proportional:
            conv1 = convertor(1.0, unit1, unit2)
            conv2 = 1.0 / convertor(1.0, unit2, unit1)
            self.assertAlmostEqual((conv1 - conv2) / conv1, 0.0)
        for unit1, unit2 in pairs_inverse:
            conv1 = convertor(1.0, unit1, unit2)
            conv2 = convertor(1.0, unit2, unit1)
            self.assertAlmostEqual((conv1 - conv2) / conv1, 0.0)

    def test_periodictable(self):
        """Does the periodic table work correctly?"""
        pt = cclib.parser.utils.PeriodicTable()
        self.assertEquals(pt.element[6], 'C')
        self.assertEquals(pt.number['C'], 6)
        self.assertEquals(pt.element[44], 'Ru')
        self.assertEquals(pt.number['Au'], 79)


if __name__ == "__main__":

    unittest.main()
