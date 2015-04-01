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

import os
import unittest

import numpy

import cclib


class ccIOTest(unittest.TestCase):
    """Unit tests for the ccio module."""

    def test_guess_filetype_fail(self):
        """Does the function fail as expected?"""
        guess = cclib.parser.ccio.guess_filetype
        self.assertIsNone(guess([]))
        self.assertIsNone(guess(None))
        self.assertIsNone(guess(os.devnull))
        self.assertIsNone(guess(['test', 'random', 'quantum chemistry']))

    def test_guess_filetype_program(self):
        """Does the function catch all expected programs?"""
        guess = cclib.parser.ccio.guess_filetype
        self.assertEqual(guess(["Amsterdam Density Functional"]), cclib.parser.ADF)
        self.assertEqual(guess(['Dalton - An Electronic Structure Program']), cclib.parser.DALTON)
        self.assertEqual(guess(['GAMESS']), cclib.parser.GAMESS)
        self.assertEqual(guess(['G A M E S S - U K']), cclib.parser.GAMESSUK)
        self.assertEqual(guess(['Gaussian, Inc.']), cclib.parser.Gaussian)
        self.assertEqual(guess(['Jaguar']), cclib.parser.Jaguar)
        self.assertEqual(guess(['PROGRAM SYSTEM MOLPRO']), cclib.parser.Molpro)
        self.assertEqual(guess(['Northwest Computational Chemistry Package']), cclib.parser.NWChem)
        self.assertEqual(guess(['O   R   C   A']), cclib.parser.ORCA)
        self.assertEqual(guess(["PSI ...Ab Initio Electronic Structure"]), cclib.parser.Psi)
        self.assertEqual(guess(['A Quantum Leap Into The Future Of Chemistry']), cclib.parser.QChem)

    def test_ccread_fail(self):
        """Does the function fail as expected?"""
        ccread = cclib.parser.ccio.ccread
        self.assertIsNone(ccread("", quiet=True))
        self.assertIsNone(ccread([], quiet=True))
        self.assertIsNone(ccread(None, quiet=True))

    def test_ccopen_fail(self):
        """Does the function fail as expected?"""
        ccopen = cclib.parser.ccio.ccread
        self.assertIsNone(ccopen("", quiet=True))
        self.assertIsNone(ccopen([], quiet=True))
        self.assertIsNone(ccopen(None, quiet=True))

    def test_fallback_fail(self):
        """Does the functin fail as expected?"""
        fallback = cclib.parser.ccio.fallback
        self.assertIsNone(fallback(None))


class DataTest(unittest.TestCase):
    """Unit tests for the data module."""

    check_array = ['atomcoords', 'scfenergies']
    check_arrlist = ['mocoeffs', 'scfvalues']
    check_arrdict = ['atomcharges', 'atomspins']

    def setUp(self):
        self.data = cclib.parser.logfileparser.ccData()

    def _set_attributes(self, names, val):
        for n in names:
            setattr(self.data, n, val)

    def test_listify_ndarray(self):
        """Does the method convert ndarrays as expected?"""
        self._set_attributes(self.check_array, numpy.array([1,2,3]))
        self.data.listify()
        for attr in self.check_array:
            self.assertIsInstance(getattr(self.data, attr), list)

    def test_listify_arraylist(self):
        """Does the method convert lists of arrays as expected?"""
        self._set_attributes(self.check_arrlist, [numpy.array([1,2,3]), numpy.array([4,5,6])])
        self.data.listify()
        for attr in self.check_arrlist:
            for a in getattr(self.data, attr):
                self.assertIsInstance(a, list)

    def test_listify_arraydict(self):
        """Does the method convert dicts of arrays as expected?"""
        self._set_attributes(self.check_arrdict, {1: numpy.array([1,2,3]), 2: numpy.array([4,5,6])})
        self.data.listify()
        for attr in self.check_arrdict:
            for a in getattr(self.data, attr).values():
                self.assertIsInstance(a, list)

    def test_arrayify_ndarray(self):
        """Does the method convert lists as expected?"""
        self._set_attributes(self.check_array, [1,2,3])
        self.data.arrayify()
        for attr in self.check_array:
            self.assertIsInstance(getattr(self.data, attr), numpy.ndarray)

    def test_arrayify_arraylist(self):
        """Does the method convert lists of lists as expected?"""
        self._set_attributes(self.check_arrdict, [[1,2,3], [4,5,6]])
        self.data.arrayify()
        for attr in self.check_arrdict:
            for a in getattr(self.data, attr):
                self.assertIsInstance(a, numpy.ndarray)

    def test_arrayify_arraydict(self):
        """Does the method convert dicts of lists as expected?"""
        self._set_attributes(self.check_arrdict, {1: [1,2,3], 2: [4,5,6]})
        self.data.arrayify()
        for attr in self.check_arrdict:
            for a in getattr(self.data, attr).values():
                self.assertIsInstance(a, numpy.ndarray)


class LogfileTest(unittest.TestCase):
    """Unit tests for the logfile module."""

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
