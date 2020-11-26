# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for utilities."""

import unittest

from cclib.parser import utils

import numpy
import scipy.spatial.transform


class FloatTest(unittest.TestCase):
    def test_float_basic(self):
        """Are floats converted from strings correctly?"""
        self.assertEqual(utils.float("0.0"), 0.0)
        self.assertEqual(utils.float("1.0"), 1.0)
        self.assertEqual(utils.float("-1.0"), -1.0)

    def test_float_numeric_format(self):
        """Does numeric formatting get converted correctly?"""
        self.assertEqual(utils.float("1.2345E+02"), 123.45)
        self.assertEqual(utils.float("1.2345D+02"), 123.45)

    def test_float_stars(self):
        """Does the function return nan for stars?"""
        self.assertTrue(numpy.isnan(utils.float("*")))
        self.assertTrue(numpy.isnan(utils.float("*****")))


class ConvertorTest(unittest.TestCase):

    def test_convertor(self):
        self.assertEqual("%.3f" % utils.convertor(8.0, "eV", "wavenumber"),
                         "64524.354")


class GetRotationTest(unittest.TestCase):
    delta = 1e-14

    def setUp(self):
        self.r = scipy.spatial.transform.Rotation.from_euler('xyz', [15, 25, 35], degrees=True)
        self.t = numpy.array([-1, 0, 2])
        self.a = numpy.array([[1., 1., 1.],
                              [0., 1., 2.],
                              [0., 0., 0.],
                              [0., 0., 4.]])
        self.b = self.r.apply(self.a + self.t)

    def test_default(self):
        """Is the rotation is correct?"""
        _r = utils.get_rotation(self.a, self.b)
        # as_dcm is renamed to from_matrix in scipy 1.4.0 and will be removed in sicpy 1.6.0
        if hasattr(self.r, "as_matrix"):
            numpy.testing.assert_allclose(self.r.as_matrix(), _r.as_matrix(), atol=self.delta)
        else:
            numpy.testing.assert_allclose(self.r.as_dcm(), _r.as_dcm(), atol=self.delta)

    def test_two_atoms(self):
        """Is the rotation is correct for 2 atoms?"""
        a2 = self.a[:2]
        b2 = self.b[:2]
        rotated_diff = self.r.apply(a2) - utils.get_rotation(a2, b2).apply(a2)
        # rotated_diff should be translation
        numpy.testing.assert_allclose(rotated_diff[0], rotated_diff[1], atol=self.delta)

    def test_one_atom(self):
        """Is the rotation is identity for 1 atom?"""
        a1 = self.a[:1]
        b1 = self.b[:1]
        if hasattr(self.r, "as_matrix"):
            numpy.testing.assert_allclose(numpy.eye(3), utils.get_rotation(a1, b1).as_matrix(), atol=self.delta)
        else:
            numpy.testing.assert_allclose(numpy.eye(3), utils.get_rotation(a1, b1).as_dcm(), atol=self.delta)


class PeriodicTableTest(unittest.TestCase):

    def setUp(self):
        self.t = utils.PeriodicTable()

    def test_periodictable(self):
        self.assertEqual(self.t.element[6], 'C')
        self.assertEqual(self.t.number['C'], 6)
        self.assertEqual(self.t.element[44], 'Ru')
        self.assertEqual(self.t.number['Au'], 79)


class WidthSplitterTest(unittest.TestCase):

    def test_default(self):
        """Does the splitter remove empty fields by default properly?"""
        fixed_splitter = utils.WidthSplitter((4, 3, 5, 6, 10, 10, 10, 10, 10, 10))
        line_full = "  60  H 10  s        0.14639   0.00000   0.00000  -0.00000  -0.00000   0.00000"
        line_truncated = "   1  C 1   s       -0.00000  -0.00000   0.00000"
        ref_full = ['60', 'H', '10', 's', '0.14639', '0.00000', '0.00000', '-0.00000', '-0.00000', '0.00000']
        ref_truncated = ['1', 'C', '1', 's', '-0.00000', '-0.00000', '0.00000']
        tokens_full = fixed_splitter.split(line_full)
        tokens_truncated = fixed_splitter.split(line_truncated)
        self.assertEqual(ref_full, tokens_full)
        self.assertEqual(ref_truncated, tokens_truncated)

    def test_no_truncation(self):
        """Does the splitter return even the empty fields when asked?"""
        fixed_splitter = utils.WidthSplitter((4, 3, 5, 6, 10, 10, 10, 10, 10, 10))
        line = "   1  C 1   s       -0.00000  -0.00000   0.00000"
        ref_not_truncated = ['1', 'C', '1', 's', '-0.00000', '-0.00000', '0.00000', '', '', '']
        tokens_not_truncated = fixed_splitter.split(line, truncate=False)
        self.assertEqual(ref_not_truncated, tokens_not_truncated)


if __name__ == "__main__":
    unittest.main()
