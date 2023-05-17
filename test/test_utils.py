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
    def test_float_basic(self) -> None:
        """Are floats converted from strings correctly?"""
        assert utils.float("0.0") == 0.0
        assert utils.float("1.0") == 1.0
        assert utils.float("-1.0") == -1.0

    def test_float_numeric_format(self) -> None:
        """Does numeric formatting get converted correctly?"""
        assert utils.float("1.2345E+02") == 123.45
        assert utils.float("1.2345D+02") == 123.45

    def test_float_stars(self) -> None:
        """Does the function return nan for stars?"""
        assert numpy.isnan(utils.float("*"))
        assert numpy.isnan(utils.float("*****"))


class ConvertorTest(unittest.TestCase):

    def test_convertor(self) -> None:
        assert f"{utils.convertor(8.0, 'eV', 'wavenumber'):.3f}" == "64524.354"


class GetRotationTest(unittest.TestCase):
    delta = 1e-14

    def setUp(self) -> None:
        self.r = scipy.spatial.transform.Rotation.from_euler('xyz', [15, 25, 35], degrees=True)
        self.t = numpy.array([-1, 0, 2])
        self.a = numpy.array([[1., 1., 1.],
                              [0., 1., 2.],
                              [0., 0., 0.],
                              [0., 0., 4.]])
        self.b = self.r.apply(self.a + self.t)

    def test_default(self) -> None:
        """Is the rotation is correct?"""
        _r = utils.get_rotation(self.a, self.b)
        # as_dcm is renamed to from_matrix in scipy 1.4.0 and will be removed in sicpy 1.6.0
        if hasattr(self.r, "as_matrix"):
            numpy.testing.assert_allclose(self.r.as_matrix(), _r.as_matrix(), atol=self.delta)
        else:
            numpy.testing.assert_allclose(self.r.as_dcm(), _r.as_dcm(), atol=self.delta)

    def test_two_atoms(self) -> None:
        """Is the rotation is correct for 2 atoms?"""
        a2 = self.a[:2]
        b2 = self.b[:2]
        rotated_diff = self.r.apply(a2) - utils.get_rotation(a2, b2).apply(a2)
        # rotated_diff should be translation
        numpy.testing.assert_allclose(rotated_diff[0], rotated_diff[1], atol=self.delta)

    def test_one_atom(self) -> None:
        """Is the rotation is identity for 1 atom?"""
        a1 = self.a[:1]
        b1 = self.b[:1]
        if hasattr(self.r, "as_matrix"):
            numpy.testing.assert_allclose(numpy.eye(3), utils.get_rotation(a1, b1).as_matrix(), atol=self.delta)
        else:
            numpy.testing.assert_allclose(numpy.eye(3), utils.get_rotation(a1, b1).as_dcm(), atol=self.delta)


class PeriodicTableTest(unittest.TestCase):

    def setUp(self) -> None:
        self.t = utils.PeriodicTable()

    def test_periodictable(self) -> None:
        assert self.t.element[6] == 'C'
        assert self.t.number['C'] == 6
        assert self.t.element[44] == 'Ru'
        assert self.t.number['Au'] == 79


class WidthSplitterTest(unittest.TestCase):

    def test_default(self) -> None:
        """Does the splitter remove empty fields by default properly?"""
        fixed_splitter = utils.WidthSplitter((4, 3, 5, 6, 10, 10, 10, 10, 10, 10))
        line_full = "  60  H 10  s        0.14639   0.00000   0.00000  -0.00000  -0.00000   0.00000"
        line_truncated = "   1  C 1   s       -0.00000  -0.00000   0.00000"
        ref_full = ['60', 'H', '10', 's', '0.14639', '0.00000', '0.00000', '-0.00000', '-0.00000', '0.00000']
        ref_truncated = ['1', 'C', '1', 's', '-0.00000', '-0.00000', '0.00000']
        tokens_full = fixed_splitter.split(line_full)
        tokens_truncated = fixed_splitter.split(line_truncated)
        assert ref_full == tokens_full
        assert ref_truncated == tokens_truncated

    def test_no_truncation(self) -> None:
        """Does the splitter return even the empty fields when asked?"""
        fixed_splitter = utils.WidthSplitter((4, 3, 5, 6, 10, 10, 10, 10, 10, 10))
        line = "   1  C 1   s       -0.00000  -0.00000   0.00000"
        ref_not_truncated = ['1', 'C', '1', 's', '-0.00000', '-0.00000', '0.00000', '', '', '']
        tokens_not_truncated = fixed_splitter.split(line, truncate=False)
        assert ref_not_truncated == tokens_not_truncated


class SymmetrizeTest(unittest.TestCase):

    def test_dim_from_tblock_size(self) -> None:
        assert utils._dim_from_tblock_size(1) == 1
        # This isn't possible until we fully move to pytest.
        # with pytest.raises(
        #     RuntimeError,
        #     match="The number of elements (2) isn't possible for a matrix triangle"
        # ):
        #     assert utils._dim_from_tblock_size(2)
        assert utils._dim_from_tblock_size(3) == 2
        assert utils._dim_from_tblock_size(6) == 3
        assert utils._dim_from_tblock_size(10) == 4

    def test_block_to_matrix(self) -> None:
        inp = numpy.array([1, 2, 3, 4, 5, 6], dtype=int)
        ref = numpy.array([[1, 2, 4], [2, 3, 5], [4, 5, 6]], dtype=int)
        numpy.testing.assert_equal(utils.block_to_matrix(inp), ref)

    def test_symmetrize(self) -> None:
        inp = numpy.array([[1, 9, 7],
                           [4, 8, 3],
                           [6, 2, 5]], dtype=int)
        ref_lower = numpy.array([[1, 4, 6],
                                 [4, 8, 2],
                                 [6, 2, 5]], dtype=int)
        ref_upper = numpy.array([[1, 9, 7],
                                 [9, 8, 3],
                                 [7, 3, 5]], dtype=int)
        numpy.testing.assert_equal(utils.symmetrize(inp, "lower"), ref_lower)
        numpy.testing.assert_equal(utils.symmetrize(inp, "upper"), ref_upper)

if __name__ == "__main__":
    unittest.main()
