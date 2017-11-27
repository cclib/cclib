# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for utilities."""

import unittest

from cclib.parser.utils import WidthSplitter


class WidthSplitterTest(unittest.TestCase):

    def test_default(self):
        """Does the splitter remove empty fields by default properly?"""
        fixed_splitter = WidthSplitter((4, 3, 5, 6, 10, 10, 10, 10, 10, 10))
        line_full = "  60  H 10  s        0.14639   0.00000   0.00000  -0.00000  -0.00000   0.00000"
        line_truncated = "   1  C 1   s       -0.00000  -0.00000   0.00000"
        ref_full = ['60', 'H', '10', 's', '0.14639', '0.00000', '0.00000', '-0.00000', '-0.00000', '0.00000']
        ref_truncated = ['1', 'C', '1', 's', '-0.00000', '-0.00000', '0.00000']
        tokens_full = fixed_splitter.split(line_full)
        tokens_truncated = fixed_splitter.split(line_truncated)
        self.assertEqual(ref_full, tokens_full)
        self.assertEqual(ref_truncated, tokens_truncated)

    def test_no_truncation(self):
        """Does the splitter return even them empty fields when asked?"""
        fixed_splitter = WidthSplitter((4, 3, 5, 6, 10, 10, 10, 10, 10, 10))
        line = "   1  C 1   s       -0.00000  -0.00000   0.00000"
        ref_not_truncated = ['1', 'C', '1', 's', '-0.00000', '-0.00000', '0.00000', '', '', '']
        tokens_not_truncated = fixed_splitter.split(line, truncate=False)
        self.assertEqual(ref_not_truncated, tokens_not_truncated)


if __name__ == "__main__":
    unittest.main()
