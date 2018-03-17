# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for main scripts (ccget, ccwrite)."""
import os
import unittest
from io import StringIO

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch


__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..", "data")


INPUT_FILE = os.path.join(
    __datadir__,
    'ADF/basicADF2007.01/dvb_gopt.adfout'
)
CJSON_OUTPUT_FILENAME = 'dvb_gopt.cjson'


class ccgetTest(unittest.TestCase):

    def setUp(self):
        try:
            from scripts import ccget
        except ImportError:
            self.fail("ccget cannot be imported")

        self.main = ccget.ccget

    @patch("scripts.ccget.sys.argv", ["ccget"])
    def test_empty_argv(self):
        """Does the script fail as expected if called without parameters?"""
        with self.assertRaises(SystemExit):
            self.main()


@patch("scripts.ccwrite.ccwrite")
class ccwriteTest(unittest.TestCase):

    def setUp(self):
        try:
            from scripts import ccwrite
        except ImportError:
            self.fail("ccwrite cannot be imported")

        self.main = ccwrite.main

    @patch('scripts.ccwrite.sys.argv', ['ccwrite'])
    def test_empty_argv(self, mock_ccwrite):
        """Does the script fail as expected if called without parameters?"""
        with self.assertRaises(SystemExit):
            self.main()

    @patch(
        "scripts.ccwrite.sys.argv",
        ["ccwrite", "cjson", INPUT_FILE]
    )
    def test_cjson_output(self, mock_ccwrite):
        """is ccwrite called with the given parameters?"""
        self.main()

        self.assertEqual(mock_ccwrite.call_count, 1)
        ccwrite_call_args, ccwrite_call_kwargs = mock_ccwrite.call_args
        self.assertEqual(ccwrite_call_args[1], 'cjson')
        self.assertEqual(ccwrite_call_args[2], CJSON_OUTPUT_FILENAME)


if __name__ == "__main__":
    unittest.main()
