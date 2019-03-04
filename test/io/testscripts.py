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

from six import add_move, MovedModule
add_move(MovedModule('mock', 'mock', 'unittest.mock'))
from six.moves import mock


__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..", "data")


INPUT_FILE = os.path.join(
    __datadir__,
    'ADF/basicADF2007.01/dvb_gopt.adfout'
)
CJSON_OUTPUT_FILENAME = 'dvb_gopt.cjson'


@mock.patch("cclib.scripts.ccget.ccread")
class ccgetTest(unittest.TestCase):

    def setUp(self):
        try:
            from cclib.scripts import ccget
        except ImportError:
            self.fail("ccget cannot be imported")

        self.main = ccget.ccget

    @mock.patch("cclib.scripts.ccget.sys.argv", ["ccget"])
    def test_empty_argv(self, mock_ccread):
        """Does the script fail as expected if called without parameters?"""
        with self.assertRaises(SystemExit):
            self.main()

    @mock.patch(
        "cclib.scripts.ccget.sys.argv",
        ["ccget", "atomcoords", INPUT_FILE]
    )
    def test_ccread_invocation(self, mock_ccread):
        self.main()

        self.assertEqual(mock_ccread.call_count, 1)
        ccread_call_args, ccread_call_kwargs = mock_ccread.call_args
        self.assertEqual(ccread_call_args[0], INPUT_FILE)


@mock.patch("cclib.scripts.ccwrite.ccwrite")
class ccwriteTest(unittest.TestCase):

    def setUp(self):
        try:
            from cclib.scripts import ccwrite
        except ImportError:
            self.fail("ccwrite cannot be imported")

        self.main = ccwrite.main

    @mock.patch('cclib.scripts.ccwrite.sys.argv', ['ccwrite'])
    def test_empty_argv(self, mock_ccwrite):
        """Does the script fail as expected if called without parameters?"""
        with self.assertRaises(SystemExit):
            self.main()

    @mock.patch(
        "cclib.scripts.ccwrite.sys.argv",
        ["ccwrite", "cjson", INPUT_FILE]
    )
    def test_ccwrite_call(self, mock_ccwrite):
        """is ccwrite called with the given parameters?"""
        self.main()

        self.assertEqual(mock_ccwrite.call_count, 1)
        ccwrite_call_args, ccwrite_call_kwargs = mock_ccwrite.call_args
        self.assertEqual(ccwrite_call_args[1], 'cjson')
        self.assertEqual(ccwrite_call_args[2], CJSON_OUTPUT_FILENAME)


@mock.patch("cclib.scripts.ccframe.ccframe")
class ccframeTest(unittest.TestCase):

    def setUp(self):
        try:
            from cclib.scripts import ccframe
        except ImportError:
            self.fail("ccframe cannot be imported")

        self.main = ccframe.main

    @mock.patch('cclib.scripts.ccframe.sys.argv', ['ccframe'])
    def test_empty_argv(self, mock_ccframe):
        """Does the script fail as expected if called without parameters?"""
        with self.assertRaises(SystemExit):
            self.main()

    @mock.patch(
        "cclib.scripts.ccframe.sys.argv",
        ["ccframe", INPUT_FILE]
    )
    @mock.patch("cclib.scripts.ccframe._has_pandas", True)
    def test_ccframe_call(self, mock_ccframe):
        """is ccframe called with the given parameters?"""
        self.main()

        self.assertEqual(mock_ccframe.call_count, 1)
        ccframe_call_args, ccframe_call_kwargs = mock_ccframe.call_args
        self.assertEqual(ccframe_call_args[0][0].filename, INPUT_FILE)

    @mock.patch(
        "cclib.scripts.ccframe.sys.argv",
        ["ccframe", INPUT_FILE]
    )
    @mock.patch("cclib.scripts.ccframe._has_pandas", False)
    def test_ccframe_call_without_pandas(self, mock_ccframe):
        """does ccframe fails cleanly if pandas can't be imported?"""
        with self.assertRaises(SystemExit):
            self.main()


if __name__ == "__main__":
    unittest.main()
