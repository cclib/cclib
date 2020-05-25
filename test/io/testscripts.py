# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Unit tests for main scripts (ccget, ccwrite)."""

from __future__ import print_function

import os
import sys
import unittest
import six
from six import add_move, MovedModule
add_move(MovedModule('mock', 'mock', 'unittest.mock'))
from six.moves import mock
import cclib


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

    @mock.patch("logging.warning")
    @mock.patch(
        "cclib.scripts.ccget.sys.argv",
        ["ccget", "atomcoord", INPUT_FILE]
    )
    def test_ccread_invocation_matching_args(self, mock_warn, mock_ccread):
        self.main()
        self.assertEqual(mock_warn.call_count, 1)
        warn_call_args, warn_call_kwargs = mock_warn.call_args
        warn_message = warn_call_args[0]
        self.assertEqual(warn_message, "Attribute 'atomcoord' not found, but attribute 'atomcoords' is close. Using 'atomcoords' instead.")
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


class ccframeTest(unittest.TestCase):

    def setUp(self):
        # It would be best to test with Pandas and not a mock!
        if not hasattr(cclib.io.ccio, "pd"):
            cclib.io.ccio.pd = mock.MagicMock()

    def test_main_empty_argv(self):
        """Does main() fail as expected if called without arguments?"""
        with self.assertRaises(SystemExit):
            cclib.scripts.ccframe.main()

    @mock.patch(
        "cclib.scripts.ccframe.sys.argv",
        ["ccframe", INPUT_FILE]
    )
    @mock.patch("cclib.io.ccio._has_pandas", False)
    def test_main_without_pandas(self):
        """Does ccframe fail if Pandas can't be imported?"""
        with six.assertRaisesRegex(
            self, ImportError, "You must install `pandas` to use this function"
        ):
            cclib.scripts.ccframe.main()

    @mock.patch(
        "cclib.scripts.ccframe.sys.argv",
        ["ccframe", INPUT_FILE]
    )
    @mock.patch("cclib.io.ccio._has_pandas", True)
    def test_main(self):
        """Is ccframe called with the given parameters?"""
        with mock.patch('sys.stdout') as mock_stdout:
            cclib.scripts.ccframe.main()
            self.assertEqual(mock_stdout.write.call_count, 2)
            df, newline = mock_stdout.write.call_args_list
            if isinstance(df[0][0], mock.MagicMock):
                self.assertEqual(df[0][0].name, 'mock.DataFrame()')
            else:
                # TODO: this is what we really should be testing
                pass
            self.assertEqual(newline[0][0], '\n')


if __name__ == "__main__":
    unittest.main()
