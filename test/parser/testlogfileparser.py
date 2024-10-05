# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for the logfileparser module."""

import io
import os
import sys
import tempfile
from unittest import mock
from urllib.request import urlopen

import cclib

import pytest
from pytest import MonkeyPatch

__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class FileWrapperTest:
    def check_seek(self, wrapper):
        """Check that a FileWrapper can seek properly"""
        wrapper.seek(0, 2)
        assert wrapper.pos == wrapper.size

        wrapper.seek(0, 0)
        assert wrapper.pos == 0

        with pytest.raises(NotImplementedError):
            wrapper.seek(0, 1)

    def test_file_seek(self):
        """Can we seek anywhere in a file object?"""
        fpath = os.path.join(__datadir__, "data/ADF/basicADF2007.01/dvb_gopt.adfout")
        with open(fpath) as fobject:
            wrapper = cclib.parser.logfileparser.FileWrapper(fobject)
            self.check_seek(wrapper)

    def test_url_seek(self):
        """Can we seek only to the end of an url stream?"""

        url = "https://raw.githubusercontent.com/cclib/cclib/master/data/ADF/basicADF2007.01/dvb_gopt.adfout"
        stream = urlopen(url)
        wrapper = cclib.parser.logfileparser.FileWrapper(stream)

        self.check_seek(wrapper)

    def test_stdin_seek(self):
        """We shouldn't be able to seek anywhere in standard input."""
        # stdin is disabled by pytest.
        # the recommended way of emulating stdin is by doing this
        monkeypatch = MonkeyPatch()
        monkeypatch.setattr("sys.stdin", io.StringIO())

        wrapper = cclib.parser.logfileparser.FileWrapper(sys.stdin)
        self.check_seek(wrapper)

    def test_data_stdin(self):
        """Check that the same attributes are parsed when a file is piped through standard input."""
        logfiles = [
            "data/ADF/basicADF2007.01/dvb_gopt.adfout",
            "data/GAMESS/basicGAMESS-US2017/C_bigbasis.out",
        ]

        def get_attributes(data):
            return [a for a in data._attrlist if hasattr(data, a)]

        for lf in logfiles:
            path = f"{__datadir__}/{lf}"
            data = cclib.io.ccread(path)
            assert data is not None
            expected_attributes = get_attributes(data)
            with open(path) as handle:
                contents = handle.read()

            # stdin emulation
            monkeypatch = MonkeyPatch()
            monkeypatch.setattr("sys.stdin", io.StringIO(contents))

            data = cclib.io.ccread(sys.stdin)
            assert data is not None
            assert get_attributes(data) == expected_attributes


class LogfileTest:
    """Unit tests for the Logfile class."""

    def test_parse_check_values(self):
        """Are custom checks performed after parsing finishes?

        The purpose of this test is not to comprehensively cover all the checks,
        but rather to make sure the call and logging works. The unit tests
        for the data class should have comprehensive coverage.
        """
        _, path = tempfile.mkstemp()
        logfileclass = cclib.parser.logfileparser.Logfile
        logfileclass.__abstractmethods__ = set()
        parser = logfileclass(path)
        parser.extract = lambda self, inputfile, line: None
        parser.logger = mock.Mock()

        parser.etenergies = [1, -1]
        parser.parse()
        parser.logger.error.assert_called_once()
