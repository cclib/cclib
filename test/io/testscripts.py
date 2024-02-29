# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Unit tests for main scripts (ccget, ccwrite)."""

import os
from test.conftest import get_program_dir, gettestdata
from test.io.testccio import BASE_URL, URL_FILES
from unittest import mock

import cclib
from cclib.io import ccread, ccwrite

import pytest

__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..", "data")


INPUT_FILE = os.path.join(__datadir__, "Gaussian/basicGaussian16/dvb_sp.out")
CJSON_OUTPUT_FILENAME = "dvb_sp.cjson"

TEST_FILES = gettestdata()

@mock.patch("cclib.scripts.ccget.ccread", wraps=ccread)
class ccgetTest:
    def setup_method(self) -> None:
        try:
            from cclib.scripts import ccget
        except ImportError:
            self.fail("ccget cannot be imported")

        self.main = ccget.ccget

    @mock.patch("cclib.scripts.ccget.sys.argv", ["ccget"])
    def test_empty_argv(self, mock_ccread) -> None:
        """Does the script fail as expected if called without parameters?"""
        with pytest.raises(SystemExit):
            self.main()

    @mock.patch("cclib.scripts.ccget.sys.argv", ["ccget", "atomcoords", INPUT_FILE])
    def test_ccread_invocation(self, mock_ccread) -> None:
        self.main()

        assert mock_ccread.call_count == 1
        ccread_call_args, ccread_call_kwargs = mock_ccread.call_args
        assert ccread_call_args[0] == INPUT_FILE

    @mock.patch("logging.Logger.warning")
    @mock.patch("cclib.scripts.ccget.sys.argv", ["ccget", "atomcoord", INPUT_FILE])
    def test_ccread_invocation_matching_args(self, mock_warn, mock_ccread):
        self.main()
        assert mock_warn.call_count == 1
        warn_call_args, warn_call_kwargs = mock_warn.call_args
        warn_message = warn_call_args[0]
        assert (
            warn_message
            == "Attribute 'atomcoord' not found, but attribute 'atomcoords' is close. Using 'atomcoords' instead."
        )
        assert mock_ccread.call_count == 1
        ccread_call_args, ccread_call_kwargs = mock_ccread.call_args
        assert ccread_call_args[0] == INPUT_FILE

    @mock.patch("cclib.scripts.ccget.sys.argv", ["ccget", "metadata", BASE_URL + URL_FILES[0]])
    def test_ccread_url(self, mock_ccread) -> None:
        self.main()

    @mock.patch("cclib.scripts.ccget.sys.argv", ["ccget", "metadata", "http://fo.bar"])
    def test_ccread_bad_url(self, mock_ccread) -> None:
        with pytest.raises(Exception):
            self.main()

    @pytest.mark.parametrize("file_path", TEST_FILES, ids=[
        "{}/{}/{}".format(file_path["parser"], file_path["subdir"], ",".join(file_path["files"]))
        for file_path in TEST_FILES
    ])
    def test_all(self, mock_ccread, file_path):
        if file_path["parser"] == "Psi3":
            pytest.skip("Psi3 is not yet supported")
        # Build a list of files.
        input_files = [
            os.path.join(
                __datadir__, get_program_dir(file_path["parser"]), file_path["subdir"], file_name
            )
            for file_name in file_path["files"]
        ]

        sig = ["ccget"]
        if len(input_files) > 1:
            sig.append("--multi")
        sig.extend(input_files)
        sig.append("metadata")

        with mock.patch("cclib.scripts.ccget.sys.argv", sig):
            self.main()


@mock.patch("cclib.scripts.ccwrite.ccwrite", wraps=ccwrite)
class ccwriteTest:
    def setup_method(self) -> None:
        try:
            from cclib.scripts import ccwrite
        except ImportError:
            self.fail("ccwrite cannot be imported")

        self.main = ccwrite.main

    @mock.patch("cclib.scripts.ccwrite.sys.argv", ["ccwrite"])
    def test_empty_argv(self, mock_ccwrite) -> None:
        """Does the script fail as expected if called without parameters?"""
        with pytest.raises(SystemExit):
            self.main()

    @mock.patch("cclib.scripts.ccwrite.sys.argv", ["ccwrite", "cjson", INPUT_FILE])
    def test_ccwrite_call(self, mock_ccwrite, tmp_path, monkeypatch) -> None:
        """is ccwrite called with the given parameters?"""
        monkeypatch.chdir(tmp_path)

        self.main()

        assert mock_ccwrite.call_count == 1
        ccwrite_call_args, ccwrite_call_kwargs = mock_ccwrite.call_args
        assert ccwrite_call_args[1] == "cjson"
        assert ccwrite_call_args[2] == CJSON_OUTPUT_FILENAME


class ccframeTest:
    def setup_method(self) -> None:
        # It would be best to test with Pandas and not a mock!
        if not hasattr(cclib.io.ccio, "pd"):
            cclib.io.ccio.pd = mock.MagicMock()

    @mock.patch("cclib.scripts.ccframe.sys.argv", ["ccframe"])
    def test_main_empty_argv(self) -> None:
        """Does main() fail as expected if called without arguments?"""
        with pytest.raises(SystemExit):
            cclib.scripts.ccframe.main()

    @mock.patch("cclib.scripts.ccframe.sys.argv", ["ccframe", INPUT_FILE])
    @mock.patch("cclib.io.ccio._has_pandas", False)
    def test_main_without_pandas(self) -> None:
        """Does ccframe fail if Pandas can't be imported?"""
        with pytest.raises(ImportError, match="You must install `pandas` to use this function"):
            cclib.scripts.ccframe.main()

    @mock.patch("cclib.scripts.ccframe.sys.argv", ["ccframe", INPUT_FILE])
    @mock.patch("cclib.io.ccio._has_pandas", True)
    def test_main(self) -> None:
        """Is ccframe called with the given parameters?"""
        with mock.patch("sys.stdout") as mock_stdout:
            cclib.scripts.ccframe.main()
            assert mock_stdout.write.call_count == 2
            df, newline = mock_stdout.write.call_args_list
            if isinstance(df[0][0], mock.MagicMock):
                assert df[0][0].name == "mock.DataFrame()"
            else:
                # TODO: this is what we really should be testing
                pass
            assert newline[0][0] == "\n"
