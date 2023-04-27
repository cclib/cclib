# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for parser ccio module."""

import os
import sys
import tempfile
from io import StringIO
import unittest
from unittest import mock

import cclib
import pytest


__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class guess_fileypeTest(unittest.TestCase):

    def setUp(self):
        self.guess = cclib.io.ccio.guess_filetype

    def test_fail(self):
        """Does the function fail as expected?"""
        assert self.guess([]) is None
        assert self.guess(None) is None
        assert self.guess(os.devnull) is None
        assert self.guess(['test', 'random', 'quantum chemistry']) is None

    def test_programs(self):
        """Does the function catch programs as expected?"""
        assert self.guess(["Amsterdam Density Functional"]) == cclib.parser.ADF
        assert self.guess(['Dalton - An Electronic Structure Program']) == cclib.parser.DALTON
        assert self.guess(['GAMESS']) == cclib.parser.GAMESS
        assert self.guess(['G A M E S S - U K']) == cclib.parser.GAMESSUK
        assert self.guess(['Gaussian, Inc.']) == cclib.parser.Gaussian
        assert self.guess(['Jaguar']) == cclib.parser.Jaguar
        assert self.guess(['PROGRAM SYSTEM MOLPRO']) == cclib.parser.Molpro
        assert self.guess(['MOPAC2016']) == cclib.parser.MOPAC
        assert self.guess(['Northwest Computational Chemistry Package']) == cclib.parser.NWChem
        assert self.guess(['O   R   C   A']) == cclib.parser.ORCA
        assert self.guess(["PSI3: An Open-Source Ab Initio Electronic Structure Package"]) == cclib.parser.Psi3
        assert self.guess(["Psi4: An Open-Source Ab Initio Electronic Structure Package"]) == cclib.parser.Psi4
        assert self.guess(['A Quantum Leap Into The Future Of Chemistry']) == cclib.parser.QChem


class ccreadTest(unittest.TestCase):

    def setUp(self):
        self.ccread = cclib.io.ccio.ccread

    def test_fail(self):
        """Does the function fail as expected?"""
        assert self.ccread("", quiet=True) is None
        assert self.ccread([], quiet=True) is None
        assert self.ccread(None, quiet=True) is None


class ccopenTest(unittest.TestCase):

    def setUp(self):
        self.ccopen = cclib.io.ccio.ccopen

    def test_ccopen_fail(self):
        """Does the function fail as expected?"""
        assert self.ccopen("", quiet=True) is None
        assert self.ccopen([], quiet=True) is None
        assert self.ccopen(None, quiet=True) is None

    def test_list_of_filenames(self):
        """Can we ccopen a list of filenames (https://github.com/cclib/cclib/issues/350)?"""
        absdir = os.path.dirname(os.path.abspath(__file__))
        rootdir = os.path.join(os.sep, *absdir.split(os.sep)[:-2])
        molprodir = os.path.join(rootdir, "data", "Molpro", "basicMolpro2012")
        filenames = ["dvb_gopt.out", "dvb_gopt.log"]
        filepaths = [os.path.join(molprodir, fn) for fn in filenames]
        assert self.ccopen(filepaths) is not None
        assert self.ccopen(filepaths).parse() is not None

    def test_cjson_empty_tempfile(self):
        """Do we get a CJSON object when the keyword argument used?"""
        with tempfile.NamedTemporaryFile() as tf:
            assert isinstance(self.ccopen(tf.name, cjson=True), cclib.io.cjsonreader.CJSON)

    def test_url_io(self):
        """Does the function works with URLs such good as with filenames?"""
        fpath = os.path.join(__datadir__, "data")
        base_url = "https://raw.githubusercontent.com/cclib/cclib/master/data/"
        filenames = ["QChem/basicQChem5.1/dvb_td.out", "Molpro/basicMolpro2012/h2o_mp2.out"]
        for fname in filenames:
            assert self.ccopen(os.path.join(fpath, fname), quiet=True).parse().getattributes(tolists=True) == \
                             self.ccopen(base_url + fname, quiet=True).parse().getattributes(tolists=True)

    def test_multi_url_io(self):
        """Does the function works with multiple URLs such good as with multiple filenames?"""
        fpath = os.path.join(__datadir__, "data")
        base_url = "https://raw.githubusercontent.com/cclib/cclib/master/data/"
        filenames = ["Molpro/basicMolpro2012/dvb_gopt.out", "Molpro/basicMolpro2012/dvb_gopt.log"]
        assert self.ccopen([os.path.join(fpath, fname) for fname in filenames], quiet=True) \
                .parse().getattributes(tolists=True) == \
            self.ccopen([base_url + fname for fname in filenames], quiet=True) \
                .parse().getattributes(tolists=True)

    @unittest.skip("This should also work if cjsonreader supported streams.")
    def test_cjson(self):
        """Do we get a CJSON object then keyword argument used?"""
        assert isinstance(self.ccopen(StringIO(""), cjson=True), cclib.io.cjsonreader.CJSON)
        
    def test_bz2_io(self):
        """Can we read from a bz2 archive?"""
        file_path = os.path.join(__filedir__, "data/dvb_gopt.out.bz2")
        # Test both single-file and multi-file parsing.
        assert self.ccopen(file_path) is not None
        assert self.ccopen([file_path, file_path]) is not None
    
    def test_gz_io(self):
        """Can we read from a gz archive?"""
        file_path = os.path.join(__filedir__, "data/dvb_gopt.out.gz")
        # Test both single-file and multi-file parsing.
        assert self.ccopen(file_path) is not None
        assert self.ccopen([file_path, file_path]) is not None
        
    def test_zip_io(self):
        """Can we read from a zip archive?"""
        file_path = os.path.join(__filedir__, "data/dvb_gopt.out.zip")
        # Test both single-file and multi-file parsing.
        assert self.ccopen(file_path) is not None
        assert self.ccopen([file_path, file_path]) is not None


class _determine_output_formatTest(unittest.TestCase):

    def setUp(self):
        self._determine_output_format = cclib.io.ccio._determine_output_format
        self.UnknownOutputFormatError = cclib.io.ccio.UnknownOutputFormatError

    def test_outputclass(self):
        """Does the function determine output class as expected."""
        outputtype = "xyz"
        outputdest = "file.xyz"
        assert self._determine_output_format(outputtype, outputdest) == \
                         cclib.io.xyzwriter.XYZ
        # Must raise a KeyError for unsuported extensions
        with pytest.raises(self.UnknownOutputFormatError):
            self._determine_output_format('ext', outputdest)
        with pytest.raises(self.UnknownOutputFormatError):
            self._determine_output_format(None, None)


class fallbackTest(unittest.TestCase):

    def setUp(self):
        self.fallback = cclib.io.ccio.fallback

    def test_fallback_fail(self):
        """Does the function fail as expected?"""
        assert self.fallback(None) is None


class ccframeTest(unittest.TestCase):

    @mock.patch("cclib.io.ccio._has_pandas", False)
    def test_ccframe_call_without_pandas(self):
        """Does ccframe fails cleanly if Pandas can't be imported?"""
        with pytest.raises(ImportError, match="You must install `pandas` to use this function"):
            cclib.io.ccio.ccframe([])


if __name__ == "__main__":
    unittest.main()
