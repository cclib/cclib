# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for parser ccio module."""

import os
import tempfile
import unittest

import cclib

from six import StringIO


__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class guess_fileypeTest(unittest.TestCase):

    def setUp(self):
        self.guess = cclib.io.ccio.guess_filetype

    def test_fail(self):
        """Does the function fail as expected?"""
        self.assertIsNone(self.guess([]))
        self.assertIsNone(self.guess(None))
        self.assertIsNone(self.guess(os.devnull))
        self.assertIsNone(self.guess(['test', 'random', 'quantum chemistry']))

    def test_programs(self):
        """Does the function catch programs as expected?"""
        self.assertEqual(self.guess(["Amsterdam Density Functional"]), cclib.parser.ADF)
        self.assertEqual(self.guess(['Dalton - An Electronic Structure Program']), cclib.parser.DALTON)
        self.assertEqual(self.guess(['GAMESS']), cclib.parser.GAMESS)
        self.assertEqual(self.guess(['G A M E S S - U K']), cclib.parser.GAMESSUK)
        self.assertEqual(self.guess(['Gaussian, Inc.']), cclib.parser.Gaussian)
        self.assertEqual(self.guess(['Jaguar']), cclib.parser.Jaguar)
        self.assertEqual(self.guess(['PROGRAM SYSTEM MOLPRO']), cclib.parser.Molpro)
        self.assertEqual(self.guess(['MOPAC2016']), cclib.parser.MOPAC)
        self.assertEqual(self.guess(['Northwest Computational Chemistry Package']), cclib.parser.NWChem)
        self.assertEqual(self.guess(['O   R   C   A']), cclib.parser.ORCA)
        self.assertEqual(self.guess(["PSI3: An Open-Source Ab Initio Electronic Structure Package"]), cclib.parser.Psi3)
        self.assertEqual(self.guess(["Psi4: An Open-Source Ab Initio Electronic Structure Package"]), cclib.parser.Psi4)
        self.assertEqual(self.guess(['A Quantum Leap Into The Future Of Chemistry']), cclib.parser.QChem)


class ccreadTest(unittest.TestCase):

    def setUp(self):
        self.ccread = cclib.io.ccio.ccread

    def test_fail(self):
        """Does the function fail as expected?"""
        self.assertIsNone(self.ccread("", quiet=True))
        self.assertIsNone(self.ccread([], quiet=True))
        self.assertIsNone(self.ccread(None, quiet=True))


class ccopenTest(unittest.TestCase):

    def setUp(self):
        self.ccopen = cclib.io.ccio.ccopen

    def test_ccopen_fail(self):
        """Does the function fail as expected?"""
        self.assertIsNone(self.ccopen("", quiet=True))
        self.assertIsNone(self.ccopen([], quiet=True))
        self.assertIsNone(self.ccopen(None, quiet=True))

    def test_list_of_filenames(self):
        """Can we ccopen a list of filenames (https://github.com/cclib/cclib/issues/350)?"""
        absdir = os.path.dirname(os.path.abspath(__file__))
        rootdir = os.path.join(os.sep, *absdir.split(os.sep)[:-2])
        molprodir = os.path.join(rootdir, "data", "Molpro", "basicMolpro2012")
        filenames = ["dvb_gopt.out", "dvb_gopt.log"]
        filepaths = [os.path.join(molprodir, fn) for fn in filenames]
        self.assertIsNotNone(self.ccopen(filepaths))
        self.assertIsNotNone(self.ccopen(filepaths).parse())

    def test_cjson_empty_tempfile(self):
        """Do we get a CJSON object when the keyword argument used?"""
        with tempfile.NamedTemporaryFile() as tf:
            self.assertIsInstance(self.ccopen(tf.name, cjson=True), cclib.io.cjsonreader.CJSON)

    def test_url_io(self):
        """Does the function works with URLs such good as with filenames?"""
        fpath = os.path.join(__datadir__, "data")
        base_url = "https://raw.githubusercontent.com/cclib/cclib/master/data/"
        filenames = ["QChem/basicQChem5.1/dvb_td.out", "Molpro/basicMolpro2012/h2o_mp2.out"]
        for fname in filenames:
            self.assertEqual(self.ccopen(os.path.join(fpath, fname), quiet=True).parse().getattributes(tolists=True),
                             self.ccopen(base_url + fname, quiet=True).parse().getattributes(tolists=True))

    def test_multi_url_io(self):
        """Does the function works with multiple URLs such good as with multiple filenames?"""
        fpath = os.path.join(__datadir__, "data")
        base_url = "https://raw.githubusercontent.com/cclib/cclib/master/data/"
        filenames = ["Molpro/basicMolpro2012/dvb_gopt.out", "Molpro/basicMolpro2012/dvb_gopt.log"]
        self.assertEqual(
            self.ccopen([os.path.join(fpath, fname) for fname in filenames], quiet=True)
                .parse().getattributes(tolists=True),
            self.ccopen([base_url + fname for fname in filenames], quiet=True)
                .parse().getattributes(tolists=True))

    @unittest.skip("This should also work if cjsonreader supported streams.")
    def test_cjson(self):
        """Do we get a CJSON object then keyword argument used?"""
        self.assertIsInstance(self.ccopen(StringIO(""), cjson=True), cclib.io.cjsonreader.CJSON)


class _determine_output_formatTest(unittest.TestCase):

    def setUp(self):
        self._determine_output_format = cclib.io.ccio._determine_output_format
        self.UnknownOutputFormatError = cclib.io.ccio.UnknownOutputFormatError

    def test_outputclass(self):
        """Does the function determine output class as expected."""
        outputtype = "xyz"
        outputdest = "file.xyz"
        self.assertEqual(self._determine_output_format(outputtype, outputdest),
                         cclib.io.xyzwriter.XYZ)
        # Must raise a KeyError for unsuported extensions
        self.assertRaises(self.UnknownOutputFormatError,
                          self._determine_output_format, 'ext', outputdest)
        self.assertRaises(self.UnknownOutputFormatError,
                          self._determine_output_format, None, None)


class fallbackTest(unittest.TestCase):

    def setUp(self):
        self.fallback = cclib.io.ccio.fallback

    def test_fallback_fail(self):
        """Does the function fail as expected?"""
        self.assertIsNone(self.fallback(None))


if __name__ == "__main__":
    unittest.main()
