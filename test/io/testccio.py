# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2015-2016, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Unit tests for parser ccio module."""

import os
import tempfile
import unittest

import cclib


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
        self.assertEqual(self.guess(['Northwest Computational Chemistry Package']), cclib.parser.NWChem)
        self.assertEqual(self.guess(['O   R   C   A']), cclib.parser.ORCA)
        self.assertEqual(self.guess(["PSI ...Ab Initio Electronic Structure"]), cclib.parser.Psi)
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

    def test_cjson_empty_tempfile(self):
        """Do we get a CJSON object when the keyword argument used?"""
        with tempfile.NamedTemporaryFile() as tf:
            self.assertIsInstance(self.ccopen(tf.name, cjson=True), cclib.io.cjsonreader.CJSON)

    # This should also work if cjsonreader supported streams.
    #def test_cjson(self):
    #    """Do we get a CJSON object then keyword argument used?"""
    #    self.assertIsInstance(self.ccopen(StringIO.StringIO(""), cjson=True), cclib.io.cjsonreader.CJSON)


class fallbackTest(unittest.TestCase):

    def setUp(self):
        self.fallback = cclib.io.ccio.fallback

    def test_fallback_fail(self):
        """Does the functin fail as expected?"""
        self.assertIsNone(self.fallback(None))


if __name__ == "__main__":
    unittest.main()
