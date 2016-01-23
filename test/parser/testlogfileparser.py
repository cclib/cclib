# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Unit tests for parser logfileparser module."""

import io
import os
import sys
import unittest

# The structure of urllib changed in Python3.
try:
    from urllib.request import urlopen
except ImportError:
    from urllib import urlopen

import numpy

import cclib


__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")

class FileWrapperTest(unittest.TestCase):

    def test_file_seek(self):
        """Can we seek anywhere in a file object?"""
        fpath = os.path.join(__datadir__,"data/ADF/basicADF2007.01/dvb_gopt.adfout")
        with open(fpath, 'r') as fobject:
            wrapper = cclib.parser.logfileparser.FileWrapper(fobject)
            wrapper.seek(0, 0)
            self.assertEqual(wrapper.pos, 0)
            wrapper.seek(10, 0)
            self.assertEqual(wrapper.pos, 10)
            wrapper.seek(0, 2)
            self.assertEqual(wrapper.pos, wrapper.size)

    def test_url_seek(self):
        """Can we seek only to the end of an url stream?"""

        url = "https://raw.githubusercontent.com/cclib/cclib/master/data/ADF/basicADF2007.01/dvb_gopt.adfout"
        stream = urlopen(url)
        wrapper = cclib.parser.logfileparser.FileWrapper(stream)

        # Unfortunately, the behavior of this wrapper differs between Python 2 and 3,
        # so we need to diverge the assertions. We should try to keep the code as
        # consistent as possible, but the Errors raised are actually different.
        wrapper.seek(0, 2)
        self.assertEqual(wrapper.pos, wrapper.size)
        if sys.version_info[0] == 2:
            self.assertRaises(AttributeError, wrapper.seek, 0, 0)
            self.assertRaises(AttributeError, wrapper.seek, 0, 1)
        else:
            self.assertRaises(io.UnsupportedOperation, wrapper.seek, 0, 0)
            self.assertRaises(io.UnsupportedOperation, wrapper.seek, 0, 1)


class LogfileTest(unittest.TestCase):
    """Unit tests for the Logfile class."""

    def test_float_basic(self):
        """Are floats converted from strings correctly?"""
        float = cclib.parser.logfileparser.Logfile('').float
        self.assertEqual(float("0.0"), 0.0)
        self.assertEqual(float("1.0"), 1.0)
        self.assertEqual(float("-1.0"), -1.0)

    def test_float_numeric_format(self):
        """Does numeric formatting get converted correctly?"""
        float = cclib.parser.logfileparser.Logfile('').float
        self.assertEqual(float("1.2345E+02"), 123.45)
        self.assertEqual(float("1.2345D+02"), 123.45)

    def test_float_stars(self):
        """Does the function return nan for stars?"""
        float = cclib.parser.logfileparser.Logfile('').float
        self.assertTrue(numpy.isnan(float("*")))
        self.assertTrue(numpy.isnan(float("*****")))

    def test_normalisesym_base_class_error(self):
        """Does this method return ERROR in base class?"""
        normalisesym = cclib.parser.logfileparser.Logfile('').normalisesym
        self.assertTrue("ERROR" in normalisesym(""))


if __name__ == "__main__":
    unittest.main()
