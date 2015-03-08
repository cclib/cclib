# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Custom unittest class for cclib"""

import unittest

import numpy

# from cclib.parser import ADF, GAMESS, GAMESSUK, Gaussian, Jaguar, Molpro
from cclib.parser import Gaussian


class TestCase(unittest.TestCase):
    """Create a class with extra 'asserts' for testing numerical data,
    and a special run() method for loading cclib test files on run-time.

    It is not possible to test equality of numpy arrays using assertEquals().
    Instead, use assertArrayEquals() as defined below. For the original solution see:
    http://mail.python.org/pipermail/python-list/2005-November/311235.html
    """

    def run(self, result=None):
        """Custom run method for cclib."""
        
        # Skip the actuall call to run() if we want to skip the test,
        #   for any reason (feature not implemented, etc.).
        # Still, increment the total, and append to a new "skipped" attribute.
        # In Python 2.4, the document string is in _TestCase__testMethodDoc.
        # In Python 2.5, the document string is in _testMethodDoc.
        try:
            doc = self._testMethodDoc
        except AttributeError:
            doc = self._TestCase__testMethodDoc
        if "PASS" in doc:
            result.stream.writeln(doc)
            result.testsRun += 1
            if not hasattr(result, "skipped"):
                result.skipped = []
            result.skipped.append(doc)
            return
        
        # If test not skipped, run the test now.
        unittest.TestCase.run(self, result=result)
