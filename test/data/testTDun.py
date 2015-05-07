# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2007,2012,2014,2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test single point, unrestricted time-dependent logfiles in cclib"""

import os
import unittest

import numpy


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericTDunTest(unittest.TestCase):
    """Generic time-dependent unrestricted HF/DFT unittest"""

    number = 24

    def testenergiesnumber(self):
        """Is the length of etenergies correct?"""
        self.assertEqual(len(self.data.etenergies), self.number)
    
    def testoscsnumber(self):
        """Is the length of eotscs correct?"""
        self.assertEqual(len(self.data.etoscs), self.number)

    def testrotatsnumber(self):
        """Is the length of etrotats correct?"""
        self.assertEqual(len(self.data.etrotats), self.number)

    def testsecsnumber(self):
        """Is the length of etsecs correct?"""
        self.assertEqual(len(self.data.etsecs), self.number)

    def testsymsnumber(self):
        """Is the length of etsyms correct?"""
        self.assertEqual(len(self.data.etsyms), self.number)

    def testsyms(self):
        """Is etsyms populated by singlets and triplets 50/50?"""
        singlets = [sym for sym in self.data.etsyms if "Singlet" in sym]
        triplets = [sym for sym in self.data.etsyms if "Triplet" in sym]
        self.assertEqual(len(singlets), self.number/2)
        self.assertEqual(len(triplets), self.number/2)


if __name__=="__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['TDun'])
    suite.testall()
