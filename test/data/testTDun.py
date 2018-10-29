# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test single point, unrestricted time-dependent logfiles in cclib"""

import os
import unittest

import numpy

from skip import skipForParser

__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericTDunTest(unittest.TestCase):
    """Generic time-dependent unrestricted HF/DFT unittest"""

    number = 24

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testenergiesnumber(self):
        """Is the length of etenergies correct?"""
        self.assertEqual(len(self.data.etenergies), self.number)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')    
    def testoscsnumber(self):
        """Is the length of eotscs correct?"""
        self.assertEqual(len(self.data.etoscs), self.number)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testrotatsnumber(self):
        """Is the length of etrotats correct?"""
        self.assertEqual(len(self.data.etrotats), self.number)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testsecsnumber(self):
        """Is the length of etsecs correct?"""
        self.assertEqual(len(self.data.etsecs), self.number)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testsymsnumber(self):
        """Is the length of etsyms correct?"""
        self.assertEqual(len(self.data.etsyms), self.number)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testsyms(self):
        """Is etsyms populated by singlets and triplets 50/50?"""
        singlets = [sym for sym in self.data.etsyms if "Singlet" in sym]
        triplets = [sym for sym in self.data.etsyms if "Triplet" in sym]
        self.assertEqual(len(singlets), self.number/2)
        self.assertEqual(len(triplets), self.number/2)


if __name__=="__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['TDun'])
    suite.testall()
