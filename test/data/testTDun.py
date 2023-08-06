# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test single point, unrestricted time-dependent logfiles in cclib"""

import os
import unittest

import numpy

from skip import skipForParser
from skip import skipForLogfile

__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericTDunTest(unittest.TestCase):
    """Generic time-dependent unrestricted HF/DFT unittest"""

    number = 24

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    def testenergiesnumber(self):
        """Is the length of etenergies correct?"""
        assert len(self.data.etenergies) == self.number

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForLogfile("Turbomole/basicTurbomole7.4/CO_cc2_TD_un", "Oscillator strengths are not available for triplets with Turbomole's ricc2")
    def testoscsnumber(self):
        """Is the length of eotscs correct?"""
        assert len(self.data.etoscs) == self.number

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForLogfile("Turbomole/basicTurbomole7.4/CO_cc2_TD", "Rotatory strengths are not currently available for ricc2")
    def testrotatsnumber(self):
        """Is the length of etrotats correct?"""
        assert len(self.data.etrotats) == self.number

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    def testsecsnumber(self):
        """Is the length of etsecs correct?"""
        assert len(self.data.etsecs) == self.number

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    def testsymsnumber(self):
        """Is the length of etsyms correct?"""
        assert len(self.data.etsyms) == self.number

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole', 'Turbomole etsyms are not available for UHF')
    def testsyms(self):
        """Is etsyms populated by singlets and triplets 50/50?"""
        singlets = [sym for sym in self.data.etsyms if "Singlet" in sym]
        triplets = [sym for sym in self.data.etsyms if "Triplet" in sym]
        assert len(singlets) == self.number/2
        assert len(triplets) == self.number/2


class TurbomoleTDunTest(GenericTDunTest):
    """Customized time-dependent unrestricted HF/DFT unittest"""
    
    number = 10



if __name__=="__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['TDun'])
    suite.testall()
