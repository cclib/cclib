# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2007,2012,2014,2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test Moller-Plesset logfiles in cclib"""

import os
import unittest

import numpy


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericMP2Test(unittest.TestCase):
    """Generic MP2 unittest"""

    level = 2
    
    def testsizeandshape(self):
        """(MP2) Are the dimensions of mpenergies correct?"""
        self.assertEqual(self.data.mpenergies.shape,
                         (len(self.data.scfenergies), self.level-1))

    def testchange(self):
        """(MP2) Are Moller-Plesset corrections negative?"""
        if self.level == 2:
            corrections = self.data.mpenergies[:,0] - self.data.scfenergies
        else:
            corrections = self.data.mpenergies[:,self.level-2] - self.data.mpenergies[:,self.level-3]
        self.failUnless(numpy.alltrue(corrections < 0.0))
        
class GenericMP3Test(GenericMP2Test):
    """Generic MP3 unittest"""
    level = 3

class GenericMP4SDQTest(GenericMP2Test):
    """Generic MP4(SDQ) unittest"""
    level = 4

class GenericMP4SDTQTest(GenericMP2Test):
    """Generic MP4(SDTQ) unittest"""
    level = 4

class GenericMP5Test(GenericMP2Test):
    """Generic MP5 unittest"""
    level = 5


class GaussianMP2Test(GenericMP2Test):
    """Customized MP2 unittest"""
        
    def testnocoeffs(self):
        """Are natural orbital coefficients the right size?"""
        self.assertEquals(self.data.nocoeffs.shape, (self.data.nmo, self.data.nbasis))

    def testnocoeffs(self):
        """Are natural orbital occupation numbers the right size?"""
        self.assertEquals(self.data.nooccnos.shape, (self.data.nmo, ))

class GaussianMP3Test(GenericMP2Test):
    """Customized MP3 unittest"""
    level = 3

class GaussianMP4SDQTest(GenericMP2Test):
    """Customized MP4-SDQ unittest"""
    level = 4

class GaussianMP4SDTQTest(GenericMP2Test):
    """Customized MP4-SDTQ unittest"""
    level = 4


class QChemMP4SDQTest(GenericMP2Test):
    """Customized MP4-SDQ unittest"""
    level = 4

class QChemMP4SDTQTest(GenericMP2Test):
    """Customized MP4-SD(T)Q unittest"""
    level = 5


if __name__=="__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['MP'])
    suite.testall()
