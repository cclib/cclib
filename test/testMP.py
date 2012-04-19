# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2007, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

__revision__ = "$Revision$"

import numpy

import bettertest


class GenericMP2Test(bettertest.TestCase):
    """Generic MP2 unittest."""

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
    """Generic MP3 unittest."""

    level = 3
    
    def testsizeandshape(self):
        """(MP3) Are the dimensions of mpenergies correct?"""
        super(GenericMP3Test,self).testsizeandshape()
        
    def testchange(self):
        """(MP3) Are Moller-Plesset corrections negative?"""
        super(GenericMP3Test,self).testchange()

class GenericMP4SDQTest(GenericMP2Test):
    """Generic MP4(SDQ) unittest."""

    level = 4
    
    def testsizeandshape(self):
        """(MP4-SDQ) Are the dimensions of mpenergies correct?"""
        super(GenericMP4SDQTest,self).testsizeandshape()
        
    def testchange(self):
        """(MP4-SDQ) Are Moller-Plesset corrections negative?"""
        super(GenericMP4SDQTest,self).testchange()

class GenericMP4SDTQTest(GenericMP2Test):
    """Generic MP4(SDTQ) unittest."""
    
    level = 4
    
    def testsizeandshape(self):
        """(MP4-SDTQ) Are the dimensions of mpenergies correct?"""
        super(GenericMP4SDTQTest,self).testsizeandshape()
        
    def testchange(self):
        """(MP4-SDTQ) Are Moller-Plesset corrections negative?"""
        super(GenericMP4SDTQTest,self).testchange()

class GenericMP5Test(GenericMP2Test):
    """Generic MP5 unittest."""

    level = 5
    
    def testsizeandshape(self):
        """(MP5) Are the dimensions of mpenergies correct?"""
        super(GenericMP5Test,self).testsizeandshape()

    def testchange(self):
        """(MP5) Are Moller-Plesset corrections negative?"""
        super(GenericMP5Test,self).testchange()

class GAMESSUKMP2Test(GenericMP2Test):
    """GAMESS-UK MP2 unittest."""

class GAMESSUKMP3Test(GenericMP3Test):
    """GAMESS-UK MP3 unittest."""

class GAMESSUSMP2Test(GenericMP2Test):
    """GAMESS-US MP2 unittest."""

class GaussianMP2Test(GenericMP2Test):
    """Gaussian MP2 unittest."""
        
    def testnocoeffs(self):
        """(MP2) Are Natural Orbital coefficients the right size?"""
        self.assertEquals(self.data.nocoeffs.shape, (self.data.nmo, self.data.nbasis))
        
class GaussianMP3Test(GenericMP3Test):
    """Gaussian MP3 unittest."""
        
    def testnocoeffs(self):
        """(MP2) Are Natural Orbital coefficients the right size?"""
        self.assertEquals(self.data.nocoeffs.shape, (self.data.nmo, self.data.nbasis))

class GaussianMP4SDQTest(GenericMP4SDQTest):
    """Gaussian MP4-SDQ unittest."""
        
    def testnocoeffs(self):
        """(MP2) Are Natural Orbital coefficients the right size?"""
        self.assertEquals(self.data.nocoeffs.shape, (self.data.nmo, self.data.nbasis))

class GaussianMP4SDTQTest(GenericMP4SDTQTest):
    """Gaussian MP4-SDTQ unittest."""

class GaussianMP5Test(GenericMP5Test):
    """Gaussian MP5 unittest."""

class JaguarLMP2Test(GenericMP2Test):
    """Jaguar LMP2 unittest."""

class MolproMP2Test(GenericMP2Test):
    """Molpro MP2 unittest."""

class MolproMP3Test(GenericMP3Test):
    """Molpro MP3 unittest."""

class MolproMP4SDTQTest(GenericMP4SDTQTest):
    """Molpro MP4-SDTQ unittest."""

class PCGAMESSMP2Test(GenericMP2Test):
    """PC-GAMESS MP2 unittest."""

class PCGAMESSMP3Test(GenericMP3Test):
    """PC-GAMESS MP3 unittest."""

class PCGAMESSMP4SDQTest(GenericMP4SDQTest):
    """PC-GAMESS MP4-SDQ unittest."""

class PCGAMESSMP4SDTQTest(GenericMP4SDTQTest):
    """PC-GAMESS MP4-SDTQ unittest."""

              
if __name__=="__main__":

    from testall import testall
    testall(modules=["MP"])
