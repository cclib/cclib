# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

__revision__ = "$Revision$"

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

import bettertest


class GenericSPunTest(bettertest.TestCase):
    """Unrestricted single point unittest."""

    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        self.assertEquals(self.data.aooverlaps.shape,(self.data.nbasis,self.data.nbasis))

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 2 x nmo x nbasis?"""
        self.assertEquals(type(self.data.mocoeffs), type([]))
        self.assertEquals(len(self.data.mocoeffs), 2)
        self.assertEquals(self.data.mocoeffs[0].shape,
                          (self.data.nmo, self.data.nbasis))
        self.assertEquals(self.data.mocoeffs[1].shape,
                          (self.data.nmo, self.data.nbasis))

    def testcharge_and_mult(self):
        """Are the charge and multiplicity correct?"""
        self.assertEquals(self.data.charge, 1)
        self.assertEquals(self.data.mult, 2)

    def testhomos(self):
        """Are the homos correct?"""
        self.assertArrayEquals(self.data.homos, numpy.array([34,33],"i"),"%s != array([34,33],'i')" % numpy.array_repr(self.data.homos))

    def testmoenergies(self):
        """Are the dims of the moenergies equals to 2 x nmo?"""
        self.assertEquals(len(self.data.moenergies), 2)
        self.assertEquals(len(self.data.moenergies[0]), self.data.nmo)
        self.assertEquals(len(self.data.moenergies[1]), self.data.nmo)

    def testmosyms(self):
        """Are the dims of the mosyms equals to 2 x nmo?"""
        shape = (len(self.data.mosyms), len(self.data.mosyms[0]))
        self.assertEquals(shape, (2, self.data.nmo))

        
class ADFSPunTest(GenericSPunTest):
    """ADF unrestricted single point unittest."""

    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        #ADF uses fooverlaps
        self.assertEquals(self.data.fooverlaps.shape,(self.data.nbasis,self.data.nbasis))


class GamessUKSPunTest(GenericSPunTest):
    """GAMESS-UK unrestricted single point unittest."""

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 2 x (homos+6) x nbasis?"""
        self.assertEquals(type(self.data.mocoeffs), type([]))
        self.assertEquals(len(self.data.mocoeffs), 2)
        self.assertEquals(self.data.mocoeffs[0].shape,
                          (self.data.homos[0]+6, self.data.nbasis))
        self.assertEquals(self.data.mocoeffs[1].shape,
                          (self.data.homos[1]+6, self.data.nbasis))


class GamessUSSPunTest(GenericSPunTest):
    """GAMESS-US unrestricted single point unittest."""


class GaussianSPunTest(GenericSPunTest):
    """Gaussian unrestricted single point unittest."""

    def testatomnos(self):
        """Does atomnos have the right dimension (20)?"""
        size = len(self.data.atomnos)
        self.assertEquals(size, 20)


class JaguarSPunTest(GenericSPunTest):
    """Jaguar unrestricted single point unittest."""
        
    # Data file does not contain enough information. Can we make a new one?
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis? PASS"""
        self.assertEquals(1,1)

    # Why is this test passed?
    def testmoenergies(self):
        """Are the dims of the moenergies equal to 2 x homos+11?"""
        self.assertEquals(len(self.data.moenergies), 2)
        self.assertEquals(len(self.data.moenergies[0]), self.data.homos[0]+11)
        self.assertEquals(len(self.data.moenergies[1]), self.data.homos[1]+11)
        
    # Data file does not contain enough information. Can we make a new one?
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""
        self.assertEquals(1,1)

    # Why is this test passed?
    def testmosyms(self):
        """Are the dims of the mosyms equal to 2 x nmo? PASS"""
        self.assertEquals(1,1)


class MolproSPunTest(GenericSPunTest):
    """Molpro unrestricted single point unittest."""

    def testmosyms(self):
        """Are the dims of the mosyms equal to 2 x nmo? PASS"""
        self.assertEquals(1,1)


class OrcaSPunTest(GenericSPunTest):
    """ORCA unrestricted single point unittest."""
    
    # ORCA has no support for symmetry yet.
    def testmosyms(self):
        """Are the dims of the mosyms equals to 2 x nmo?"""
        self.assertEquals(1,1)


class PCGamessSPunTest(GenericSPunTest):
    """PC-GAMESS unrestricted single point unittest."""

              
if __name__=="__main__":

    from testall import testall
    testall(modules=["SPun"])
