# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006,2007,2012,2014,2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test unrestrictied single point logfiles in cclib"""

import os
import unittest

import numpy

from skip import skipForParser
from skip import skipForLogfile


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericSPunTest(unittest.TestCase):
    """Generic unrestricted single point unittest"""

    def testnatom(self):
        """Is the number of atoms equal to 20?"""
        self.assertEquals(self.data.natom, 20)

    def testatomnos(self):
        """Are the atomnos correct?"""
        self.failUnless(numpy.alltrue([numpy.issubdtype(atomno,int) for atomno in self.data.atomnos]))
        self.assertEquals(self.data.atomnos.shape, (20,) )
        self.assertEquals(sum(self.data.atomnos==6) + sum(self.data.atomnos==1), 20)

    def testatomcoords(self):
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        self.assertEquals(self.data.atomcoords.shape,(1,self.data.natom,3))

    @skipForParser('Jaguar', 'Data file does not contain enough information')
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
        msg = "%s != array([34,33],'i')" % numpy.array_repr(self.data.homos)
        numpy.testing.assert_array_equal(self.data.homos, numpy.array([34,33],"i"), msg)

    def testmoenergies(self):
        """Are the dims of the moenergies equals to 2 x nmo?"""
        self.assertEquals(len(self.data.moenergies), 2)
        self.assertEquals(len(self.data.moenergies[0]), self.data.nmo)
        self.assertEquals(len(self.data.moenergies[1]), self.data.nmo)

    @skipForParser('Molpro', '?')
    @skipForParser('ORCA', 'ORCA has no support for symmetry yet')
    def testmosyms(self):
        """Are the dims of the mosyms equals to 2 x nmo?"""
        shape = (len(self.data.mosyms), len(self.data.mosyms[0]))
        self.assertEquals(shape, (2, self.data.nmo))


class GenericROSPTest(GenericSPunTest):
    """Customized restricted open-shell single point unittest"""

    @skipForParser('DALTON', 'mocoeffs not implemented yet')
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        self.assertEquals(type(self.data.mocoeffs), type([]))
        self.assertEquals(len(self.data.mocoeffs), 1)
        self.assertEquals(self.data.mocoeffs[0].shape,
                          (self.data.nmo, self.data.nbasis))

    def testhomos(self):
        """Are the HOMO indices equal to 34 and 33 (one more alpha electron
        than beta electron)?
        """
        msg = "%s != array([34, 33], 'i')" % numpy.array_repr(self.data.homos)
        numpy.testing.assert_array_equal(self.data.homos, numpy.array([34, 33], "i"), msg)

    @skipForParser('QChem', 'prints 2 sets of different MO energies?')
    def testmoenergies(self):
        """Are the dims of the moenergies equals to 1 x nmo?"""
        self.assertEquals(len(self.data.moenergies), 1)
        self.assertEquals(len(self.data.moenergies[0]), self.data.nmo)

    def testmosyms(self):
        """Are the dims of the mosyms equals to 1 x nmo?"""
        shape = (len(self.data.mosyms), len(self.data.mosyms[0]))
        self.assertEquals(shape, (1, self.data.nmo))


class GamessUK70SPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 2 x (homos+6) x nbasis?"""

        self.assertEquals(type(self.data.mocoeffs), type([]))
        self.assertEquals(len(self.data.mocoeffs), 2)

        # This is only an issue in version 7.0 (and before?), since in the version 8.0
        # logfile all eigenvectors are happily printed.
        shape_alpha = (self.data.homos[0]+6, self.data.nbasis)
        shape_beta = (self.data.homos[1]+6, self.data.nbasis)
        self.assertEquals(self.data.mocoeffs[0].shape, shape_alpha)
        self.assertEquals(self.data.mocoeffs[1].shape, shape_beta)

    def testnooccnos(self):
        """Are natural orbital occupation numbers the right size?"""
        self.assertEquals(self.data.nooccnos.shape, (self.data.nmo, ))


class GamessUK80SPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testnooccnos(self):
        """Are natural orbital occupation numbers the right size?"""
        self.assertEquals(self.data.nooccnos.shape, (self.data.nmo, ))


class GaussianSPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testatomnos(self):
        """Does atomnos have the right dimension (20)?"""
        size = len(self.data.atomnos)
        self.assertEquals(size, 20)


class JaguarSPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testmoenergies(self):
        """Are the dims of the moenergies equal to 2 x homos+11?"""
        self.assertEquals(len(self.data.moenergies), 2)
        self.assertEquals(len(self.data.moenergies[0]), self.data.homos[0]+11)
        self.assertEquals(len(self.data.moenergies[1]), self.data.homos[1]+11)

    def testmosyms(self):
        """Are the dims of the mosyms equals to 2 x nmo?"""
        shape0 = (len(self.data.mosyms), len(self.data.mosyms[0]))
        shape1 = (len(self.data.mosyms), len(self.data.mosyms[1]))
        self.assertEquals(shape0, (2, self.data.homos[0]+11))
        self.assertEquals(shape1, (2, self.data.homos[1]+11))


if __name__=="__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['SPun'])
    suite.testall()
