# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test unrestrictied single point logfiles in cclib"""

import os
import unittest

import numpy

from skip import skipForParser
from skip import skipForLogfile


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericSPunTest(unittest.TestCase):
    """Generic unrestricted single point unittest"""

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testnatom(self):
        """Is the number of atoms equal to 20?"""
        self.assertEqual(self.data.natom, 20)

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testatomnos(self):
        """Are the atomnos correct?"""
        self.assertTrue(numpy.alltrue([numpy.issubdtype(atomno, numpy.signedinteger)
                                       for atomno in self.data.atomnos]))
        self.assertEqual(self.data.atomnos.shape, (20,) )
        self.assertEqual(sum(self.data.atomnos==6) + sum(self.data.atomnos==1), 20)

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testatomcoords(self):
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        self.assertEqual(self.data.atomcoords.shape,(1,self.data.natom,3))

    @skipForParser('Jaguar', 'Data file does not contain enough information')
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 2 x nmo x nbasis?"""
        if hasattr(self.data, "mocoeffs"):
            self.assertIsInstance(self.data.mocoeffs, list)
            self.assertEqual(len(self.data.mocoeffs), 2)
            self.assertEqual(self.data.mocoeffs[0].shape,
                             (self.data.nmo, self.data.nbasis))
            self.assertEqual(self.data.mocoeffs[1].shape,
                             (self.data.nmo, self.data.nbasis))

    @skipForParser('Jaguar', 'Data file does not contain enough information')
    @skipForParser('DALTON', 'mocoeffs not implemented yet')
    def testfornsoormo(self):
        """Do we have NSOs or MOs?"""
        self.assertEquals(hasattr(self.data, "nsocoeffs") or hasattr(self.data, "mocoeffs"), True)

    def testdimnsoccnos(self):
        """Are the dimensions of nsooccnos equal to 2 x nmo?"""
        if hasattr(self.data, "nsooccnos"):
            self.assertIsInstance(self.data.nsooccnos, list)
            self.assertIsInstance(self.data.nsooccnos[0], list)
            self.assertIsInstance(self.data.nsooccnos[1], list)
            self.assertEquals(len(self.data.nsooccnos), 2)
            self.assertEquals(len(self.data.nsooccnos[0]), self.data.nmo)
            self.assertEquals(len(self.data.nsooccnos[1]), self.data.nmo)

    def testdimnsocoeffs(self):
        """Are the dimensions of nsocoeffs equal to 2 x nmo x nmo?"""
        if hasattr(self.data, "nsocoeffs"):
            self.assertIsInstance(self.data.nsocoeffs, list)
            self.assertIsInstance(self.data.nsocoeffs[0], numpy.ndarray)
            self.assertIsInstance(self.data.nsocoeffs[1], numpy.ndarray)
            self.assertEquals(len(self.data.nsocoeffs), 2)
            self.assertEquals(self.data.nsocoeffs[0].shape, (self.data.nmo, self.data.nmo))
            self.assertEquals(self.data.nsocoeffs[1].shape, (self.data.nmo, self.data.nmo))

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testcharge_and_mult(self):
        """Are the charge and multiplicity correct?"""
        self.assertEqual(self.data.charge, 1)
        self.assertEqual(self.data.mult, 2)

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testhomos(self):
        """Are the homos correct?"""
        msg = "%s != array([34,33],'i')" % numpy.array_repr(self.data.homos)
        numpy.testing.assert_array_equal(self.data.homos, numpy.array([34,33],"i"), msg)

    def testmoenergies(self):
        """Are the dims of the moenergies equals to 2 x nmo?"""
        if hasattr(self.data, "moenergies"):
            self.assertEqual(len(self.data.moenergies), 2)
            self.assertEqual(len(self.data.moenergies[0]), self.data.nmo)
            self.assertEqual(len(self.data.moenergies[1]), self.data.nmo)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Molpro', '?')
    @skipForParser('ORCA', 'ORCA has no support for symmetry yet')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testmosyms(self):
        """Are the dims of the mosyms equals to 2 x nmo?"""
        shape = (len(self.data.mosyms), len(self.data.mosyms[0]))
        self.assertEqual(shape, (2, self.data.nmo))


class GenericROSPTest(GenericSPunTest):
    """Customized restricted open-shell single point unittest"""

    @skipForParser('DALTON', 'mocoeffs not implemented yet')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        self.assertEqual(type(self.data.mocoeffs), type([]))
        self.assertEqual(len(self.data.mocoeffs), 1)
        self.assertEqual(self.data.mocoeffs[0].shape,
                          (self.data.nmo, self.data.nbasis))

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testhomos(self):
        """Are the HOMO indices equal to 34 and 33 (one more alpha electron
        than beta electron)?
        """
        msg = "%s != array([34, 33], 'i')" % numpy.array_repr(self.data.homos)
        numpy.testing.assert_array_equal(self.data.homos, numpy.array([34, 33], "i"), msg)

    @skipForParser('QChem', 'prints 2 sets of different MO energies?')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testmoenergies(self):
        """Are the dims of the moenergies equals to 1 x nmo?"""
        self.assertEqual(len(self.data.moenergies), 1)
        self.assertEqual(len(self.data.moenergies[0]), self.data.nmo)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testmosyms(self):
        """Are the dims of the mosyms equals to 1 x nmo?"""
        shape = (len(self.data.mosyms), len(self.data.mosyms[0]))
        self.assertEqual(shape, (1, self.data.nmo))


class GamessUK70SPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 2 x (homos+6) x nbasis?"""

        self.assertEqual(type(self.data.mocoeffs), type([]))
        self.assertEqual(len(self.data.mocoeffs), 2)

        # This is only an issue in version 7.0 (and before?), since in the version 8.0
        # logfile all eigenvectors are happily printed.
        shape_alpha = (self.data.homos[0]+6, self.data.nbasis)
        shape_beta = (self.data.homos[1]+6, self.data.nbasis)
        self.assertEqual(self.data.mocoeffs[0].shape, shape_alpha)
        self.assertEqual(self.data.mocoeffs[1].shape, shape_beta)

    def testnooccnos(self):
        """Are natural orbital occupation numbers the right size?"""
        self.assertEqual(self.data.nooccnos.shape, (self.data.nmo, ))


class GamessUK80SPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testnooccnos(self):
        """Are natural orbital occupation numbers the right size?"""
        self.assertEqual(self.data.nooccnos.shape, (self.data.nmo, ))


class GaussianSPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testatomnos(self):
        """Does atomnos have the right dimension (20)?"""
        size = len(self.data.atomnos)
        self.assertEqual(size, 20)


class JaguarSPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testmoenergies(self):
        """Are the dims of the moenergies equal to 2 x homos+11?"""
        self.assertEqual(len(self.data.moenergies), 2)
        self.assertEqual(len(self.data.moenergies[0]), self.data.homos[0]+11)
        self.assertEqual(len(self.data.moenergies[1]), self.data.homos[1]+11)

    def testmosyms(self):
        """Are the dims of the mosyms equals to 2 x nmo?"""
        shape0 = (len(self.data.mosyms), len(self.data.mosyms[0]))
        shape1 = (len(self.data.mosyms), len(self.data.mosyms[1]))
        self.assertEqual(shape0, (2, self.data.homos[0]+11))
        self.assertEqual(shape1, (2, self.data.homos[1]+11))


if __name__=="__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['SPun'])
    suite.testall()
