# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
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

    def testnatom(self):
        """Is the number of atoms equal to 20?"""
        assert self.data.natom == 20

    def testatomnos(self):
        """Are the atomnos correct?"""
        assert numpy.alltrue([numpy.issubdtype(atomno, numpy.signedinteger)
                              for atomno in self.data.atomnos])
        assert self.data.atomnos.shape == (20,)
        assert sum(self.data.atomnos == 6) + sum(self.data.atomnos == 1) == 20

    @skipForParser('ADF', '???')
    @skipForParser('DALTON', 'DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now')
    @skipForParser('FChk', 'The parser is still being developed so we skip this test')
    @skipForParser('Jaguar', '???')
    @skipForParser('Molcas', 'Length is zero for some reason')
    @skipForParser('Molpro', '???')
    @skipForParser('Turbomole', '???')
    def testatomcharges(self):
        """Are atomic charges consistent with natom?"""
        for atomcharge_type in self.data.atomcharges:
            charges = self.data.atomcharges[atomcharge_type]
            natom = self.data.natom
            assert len(charges) == natom, f"len(atomcharges['{atomcharge_type}']) = {len(charges)}, natom = {natom}"

    @skipForParser('ADF', '???')
    @skipForParser('DALTON', 'DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now')
    @skipForParser('FChk', 'The parser is still being developed so we skip this test')
    @skipForParser('Jaguar', '???')
    @skipForParser('Molcas', 'Length is zero for some reason')
    @skipForParser('Molpro', '???')
    @skipForParser('Turbomole', '???')
    def testatomcharges_mulliken(self):
        """Do Mulliken atomic charges sum to positive one?"""
        charges = self.data.atomcharges["mulliken"]
        assert abs(sum(charges)-1.0) < 1.0e-2

    def testatomcoords(self):
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        assert self.data.atomcoords.shape == (1,self.data.natom,3)

    @skipForParser('Jaguar', 'Data file does not contain enough information')
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 2 x nmo x nbasis?"""
        if hasattr(self.data, "mocoeffs"):
            assert isinstance(self.data.mocoeffs, list)
            assert len(self.data.mocoeffs) == 2
            assert self.data.mocoeffs[0].shape == (self.data.nmo, self.data.nbasis)
            assert self.data.mocoeffs[1].shape == (self.data.nmo, self.data.nbasis)

    @skipForParser('Jaguar', 'Data file does not contain enough information')
    @skipForParser('DALTON', 'mocoeffs not implemented yet')
    def testfornsoormo(self):
        """Do we have NSOs or MOs?"""
        assert hasattr(self.data, "nsocoeffs") or hasattr(self.data, "mocoeffs")

    def testdimnsoccnos(self):
        """Are the dimensions of nsooccnos equal to 2 x nmo?"""
        if hasattr(self.data, "nsooccnos"):
            assert isinstance(self.data.nsooccnos, list)
            assert isinstance(self.data.nsooccnos[0], list)
            assert isinstance(self.data.nsooccnos[1], list)
            assert len(self.data.nsooccnos) == 2
            assert len(self.data.nsooccnos[0]) == self.data.nmo
            assert len(self.data.nsooccnos[1]) == self.data.nmo

    def testdimnsocoeffs(self):
        """Are the dimensions of nsocoeffs equal to 2 x nmo x nmo?"""
        if hasattr(self.data, "nsocoeffs"):
            assert isinstance(self.data.nsocoeffs, list)
            assert isinstance(self.data.nsocoeffs[0], numpy.ndarray)
            assert isinstance(self.data.nsocoeffs[1], numpy.ndarray)
            assert len(self.data.nsocoeffs) == 2
            assert self.data.nsocoeffs[0].shape == (self.data.nmo, self.data.nmo)
            assert self.data.nsocoeffs[1].shape == (self.data.nmo, self.data.nmo)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    def testcharge_and_mult(self):
        """Are the charge and multiplicity correct?"""
        assert self.data.charge == 1
        assert self.data.mult == 2

    def testhomos(self):
        """Are the homos correct?"""
        msg = f"{numpy.array_repr(self.data.homos)} != array([34,33],'i')"
        numpy.testing.assert_array_equal(self.data.homos, numpy.array([34,33],"i"), msg)

    def testmoenergies(self):
        """Are the dims of the moenergies equals to 2 x nmo?"""
        if hasattr(self.data, "moenergies"):
            assert len(self.data.moenergies) == 2
            assert len(self.data.moenergies[0]) == self.data.nmo
            assert len(self.data.moenergies[1]) == self.data.nmo

    @skipForParser('FChk', 'Fchk files do not have a section for symmetry')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Molpro', '?')
    @skipForParser('ORCA', 'ORCA has no support for symmetry yet')
    def testmosyms(self):
        """Are the dims of the mosyms equals to 2 x nmo?"""
        shape = (len(self.data.mosyms), len(self.data.mosyms[0]))
        assert shape == (2, self.data.nmo)


class GenericROSPTest(GenericSPunTest):
    """Customized restricted open-shell single point unittest"""

    @skipForParser('DALTON', 'mocoeffs not implemented yet')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        assert isinstance(self.data.mocoeffs, list)
        assert len(self.data.mocoeffs) == 1
        assert self.data.mocoeffs[0].shape == (self.data.nmo, self.data.nbasis)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testhomos(self):
        """Are the HOMO indices equal to 34 and 33 (one more alpha electron
        than beta electron)?
        """
        msg = f"{numpy.array_repr(self.data.homos)} != array([34, 33], 'i')"
        numpy.testing.assert_array_equal(self.data.homos, numpy.array([34, 33], "i"), msg)

    @skipForParser('QChem', 'prints 2 sets of different MO energies?')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testmoenergies(self):
        """Are the dims of the moenergies equals to 1 x nmo?"""
        assert len(self.data.moenergies) == 1
        assert len(self.data.moenergies[0]) == self.data.nmo

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testmosyms(self):
        """Are the dims of the mosyms equals to 1 x nmo?"""
        shape = (len(self.data.mosyms), len(self.data.mosyms[0]))
        assert shape == (1, self.data.nmo)


class GamessUK70SPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 2 x (homos+6) x nbasis?"""

        assert isinstance(self.data.mocoeffs, list)
        assert len(self.data.mocoeffs) == 2

        # This is only an issue in version 7.0 (and before?), since in the version 8.0
        # logfile all eigenvectors are happily printed.
        shape_alpha = (self.data.homos[0]+6, self.data.nbasis)
        shape_beta = (self.data.homos[1]+6, self.data.nbasis)
        assert self.data.mocoeffs[0].shape == shape_alpha
        assert self.data.mocoeffs[1].shape == shape_beta

    def testnooccnos(self):
        """Are natural orbital occupation numbers the right size?"""
        assert self.data.nooccnos.shape == (self.data.nmo, )


class GamessUK80SPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testnooccnos(self):
        """Are natural orbital occupation numbers the right size?"""
        assert self.data.nooccnos.shape == (self.data.nmo, )


class GaussianSPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testatomspins(self):
        """Are atomic spins from Mulliken population analysis consistent with
        natom and sum to one (doublet)?
        """
        spins = self.data.atomspins['mulliken']
        assert len(spins) == self.data.natom
        assert abs(sum(spins)-1.0) < 0.001

            
class JaguarSPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testmoenergies(self):
        """Are the dims of the moenergies equal to 2 x homos+11?"""
        assert len(self.data.moenergies) == 2
        assert len(self.data.moenergies[0]) == self.data.homos[0]+11
        assert len(self.data.moenergies[1]) == self.data.homos[1]+11

    def testmosyms(self):
        """Are the dims of the mosyms equals to 2 x nmo?"""
        shape0 = (len(self.data.mosyms), len(self.data.mosyms[0]))
        shape1 = (len(self.data.mosyms), len(self.data.mosyms[1]))
        assert shape0 == (2, self.data.homos[0]+11)
        assert shape1 == (2, self.data.homos[1]+11)


if __name__=="__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['SPun'])
    suite.testall()
