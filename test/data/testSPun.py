# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test unrestrictied single point logfiles in cclib"""

import numpy
from skip import skipForParser


class GenericSPunTest:
    """Generic unrestricted single point unittest"""

    def testnatom(self, data) -> None:
        """Is the number of atoms equal to 20?"""
        assert data.natom == 20

    def testatomnos(self, data) -> None:
        """Are the atomnos correct?"""
        assert numpy.all([numpy.issubdtype(atomno, numpy.signedinteger) for atomno in data.atomnos])
        assert data.atomnos.shape == (20,)
        assert sum(data.atomnos == 6) + sum(data.atomnos == 1) == 20

    @skipForParser("ADF", "???")
    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser(
        "DALTON",
        "DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now",
    )
    @skipForParser("FChk", "The parser is still being developed so we skip this test")
    @skipForParser("Jaguar", "???")
    @skipForParser("Molcas", "Length is zero for some reason")
    @skipForParser("Molpro", "???")
    @skipForParser("Turbomole", "???")
    @skipForParser("Serenity", "not implemented yet")
    def testatomcharges(self, data) -> None:
        """Are atomic charges consistent with natom?"""
        for atomcharge_type in data.atomcharges:
            charges = data.atomcharges[atomcharge_type]
            natom = data.natom
            assert len(charges) == natom, (
                f"len(atomcharges['{atomcharge_type}']) = {len(charges)}, natom = {natom}"
            )

    @skipForParser("ADF", "???")
    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser(
        "DALTON",
        "DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now",
    )
    @skipForParser("FChk", "The parser is still being developed so we skip this test")
    @skipForParser("Jaguar", "???")
    @skipForParser("Molcas", "Length is zero for some reason")
    @skipForParser("Molpro", "???")
    @skipForParser("Turbomole", "???")
    @skipForParser("Serenity", "not included in testfile yet")
    def testatomcharges_mulliken(self, data) -> None:
        """Do Mulliken atomic charges sum to positive one?"""
        charges = data.atomcharges["mulliken"]
        assert abs(sum(charges) - 1.0) < 1.0e-2

    def testatomcoords(self, data) -> None:
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        assert data.atomcoords.shape == (1, data.natom, 3)

    @skipForParser("Jaguar", "Data file does not contain enough information")
    def testdimmocoeffs(self, data) -> None:
        """Are the dimensions of mocoeffs equal to 2 x nmo x nbasis?"""
        if hasattr(data, "mocoeffs"):
            assert isinstance(data.mocoeffs, list)
            assert len(data.mocoeffs) == 2
            assert data.mocoeffs[0].shape == (data.nmo, data.nbasis)
            assert data.mocoeffs[1].shape == (data.nmo, data.nbasis)

    @skipForParser("Jaguar", "Data file does not contain enough information")
    @skipForParser("DALTON", "mocoeffs not implemented yet")
    def testfornoormo(self, data) -> None:
        """Do we have NOs or MOs?"""
        assert hasattr(data, "nocoeffs") or hasattr(data, "mocoeffs")

    @skipForParser("Serenity", "no NOs in Serenity")
    def testdimnoccnos(self, data) -> None:
        """Is the length of nooccnos equal to nmo?"""
        if hasattr(data, "nooccnos"):
            assert isinstance(data.nooccnos, numpy.ndarray)
            # FIXME
            assert data.nooccnos.shape in [(data.nmo,), (2, data.nmo)]

    @skipForParser("Serenity", "no NOs in Serenity")
    def testdimnocoeffs(self, data) -> None:
        """Are the dimensions of nocoeffs equal to 2 x nmo x nmo?"""
        if hasattr(data, "nocoeffs"):
            assert isinstance(data.nocoeffs, numpy.ndarray)
            assert data.nocoeffs.shape == (2, data.nmo, data.nmo)

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    def testcharge_and_mult(self, data) -> None:
        """Are the charge and multiplicity correct?"""
        assert data.charge == 1
        assert data.mult == 2

    def testhomos(self, data) -> None:
        """Are the homos correct?"""
        msg = f"{numpy.array_repr(data.homos)} != array([34,33],'i')"
        numpy.testing.assert_array_equal(data.homos, numpy.array([34, 33], "i"), msg)

    def testmoenergies(self, data) -> None:
        """Are the dims of the moenergies equals to 2 x nmo?"""
        if hasattr(data, "moenergies"):
            assert len(data.moenergies) == 2
            assert len(data.moenergies[0]) == data.nmo
            assert len(data.moenergies[1]) == data.nmo

    @skipForParser("FChk", "Fchk files do not have a section for symmetry")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Molpro", "?")
    @skipForParser("ORCA", "ORCA has no support for symmetry yet")
    @skipForParser("Serenity", "Serenity does not use symmetry.")
    def testmosyms(self, data) -> None:
        """Are the dims of the mosyms equals to 2 x nmo?"""
        shape = (len(data.mosyms), len(data.mosyms[0]))
        assert shape == (2, data.nmo)


class GenericROSPTest(GenericSPunTest):
    """Customized restricted open-shell single point unittest"""

    @skipForParser("DALTON", "mocoeffs not implemented yet")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testdimmocoeffs(self, data) -> None:
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        assert isinstance(data.mocoeffs, list)
        assert len(data.mocoeffs) == 1
        assert data.mocoeffs[0].shape == (data.nmo, data.nbasis)

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testhomos(self, data) -> None:
        """Are the HOMO indices equal to 34 and 33 (one more alpha electron
        than beta electron)?
        """
        msg = f"{numpy.array_repr(data.homos)} != array([34, 33], 'i')"
        numpy.testing.assert_array_equal(data.homos, numpy.array([34, 33], "i"), msg)

    @skipForParser("QChem", "prints 2 sets of different MO energies?")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testmoenergies(self, data) -> None:
        """Are the dims of the moenergies equals to 1 x nmo?"""
        assert len(data.moenergies) == 1
        assert len(data.moenergies[0]) == data.nmo

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    @skipForParser("Serenity", "Serenity does not use symmetry.")
    def testmosyms(self, data) -> None:
        """Are the dims of the mosyms equals to 1 x nmo?"""
        shape = (len(data.mosyms), len(data.mosyms[0]))
        assert shape == (1, data.nmo)


class GamessUK70SPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testdimmocoeffs(self, data) -> None:
        """Are the dimensions of mocoeffs equal to 2 x (homos+6) x nbasis?"""

        assert isinstance(data.mocoeffs, list)
        assert len(data.mocoeffs) == 2

        # This is only an issue in version 7.0 (and before?), since in the version 8.0
        # logfile all eigenvectors are happily printed.
        shape_alpha = (data.homos[0] + 6, data.nbasis)
        shape_beta = (data.homos[1] + 6, data.nbasis)
        assert data.mocoeffs[0].shape == shape_alpha
        assert data.mocoeffs[1].shape == shape_beta

    def testnooccnos(self, data) -> None:
        """Are natural orbital occupation numbers the right size?"""
        assert data.nooccnos.shape == (data.nmo,)


class GamessUK80SPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testnooccnos(self, data) -> None:
        """Are natural orbital occupation numbers the right size?"""
        assert data.nooccnos.shape == (data.nmo,)


class GaussianSPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testatomspins(self, data) -> None:
        """Are atomic spins from Mulliken population analysis consistent with
        natom and sum to one (doublet)?
        """
        spins = data.atomspins["mulliken"]
        assert len(spins) == data.natom
        assert abs(sum(spins) - 1.0) < 0.001


class JaguarSPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testmoenergies(self, data) -> None:
        """Are the dims of the moenergies equal to 2 x homos+11?"""
        assert len(data.moenergies) == 2
        assert len(data.moenergies[0]) == data.homos[0] + 11
        assert len(data.moenergies[1]) == data.homos[1] + 11

    def testmosyms(self, data) -> None:
        """Are the dims of the mosyms equals to 2 x nmo?"""
        shape0 = (len(data.mosyms), len(data.mosyms[0]))
        shape1 = (len(data.mosyms), len(data.mosyms[1]))
        assert shape0 == (2, data.homos[0] + 11)
        assert shape1 == (2, data.homos[1] + 11)
