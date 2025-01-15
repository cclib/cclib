# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test single point logfiles in cclib."""

import datetime

from cclib.parser import utils

import numpy
import packaging
from common import get_minimum_carbon_separation
from skip import skipForLogfile, skipForParser


class GenericSPTest:
    """Generic restricted single point unittest"""

    # Molecular mass of DVB in mD, and expected precision.
    molecularmass = 130078.25
    mass_precision = 0.10

    # In STO-3G, H has 1, C has 5 (1 S and 4 SP).
    nbasisdict = {1: 1, 6: 5}

    # Approximate B3LYP energy of dvb after SCF in STO-3G (Gaussian 16).
    scfenergy = -382.308266602
    scfenergy_delta = 3.0e-1

    # Approximate energy of the innermost molecular orbital of DVB with
    # B3LYP/STO-3G (from Q-Chem 5.4 fchk).
    moenergy = -10.0179353
    moenergy_delta = 3.0e-3

    # Overlap first two atomic orbitals.
    overlap01 = 0.24

    # Generally, one criteria for SCF energy convergence.
    num_scf_criteria = 1

    def testnatom(self, data) -> None:
        """Is the number of atoms equal to 20?"""
        assert data.natom == 20

    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testatomnos(self, data) -> None:
        """Are the atomnos correct?"""

        # The nuclear charges should be integer values in a NumPy array.
        assert numpy.all([numpy.issubdtype(atomno, numpy.signedinteger) for atomno in data.atomnos])
        assert data.atomnos.dtype.char == "i"

        assert data.atomnos.shape == (20,)
        assert sum(data.atomnos == 6) + sum(data.atomnos == 1) == 20

    @skipForParser(
        "DALTON",
        "DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now",
    )
    @skipForParser("FChk", "The parser is still being developed so we skip this test")
    @skipForParser("GAMESSDAT", "We are not going to implement atom charges in this version.")
    @skipForLogfile(
        "Jaguar/basicJaguar7",
        "We did not print the atomic partial charges in the unit tests for this version",
    )
    @skipForLogfile(
        "Molpro/basicMolpro2006",
        "These tests were run a long time ago and since we don't have access to Molpro 2006 anymore, we can skip this test (it is tested in 2012)",
    )
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    @skipForParser("xTB", "not implemented yet")
    def testatomcharges(self, data) -> None:
        """Are atomic charges consistent with natom?"""
        for atomcharge_type in data.atomcharges:
            charges = data.atomcharges[atomcharge_type]
            natom = data.natom
            assert len(charges) == natom, (
                f"len(atomcharges['{atomcharge_type}']) = {len(charges)}, natom = {natom}"
            )

    @skipForParser(
        "DALTON",
        "DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now",
    )
    @skipForLogfile("FChk/basicQChem5.2", "not printed for Q-Chem")
    @skipForLogfile("FChk/basicQChem5.4", "not printed for Q-Chem")
    @skipForParser(
        "GAMESSDAT",
        "We are not sure about the specific type of atom charges, it is best to skip the test for now.",
    )
    @skipForLogfile(
        "Jaguar/basicJaguar7",
        "We did not print the atomic partial charges in the unit tests for this version",
    )
    @skipForLogfile(
        "Molpro/basicMolpro2006",
        "These tests were run a long time ago and since we don't have access to Molpro 2006 anymore, we can skip this test (it is tested in 2012)",
    )
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testatomcharges_mulliken(self, data) -> None:
        """Do Mulliken atomic charges sum to zero?"""
        charges = data.atomcharges["mulliken"]
        assert abs(sum(charges)) < 1.0e-2

    @skipForParser("ADF", "Lowdin charges not present by default")
    @skipForParser(
        "DALTON",
        "DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now",
    )
    @skipForParser("FChk", "The parser is still being developed so we skip this test")
    @skipForParser(
        "GAMESSDAT",
        "We are not sure about the specific type of atom charges, it is best to skip the test for now.",
    )
    @skipForParser("Gaussian", "Lowdin charges not present by default")
    @skipForParser("Jaguar", "Lowdin charges not present by default")
    @skipForParser("NWChem", "Lowdin charges not present by default")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("Molcas", "Lowdin charges not present by default")
    @skipForParser("Molpro", "Lowdin charges not present by default")
    @skipForParser("QChem", "Lowdin charges not present by default")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    @skipForParser("xTB", "not implemented yet")
    def testatomcharges_lowdin(self, data) -> None:
        """Do Lowdin atomic charges sum to zero?"""
        charges = data.atomcharges["lowdin"]
        assert abs(sum(charges)) < 1.0e-2

    @skipForParser("ADF", "Hirshfeld charges not implemented")
    @skipForParser(
        "DALTON",
        "DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now",
    )
    @skipForParser("FChk", "The parser is still being developed so we skip this test")
    @skipForParser("Gaussian", "Hirshfeld charges not implemented")
    @skipForParser("GAMESS", "Hirshfeld charges not implemented")
    @skipForParser("GAMESSUK", "Hirshfeld charges not implemented")
    @skipForParser(
        "GAMESSDAT",
        "We are not sure about the specific type of atom charges, it is best to skip the test for now.",
    )
    @skipForParser("Jaguar", "Hirshfeld charges not implemented")
    @skipForParser("NWChem", "Hirshfeld charges not implemented")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("Molcas", "Hirshfeld charges not implemented")
    @skipForParser("Molpro", "Hirshfeld charges not implemented")
    @skipForLogfile("ORCA/basicORCA4.1", "This needs to be moved to regressions")
    @skipForParser("Psi4", "Hirshfeld charges not implemented")
    @skipForParser("QChem", "Hirshfeld charges not implemented")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    @skipForParser("xTB", "not implemented yet")
    def testatomcharges_hirshfeld(self, data) -> None:
        """Do Hirshfeld atomic charges sum to roughly zero?"""
        charges = data.atomcharges["hirshfeld"]
        assert abs(sum(charges)) < 4.0e-3

    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testatomcoords(self, data) -> None:
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        expected_shape = (1, data.natom, 3)
        assert data.atomcoords.shape == expected_shape

    @skipForParser(
        "GAMESSDAT", "Vectors need some calculations to transform them. Current mm value is 2.54"
    )
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testatomcoords_units(self, data) -> None:
        """Are atomcoords consistent with Angstroms?"""
        min_carbon_dist = get_minimum_carbon_separation(data)
        dev = abs(min_carbon_dist - 1.34)
        assert dev < 0.03, f"Minimum carbon dist is {min_carbon_dist:.2f} (not 1.34)"

    @skipForParser("GAMESSDAT", "Neither charge nor mult exists in the files.")
    @skipForParser("Molcas", "missing mult")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testcharge_and_mult(self, data) -> None:
        """Are the charge and multiplicity correct?"""
        assert data.charge == 0
        assert data.mult == 1

    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testnbasis(self, data) -> None:
        """Is the number of basis set functions correct?"""
        count = sum([self.nbasisdict[n] for n in data.atomnos])
        assert data.nbasis == count

    @skipForParser("ADF", "ADF parser does not extract atombasis")
    @skipForLogfile(
        "Jaguar/basicJaguar7",
        "Data file does not contain enough information. Can we make a new one?",
    )
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    @skipForParser("xTB", "not implemented yet")
    def testatombasis(self, data) -> None:
        """Are the indices in atombasis the right amount and unique?"""
        all = []
        for i, atom in enumerate(data.atombasis):
            assert len(atom) == self.nbasisdict[data.atomnos[i]]
            all += atom
        # Test if there are as many indices as atomic orbitals.
        assert len(all) == data.nbasis
        # Check if all are different (every orbital indexed once).
        assert len(set(all)) == len(all)

    @skipForLogfile("FChk/basicQChem5.2", "Q-Chem doesn't print SCF energy to fchk")
    @skipForLogfile("FChk/basicQChem5.4", "Q-Chem doesn't print SCF energy to fchk")
    @skipForParser("GAMESS", "atommasses not implemented yet")
    @skipForParser("GAMESSUK", "atommasses not implemented yet")
    @skipForParser("GAMESSDAT", "Atommasses implemented, but it does not pass the test.")
    @skipForParser("Jaguar", "atommasses not implemented yet")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Molpro", "atommasses not implemented yet")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForLogfile("Psi4/basicPsi4.0b5", "atommasses not implemented yet")
    @skipForParser("QChem", "atommasses not implemented yet")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    @skipForParser("xTB", "not implemented yet")
    def testatommasses(self, data) -> None:
        """Do the atom masses sum up to the molecular mass?"""
        mm = 1000 * sum(data.atommasses)
        msg = f"Molecule mass: {mm:f} not {self.molecularmass:f} +- {self.mass_precision:f}mD"
        assert abs(mm - self.molecularmass) < self.mass_precision, msg

    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testcoreelectrons(self, data) -> None:
        """Are the coreelectrons all 0?"""
        ans = numpy.zeros(data.natom, "i")
        numpy.testing.assert_array_equal(data.coreelectrons, ans)

    @skipForParser("FChk", "Formatted checkpoint files do not have a section for symmetry")
    @skipForParser("GAMESSDAT", "Mosyms do not exist in the file")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Molpro", "?")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testsymlabels(self, data) -> None:
        """Are all the symmetry labels either Ag/u or Bg/u?"""
        sumwronglabels = sum([x not in ["Ag", "Bu", "Au", "Bg"] for x in data.mosyms[0]])
        assert sumwronglabels == 0

    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "xTB does not print them all")
    def testhomos(self, data) -> None:
        """Is the index of the HOMO equal to 34?"""
        numpy.testing.assert_array_equal(
            data.homos, numpy.array([34], "i"), f"{numpy.array_repr(data.homos)} != array([34],'i')"
        )

    @skipForParser("FChk", "Formatted Checkpoint files do not have a section for SCF energy")
    @skipForParser("GAMESSDAT", "Scfvalues probably do not exist in the file")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testscfvaluetype(self, data) -> None:
        """Are scfvalues and its elements the right type??"""
        assert isinstance(data.scfvalues, list)
        assert isinstance(data.scfvalues[0], numpy.ndarray)

    @skipForLogfile("FChk/basicQChem5.2", "Q-Chem doesn't print SCF energy to fchk")
    @skipForLogfile("FChk/basicQChem5.4", "Q-Chem doesn't print SCF energy to fchk")
    @skipForParser("GAMESSDAT", "Scfenergies probably do not exist in the file")
    @skipForParser("NBO", "attribute not implemented in this version")
    def testscfenergy(self, data) -> None:
        """Is the SCF energy within the target?"""
        assert abs(
            data.scfenergies[-1] - utils.convertor(self.scfenergy, "hartree", "eV")
        ) < utils.convertor(self.scfenergy_delta, "hartree", "eV")

    @skipForParser("FChk", "Formatted Checkpoint files do not have a section for SCF convergence")
    @skipForParser("GAMESSDAT", "Scftargets probably do not exist in the file")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testscftargetdim(self, data) -> None:
        """Do the scf targets have the right dimensions?"""
        assert data.scftargets.shape == (len(data.scfvalues), len(data.scfvalues[0][0]))

    @skipForParser("FChk", "Formatted Checkpoint files do not have a section for SCF convergence")
    @skipForParser("GAMESSDAT", "Scftargets probably do not exist in the file")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testscftargets(self, data) -> None:
        """Are correct number of SCF convergence criteria being parsed?"""
        assert len(data.scftargets[0]) == self.num_scf_criteria

    @skipForParser("xTB", "not implemented yet")
    def testlengthmoenergies(self, data) -> None:
        """Is the number of evalues equal to nmo?"""
        if hasattr(data, "moenergies"):
            assert len(data.moenergies[0]) == data.nmo

    def testtypemoenergies(self, data) -> None:
        """Is moenergies a list containing one numpy array?"""
        if hasattr(data, "moenergies"):
            assert isinstance(data.moenergies, list)
            assert isinstance(data.moenergies[0], numpy.ndarray)

    @skipForParser("GAMESSDAT", "Moenergies probably do not exist in the file")
    @skipForLogfile("Gaussian/basicGaussian16/dvb_sp_no.out", "no energies for natural orbitals")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForLogfile("Turbomole/basicTurbomole5.9/dvb_sp_symm", "delta of 7.4, everything else ok")
    @skipForParser("xTB", "not implemented yet")
    def testfirstmoenergy(self, data) -> None:
        """Is the lowest energy molecular orbital within the target?"""
        assert abs(
            data.moenergies[0][0] - utils.convertor(self.moenergy, "hartree", "eV")
        ) < utils.convertor(self.moenergy_delta, "hartree", "eV")

    @skipForParser("DALTON", "mocoeffs not implemented yet")
    @skipForLogfile(
        "Jaguar/basicJaguar7",
        "Data file does not contain enough information. Can we make a new one?",
    )
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("Turbomole", "Use of symmetry has reduced the number of mo coeffs")
    def testdimmocoeffs(self, data) -> None:
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        if hasattr(data, "mocoeffs"):
            assert isinstance(data.mocoeffs, list)
            assert len(data.mocoeffs) == 1
            assert data.mocoeffs[0].shape == (data.nmo, data.nbasis)

    @skipForParser("DALTON", "mocoeffs not implemented yet")
    @skipForLogfile(
        "Jaguar/basicJaguar7",
        "Data file does not contain enough information. Can we make a new one?",
    )
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testfornoormo(self, data) -> None:
        """Do we have NOs or MOs?"""
        assert hasattr(data, "nocoeffs") or hasattr(data, "mocoeffs")

    @skipForParser("NBO", "attribute not implemented in this version")
    def testdimnoccnos(self, data) -> None:
        """Is the length of nooccnos equal to nmo?"""
        if hasattr(data, "nooccnos"):
            assert isinstance(data.nooccnos, numpy.ndarray)
            assert len(data.nooccnos) == data.nmo

    @skipForParser("NBO", "attribute not implemented in this version")
    def testdimnocoeffs(self, data) -> None:
        """Are the dimensions of nocoeffs equal to nmo x nmo?"""
        if hasattr(data, "nocoeffs"):
            assert isinstance(data.nocoeffs, numpy.ndarray)
            assert data.nocoeffs.shape == (data.nmo, data.nmo)

    @skipForParser("DALTON", "To print: **INTEGRALS\n.PROPRI")
    @skipForLogfile("FChk/basicGaussian09", "Only available in QChem")
    @skipForLogfile("FChk/basicGaussian16", "Only available in QChem")
    @skipForParser("GAMESSDAT", "Aooverlaps probably do not exist in the file.")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("Psi4", "Psi4 does not currently have the option to print the overlap matrix")
    @skipForParser("QChem", "QChem cannot print the overlap matrix")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    @skipForParser("xTB", "not implemented yet")
    def testaooverlaps(self, data) -> None:
        """Are the dims and values of the overlap matrix correct?"""

        assert data.aooverlaps.shape == (data.nbasis, data.nbasis)

        # The matrix is symmetric.
        row = data.aooverlaps[0, :]
        col = data.aooverlaps[:, 0]
        assert sum(col - row) == 0.0

        # All values on diagonal should be exactly one.
        for i in range(data.nbasis):
            assert data.aooverlaps[i, i] == 1.0

        # Check some additional values that don't seem to move around between programs.
        assert abs(data.aooverlaps[0, 1] - self.overlap01) < 0.01
        assert abs(data.aooverlaps[1, 0] - self.overlap01) < 0.01
        assert round(abs(data.aooverlaps[3, 0]), 7) == 0
        assert round(abs(data.aooverlaps[0, 3]), 7) == 0

    def testoptdone(self, data) -> None:
        """There should be no optdone attribute set."""
        assert not hasattr(data, "optdone")

    @skipForParser("ADF", "Not implemented yes")
    @skipForParser("DALTON", "Not implemented yes")
    @skipForParser("FChk", "Rotational constants are never written to fchk files")
    @skipForParser("GAMESS", "Not implemented yes")
    @skipForParser("GAMESSUK", "Not implemented yet")
    @skipForParser("GAMESSDAT", "Not implemented yet")
    @skipForParser("Jaguar", "Not implemented yet")
    @skipForParser("Molcas", "Not implemented yes")
    @skipForParser("Molpro", "Not implemented yes")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("NWChem", "Not implemented yes")
    @skipForParser("ORCA", "Not implemented yes")
    @skipForParser("Psi4", "Not implemented yes")
    @skipForParser("QChem", "Not implemented yes")
    @skipForParser("Turbomole", "Not implemented yes")
    @skipForParser("xTB", "not implemented yet")
    def testrotconsts(self, data) -> None:
        """A single geometry leads to single set of rotational constants."""
        assert data.rotconsts.shape == (1, 3)
        # taken from Gaussian16/dvb_sp.out
        ref = [4.6266363, 0.6849065, 0.5965900]
        numpy.testing.assert_allclose(data.rotconsts[0], ref, rtol=0, atol=1.0e-3)

    @skipForParser("FChk", "The parser is still being developed so we skip this test")
    @skipForParser("Gaussian", "Logfile needs to be updated")
    @skipForParser("Jaguar", "No dipole moments in the logfile")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testmoments(self, data) -> None:
        """Does the dipole and possible higher molecular moments look reasonable?"""

        # The reference point is always a vector, but not necessarily the
        # origin or center of mass. In this case, however, the center of mass
        # is at the origin, so we now what to expect.
        reference = data.moments[0]
        assert len(reference) == 3
        for x in reference:
            assert x == 0.0

        # Length and value of dipole moment should always be correct (zero for this test).
        dipole = data.moments[1]
        assert len(dipole) == 3
        for d in dipole:
            assert round(abs(d), 7) == 0

        # If the quadrupole is there, we can expect roughly -50B for the XX moment,
        # -50B for the YY moment and and -60B for the ZZ moment.
        if len(data.moments) > 2:
            quadrupole = data.moments[2]
            assert len(quadrupole) == 6
            assert abs(quadrupole[0] - -50) < 2.5
            assert abs(quadrupole[3] - -50) < 2.5
            assert abs(quadrupole[5] - -60) < 3

        # If the octupole is there, it should have 10 components and be zero.
        if len(data.moments) > 3:
            octupole = data.moments[3]
            assert len(octupole) == 10
            for m in octupole:
                assert abs(m) < 0.001

        # The hexadecapole should have 15 elements, an XXXX component of around -1900 Debye*ang^2,
        # a YYYY component of -330B and a ZZZZ component of -50B.
        if len(data.moments) > 4:
            hexadecapole = data.moments[4]
            assert len(hexadecapole) == 15
            assert abs(hexadecapole[0] - -1900) < 90
            assert abs(hexadecapole[10] - -330) < 11
            assert abs(hexadecapole[14] - -50) < 2.5

        # The are 21 unique 32-pole moments, and all are zero in this test case.
        if len(data.moments) > 5:
            moment32 = data.moments[5]
            assert len(moment32) == 21
            for m in moment32:
                assert m == 0.0

    @skipForParser("ADF", "reading basis set names is not implemented")
    @skipForParser("GAMESSDAT", "Basis set not implemented in this version.")
    @skipForParser("GAMESSUK", "reading basis set names is not implemented")
    @skipForParser("Molcas", "reading basis set names is not implemented")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("Psi4", "reading basis set names is not implemented")
    @skipForParser("xTB", "not implemented yet")
    def testmetadata_basis_set(self, data) -> None:
        """Does metadata have expected keys and values?"""
        assert data.metadata["basis_set"].lower() == "sto-3g"

    @skipForParser("ADF", "reading input file contents and name is not implemented")
    @skipForParser("DALTON", "reading input file contents and name is not implemented")
    @skipForParser("FChk", "Formatted checkpoint files do not have an input file section")
    @skipForParser("GAMESS", "reading input file contents and name is not implemented")
    @skipForParser("GAMESSUK", "reading input file contents and name is not implemented")
    @skipForParser("GAMESSDAT", "reading input file contents and name is not implemented")
    @skipForParser("Gaussian", "reading input file contents and name is not implemented")
    @skipForParser("Jaguar", "reading input file contents and name is not implemented")
    @skipForParser("Molcas", "reading input file contents and name is not implemented")
    @skipForParser("Molpro", "reading input file contents and name is not implemented")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("NWChem", "reading input file contents and name is not implemented")
    @skipForParser("Psi4", "reading input file contents and name is not implemented")
    @skipForParser("QChem", "reading input file contents and name is not implemented")
    @skipForParser("Turbomole", "reading input file contents and name is not implemented")
    @skipForParser("xTB", "not implemented yet")
    def testmetadata_input_file(self, data) -> None:
        """Does metadata have expected keys and values?"""
        assert "input_file_contents" in data.metadata
        # TODO make input file names consistent where possible, though some
        # programs do not allow arbitrary file extensions; for example, DALTON
        # must end in `dal`.
        assert "dvb_sp.in" in data.metadata["input_file_name"]

    @skipForParser("NBO", "attribute not implemented in this version")
    def testmetadata_methods(self, data) -> None:
        """Does metadata have expected keys and values?"""
        # TODO implement and unify across parsers; current values are [],
        # ["HF"], ["RHF"], and ["DFT"]
        assert "methods" in data.metadata

    @skipForParser("NBO", "attribute not implemented in this version")
    def testmetadata_package(self, data) -> None:
        """Does metadata have expected keys and values?"""
        # TODO How can the value be tested when the package name comes from
        # the parser and isn't stored on ccData?
        assert "package" in data.metadata

    @skipForParser(
        "FChk", "Formatted Checkpoint files do not have section for legacy package version"
    )
    @skipForParser("GAMESSDAT", "Files do not contain information about the legacy package version")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("xTB", "not implemented yet")
    def testmetadata_legacy_package_version(self, data) -> None:
        """Does metadata have expected keys and values?"""
        # TODO Test specific values for each unit test.
        assert "legacy_package_version" in data.metadata

    @skipForParser("FChk", "Formatted Checkpoint files do not have section for package version")
    @skipForParser("GAMESSDAT", "Files do not contain information about the package version")
    @skipForParser("NBO", "attribute not implemented in this version")
    def testmetadata_package_version(self, data) -> None:
        """Does metadata have expected keys and values?"""
        # TODO Test specific values for each unit test.
        assert isinstance(
            packaging.version.parse(data.metadata["package_version"]), packaging.version.Version
        )

    @skipForLogfile("NBO/basicNBO7.0/basicORCA5.0/dvb_sp.nbo.out", "TODO impossible to determine?")
    @skipForLogfile("FChk/basicGaussian09/dvb_sp.fchk", "impossible to determine")
    @skipForLogfile("FChk/basicQChem5.2/dvb_sp_modified.fchk", "impossible to determine")
    @skipForLogfile("FChk/basicQChem5.4/dvb_sp.fchk", "impossible to determine")
    @skipForLogfile("GAMESSDAT/basicGAMESS-US2018/dvb_sp.dat", "TODO impossible to determine?")
    def testmetadata_success(self, data) -> None:
        """Does metadata have expected keys and values?"""
        assert "success" in data.metadata
        assert data.metadata["success"]

    @skipForParser("FChk", "point group symmetry cannot be printed")
    @skipForParser("GAMESSDAT", "Files probably do not contain information about symmetry_detected")
    @skipForParser("Molcas", "reading point group symmetry and name is not implemented")
    @skipForParser("Molpro", "reading point group symmetry and name is not implemented")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("Turbomole", "reading point group symmetry and name is not implemented")
    @skipForParser("xTB", "not implemented yet")
    def testmetadata_symmetry_detected(self, data) -> None:
        """Does metadata have expected keys and values?"""
        assert data.metadata["symmetry_detected"] == "c2h"

    @skipForParser("FChk", "point group symmetry cannot be printed")
    @skipForParser("GAMESSDAT", "Files probably do not contain information about symmetry_used")
    @skipForParser("Molcas", "reading point group symmetry and name is not implemented")
    @skipForParser("Molpro", "reading point group symmetry and name is not implemented")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("Turbomole", "reading point group symmetry and name is not implemented")
    @skipForParser("xTB", "not implemented yet")
    def testmetadata_symmetry_used(self, data) -> None:
        """Does metadata have expected keys and values?"""
        assert data.metadata["symmetry_used"] == "c2h"

    @skipForParser("ADF", "reading cpu/wall time is not implemented for this parser")
    @skipForParser("DALTON", "reading cpu/wall time is not implemented for this parser")
    @skipForParser("FChk", "reading cpu/wall time is not implemented for this parser")
    @skipForParser("GAMESS", "reading cpu/wall time is not implemented for this parser")
    @skipForParser("GAMESSUK", "reading cpu/wall time is not implemented for this parser")
    @skipForParser("GAMESSUS", "reading cpu/wall time is not implemented for this parser")
    @skipForParser("GAMESSDAT", "reading cpu/wall time is not implemented for this parser")
    @skipForParser("Jaguar", "reading cpu/wall time is not implemented for this parser")
    @skipForParser("Molcas", " reading cpu/wall time is not implemented for this parser")
    @skipForParser("Molpro", "reading cpu/wall time is not implemented for this parser")
    @skipForParser("NBO", "attribute not implemented in this version")
    @skipForParser("NWChem", "reading cpu/wall time is not implemented for this parser")
    @skipForParser("Psi4", "reading cpu/wall time is not implemented for this parser")
    def testmetadata_times(self, data) -> None:
        """Does metadata have expected keys and values of correct types?"""
        if "wall_time" in data.metadata:
            assert data.metadata["wall_time"]
            assert all(
                isinstance(wall_time, datetime.timedelta)
                for wall_time in data.metadata["wall_time"]
            )
        if "cpu_time" in data.metadata:
            assert data.metadata["cpu_time"]
            assert all(
                isinstance(cpu_time, datetime.timedelta) for cpu_time in data.metadata["cpu_time"]
            )


class GenericHFSPTest(GenericSPTest):
    # Approximate HF energy of dvb after SCF in STO-3G (from DALTON 2015).
    scfenergy = -379.7689629312
    scfenergy_delta = 6.5e-1

    # Approximate energy of the innermost molecular orbital of DVB with
    # HF/STO-3G (from Psi4 1.3.1).
    moenergy = -11.0407466
    moenergy_delta = 1.6e-1


class ADFSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    # ADF only prints up to 0.1mD per atom, so the precision here is worse than 0.1mD.
    mass_precision = 0.3

    foverlap00 = 1.00003
    foverlap11 = 1.02672
    foverlap22 = 1.03585
    num_scf_criteria = 2
    # 2013.1/dvb_sp_b.adfout
    scfenergy = -5.162850967929650
    moenergy = -9.9079095713775

    def testfoverlaps(self, data) -> None:
        """Are the dims and values of the fragment orbital overlap matrix correct?"""

        assert data.fooverlaps.shape == (data.nbasis, data.nbasis)

        # The matrix is symmetric.
        row = data.fooverlaps[0, :]
        col = data.fooverlaps[:, 0]
        assert sum(col - row) == 0.0

        # Although the diagonal elements are close to zero, the SFOs
        # are generally not normalized, so test for a few specific values.
        assert abs(data.fooverlaps[0, 0] - self.foverlap00) < 0.0001
        assert abs(data.fooverlaps[1, 1] - self.foverlap11) < 0.0001
        assert abs(data.fooverlaps[2, 2] - self.foverlap22) < 0.0001


class GaussianSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    num_scf_criteria = 3


class JaguarSPTest(GenericSPTest):
    """Customized restricted single point KS unittest"""

    num_scf_criteria = 2


class JaguarHFSPTest(JaguarSPTest, GenericHFSPTest):
    """Customized restricted single point KS unittest"""


class Jaguar7SPTest(JaguarSPTest):
    """Customized restricted single point unittest"""

    # Jaguar prints only 10 virtual MOs by default. Can we re-run with full output?
    def testlengthmoenergies(self, data) -> None:
        """Is the number of evalues equal to the number of occ. MOs + 10?"""
        assert len(data.moenergies[0]) == data.homos[0] + 11


class MolcasSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    num_scf_criteria = 4


class MolproSPTest(GenericHFSPTest):
    """Customized restricted single point HF unittest"""

    num_scf_criteria = 2


class NWChemKSSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    num_scf_criteria = 3


class PsiSPTest(GenericSPTest):
    """Customized restricted single point KS unittest"""

    num_scf_criteria = 2


class PsiHFSPTest(PsiSPTest, GenericHFSPTest):
    """Customized restricted single point HF unittest"""


class OrcaSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    # Orca has different weights for the masses
    molecularmass = 130190

    moenergy_delta = 1.2e-1

    num_scf_criteria = 3


class OrcaHFSPTest(OrcaSPTest, GenericHFSPTest):
    """Customized restricted single point unittest"""

    def testmetadata_input_file(self, data) -> None:
        """Does metadata have expected keys and values?"""
        assert "input_file_contents" in data.metadata
        # TODO make input file names consistent where possible, though some
        # programs do not allow arbitrary file extensions; for example, DALTON
        # must end in `dal`.
        assert "dvb_sp_hf.in" in data.metadata["input_file_name"]


class NBOSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    def testpopulations(self, data) -> None:
        if hasattr(self, "populations"):
            population_key = "npa"

            expected_types = {
                "nao": int,
                "atom": str,
                "no": int,
                "lang": str,
                "type": str,
                "occupancy": float,
                "energy": float,
            }
            assert data.populations[population_key].keys == list(expected_types.keys())
            for key, exp_type in expected_types.items():
                assert isinstance(data.populations[population_key][key], list)
                assert isinstance(data.populations[population_key][key][0], exp_type)


class TurbomoleSPTest(GenericSPTest):
    """Customized restricted single point KS unittest"""

    num_scf_criteria = 2

    def testmetadata_basis_set(self, data) -> None:
        """Does metadata have expected keys and values?"""
        # One of our test cases used sto-3g hondo
        valid_basis = data.metadata["basis_set"].lower() in ("sto-3g", "sto-3g hondo")
        assert valid_basis


class TurbomoleHFSPTest(TurbomoleSPTest, GenericHFSPTest):
    """Customized restricted single point HF unittest"""


class XTBSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    scfenergy = -26.425939358406
    scfenergy_delta = 1.0e-6


class GenericDispersionTest:
    """Generic single-geometry dispersion correction unittest"""

    # Q-Chem 5.4
    dispersionenergy = -0.0147199319
    dispersionenergy_delta = 2.0e-7

    def testdispersionenergies(self, data) -> None:
        """Is the dispersion energy parsed correctly?"""
        assert len(data.dispersionenergies) == 1
        assert abs(
            data.dispersionenergies[0] - utils.convertor(self.dispersionenergy, "hartree", "eV")
        ) < utils.convertor(self.dispersionenergy_delta, "hartree", "eV")


class FireflyDispersionTest(GenericDispersionTest):
    """Customized single-geometry dispersion correction unittest"""

    # Firefly 8.1
    dispersionenergy = -0.015801551434377520


class SolventMetadataTest:
    """Check we can parse implicit solvent data."""

    model = ""
    # Toluene
    static_dielectric_constant = 2.3741

    def test_solvent_model(self, data) -> None:
        """Check solvent model was parsed correctly"""
        assert data.metadata["solvent_model"] == self.model

    @skipForLogfile(
        "basicQChem6.0/water_hf_solvent_smd_iefpcm.out",
        "the internally-used dielectric constant isn't printed, only solvent name",
    )
    def test_solvent_dielectric(self, data) -> None:
        """Check solvent dielectric was parsed correctly"""
        assert (
            abs(data.metadata["solvent_params"]["epsilon"] - self.static_dielectric_constant)
            < 1.0e-4
        )


class QChemSolventMetadataTest(SolventMetadataTest):
    static_dielectric_constant = 2.370


class IEFPCMMetadataTest(SolventMetadataTest):
    """Check we can parse implicit solvent data."""

    model = "IEFPCM"


class SCIPCMMetadataTest(SolventMetadataTest):
    """Check we can parse implicit solvent data."""

    model = "SCIPCM"


class IPCMMetadataTest(SolventMetadataTest):
    """Check we can parse implicit solvent data."""

    model = "IPCM"
    static_dielectric_constant = 78.3


class COSMOMetadataTest(SolventMetadataTest):
    """Check we can parse implicit solvent data."""

    model = "COSMO"


class CPCMMetadataTest(SolventMetadataTest):
    """Check we can parse implicit solvent data."""

    model = "CPCM"


class CPCMCOSMOMetadataTest(SolventMetadataTest):
    """Check we can parse implicit solvent data."""

    model = "CPCM-COSMO"


class SMDIEFPCMMetadataTest(SolventMetadataTest):
    """Check we can parse implicit solvent data."""

    model = "SMD-IEFPCM"


class SMDCPCMMetadataTest(SolventMetadataTest):
    """Check we can parse implicit solvent data."""

    model = "SMD-CPCM"


class QChemSMDIEFPCMMetadataTest(QChemSolventMetadataTest, SMDIEFPCMMetadataTest):
    """Check we can parse implicit solvent data."""


class QChemSMDCPCMMetadataTest(QChemSolventMetadataTest, SMDCPCMMetadataTest):
    """Check we can parse implicit solvent data."""


class GaussianPerformanceMetadataTest:
    """Check we can parse CPU/memory metadata."""

    def testmetadata_cpu(self, data) -> None:
        """Does metadata have the expected number of CPUs used?"""
        assert data.metadata["num_cpu"] == 1

    def testmetadata_memory_available(self, data) -> None:
        """Does metadata have the expected amount of memory?"""
        # 400 MB
        assert data.metadata["memory_available"] == 400000000

    def testmetadata_memory_used(self, data) -> None:
        """Does metadata have the expected amount of memory?"""
        assert data.metadata["memory_used"] == 52428800


class ORCAPerformanceMetadataTest:
    """Check we can parse CPU/memory metadata."""

    def testmetadata_cpu(self, data) -> None:
        """Does metadata have the expected number of CPUs used?"""
        assert data.metadata["num_cpu"] == 2

    def testmetadata_memory_available(self, data) -> None:
        """Does metadata have the expected amount of memory?"""
        # 400 MB
        assert data.metadata["memory_available"] == 1000000000

    def testmetadata_memory_used(self, data) -> None:
        """Does metadata have the expected amount of memory?"""
        assert data.metadata["memory_used"] == 463000000


class TurbomolePerformanceMetadataTest:
    """Check we can parse CPU/memory metadata."""

    def testmetadata_cpu(self, data) -> None:
        """Does metadata have the expected number of CPUs used?"""
        assert data.metadata["num_cpu"] == 1

    def testmetadata_memory_available(self, data) -> None:
        """Does metadata have the expected amount of memory?"""
        assert data.metadata["memory_available"] == 524288000

    @skipForParser("Turbomole", "memory used is not available for Turbomole")
    def testmetadata_memory_used(self, data) -> None:
        """Does metadata have the expected amount of memory?"""
        assert data.metadata["memory_used"] == 0
