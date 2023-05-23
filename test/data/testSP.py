# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test single point logfiles in cclib."""

import datetime
import os
import unittest

import numpy
import packaging

from common import get_minimum_carbon_separation

from skip import skipForParser
from skip import skipForLogfile


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericSPTest(unittest.TestCase):
    """Generic restricted single point unittest"""

    # Molecular mass of DVB in mD, and expected precision.
    molecularmass = 130078.25
    mass_precision = 0.10

    # In STO-3G, H has 1, C has 5 (1 S and 4 SP).
    nbasisdict = {1:1, 6:5}

    # Approximate B3LYP energy of dvb after SCF in STO-3G.
    b3lyp_energy = -10365

    # Approximate energy of the innermost molecular orbital of DVB with
    # B3LYP/STO-3G (from Q-Chem 5.4).
    b3lyp_moenergy = -272.60365543
    b3lyp_moenergy_delta = 7.55e-2

    # Overlap first two atomic orbitals.
    overlap01 = 0.24

    # Generally, one criteria for SCF energy convergence.
    num_scf_criteria = 1

    def testnatom(self):
        """Is the number of atoms equal to 20?"""
        assert self.data.natom == 20

    def testatomnos(self):
        """Are the atomnos correct?"""

        # The nuclear charges should be integer values in a NumPy array.
        assert numpy.alltrue([numpy.issubdtype(atomno, numpy.signedinteger)
                                       for atomno in self.data.atomnos])
        assert self.data.atomnos.dtype.char == 'i'

        assert self.data.atomnos.shape == (20,)
        assert sum(self.data.atomnos == 6) + sum(self.data.atomnos == 1) == 20

    @skipForParser('DALTON', 'DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now')
    @skipForParser('FChk', 'The parser is still being developed so we skip this test')
    @skipForLogfile('Jaguar/basicJaguar7', 'We did not print the atomic partial charges in the unit tests for this version')
    @skipForLogfile('Molpro/basicMolpro2006', "These tests were run a long time ago and since we don't have access to Molpro 2006 anymore, we can skip this test (it is tested in 2012)")
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testatomcharges(self):
        """Are atomic charges consistent with natom?"""
        for atomcharge_type in self.data.atomcharges:
            charges = self.data.atomcharges[atomcharge_type]
            natom = self.data.natom
            assert len(charges) == natom, f"len(atomcharges['{atomcharge_type}']) = {len(charges)}, natom = {natom}"

    @skipForParser('DALTON', 'DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now')
    @skipForParser('FChk', 'The parser is still being developed so we skip this test')
    @skipForLogfile('Jaguar/basicJaguar7', 'We did not print the atomic partial charges in the unit tests for this version')
    @skipForLogfile('Molpro/basicMolpro2006', "These tests were run a long time ago and since we don't have access to Molpro 2006 anymore, we can skip this test (it is tested in 2012)")
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testatomcharges_mulliken(self):
        """Do Mulliken atomic charges sum to zero?"""
        charges = self.data.atomcharges["mulliken"]
        assert abs(sum(charges)) < 1.0e-2

    @skipForParser('ADF', 'Lowdin charges not present by default')
    @skipForParser('DALTON', 'DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now')
    @skipForParser('FChk', 'The parser is still being developed so we skip this test')
    @skipForParser('Gaussian', 'Lowdin charges not present by default')
    @skipForParser('Jaguar', 'Lowdin charges not present by default')
    @skipForParser('NWChem', 'Lowdin charges not present by default')
    @skipForParser('Molcas', 'Lowdin charges not present by default')
    @skipForParser('Molpro', 'Lowdin charges not present by default')
    @skipForParser('QChem', 'Lowdin charges not present by default')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testatomcharges_lowdin(self):
        """Do Lowdin atomic charges sum to zero?"""
        charges = self.data.atomcharges["lowdin"]
        assert abs(sum(charges)) < 1.0e-2

    @skipForParser('ADF', 'Hirshfeld charges not implemented')
    @skipForParser('DALTON', 'DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now')
    @skipForParser('FChk', 'The parser is still being developed so we skip this test')
    @skipForParser('Gaussian', 'Hirshfeld charges not implemented')
    @skipForParser('GAMESS', 'Hirshfeld charges not implemented')
    @skipForParser('GAMESSUK', 'Hirshfeld charges not implemented')
    @skipForParser('Jaguar', 'Hirshfeld charges not implemented')
    @skipForParser('NWChem', 'Hirshfeld charges not implemented')
    @skipForParser('Molcas', 'Hirshfeld charges not implemented')
    @skipForParser('Molpro', 'Hirshfeld charges not implemented')
    @skipForLogfile('ORCA/basicORCA4.1', 'This needs to be moved to regressions')
    @skipForParser('Psi4', 'Hirshfeld charges not implemented')
    @skipForParser('QChem', 'Hirshfeld charges not implemented')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testatomcharges_hirshfeld(self):
        """Do Hirshfeld atomic charges sum to roughly zero?"""
        charges = self.data.atomcharges["hirshfeld"]
        assert abs(sum(charges)) < 5.0e-2

    def testatomcoords(self):
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        expected_shape = (1, self.data.natom, 3)
        assert self.data.atomcoords.shape == expected_shape

    def testatomcoords_units(self):
        """Are atomcoords consistent with Angstroms?"""
        min_carbon_dist = get_minimum_carbon_separation(self.data)
        dev = abs(min_carbon_dist - 1.34)
        assert dev < 0.03, f"Minimum carbon dist is {min_carbon_dist:.2f} (not 1.34)"

    @skipForParser('Molcas', 'missing mult')
    def testcharge_and_mult(self):
        """Are the charge and multiplicity correct?"""
        assert self.data.charge == 0
        assert self.data.mult == 1

    def testnbasis(self):
        """Is the number of basis set functions correct?"""
        count = sum([self.nbasisdict[n] for n in self.data.atomnos])
        assert self.data.nbasis == count

    @skipForParser('ADF', 'ADF parser does not extract atombasis')
    @skipForLogfile('Jaguar/basicJaguar7', 'Data file does not contain enough information. Can we make a new one?')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique?"""
        all = []
        for i, atom in enumerate(self.data.atombasis):
            assert len(atom) == self.nbasisdict[self.data.atomnos[i]]
            all += atom
        # Test if there are as many indices as atomic orbitals.
        assert len(all) == self.data.nbasis
        # Check if all are different (every orbital indexed once).
        assert len(set(all)) == len(all)

    @skipForParser('FChk', 'Formatted checkpoint files do not have a section for atommasses')
    @skipForParser('GAMESS', 'atommasses not implemented yet')
    @skipForParser('GAMESSUK', 'atommasses not implemented yet')
    @skipForParser('Jaguar', 'atommasses not implemented yet')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Molpro', 'atommasses not implemented yet')
    @skipForLogfile('Psi4/basicPsi4.0b5', 'atommasses not implemented yet')
    @skipForParser('QChem', 'atommasses not implemented yet')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testatommasses(self):
        """Do the atom masses sum up to the molecular mass?"""
        mm = 1000*sum(self.data.atommasses)
        msg = f"Molecule mass: {mm:f} not {self.molecularmass:f} +- {self.mass_precision:f}mD"
        assert abs(mm-self.molecularmass) < self.mass_precision, msg

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testcoreelectrons(self):
        """Are the coreelectrons all 0?"""
        ans = numpy.zeros(self.data.natom, 'i')
        numpy.testing.assert_array_equal(self.data.coreelectrons, ans)

    @skipForParser('FChk', 'Formatted checkpoint files do not have a section for symmetry')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Molpro', '?')
    def testsymlabels(self):
        """Are all the symmetry labels either Ag/u or Bg/u?"""
        sumwronglabels = sum([x not in ['Ag', 'Bu', 'Au', 'Bg'] for x in self.data.mosyms[0]])
        assert sumwronglabels == 0

    def testhomos(self):
        """Is the index of the HOMO equal to 34?"""
        numpy.testing.assert_array_equal(
            self.data.homos,
            numpy.array([34], "i"),
            f"{numpy.array_repr(self.data.homos)} != array([34],'i')",
        )

    @skipForParser('FChk', 'Formatted Checkpoint files do not have a section for SCF energy')
    def testscfvaluetype(self):
        """Are scfvalues and its elements the right type??"""
        assert isinstance(self.data.scfvalues, list)
        assert isinstance(self.data.scfvalues[0], numpy.ndarray)

    @skipForParser('FChk', 'Formatted Checkpoint files do not have a section for SCF energy')
    def testscfenergy(self):
        """Is the SCF energy within the target?"""
        assert abs(self.data.scfenergies[-1]-self.b3lyp_energy) < 40

    @skipForParser('FChk', 'Formatted Checkpoint files do not have a section for SCF convergence')
    def testscftargetdim(self):
        """Do the scf targets have the right dimensions?"""
        assert self.data.scftargets.shape == (len(self.data.scfvalues), len(self.data.scfvalues[0][0]))

    @skipForParser('FChk', 'Formatted Checkpoint files do not have a section for SCF convergence')
    def testscftargets(self):
        """Are correct number of SCF convergence criteria being parsed?"""
        assert len(self.data.scftargets[0]) == self.num_scf_criteria

    def testlengthmoenergies(self):
        """Is the number of evalues equal to nmo?"""
        if hasattr(self.data, "moenergies"):
            assert len(self.data.moenergies[0]) == self.data.nmo

    def testtypemoenergies(self):
        """Is moenergies a list containing one numpy array?"""
        if hasattr(self.data, "moenergies"):
            assert isinstance(self.data.moenergies, list)
            assert isinstance(self.data.moenergies[0], numpy.ndarray)

    @skipForLogfile('Gaussian/basicGaussian16/dvb_sp_no.out', 'no energies for natural orbitals')
    @skipForLogfile('Turbomole/basicTurbomole5.9/dvb_sp_symm', 'delta of 7.4, everything else ok')
    def testfirstmoenergy(self):
        """Is the lowest energy molecular orbital within the target?"""
        assert abs(self.data.moenergies[0][0]-self.b3lyp_moenergy) < self.b3lyp_moenergy_delta

    @skipForParser('DALTON', 'mocoeffs not implemented yet')
    @skipForLogfile('Jaguar/basicJaguar7', 'Data file does not contain enough information. Can we make a new one?')
    @skipForParser('Turbomole', 'Use of symmetry has reduced the number of mo coeffs')
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        if hasattr(self.data, "mocoeffs"):
            assert isinstance(self.data.mocoeffs, list)
            assert len(self.data.mocoeffs) == 1
            assert self.data.mocoeffs[0].shape == (self.data.nmo, self.data.nbasis)

    @skipForParser('DALTON', 'mocoeffs not implemented yet')
    @skipForLogfile('Jaguar/basicJaguar7', 'Data file does not contain enough information. Can we make a new one?')
    def testfornoormo(self):
        """Do we have NOs or MOs?"""
        assert hasattr(self.data, "nocoeffs") or hasattr(self.data, "mocoeffs")

    def testdimnoccnos(self):
        """Is the length of nooccnos equal to nmo?"""
        if hasattr(self.data, "nooccnos"):
            assert isinstance(self.data.nooccnos, numpy.ndarray)
            assert len(self.data.nooccnos) == self.data.nmo

    def testdimnocoeffs(self):
        """Are the dimensions of nocoeffs equal to nmo x nmo?"""
        if hasattr(self.data, "nocoeffs"):
            assert isinstance(self.data.nocoeffs, numpy.ndarray)
            assert self.data.nocoeffs.shape == (self.data.nmo, self.data.nmo)

    @skipForParser('DALTON', 'To print: **INTEGRALS\n.PROPRI')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Psi4', 'Psi4 does not currently have the option to print the overlap matrix')
    @skipForParser('QChem', 'QChem cannot print the overlap matrix')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testaooverlaps(self):
        """Are the dims and values of the overlap matrix correct?"""

        assert self.data.aooverlaps.shape == (self.data.nbasis, self.data.nbasis)

        # The matrix is symmetric.
        row = self.data.aooverlaps[0,:]
        col = self.data.aooverlaps[:,0]
        assert sum(col - row) == 0.0

        # All values on diagonal should be exactly one.
        for i in range(self.data.nbasis):
            assert self.data.aooverlaps[i,i] == 1.0

        # Check some additional values that don't seem to move around between programs.
        assert abs(self.data.aooverlaps[0, 1] - self.overlap01) < 0.01
        assert abs(self.data.aooverlaps[1, 0] - self.overlap01) < 0.01
        assert round(abs(self.data.aooverlaps[3, 0]), 7) == 0
        assert round(abs(self.data.aooverlaps[0, 3]), 7) == 0

    def testoptdone(self):
        """There should be no optdone attribute set."""
        assert not hasattr(self.data, 'optdone')

    @skipForParser('ADF', 'Not implemented yes')
    @skipForParser('DALTON', 'Not implemented yes')
    @skipForParser('FChk', 'Rotational constants are never written to fchk files')
    @skipForParser('GAMESS', 'Not implemented yes')
    @skipForParser('GAMESSUK', 'Not implemented yet')
    @skipForParser('Jaguar', 'Not implemented yet')
    @skipForParser('Molcas', 'Not implemented yes')
    @skipForParser('Molpro', 'Not implemented yes')
    @skipForParser('NWChem', 'Not implemented yes')
    @skipForParser('ORCA', 'Not implemented yes')
    @skipForParser('Psi4', 'Not implemented yes')
    @skipForParser('QChem', 'Not implemented yes')
    @skipForParser('Turbomole', 'Not implemented yes')
    def testrotconsts(self):
        """A single geometry leads to single set of rotational constants."""
        assert self.data.rotconsts.shape == (1, 3)
        # taken from Gaussian16/dvb_sp.out
        ref = [4.6266363, 0.6849065, 0.5965900]
        numpy.testing.assert_allclose(self.data.rotconsts[0], ref, rtol=0, atol=1.0e-3)

    @skipForParser('FChk', 'The parser is still being developed so we skip this test')
    @skipForParser('Gaussian', 'Logfile needs to be updated')
    @skipForParser('Jaguar', 'No dipole moments in the logfile')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    def testmoments(self):
        """Does the dipole and possible higher molecular moments look reasonable?"""

        # The reference point is always a vector, but not necessarily the
        # origin or center of mass. In this case, however, the center of mass
        # is at the origin, so we now what to expect.
        reference = self.data.moments[0]
        assert len(reference) == 3
        for x in reference:
            assert x == 0.0

        # Length and value of dipole moment should always be correct (zero for this test).
        dipole = self.data.moments[1]
        assert len(dipole) == 3
        for d in dipole:
            assert round(abs(d), 7) == 0

        # If the quadrupole is there, we can expect roughly -50B for the XX moment,
        # -50B for the YY moment and and -60B for the ZZ moment.
        if len(self.data.moments) > 2:
            quadrupole = self.data.moments[2]
            assert len(quadrupole) == 6
            assert abs(quadrupole[0] - -50) < 2.5
            assert abs(quadrupole[3] - -50) < 2.5
            assert abs(quadrupole[5] - -60) < 3

        # If the octupole is there, it should have 10 components and be zero.
        if len(self.data.moments) > 3:
            octupole = self.data.moments[3]
            assert len(octupole) == 10
            for m in octupole:
                assert abs(m) < 0.001

        # The hexadecapole should have 15 elements, an XXXX component of around -1900 Debye*ang^2,
        # a YYYY component of -330B and a ZZZZ component of -50B.
        if len(self.data.moments) > 4:
            hexadecapole = self.data.moments[4]
            assert len(hexadecapole) == 15
            assert abs(hexadecapole[0] - -1900) < 90
            assert abs(hexadecapole[10] - -330) < 11
            assert abs(hexadecapole[14] - -50) < 2.5

        # The are 21 unique 32-pole moments, and all are zero in this test case.
        if len(self.data.moments) > 5:
            moment32 = self.data.moments[5]
            assert len(moment32) == 21
            for m in moment32:
                assert m == 0.0

    @skipForParser('ADF', 'reading basis set names is not implemented')
    @skipForParser('GAMESSUK', 'reading basis set names is not implemented')
    @skipForParser('Molcas', 'reading basis set names is not implemented')
    @skipForParser('Psi4', 'reading basis set names is not implemented')
    def testmetadata_basis_set(self):
        """Does metadata have expected keys and values?"""
        assert self.data.metadata["basis_set"].lower() == "sto-3g"

    @skipForParser('ADF', 'reading input file contents and name is not implemented')
    @skipForParser('DALTON', 'reading input file contents and name is not implemented')
    @skipForParser('FChk', 'Formatted checkpoint files do not have an input file section')
    @skipForParser('GAMESS', 'reading input file contents and name is not implemented')
    @skipForParser('GAMESSUK', 'reading input file contents and name is not implemented')
    @skipForParser('Gaussian', 'reading input file contents and name is not implemented')
    @skipForParser('Jaguar', 'reading input file contents and name is not implemented')
    @skipForParser('Molcas', 'reading input file contents and name is not implemented')
    @skipForParser('Molpro', 'reading input file contents and name is not implemented')
    @skipForParser('NWChem', 'reading input file contents and name is not implemented')
    @skipForParser('Psi4', 'reading input file contents and name is not implemented')
    @skipForParser('QChem', 'reading input file contents and name is not implemented')
    @skipForParser('Turbomole', 'reading input file contents and name is not implemented')
    def testmetadata_input_file(self):
        """Does metadata have expected keys and values?"""
        assert "input_file_contents" in self.data.metadata
        # TODO make input file names consistent where possible, though some
        # programs do not allow arbitrary file extensions; for example, DALTON
        # must end in `dal`.
        assert "dvb_sp.in" in self.data.metadata["input_file_name"]

    def testmetadata_methods(self):
        """Does metadata have expected keys and values?"""
        # TODO implement and unify across parsers; current values are [],
        # ["HF"], ["RHF"], and ["DFT"]
        assert "methods" in self.data.metadata

    def testmetadata_package(self):
        """Does metadata have expected keys and values?"""
        # TODO How can the value be tested when the package name comes from
        # the parser and isn't stored on ccData?
        assert "package" in self.data.metadata

    @skipForParser('FChk', 'Formatted Checkpoint files do not have section for legacy package version')
    def testmetadata_legacy_package_version(self):
        """Does metadata have expected keys and values?"""
        # TODO Test specific values for each unit test.
        assert "legacy_package_version" in self.data.metadata

    @skipForParser('FChk', 'Formatted Checkpoint files do not have section for package version')
    def testmetadata_package_version(self):
        """Does metadata have expected keys and values?"""
        # TODO Test specific values for each unit test.
        assert isinstance(packaging.version.parse(self.data.metadata["package_version"]), packaging.version.Version)

    @skipForParser('FChk', 'point group symmetry cannot be printed')
    @skipForParser('Molcas', 'reading point group symmetry and name is not implemented')
    @skipForParser('Molpro', 'reading point group symmetry and name is not implemented')
    @skipForParser('Turbomole', 'reading point group symmetry and name is not implemented')
    def testmetadata_symmetry_detected(self):
        """Does metadata have expected keys and values?"""
        assert self.data.metadata["symmetry_detected"] == "c2h"

    @skipForParser('FChk', 'point group symmetry cannot be printed')
    @skipForParser('Molcas', 'reading point group symmetry and name is not implemented')
    @skipForParser('Molpro', 'reading point group symmetry and name is not implemented')
    @skipForParser('Turbomole', 'reading point group symmetry and name is not implemented')
    def testmetadata_symmetry_used(self):
        """Does metadata have expected keys and values?"""
        assert self.data.metadata["symmetry_used"] == "c2h"

    @skipForParser('ADF', 'reading cpu/wall time is not implemented for this parser')
    @skipForParser('DALTON', 'reading cpu/wall time is not implemented for this parser') 
    @skipForParser('FChk', 'reading cpu/wall time is not implemented for this parser') 
    @skipForParser('GAMESS', 'reading cpu/wall time is not implemented for this parser') 
    @skipForParser('GAMESSUK', 'reading cpu/wall time is not implemented for this parser') 
    @skipForParser('GAMESSUS', 'reading cpu/wall time is not implemented for this parser') 
    @skipForParser('Jaguar', 'reading cpu/wall time is not implemented for this parser') 
    @skipForParser('Molcas', ' reading cpu/wall time is not implemented for this parser') 
    @skipForParser('Molpro', 'reading cpu/wall time is not implemented for this parser') 
    @skipForParser('NWChem', 'reading cpu/wall time is not implemented for this parser')  
    @skipForParser('Psi3', 'reading cpu/wall time is not implemented for this parser') 
    @skipForParser('Psi4', 'reading cpu/wall time is not implemented for this parser')
    def testmetadata_times(self):
        """Does metadata have expected keys and values of correct types?"""
        if "wall_time" in self.data.metadata:
            assert self.data.metadata["wall_time"]
            assert all(isinstance(wall_time, datetime.timedelta)
                       for wall_time in self.data.metadata["wall_time"])
        if "cpu_time" in self.data.metadata:
            assert self.data.metadata["cpu_time"]
            assert all(isinstance(cpu_time, datetime.timedelta)
                       for cpu_time in self.data.metadata["cpu_time"])


class GenericHFSPTest(GenericSPTest):

    # Approximate HF energy of dvb after SCF in STO-3G (from DALTON 2015).
    hf_scfenergy = -10334.03948035995

    # Approximate energy of the innermost molecular orbital of DVB with
    # HF/STO-3G (from Psi4 1.3.1).
    hf_moenergy = -300.43401785663235

    @skipForParser('FChk', 'Formatted Checkpoint files do not have a section for SCF energy')
    def testscfenergy(self):
        """Is the SCF energy within the target?"""
        assert abs(self.data.scfenergies[-1]-self.hf_scfenergy) < 6.5e-1

    def testfirstmoenergy(self):
        """Is the lowest energy molecular orbital within the target?"""
        assert abs(self.data.moenergies[0][0]-self.hf_moenergy) < 1.6e-1


class ADFSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    # ADF only prints up to 0.1mD per atom, so the precision here is worse than 0.1mD.
    mass_precision = 0.3

    foverlap00 = 1.00003
    foverlap11 = 1.02672
    foverlap22 = 1.03585
    num_scf_criteria = 2
    b3lyp_energy = -140
    b3lyp_moenergy = -269.6079423873336

    def testfoverlaps(self):
        """Are the dims and values of the fragment orbital overlap matrix correct?"""

        assert self.data.fooverlaps.shape == (self.data.nbasis, self.data.nbasis)

        # The matrix is symmetric.
        row = self.data.fooverlaps[0,:]
        col = self.data.fooverlaps[:,0]
        assert sum(col - row) == 0.0

        # Although the diagonal elements are close to zero, the SFOs
        # are generally not normalized, so test for a few specific values.
        assert abs(self.data.fooverlaps[0, 0]-self.foverlap00) < 0.0001
        assert abs(self.data.fooverlaps[1, 1]-self.foverlap11) < 0.0001
        assert abs(self.data.fooverlaps[2, 2]-self.foverlap22) < 0.0001


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
    def testlengthmoenergies(self):
        """Is the number of evalues equal to the number of occ. MOs + 10?"""
        assert len(self.data.moenergies[0]) == self.data.homos[0]+11

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

    b3lyp_moenergy_delta = 1.2e-1

    num_scf_criteria = 3


class TurbomoleSPTest(GenericSPTest):
    """Customized restricted single point KS unittest"""

    num_scf_criteria = 2
    
    def testmetadata_basis_set(self):
        """Does metadata have expected keys and values?"""
        # One of our test cases used sto-3g hondo
        valid_basis = self.data.metadata["basis_set"].lower() in ("sto-3g", "sto-3g hondo")
        assert valid_basis


class TurbomoleHFSPTest(TurbomoleSPTest, GenericHFSPTest):
    """Customized restricted single point HF unittest"""


class GenericDispersionTest(unittest.TestCase):
    """Generic single-geometry dispersion correction unittest"""

    dispersionenergy = -0.4005496

    def testdispersionenergies(self):
        """Is the dispersion energy parsed correctly?"""
        assert len(self.data.dispersionenergies) == 1
        assert abs(self.data.dispersionenergies[0]-self.dispersionenergy) < 2.0e-7


class FireflyDispersionTest(GenericDispersionTest):
    """Customized single-geometry dispersion correction unittest"""
    dispersionenergy = -0.4299821


if __name__ == "__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['SP'])
    suite.testall()
