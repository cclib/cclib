# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test single point logfiles in cclib."""

import os
import unittest

import numpy

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

    # Overlap first two atomic orbitals.
    overlap01 = 0.24

    # Generally, one criteria for SCF energy convergence. 
    num_scf_criteria = 1

    def testnatom(self):
        """Is the number of atoms equal to 20?"""
        self.assertEqual(self.data.natom, 20)

    def testatomnos(self):
        """Are the atomnos correct?"""

        # The nuclear charges should be integer values in a NumPy array.
        self.assertTrue(numpy.alltrue([numpy.issubdtype(atomno, numpy.signedinteger)
                                       for atomno in self.data.atomnos]))
        self.assertEqual(self.data.atomnos.dtype.char, 'i')

        self.assertEqual(self.data.atomnos.shape, (20,) )
        self.assertEqual(sum(self.data.atomnos == 6) + sum(self.data.atomnos == 1), 20)

    @skipForParser('DALTON', 'DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now')
    @skipForLogfile('Jaguar/basicJaguar7', 'We did not print the atomic partial charges in the unit tests for this version')
    @skipForLogfile('Molpro/basicMolpro2006', "These tests were run a long time ago and since we don't have access to Molpro 2006 anymore, we can skip this test (it is tested in 2012)")
    @skipForLogfile('Psi3/basicPsi3', 'Psi3 did not print partial atomic charges')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testatomcharges(self):
        """Are atomcharges (at least Mulliken) consistent with natom and sum to zero?"""
        for type in set(['mulliken'] + list(self.data.atomcharges.keys())):
            charges = self.data.atomcharges[type]
            self.assertEqual(len(charges), self.data.natom)
            self.assertAlmostEqual(sum(charges), 0.0, delta=0.001)

    def testatomcoords(self):
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        expected_shape = (1, self.data.natom, 3)
        self.assertEqual(self.data.atomcoords.shape, expected_shape)

    def testatomcoords_units(self):
        """Are atomcoords consistent with Angstroms?"""
        min_carbon_dist = get_minimum_carbon_separation(self.data)
        dev = abs(min_carbon_dist - 1.34)
        self.assertTrue(dev < 0.03, "Minimum carbon dist is %.2f (not 1.34)" % min_carbon_dist)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testcharge_and_mult(self):
        """Are the charge and multiplicity correct?"""
        self.assertEqual(self.data.charge, 0)
        self.assertEqual(self.data.mult, 1)

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testnbasis(self):
        """Is the number of basis set functions correct?"""
        count = sum([self.nbasisdict[n] for n in self.data.atomnos])
        self.assertEqual(self.data.nbasis, count)

    @skipForParser('ADF', 'ADF parser does not extract atombasis')
    @skipForLogfile('Jaguar/basicJaguar7', 'Data file does not contain enough information. Can we make a new one?')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique?"""
        all = []
        for i, atom in enumerate(self.data.atombasis):
            self.assertEqual(len(atom), self.nbasisdict[self.data.atomnos[i]])
            all += atom
        # Test if there are as many indices as atomic orbitals.
        self.assertEqual(len(all), self.data.nbasis)
        # Check if all are different (every orbital indexed once).
        self.assertEqual(len(set(all)), len(all))

    @skipForParser('GAMESS', 'atommasses not implemented yet')
    @skipForParser('GAMESSUK', 'atommasses not implemented yet')
    @skipForParser('Jaguar', 'atommasses not implemented yet')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Molpro', 'atommasses not implemented yet')
    @skipForParser('NWChem', 'atommasses not implemented yet')
    @skipForLogfile('Psi3/basicPsi3', 'atommasses not implemented yet')
    @skipForLogfile('Psi4/basicPsi4.0b5', 'atommasses not implemented yet')
    @skipForParser('QChem', 'atommasses not implemented yet')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testatommasses(self):
        """Do the atom masses sum up to the molecular mass?"""
        mm = 1000*sum(self.data.atommasses)
        msg = "Molecule mass: %f not %f +- %fmD" % (mm, self.molecularmass, self.mass_precision)
        self.assertAlmostEqual(mm, self.molecularmass, delta=self.mass_precision, msg=msg)

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testcoreelectrons(self):
        """Are the coreelectrons all 0?"""
        ans = numpy.zeros(self.data.natom, 'i')
        numpy.testing.assert_array_equal(self.data.coreelectrons, ans)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testnormalisesym(self):
        """Did this subclass overwrite normalisesym?"""
        # https://stackoverflow.com/a/8747890
        self.logfile.normalisesym("A")

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Molpro', '?')
    @skipForParser('ORCA', 'ORCA has no support for symmetry yet')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testsymlabels(self):
        """Are all the symmetry labels either Ag/u or Bg/u?"""
        sumwronglabels = sum([x not in ['Ag', 'Bu', 'Au', 'Bg'] for x in self.data.mosyms[0]])
        self.assertEqual(sumwronglabels, 0)

    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testhomos(self):
        """Is the index of the HOMO equal to 34?"""
        numpy.testing.assert_array_equal(self.data.homos, numpy.array([34],"i"), "%s != array([34],'i')" % numpy.array_repr(self.data.homos))

    def testscfvaluetype(self):
        """Are scfvalues and its elements the right type??"""
        self.assertEqual(type(self.data.scfvalues),type([]))
        self.assertEqual(type(self.data.scfvalues[0]),type(numpy.array([])))

    def testscfenergy(self):
        """Is the SCF energy within the target?"""
        self.assertAlmostEqual(self.data.scfenergies[-1], self.b3lyp_energy, delta=40, msg="Final scf energy: %f not %i +- 40eV" %(self.data.scfenergies[-1], self.b3lyp_energy))

    def testscftargetdim(self):
        """Do the scf targets have the right dimensions?"""
        self.assertEqual(self.data.scftargets.shape, (len(self.data.scfvalues), len(self.data.scfvalues[0][0])))

    def testscftargets(self):
        """Are correct number of SCF convergence criteria being parsed?"""
        self.assertEqual(len(self.data.scftargets[0]), self.num_scf_criteria)

    def testlengthmoenergies(self):
        """Is the number of evalues equal to nmo?"""
        if hasattr(self.data, "moenergies"):
            self.assertEqual(len(self.data.moenergies[0]), self.data.nmo)

    def testtypemoenergies(self):
        """Is moenergies a list containing one numpy array?"""
        if hasattr(self.data, "moenergies"):
            self.assertIsInstance(self.data.moenergies, list)
            self.assertIsInstance(self.data.moenergies[0], numpy.ndarray)

    @skipForParser('DALTON', 'mocoeffs not implemented yet')
    @skipForLogfile('Jaguar/basicJaguar7', 'Data file does not contain enough information. Can we make a new one?')
    @skipForLogfile('Psi3/basicPsi3', 'MO coefficients are printed separately for each SALC')
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        if hasattr(self.data, "mocoeffs"):
            self.assertIsInstance(self.data.mocoeffs, list)
            self.assertEqual(len(self.data.mocoeffs), 1)
            self.assertEqual(self.data.mocoeffs[0].shape,
                             (self.data.nmo, self.data.nbasis))
    
    @skipForParser('DALTON', 'mocoeffs not implemented yet')
    @skipForLogfile('Jaguar/basicJaguar7', 'Data file does not contain enough information. Can we make a new one?')
    @skipForLogfile('Psi3/basicPsi3', 'MO coefficients are printed separately for each SALC')
    def testfornoormo(self):
        """Do we have NOs or MOs?"""
        self.assertEquals(hasattr(self.data, "nocoeffs") or hasattr(self.data, "mocoeffs"), True)

    def testdimnoccnos(self):
        """Is the length of nooccnos equal to nmo?"""
        if hasattr(self.data, "nooccnos"):
            self.assertIsInstance(self.data.nooccnos, numpy.ndarray)
            self.assertEquals(len(self.data.nooccnos), self.data.nmo)

    def testdimnocoeffs(self):
        """Are the dimensions of nocoeffs equal to nmo x nmo?"""
        if hasattr(self.data, "nocoeffs"):
            self.assertIsInstance(self.data.nocoeffs, numpy.ndarray)
            self.assertEquals(self.data.nocoeffs.shape, (self.data.nmo, self.data.nmo))

    @skipForParser('DALTON', 'To print: **INTEGRALS\n.PROPRI')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Psi3', 'Psi3 does not currently have the option to print the overlap matrix')
    @skipForParser('Psi4', 'Psi4 does not currently have the option to print the overlap matrix')
    @skipForParser('QChem', 'QChem cannot print the overlap matrix')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testaooverlaps(self):
        """Are the dims and values of the overlap matrix correct?"""

        self.assertEqual(self.data.aooverlaps.shape, (self.data.nbasis, self.data.nbasis))

        # The matrix is symmetric.
        row = self.data.aooverlaps[0,:]
        col = self.data.aooverlaps[:,0]
        self.assertEqual(sum(col - row), 0.0)

        # All values on diagonal should be exactly zero.
        for i in range(self.data.nbasis):
            self.assertEqual(self.data.aooverlaps[i,i], 1.0)

        # Check some additional values that don't seem to move around between programs.
        self.assertAlmostEqual(self.data.aooverlaps[0, 1], self.overlap01, delta=0.01)
        self.assertAlmostEqual(self.data.aooverlaps[1, 0], self.overlap01, delta=0.01)
        self.assertEqual(self.data.aooverlaps[3,0], 0.0)
        self.assertEqual(self.data.aooverlaps[0,3], 0.0)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testoptdone(self):
        """There should be no optdone attribute set."""
        self.assertFalse(hasattr(self.data, 'optdone'))

    @skipForParser('Gaussian', 'Logfile needs to be updated')
    @skipForParser('Jaguar', 'No dipole moments in the logfile')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testmoments(self):
        """Does the dipole and possible higher molecular moments look reasonable?"""

        # The reference point is always a vector, but not necessarily the
        # origin or center of mass. In this case, however, the center of mass
        # is at the origin, so we now what to expect.
        reference = self.data.moments[0]
        self.assertEqual(len(reference), 3)
        for x in reference:
            self.assertEqual(x, 0.0)

        # Length and value of dipole moment should always be correct (zero for this test).
        dipole = self.data.moments[1]
        self.assertEqual(len(dipole), 3)
        for d in dipole:
            self.assertAlmostEqual(d, 0.0, places=7)

        # If the quadrupole is there, we can expect roughly -50B for the XX moment,
        # -50B for the YY moment and and -60B for the ZZ moment.
        if len(self.data.moments) > 2:
            quadrupole = self.data.moments[2]
            self.assertEqual(len(quadrupole), 6)
            self.assertAlmostEqual(quadrupole[0], -50, delta=2.5)
            self.assertAlmostEqual(quadrupole[3], -50, delta=2.5)
            self.assertAlmostEqual(quadrupole[5], -60, delta=3)

        # If the octupole is there, it should have 10 components and be zero.
        if len(self.data.moments) > 3:
            octupole = self.data.moments[3]
            self.assertEqual(len(octupole), 10)
            for m in octupole:
                self.assertAlmostEqual(m, 0.0, delta=0.001)

        # The hexadecapole should have 15 elements, an XXXX component of around -1900 Debye*ang^2,
        # a YYYY component of -330B and a ZZZZ component of -50B.
        if len(self.data.moments) > 4:
            hexadecapole = self.data.moments[4]
            self.assertEqual(len(hexadecapole), 15)
            self.assertAlmostEqual(hexadecapole[0], -1900, delta=90)
            self.assertAlmostEqual(hexadecapole[10], -330, delta=11)
            self.assertAlmostEqual(hexadecapole[14], -50, delta=2.5)

        # The are 21 unique 32-pole moments, and all are zero in this test case.
        if len(self.data.moments) > 5:
            moment32 = self.data.moments[5]
            self.assertEqual(len(moment32), 21)
            for m in moment32:
                self.assertEqual(m, 0.0)

    @skipForParser('ADF', 'Does not support metadata yet')
    @skipForParser('GAMESSUK', 'Does not support metadata yet')
    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Molpro', 'Does not support metadata yet')
    @skipForParser('NWChem', 'Does not support metadata yet')
    @skipForParser('ORCA', 'Does not support metadata yet')
    @skipForParser('Psi3', 'Does not support metadata yet')
    @skipForParser('Psi4', 'Does not support metadata yet')
    @skipForParser('QChem', 'Does not support metadata yet')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testmetadata(self):
        """Does metadata have expected keys and values?"""
        self.assertTrue(hasattr(self.data, "metadata"))
        if self.logfile.logname not in ['ORCA', 'Psi']:
            self.assertIn("basis_set", self.data.metadata)
        if self.logfile.logname == 'ORCA':
            self.assertIn("input_file_name", self.data.metadata)
            self.assertIn("input_file_contents", self.data.metadata)
        self.assertIn("methods", self.data.metadata)
        self.assertIn("package", self.data.metadata)
        self.assertIn("package_version", self.data.metadata)

class ADFSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    # ADF only prints up to 0.1mD per atom, so the precision here is worse than 0.1mD.
    mass_precision = 0.3

    foverlap00 = 1.00003
    foverlap11 = 1.02672
    foverlap22 = 1.03585
    num_scf_criteria = 2
    b3lyp_energy = -140

    def testfoverlaps(self):
        """Are the dims and values of the fragment orbital overlap matrix correct?"""

        self.assertEqual(self.data.fooverlaps.shape, (self.data.nbasis, self.data.nbasis))

        # The matrix is symmetric.
        row = self.data.fooverlaps[0,:]
        col = self.data.fooverlaps[:,0]
        self.assertEqual(sum(col - row), 0.0)

        # Although the diagonal elements are close to zero, the SFOs
        # are generally not normalized, so test for a few specific values.
        self.assertAlmostEqual(self.data.fooverlaps[0, 0], self.foverlap00, delta=0.0001)
        self.assertAlmostEqual(self.data.fooverlaps[1, 1], self.foverlap11, delta=0.0001)
        self.assertAlmostEqual(self.data.fooverlaps[2, 2], self.foverlap22, delta=0.0001)

class GaussianSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    num_scf_criteria = 3

class JaguarSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    num_scf_criteria = 2

class Jaguar7SPTest(JaguarSPTest):
    """Customized restricted single point unittest"""

    # Jaguar prints only 10 virtual MOs by default. Can we re-run with full output?
    def testlengthmoenergies(self):
        """Is the number of evalues equal to the number of occ. MOs + 10?"""
        self.assertEqual(len(self.data.moenergies[0]), self.data.homos[0]+11)

class MolcasSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    num_scf_criteria = 4

class MolproSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    num_scf_criteria = 2

class NWChemKSSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    num_scf_criteria = 3

class Psi3SPTest(GenericSPTest):
    """Customized restricted single point HF/KS unittest"""

    # The final energy is also a bit higher here, I think due to the fact
    # that a SALC calculation is done instead of a full LCAO.
    b3lyp_energy = -10300

class PsiSPTest(GenericSPTest):
    """Customized restricted single point HF/KS unittest"""

    num_scf_criteria = 2

class OrcaSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    # Orca has different weights for the masses
    molecularmass = 130190

    num_scf_criteria = 3

class TurbomoleSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    num_scf_criteria = 2


if __name__=="__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['SP'])
    suite.testall()
