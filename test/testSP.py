# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test single point logfiles in cclib"""

import sys
import numpy

import bettertest


class GenericSPTest(bettertest.TestCase):
    """Generic restricted single point unittest"""

    # In STO-3G, H has 1, C has 5 (1 S and 4 SP).
    nbasisdict = {1:1, 6:5}
    
    # Approximate B3LYP energy of dvb after SCF in STO-3G.
    b3lyp_energy = -10365

    def testnatom(self):
        """Is the number of atoms equal to 20?"""
        self.assertEquals(self.data.natom, 20)

    def testatomnos(self):
        """Are the atomnos correct?"""
        self.failUnless(numpy.alltrue([numpy.issubdtype(atomno,int) for atomno in self.data.atomnos]))
        # This will work only for numpy
        #self.assertEquals(self.data.atomnos.dtype.char, 'i')
        self.assertEquals(self.data.atomnos.shape, (20,) )
        self.assertEquals(sum(self.data.atomnos==6) + sum(self.data.atomnos==1), 20)

    def testatomcharges(self):
        """Are atomcharges (at least Mulliken) consistent with natom and sum to zero?"""
        for type in set(['mulliken'] + list(self.data.atomcharges.keys())):
            charges = self.data.atomcharges[type]
            self.assertEquals(len(charges), self.data.natom)
            self.assertInside(sum(charges), 0.0, 0.001)

    def testatomcoords(self):
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        expected_shape = (1, self.data.natom, 3)
        self.assertEquals(self.data.atomcoords.shape, expected_shape)
    
    def testcharge_and_mult(self):
        """Are the charge and multiplicity correct?"""
        self.assertEquals(self.data.charge, 0)
        self.assertEquals(self.data.mult, 1)

    def testnbasis(self):
        """Is the number of basis set functions correct?"""
        count = sum([self.nbasisdict[n] for n in self.data.atomnos])
        self.assertEquals(self.data.nbasis, count)

    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique?"""
        all = []
        for i,atom in enumerate(self.data.atombasis):
            self.assertEquals(len(atom), self.nbasisdict[self.data.atomnos[i]])
            all += atom
        # Test if there are as many indices as atomic orbitals.
        self.assertEquals(len(all), self.data.nbasis)
        # Check if all are different (every orbital indexed once).
        self.assertEquals(len(set(all)), len(all))

    def testcoreelectrons(self):
        """Are the coreelectrons all 0?"""
        ans = numpy.zeros(self.data.natom, 'i')
        self.assertArrayEquals(self.data.coreelectrons, ans)

    def testnormalisesym(self):
        """Did this subclass overwrite normalisesym?"""
        self.assertNotEquals(self.logfile.normalisesym("A"), "ERROR: This should be overwritten by this subclass")

    def testsymlabels(self):
        """Are all the symmetry labels either Ag/u or Bg/u?"""
        sumwronglabels = sum([x not in ['Ag','Bu','Au','Bg'] for x in self.data.mosyms[0]])
        self.assertEquals(sumwronglabels,0)

    def testhomos(self):
        """Is the index of the HOMO equal to 34?"""
        self.assertArrayEquals(self.data.homos, numpy.array([34],"i"),"%s != array([34],'i')" % numpy.array_repr(self.data.homos))

    def testscfvaluetype(self):
        """Are scfvalues and its elements the right type??"""
        self.assertEquals(type(self.data.scfvalues),type([]))
        self.assertEquals(type(self.data.scfvalues[0]),type(numpy.array([])))

    def testscfenergy(self):
        """Is the SCF energy within 40eV of target?"""
        self.assertInside(self.data.scfenergies[-1], self.b3lyp_energy, 40, "Final scf energy: %f not %i +- 40eV" %(self.data.scfenergies[-1], self.b3lyp_energy))

    def testscftargetdim(self):
        """Do the scf targets have the right dimensions?"""
        self.assertEquals(self.data.scftargets.shape, (len(self.data.scfvalues),len(self.data.scfvalues[0][0])))

    def testlengthmoenergies(self):
        """Is the number of evalues equal to nmo?"""
        self.assertEquals(len(self.data.moenergies[0]), self.data.nmo)

    def testtypemoenergies(self):
        """Is moenergies a list containing one numpy array?"""
        self.assertEquals(type(self.data.moenergies), type([]))
        self.assertEquals(type(self.data.moenergies[0]), type(numpy.array([])))

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        self.assertEquals(type(self.data.mocoeffs), type([]))
        self.assertEquals(len(self.data.mocoeffs), 1)
        self.assertEquals(self.data.mocoeffs[0].shape,
                          (self.data.nmo, self.data.nbasis))

    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis?"""
        # ADF has the attribute fooverlaps instead of aooverlaps.
        if not hasattr(self.data, "aooverlaps") and hasattr(self.data, "fooverlaps"):
            self.data.aooverlaps = self.data.fooverlaps
        self.assertEquals(self.data.aooverlaps.shape, (self.data.nbasis, self.data.nbasis))

    def testaooverlaps(self):
        """Are the first row and column of the overlap matrix identical?"""
        # ADF has the attribute fooverlaps instead of aooverlaps.
        if not hasattr(self.data, "aooverlaps") and hasattr(self.data, "fooverlaps"):
            self.data.aooverlaps = self.data.fooverlaps
        self.assertEquals(sum(self.data.aooverlaps[0,:] -
                              self.data.aooverlaps[:,0]),
                          0)

    def testoptdone(self):
        """There should be no optdone attribute set."""
        self.assertFalse(hasattr(self.data, 'optdone'))

    def testmoments(self):
        """Does the dipole and possible higher molecular moments look reasonable?"""

        # The reference point is always a vector, but not necessarily the
        # origin or center of mass. In this case, however, the center of mass
        # is at the origin, so we now what to expect.
        reference = self.data.moments[0]
        self.assertEquals(len(reference), 3)
        for x in reference:
            self.assertInside(x, 0.0, 0.001)

        # Length and value of dipole moment should always be correct (zero for this test).
        dipole = self.data.moments[1]
        self.assertEquals(len(dipole), 3)
        for d in dipole:
            self.assertInside(d, 0.0, 0.001)

        # If the quadrupole is there, we can expect roughly -50B for the XX moment,
        # -50B for the YY moment and and -60B for the ZZ moment.
        if len(self.data.moments) > 2:
            quadrupole = self.data.moments[2]
            self.assertEquals(len(quadrupole), 6)
            self.assertInside(quadrupole[0], -50, 5)
            self.assertInside(quadrupole[3], -50, 5)
            self.assertInside(quadrupole[5], -60, 5)

        # If the octupole is there, it should have 10 components and be zero.
        if len(self.data.moments) > 3:
            octupole = self.data.moments[3]
            self.assertEquals(len(octupole), 10)
            for m in octupole:
                self.assertInside(m, 0.0, 0.001)

        # Hexadecapole should have 15 elements and an XXXX componenet of around -1900 Debye*ang^2,
        # YYYY component of -330B and ZZZZ component of -50B.
        if len(self.data.moments) > 4:
            hexadecapole = self.data.moments[4]
            self.assertEquals(len(hexadecapole), 15)
            self.assertInside(hexadecapole[0], -1900, 100)
            self.assertInside(hexadecapole[10], -330, 10)
            self.assertInside(hexadecapole[14], -50, 5)

        # The are 21 unique 32-pole moments, and all are zero in this test case.
        if len(self.data.moments) > 5:
            moment32 = self.data.moments[5]
            self.assertEquals(len(moment32), 21)
            for m in moment32:
                self.assertInside(m, 0.0, 0.001)



class ADFSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    # ADF parser does not extract atombasis.
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique? PASS"""
        self.assertEquals(1, 1)
       
    def testscfenergy(self):
        """Is the SCF energy within 1eV of -140eV?"""
        self.assertInside(self.data.scfenergies[-1],-140,1,"Final scf energy: %f not -140+-1eV" % self.data.scfenergies[-1])


class GaussianSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    # Molecular mass of DVB in mD.
    molecularmass = 130078.25

    def testatommasses(self):
        """Do the atom masses sum up to the molecular mass (130078.25+-0.1mD)?"""
        mm = 1000*sum(self.data.atommasses)
        self.assertInside(mm, 130078.25, 0.1, "Molecule mass: %f not 130078 +- 0.1mD" %mm)

    def testmoments(self):
        """We don't yet have dipole moments printed in the unit test logfile. PASS"""


class Jaguar7SPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    # Data file does not contain enough information. Can we make a new one?
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique? PASS"""
        self.assertEquals(1, 1)

    # We did not print the atomic partial charges in the unit tests for this version.
    def testatomcharges(self):
        """Are all atomcharges consistent with natom and do they sum to zero? PASS"""
        self.assertEquals(1, 1)

    # Jaguar prints only 10 virtual MOs by default. Can we re-run with full output?
    def testlengthmoenergies(self):
        """Is the number of evalues equal to the number of occ. MOs + 10?"""
        self.assertEquals(len(self.data.moenergies[0]), self.data.homos[0]+11)

    # Data file does not contain enough information. Can we make a new one?
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""
        self.assertEquals(1, 1)

    def testmoments(self):
        """No dipole moments in the logfile. PASS"""

class JaguarSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    def testmoments(self):
        """No dipole moments in the logfile. PASS"""


class MolproSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    def testsymlabels(self):
        """Are all the symmetry labels either Ag/u or Bg/u? PASS"""
        self.assertEquals(1,1)


class MolproSPTest2006(MolproSPTest):
    """Customized restricted single point unittest"""

    # These tests were run a long time ago and since we don't have access
    # to Molpro 2006 anymore, we can skip this test (it is tested in 2012).
    def testatomcharges(self):
        """Are atomcharges (at least Mulliken) consistent with natom and sum to zero? PASS"""
        self.assertEquals(1,1)


class OrcaSPTest(GenericSPTest):
    """Customized restricted single point unittest"""
    
    # ORCA has no support for symmetry yet.
    def testsymlabels(self):
        """Are all the symmetry labels either Ag/u or Bg/u? PASS"""
        self.assertEquals(1,1)


class PsiSPTest(GenericSPTest):
    """Customized restricted single point HF/KS unittest"""

    # Psi does not currently have the option to print the overlap matrix.
    def testaooverlaps(self):
        """Are the first row and colm of the overlap matrix identical? PASS"""
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis? PASS"""


class Psi3SPTest(PsiSPTest):
    """Customized restricted single point HF/KS unittest"""

    # The final energy is also a bit higher here, I think due to the fact
    # that a SALC calculation is done instead of a full LCAO.
    b3lyp_energy = -10300

    # Psi3 did not print partial atomic charges.
    def testatomcharges(self):
        """Are atomcharges (at least Mulliken) consistent with natom and sum to zero? PASS"""

    # The molecular orbitals in Psi3 are printed within each irreducible representation,
    # but I don't know if that means there is no mixing between them (SALC calculation).
    # In any case, Psi4 prints a standard LCAO, with coefficients between all basis functions
    # and molecular orbitals, so we do not parse mocoeffs in Psi3 at all.
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""


class QChemSPTest(GenericSPTest):
    """Customized restricted single point unittest"""

    # Q-Chem cannot print the overlap matrix.
    def testdimaooverlaps(self):
        """Are the dims of the overlap matrix consistent with nbasis? PASS"""

    # Q-Chem cannot print the overlap matrix.
    def testaooverlaps(self):
        """Are the first row and column of the overlap matrix identical? PASS"""

    # `mocoeffs` not implemented yet.
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""

    # `atombasis` not implemented yet.
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique?"""

if __name__=="__main__":

    from testall import testall
    testall(modules=["SP"])
