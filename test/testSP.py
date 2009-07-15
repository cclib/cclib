__revision__ = "$Revision$"

import numpy

import bettertest


class GenericSPTest(bettertest.TestCase):
    """Restricted single point unittest."""

    # In STO-3G, H has 1, C has 3.
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

    def testatomcoords(self):
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        self.assertEquals(self.data.atomcoords.shape,(1,self.data.natom,3))
    
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
        self.assertNotEquals(self.logfile.normalisesym("A"),"ERROR: This should be overwritten by this subclass")

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
        self.assertEquals(self.data.scftargets.shape,(len(self.data.scfvalues),len(self.data.scfvalues[0][0])))

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
        self.assertEquals(self.data.aooverlaps.shape,(self.data.nbasis,self.data.nbasis))

    def testaooverlaps(self):
        """Are the first row and colm of the overlap matrix identical?"""
        # ADF has the attribute fooverlaps instead of aooverlaps.
        if not hasattr(self.data, "aooverlaps") and hasattr(self.data, "fooverlaps"):
            self.data.aooverlaps = self.data.fooverlaps
        self.assertEquals(sum(self.data.aooverlaps[0,:] -
                              self.data.aooverlaps[:,0]),
                          0)

class ADFSPTest(GenericSPTest):
    """ADF restricted single point unittest."""

    # ADF parser does not extract atombasis.
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique? PASS"""
        self.assertEquals(1, 1)
       
    def testscfenergy(self):
        """Is the SCF energy within 1eV of -140eV?"""
        self.assertInside(self.data.scfenergies[-1],-140,1,"Final scf energy: %f not -140+-1eV" % self.data.scfenergies[-1])

class GamessUKSPTest(GenericSPTest):
    """GAMESS-UK restricted single point unittest."""

class GamessUSSPTest(GenericSPTest):
    """GAMESS-US restricted single point unittest."""

class GaussianSPTest(GenericSPTest):
    """Gaussian restricted single point unittest."""

class JaguarSPTest(GenericSPTest):
    """Jaguar restricted single point unittest."""

    # Data file does not contain enough information. Can we make a new one?
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique? PASS"""
        self.assertEquals(1, 1)
       
    # Jaguar prints only 10 virtual MOs by default. Can we re-run with full output?
    def testlengthmoenergies(self):
        """Is the number of evalues equal to the number of occ. MOs + 10?"""
        self.assertEquals(len(self.data.moenergies[0]), self.data.homos[0]+11)

    # Data file does not contain enough information. Can we make a new one?
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""
        self.assertEquals(1, 1)

class MolproSPTest(GenericSPTest):
    """Molpro restricted single point unittest."""

    def testsymlabels(self):
        """Are all the symmetry labels either Ag/u or Bg/u? PASS"""
        self.assertEquals(1,1)

class OrcaSPTest(GenericSPTest):
    """ORCA restricted single point unittest."""

    # This was run in 3-21G; H has 2, C has 9.
    nbasisdict = {1:2, 6:9}
    
    # Approximate B3LYP energy of dvb after SCF in 3-21G.
    b3lyp_energy = -10470

    # ORCA has no support for symmetry yet.
    def testsymlabels(self):
        """Are all the symmetry labels either Ag/u or Bg/u? PASS"""
        self.assertEquals(1,1)

class PCGamessSPTest(GenericSPTest):
    """PC-GAMESS restricted single point unittest."""


if __name__=="__main__":

    from testall import testall
    testall(modules=["SP"])
