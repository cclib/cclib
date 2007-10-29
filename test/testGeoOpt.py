__revision__ = "$Revision$"

import math

import numpy

import bettertest
import testSP


class GenericGeoOptTest(bettertest.TestCase):
    """Geometry optimization unittest."""

    # In STO-3G, H has 1, C has 3.
    nbasisdict = {1:1, 6:5}
    
    # Some programs print surplus atom coordinates by default.
    extracoords = 0
    
    # Some programs do surplus SCF cycles by default.
    extrascfs = 0

    # Approximate B3LYP energy of dvb after SCF in STO-3G.
    b3lyp_energy = -10365

    def testnatom(self):
        """Is the number of atoms equal to 20?"""
        self.assertEquals(self.data.natom, 20)

    def testatomnos(self):
        """Are the atomnos correct?"""
        self.failUnless(numpy.alltrue([isinstance(atomno,int) for atomno in self.data.atomnos]))
        # This will work only for numpy
        #self.assertEquals(self.data.atomnos.dtype.char, 'i')
        self.assertEquals(self.data.atomnos.shape, (20,) )
        self.assertEquals(sum(self.data.atomnos==6) + sum(self.data.atomnos==1), 20)        

    def testatomcoords(self):
        """Are atomcoords consistent with natom and Angstroms?"""
        coords = self.data.atomcoords
        self.assertEquals(self.data.natom,len(coords[0]),"natom is %d but len(atomcoords[0]) is %d" % (self.data.natom,len(coords[0])))

        # Find the minimum distance between two C atoms
        mindist = 999
        for i in range(self.data.natom-1):
            if self.data.atomnos[i]==6:
                for j in range(i+1,self.data.natom):
                    if self.data.atomnos[j]==6:
                        # Find the distance in the final iteration
                        dist = math.sqrt(sum((coords[-1][i]-coords[-1][j])**2))
                        mindist = min(mindist,dist)
        self.assert_(abs(mindist-1.34)<0.03,"Mindist is %f (not 1.34)" % mindist)

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
        """Are scfvalues and its elements the right type?"""
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

    def testgeovalues_atomcoords(self):
        """Are atomcoords consistent with geovalues?"""
        coords = self.data.atomcoords
        self.assertEquals(len(self.data.geovalues),len(coords)-self.extracoords,"len(atomcoords) is %d but len(geovalues) is %d" % (len(coords),len(self.data.geovalues)))
        
    def testgeovalues_scfvalues(self):
        """Are scfvalues consistent with geovalues?"""
        self.assertEquals(len(self.data.scfvalues)-self.extrascfs,len(self.data.geovalues))

    def testgeotargets(self):
        """Do the geo targets have the right dimensions?"""
        self.assertEquals(self.data.geotargets.shape,(len(self.data.geovalues[0]),))
    
class ADFGeoOptTest(GenericGeoOptTest):
    """ADF geometry optimization unittest."""

    extracoords = 1
    extrascfs = 1

    def testscfenergy(self):
        """Is the SCF energy within 1eV of -140eV?"""
        self.assertInside(self.data.scfenergies[-1],-140,1,"Final scf energy: %f not -140+-1eV" % self.data.scfenergies[-1])

class GamessUKGeoOptTest(GenericGeoOptTest):
    """GAMESS-UK geometry optimization unittest."""

    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x (homo+5) x nbasis?"""
        self.assertEquals(type(self.data.mocoeffs), type([]))
        self.assertEquals(len(self.data.mocoeffs), 1)
        self.assertEquals(self.data.mocoeffs[0].shape,
                          (self.data.homos[0]+1+5, self.data.nbasis))
        
class GamessUSGeoOptTest(GenericGeoOptTest):
    """GAMESS-US geometry optimization unittest."""

class GaussianGeoOptTest(GenericGeoOptTest):
    """Gaussian geometry optimization unittest."""

    # Data file does not contain enough information. Can we make a new one?
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique? PASS"""
        self.assertEquals(1, 1)
       
    # Data file does not contain enough information. Can we make a new one?
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""
        self.assertEquals(1, 1)

class JaguarGeoOptTest(GenericGeoOptTest):
    """Jaguar geometry optimization unittest."""

    # Data file does not contain enough information. Can we make a new one?
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique? PASS"""
        self.assertEquals(1, 1)
       
    # Data file does not contain enough information. Can we make a new one?
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""
        self.assertEquals(1, 1)

class MolproGeoOptTest(GenericGeoOptTest):
    """Molpro geometry optimization unittest."""
    
    extracoords = 1
    extrascfs = 2

    def testsymlabels(self):
        """Are all the symmetry labels either Ag/u or Bg/u? PASS"""
        self.assertEquals(1,1)

class OrcaGeoOptTest(GenericGeoOptTest):
    """ORCA geometry optimization unittest."""
    def testnbasis(self):
        """Is the number of basis set functions equal to 110 (3-21G)?"""
        self.assertEquals(self.data.nbasis,110)
    def testlengthmoenergies(self):
        """Is the number of evalues equal to 110 (3-21G)?"""
        self.assertEquals(110,len(self.data.moenergies[0]))

class PCGamessGeoOptTest(GenericGeoOptTest):
    """PC-GAMESS geometry optimization unittest."""


if __name__=="__main__":

    from testall import testall
    testall(modules=["GeoOpt"])
