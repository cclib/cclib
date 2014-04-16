# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

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
        # This will work only for numpy
        #self.assertEquals(self.data.atomnos.dtype.char, 'i')

        atomnos_types = [numpy.issubdtype(atomno,int) for atomno in self.data.atomnos]
        self.failUnless(numpy.alltrue(atomnos_types))

        self.assertEquals(self.data.atomnos.shape, (20,) )

        count_C = sum(self.data.atomnos == 6)
        count_H = sum(self.data.atomnos == 1)
        self.assertEquals(count_C + count_H, 20)        

    def testatomcharges(self):
        """Are all atomcharges consistent with natom and do they sum to zero?"""
        for type,charges in self.data.atomcharges.items():
            self.assertEquals(self.data.natom, len(charges))
            self.assertInside(sum(charges), 0, 0.001)

    def testatomcoords(self):
        """Are atomcoords consistent with natom and Angstroms?"""
        natom = len(self.data.atomcoords[0])
        ref = self.data.natom
        msg = "natom is %d but len(atomcoords[0]) is %d" % (ref, natom)
        self.assertEquals(natom, ref, msg)

        # Find the minimum distance between two C atoms.
        mindist = 999
        for i in range(self.data.natom-1):
            if self.data.atomnos[i]==6:
                for j in range(i+1,self.data.natom):
                    if self.data.atomnos[j]==6:
                        # Find the distance in the final iteration
                        final_x = self.data.atomcoords[-1][i]
                        final_y = self.data.atomcoords[-1][j]
                        dist = math.sqrt(sum((final_x - final_y)**2))
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
        ref = numpy.array([34], "i")
        msg = "%s != array([34], 'i')" % numpy.array_repr(self.data.homos)
        self.assertArrayEquals(self.data.homos, ref, msg)

    def testscfvaluetype(self):
        """Are scfvalues and its elements the right type?"""
        self.assertEquals(type(self.data.scfvalues),type([]))
        self.assertEquals(type(self.data.scfvalues[0]),type(numpy.array([])))

    def testscfenergy(self):
        """Is the SCF energy within 40eV of target?"""
        scf = self.data.scfenergies[-1]
        ref = self.b3lyp_energy
        msg = "Final scf energy: %f not %i +- 40eV" %(scf, ref)
        self.assertInside(scf, ref, 40, msg)

    def testscfenergydim(self):
        """Is the number of SCF energies consistent with atomcoords?"""
        count_scfenergies = self.data.scfenergies.shape[0] - self.extrascfs
        count_atomcoords = self.data.atomcoords.shape[0] - self.extracoords
        self.assertEquals(count_scfenergies, count_atomcoords)

    def testscftargetdim(self):
        """Do the scf targets have the right dimensions?"""
        dim_scftargets = self.data.scftargets.shape
        dim_scfvalues = (len(self.data.scfvalues),len(self.data.scfvalues[0][0]))
        self.assertEquals(dim_scftargets, dim_scfvalues)

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
        count_geovalues = len(self.data.geovalues)
        count_coords = len(self.data.atomcoords) - self.extracoords
        msg = "len(atomcoords) is %d but len(geovalues) is %d" % (count_coords, count_geovalues)
        self.assertEquals(count_geovalues, count_coords, msg)
        
    def testgeovalues_scfvalues(self):
        """Are scfvalues consistent with geovalues?"""
        count_scfvalues = len(self.data.scfvalues) - self.extrascfs
        count_geovalues = len(self.data.geovalues)
        self.assertEquals(count_scfvalues, count_geovalues)

    def testgeotargets(self):
        """Do the geo targets have the right dimensions?"""
        dim_geotargets = self.data.geotargets.shape
        dim_geovalues = (len(self.data.geovalues[0]), )
        self.assertEquals(dim_geotargets, dim_geovalues)

    def testoptdone(self):
        """Has the geometry converged and set optdone to True?"""
        self.assertTrue(self.data.optdone)
        temp = numpy.all(numpy.abs(self.data.geovalues) <= self.data.geotargets, axis=1)
        self.assertTrue(temp[-1])

class ADFGeoOptTest(GenericGeoOptTest):
    """ADF geometry optimization unittest."""

    extracoords = 1
    extrascfs = 1

    # ADF parser does not extract atombasis.
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique? PASS"""
        self.assertEquals(1, 1)
       
    def testscfenergy(self):
        """Is the SCF energy within 1eV of -140eV?"""
        scf = self.data.scfenergies[-1]
        ref = -140
        msg = "Final scf energy: %f not -140+-1eV" % scf
        self.assertInside(scf, ref, 1, msg)


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

    def testgrads(self):
        """Do the grads have the right dimensions?"""
        self.assertEquals(self.data.grads.shape,(len(self.data.geovalues),self.data.natom,3))


class JaguarGeoOptTest(GenericGeoOptTest):
    """Jaguar geometry optimization unittest."""

    # Data file does not contain enough information. Can we make a new one?
    def testatombasis(self):
        """Are the indices in atombasis the right amount and unique? PASS"""
        self.assertEquals(1, 1)

    # We did not print the atomic partial charges in the unit tests for this version.
    def testatomcharges(self):
        """Are all atomcharges consistent with natom and do they sum to zero? PASS"""
        self.assertEquals(1, 1)

    # Data file does not contain enough information. Can we make a new one?
    def testdimmocoeffs(self):
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis? PASS"""
        self.assertEquals(1, 1)


class MolproGeoOptTest(GenericGeoOptTest):
    """Molpro geometry optimization unittest."""

    # Note that these extra coordianates and energies will be available only
    # if the appropriate output is parsed, and Molpro often saves the initial
    # SCF run and subsequent geometry optimization to separate files, which
    # both need to be given to the cclib parser (as a list).
    extracoords = 1
    extrascfs = 2

    def testsymlabels(self):
        """Are all the symmetry labels either Ag/u or Bg/u? PASS"""
        self.assertEquals(1,1)

    # Here is what the manual has to say about convergence:
    # The standard MOLPRO convergency criterion requires the maximum component of the gradient
    # to be less then $3 \cdot 10^{-4}$ [a.u.] and the maximum energy change to be less than
    # $1 \cdot 10^{-6}$ [H] or the maximum component of the gradient to be less then
    # $3 \cdot 10^{-4}$ [a.u.] and the maximum component of the step to be less then
    # $3 \cdot 10^{-4}$ [a.u.].
    #
    # It is also possible to use the convergency criterion of (...)
    #
    # Source: https://www.molpro.net/info/2012.1/doc/manual/node592.html
    def testoptdone(self):
        """Has the geometry converged and set optdone to True?"""
        self.assertTrue(self.data.optdone)
        target_e, target_g, target_s = self.data.geotargets
        value_e, value_g, value_s = self.data.geovalues[-1]
        converged = (value_e < target_e and value_g < target_g) or (value_g < target_g and value_s < target_s)
        self.assertTrue(converged)


class MolproGeoOptTest2006(MolproGeoOptTest):
    """Molpro 2006 geometry optimization unittest."""

    # Same situation as SP -- this is tested for in the 2012 logfiles, but
    # the 2006 logfiles were created before atomcharges was an attribute and
    # we don't have access to Molpro 2006 anymore.
    def testatomcharges(self):
        """Are all atomcharges consistent with natom and do they sum to zero? PASS"""
        self.assertEquals(1,1)


class OrcaGeoOptTest(GenericGeoOptTest):
    """ORCA geometry optimization unittest."""

    extracoords = 1
    extrascfs = 1

    # ORCA has no support for symmetry yet.
    def testsymlabels(self):
        """Are all the symmetry labels either Ag/u or Bg/u? PASS"""
        self.assertEquals(1,1)

    # Besides all the geovalues being below their tolerances, ORCA also considers
    # an optimization finished in some extra cases. These are:
    #   1) everything converged except the energy (within 25 x tolerance)
    #   2) gradient is overachieved and displacement is reasonable (3 x tolerance)
    #   3) displacement is overachieved and gradient is reasonable (3 x tolerance)
    #   4) energy, gradients and angles are converged (displacements not considered)
    # All these exceptions are signaleld in the output with some comments, and here
    # we include the first three exceptions for the pruposes of the unit test.
    def testoptdone(self):
        """Has the geometry converged and set optdone to True?"""

        self.assertTrue(self.data.optdone)

        targets = self.data.geotargets
        values = numpy.abs(self.data.geovalues[-1])

        target_e = targets[0]
        target_g = targets[1:3]
        target_x = targets[3:]
        value_e = values[0]
        value_g = values[1:3]
        value_x = values[3:]

        conv_all = all(values < targets)
        conv_e = value_e < 25*target_e and all(value_g < target_g) and all(value_x < target_x)
        conv_g = value_e < target_e and all(value_g < target_g/3.0) and all(value_x < target_x*3.0)
        conv_x = value_e < target_e and all(value_g < target_g*3.0) and all(value_x < target_x/3.0)
        converged = conv_all or conv_e or conv_g or conv_x
        self.assertTrue(converged)


class PCGamessGeoOptTest(GenericGeoOptTest):
    """PC-GAMESS geometry optimization unittest."""


if __name__=="__main__":

    from testall import testall
    testall(modules=["GeoOpt"])
