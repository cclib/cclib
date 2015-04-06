# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006,2007,2009,2012-2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test geometry optimization logfiles in cclib"""

import os
import unittest

import numpy


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericGeoOptTest(unittest.TestCase):
    """Generic geometry optimization unittest"""

    # In STO-3G, H has 1, C has 3.
    nbasisdict = {1:1, 6:5}
    
    # Some programs print surplus atom coordinates by default.
    extracoords = 0
    
    # Some programs do surplus SCF cycles by default.
    extrascfs = 0

    # Approximate B3LYP energy of dvb after SCF in STO-3G.
    b3lyp_energy = -10365
    b3lyp_tolerance = 40

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
                        dist = numpy.sqrt(sum((final_x - final_y)**2))
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

    def testcoreelectrons(self):
        """Are the coreelectrons all 0?"""
        ans = numpy.zeros(self.data.natom, 'i')
        numpy.testing.assert_array_equal(self.data.coreelectrons, ans)

    def testnormalisesym(self):
        """Did this subclass overwrite normalisesym?"""
        self.assertNotEquals(self.logfile.normalisesym("A"),"ERROR: This should be overwritten by this subclass")

    def testhomos(self):
        """Is the index of the HOMO equal to 34?"""
        ref = numpy.array([34], "i")
        msg = "%s != array([34], 'i')" % numpy.array_repr(self.data.homos)
        numpy.testing.assert_array_equal(self.data.homos, ref, msg)

    def testscfvaluetype(self):
        """Are scfvalues and its elements the right type?"""
        self.assertEquals(type(self.data.scfvalues),type([]))
        self.assertEquals(type(self.data.scfvalues[0]),type(numpy.array([])))

    def testscfenergy(self):
        """Is the SCF energy close to target?"""
        scf = self.data.scfenergies[-1]
        ref = self.b3lyp_energy
        tol = self.b3lyp_tolerance
        msg = "Final scf energy: %f not %i +- %ieV" %(scf, ref, tol)
        self.assertAlmostEquals(scf, ref, delta=40, msg=msg)

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

    def testmoenergies(self):
        """Are only the final MOs parsed?"""
        self.assertEquals(len(self.data.moenergies), 1)
        if hasattr(self.data, "mocoeffs"):
            self.assertEquals(len(self.data.mocoeffs), 1)


class ADFGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    extracoords = 1
    extrascfs = 1

    b3lyp_energy = -140
    b3lyp_tolerance = 1


class GaussianGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    def testgrads(self):
        """Do the grads have the right dimensions?"""
        self.assertEquals(self.data.grads.shape,(len(self.data.geovalues),self.data.natom,3))


class MolproGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    # Note that these extra coordinates and energies will be available only
    # if the appropriate output is parsed, and Molpro often saves the initial
    # SCF run and subsequent geometry optimization to separate files, which
    # both need to be given to the cclib parser (as a list).
    extracoords = 1
    extrascfs = 2

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

class NWChemGeoOptTest(GenericGeoOptTest):
    """Customized restricted single point HF unittest"""

    # NWChem typically prints the coordinates in the input module, at the
    # beginning of each geometry optimization step, and then again after
    # the optimziation is finished, so the first and last coordinates
    # are repeated. On the other hand, each optimization step often
    # involves a line search which we don't parse (see parse code for details).
    extracoords = 2
    extrascfs = 0


class OrcaGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    extracoords = 1
    extrascfs = 1

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


class PsiGeoOptTest(GenericGeoOptTest):
    """Customized geometry optimization unittest"""

    # Psi has a number of different convergence strategies to choose from, as described here:
    #     http://sirius.chem.vt.edu/psi4manual/latest/optking.html
    # and the default is to check that the max. force is converged and if the max energy change
    # or dispalcement is converged. This is in fact what is tested below.
    def testoptdone(self):
        """Has the geometry converged and set optdone to True?"""

        self.assertTrue(self.data.optdone)

        targets = self.data.geotargets
        values = numpy.abs(self.data.geovalues[-1])

        # Since the other criteria are not used and are not printed in this case, they should
        # be parsed as numpy.inf, for which we can check.
        self.assertTrue(numpy.isinf(targets[2]))
        self.assertTrue(numpy.isinf(targets[4]))

        conv = values[1] < targets[1] and (values[0] < targets[0] or values[3] < targets[3])
        self.assertTrue(conv)


if __name__=="__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['GeoOpt'])
    suite.testall()
