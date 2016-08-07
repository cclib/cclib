# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2016 the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Test logfiles with (non)linear response output in cclib"""

import os
import unittest

import numpy


__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericPolarTest(unittest.TestCase):
    """Generic static polarizability unittest"""

    # Reference values are from DALTON2015/Trp_polar_abalnr.out
    isotropic = 74.12424
    principal_components = [30.29431523, 91.5361917, 100.54220307]
    isotropic_delta = 0.01
    principal_components_delta = 0.01

    def testshape(self):
        """Is the dimension of the polarizability tensor 3 x 3?"""
        self.assertEqual(len(self.data.polarizabilities), 1)
        self.assertEqual(self.data.polarizabilities[0].shape, (3, 3))

    def testisotropic(self):
        """Is the isotropic polarizability (average of the diagonal elements)
        +/- 0.01 from a reference?
        """
        isotropic = numpy.average(numpy.diag(self.data.polarizabilities[0]))
        self.assertAlmostEqual(isotropic, self.isotropic, delta=self.isotropic_delta)

    def testprincomponents(self):
        """Are each of the principal components (eigenvalues) of the
        polarizability tensor +/- 0.01 from a reference?
        """
        principal_components = numpy.linalg.eigvalsh(self.data.polarizabilities[0])
        for c in range(3):
            self.assertAlmostEqual(principal_components[c],
                                   self.principal_components[c],
                                   delta=self.principal_components_delta)


class GaussianPolarTest(GenericPolarTest):
    """Customized static polarizability unittest"""

    # Reference values are from Q-Chem 4.2/trithiolane_freq.out, since
    # with symmetry enabled Q-Chem reorients molecules similarly to
    # Gaussian.
    isotropic = 66.0955766
    principal_components = [46.71020322, 75.50778705, 76.06873953]
    # Make the thresholds looser because these test jobs use symmetry,
    # and the polarizability is orientation dependent.
    isotropic_delta = 2.0
    principal_components_delta = 0.7


if __name__=="__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['Polar'])
    suite.testall()
