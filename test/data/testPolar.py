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

    # Reference values are from DALTON 2015/Trp_polar_abalnr.out
    min_component = -22.538044
    max_component = 92.329638
    principal_components = [30.29430402, 91.53628235, 100.54212364]

    def testshape(self):
        """Is the dimension of each polarizability tensor 3 x 3?"""
        for polarizability_tensor in self.data.polarizabilities:
            self.assertEqual(polarizability_tensor.shape, (3, 3))

    # def testmaxcomponent(self):
    #     """Is the max component of the polarizbality +/- 1 from a reference?"""
    #     self.assertAlmostEqual(np.max(self.data.polarizabilities[0]),
    #                            self.max_component,
    #                            delta=1.0)

    # def testmincomponent(self):
    #     """Is the min component of the polarizbality +/- 1 from a reference?"""
    #     self.assertAlmostEqual(np.min(self.data.polarizabilities[0]),
    #                            self.min_component,
    #                            delta=1.0)

    def testprincomponents(self):
        """Are each of the principal components (eigenvalues) of the
        polarizability tensor +/- 0.01 from a reference?
        """
        # It is much easier to compare the the eigenvalues of the
        # matrix rather than its individual components.
        principal_components = numpy.linalg.eigvalsh(self.data.polarizabilities[0])
        for c in range(3):
            self.assertAlmostEqual(principal_components[c],
                                   self.principal_components[c],
                                   delta=0.01)

if __name__=="__main__":

    import sys
    sys.path.append(os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['Polar'])
    suite.testall()
