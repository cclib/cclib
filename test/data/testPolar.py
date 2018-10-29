# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test logfiles with (non)linear response output in cclib"""

import os
import unittest

import numpy

from skip import skipForParser

__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericPolarTest(unittest.TestCase):
    """Generic static polarizability unittest"""

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testshape(self):
        """Is the dimension of the polarizability tensor 3 x 3?"""
        self.assertEqual(len(self.data.polarizabilities), 1)
        self.assertEqual(self.data.polarizabilities[0].shape, (3, 3))


class ReferencePolarTest(GenericPolarTest):
    """Customized static polarizability unittest"""

    # Reference values are from DALTON2015/Trp_polar_abalnr.out
    isotropic = 74.12424
    principal_components = [30.29431523, 91.5361917, 100.54220307]
    isotropic_delta = 0.01
    principal_components_delta = 0.01

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testisotropic(self):
        """Is the isotropic polarizability (average of the diagonal elements)
        +/- 0.01 from a reference?
        """
        isotropic = numpy.average(numpy.diag(self.data.polarizabilities[0]))
        self.assertAlmostEqual(isotropic, self.isotropic, delta=self.isotropic_delta)

    @skipForParser('Molcas','The parser is still being developed so we skip this test')
    @skipForParser('Turbomole','The parser is still being developed so we skip this test')
    def testprincomponents(self):
        """Are each of the principal components (eigenvalues) of the
        polarizability tensor +/- 0.01 from a reference?
        """
        principal_components = numpy.linalg.eigvalsh(self.data.polarizabilities[0])
        for c in range(3):
            self.assertAlmostEqual(principal_components[c],
                                   self.principal_components[c],
                                   delta=self.principal_components_delta)


if __name__=="__main__":

    import sys
    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite
    suite = DataSuite(['Polar'])
    suite.testall()
