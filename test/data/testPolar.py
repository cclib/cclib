# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test logfiles with (non)linear response output in cclib"""

import numpy
from skip import skipForParser


class GenericPolarTest:
    """Generic static polarizability unittest"""

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testshape(self, data) -> None:
        """Is the dimension of the polarizability tensor 3 x 3?"""
        assert len(data.polarizabilities) == 1
        assert data.polarizabilities[0].shape == (3, 3)


class ReferencePolarTest(GenericPolarTest):
    """Customized static polarizability unittest"""

    # Reference values are from DALTON2015/Trp_polar_abalnr.out
    isotropic = 74.12424
    principal_components = [30.29431523, 91.5361917, 100.54220307]
    isotropic_delta = 0.01
    principal_components_delta = 0.01

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testisotropic(self, data) -> None:
        """Is the isotropic polarizability (average of the diagonal elements)
        +/- 0.01 from a reference?
        """
        isotropic = numpy.average(numpy.diag(data.polarizabilities[0]))
        assert abs(isotropic - self.isotropic) < self.isotropic_delta

    @skipForParser("CFOUR", "The parser is still being developed so we skip this test")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testprincomponents(self, data) -> None:
        """Are each of the principal components (eigenvalues) of the
        polarizability tensor +/- 0.01 from a reference?
        """
        principal_components = numpy.linalg.eigvalsh(data.polarizabilities[0])
        for c in range(3):
            assert (
                abs(principal_components[c] - self.principal_components[c])
                < self.principal_components_delta
            )
