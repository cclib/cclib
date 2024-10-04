# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test Moller-Plesset logfiles in cclib"""

import numpy


class GenericMP2Test:
    """Generic MP2 unittest"""

    level = 2

    def testsizeandshape(self, data) -> None:
        """(MP2) Are the dimensions of mpenergies correct?"""
        assert data.mpenergies.shape == (len(data.scfenergies), self.level - 1)

    def testsign(self, data) -> None:
        """Are the Moller-Plesset corrections negative?"""
        if self.level == 2:
            corrections = data.mpenergies[:, 0] - data.scfenergies
        else:
            corrections = data.mpenergies[:, self.level - 2] - data.mpenergies[:, self.level - 3]
        assert numpy.all(corrections < 0.0)


class GenericMP3Test(GenericMP2Test):
    """Generic MP3 unittest"""

    level = 3


class GenericMP4SDQTest(GenericMP2Test):
    """Generic MP4(SDQ) unittest"""

    level = 4


class GenericMP4SDTQTest(GenericMP2Test):
    """Generic MP4(SDTQ) unittest"""

    level = 4


class GenericMP5Test(GenericMP2Test):
    """Generic MP5 unittest"""

    level = 5


class GaussianMP2Test(GenericMP2Test):
    """Customized MP2 unittest"""

    def testnocoeffs(self, data) -> None:
        """Are natural orbital coefficients the right size?"""
        assert data.nocoeffs.shape == (data.nmo, data.nbasis)

    def testnooccnos(self, data) -> None:
        """Are natural orbital occupation numbers the right size?"""
        assert data.nooccnos.shape == (data.nmo,)


class GaussianMP3Test(GenericMP2Test):
    """Customized MP3 unittest"""

    level = 3


class GaussianMP4SDQTest(GenericMP2Test):
    """Customized MP4-SDQ unittest"""

    level = 4


class GaussianMP4SDTQTest(GenericMP2Test):
    """Customized MP4-SDTQ unittest"""

    level = 4


class QChemMP4SDQTest(GenericMP2Test):
    """Customized MP4-SDQ unittest"""

    level = 4


class QChemMP4SDTQTest(GenericMP2Test):
    """Customized MP4-SD(T)Q unittest"""

    level = 5
