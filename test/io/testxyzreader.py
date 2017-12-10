# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for XYZ reader."""

import os
import unittest

import numpy as np

import cclib


__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class XYZReaderTest(unittest.TestCase):

    contents_single = """12
uracil
N   5.435   3.248   -0.916
C   4.091   3.266   -0.635
H   6.109   3.347   -0.144
O   3.647   3.394   0.495
N   3.253   3.118   -1.717
C   3.636   2.955   -3.039
O   2.790   2.842   -3.923
C   5.065   2.953   -3.239
C   5.906   3.101   -2.202
H   2.243   3.138   -1.514
H   5.448   2.833   -4.258
H   7.003   3.119   -2.345
"""

    def test_attributes_single(self):
        """Is an XYZ file with a single geometry read into a ccData properly?
        """
        fpath = os.path.join(__datadir__, "test/bridge/uracil.xyz")
        assert os.path.exists(fpath)
        xyz = cclib.io.xyzreader.XYZ(fpath)
        xyz.read()
        data = xyz.generate_repr()

        attrs = ('natom', 'atomnos', 'atomcoords')
        for attr in attrs:
            self.assertTrue(hasattr(data, attr))

        self.assertEqual(data.natom, 12)

        atomnos = np.array([7, 6, 1, 8, 7, 6, 8, 6, 6, 1, 1, 1], dtype=int)

        self.assertEqual(data.atomnos.shape, atomnos.shape)
        self.assertTrue(np.all(np.equal(data.atomnos, atomnos)))

        atomcoords = np.array([[
            [5.435, 3.248, -0.916],
            [4.091, 3.266, -0.635],
            [6.109, 3.347, -0.144],
            [3.647, 3.394,  0.495],
            [3.253, 3.118, -1.717],
            [3.636, 2.955, -3.039],
            [2.790, 2.842, -3.923],
            [5.065, 2.953, -3.239],
            [5.906, 3.101, -2.202],
            [2.243, 3.138, -1.514],
            [5.448, 2.833, -4.258],
            [7.003, 3.119, -2.345]
        ]])

        self.assertEqual(data.atomcoords.shape, atomcoords.shape)
        self.assertTrue(np.all(np.equal(data.atomcoords, atomcoords)))

if __name__ == "__main__":
    unittest.main()
