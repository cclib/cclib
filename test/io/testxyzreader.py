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

    def test_attributes_single(self):
        """Is an XYZ file with a single geometry read into a ccData properly?
        """
        fpath = os.path.join(__datadir__, "test/bridge/uracil.xyz")
        assert os.path.exists(fpath)
        xyz = cclib.io.xyzreader.XYZ(fpath)
        data = xyz.parse()

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

    def test_attributes_two(self):
        """Is an XYZ file with a two geometries read into a ccData properly?
        """
        fpath = os.path.join(__filedir__, "data/uracil_two.xyz")
        assert os.path.exists(fpath)
        xyz = cclib.io.xyzreader.XYZ(fpath)
        data = xyz.parse()

        attrs = ('natom', 'atomnos', 'atomcoords')
        for attr in attrs:
            self.assertTrue(hasattr(data, attr))

        self.assertEqual(data.natom, 12)

        atomnos = np.array([7, 6, 1, 8, 7, 6, 8, 6, 6, 1, 1, 1], dtype=int)

        self.assertEqual(data.atomnos.shape, atomnos.shape)
        self.assertTrue(np.all(np.equal(data.atomnos, atomnos)))

        atomcoords = np.array([
            [
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
            ],
            [
                [5.42432, 3.24732, -0.94320],
                [4.08846, 3.26356, -0.64127],
                [6.05396, 3.35400, -0.15563],
                [3.68445, 3.39863,  0.50851],
                [3.24195, 3.11963, -1.70419],
                [3.62987, 2.96568, -3.00540],
                [2.83542, 2.83961, -3.92970],
                [5.08640, 2.95832, -3.25044],
                [5.91213, 3.09904, -2.21101],
                [2.25340, 3.12801, -1.51035],
                [5.42415, 2.83814, -4.27067],
                [6.99150, 3.10206, -2.32366]
            ]
        ])

        self.assertEqual(data.atomcoords.shape, atomcoords.shape)
        self.assertTrue(np.all(np.equal(data.atomcoords, atomcoords)))


if __name__ == "__main__":
    unittest.main()
