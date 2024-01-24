# -*- coding: utf-8 -*-
#
# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for XYZ reader."""

import os

import cclib

import numpy as np

__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class XYZReaderTest:
    def test_attributes_one(self):
        """Is an XYZ file with a single geometry read into a ccData properly?"""
        fpath = os.path.join(__datadir__, "test/bridge/uracil.xyz")
        xyz = cclib.io.xyzreader.XYZ(fpath)
        data = xyz.parse()

        assert data.natom == 12

        atomnos = np.array([7, 6, 1, 8, 7, 6, 8, 6, 6, 1, 1, 1], dtype=int)
        np.testing.assert_equal(data.atomnos, atomnos)

        assert data.atomcoords.shape == (1, 12, 3)

    def test_attributes_two(self):
        """Is an XYZ file with a two geometries read into a ccData properly?"""
        fpath = os.path.join(__filedir__, "data/uracil_two.xyz")
        xyz = cclib.io.xyzreader.XYZ(fpath)
        data = xyz.parse()

        assert data.natom == 12

        atomnos = np.array([7, 6, 1, 8, 7, 6, 8, 6, 6, 1, 1, 1], dtype=int)
        np.testing.assert_equal(data.atomnos, atomnos)

        assert data.atomcoords.shape == (2, 12, 3)
        assert data.metadata["comments"] == ["uracil", "	Energy:     -97.0646597"]
