# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for the scfenergies attribute."""

import inspect
import os
import unittest

import cclib

import utils


class AtomCoordsTest(utils.DataTestCase):

    @utils.for_job_type('BOMD', nenergies=36, natoms=20)
    @utils.for_program(
        'Gaussian', nenergies=35,
        msg='This may have something to do with our corrections for extrapolation step rejections.')
    def test_shape(self, data, nenergies, natoms):
        self.assertEqual(data.atomcoords.shape, (nenergies, natoms, 3))


if __name__ == '__main__':
    unittest.main()
