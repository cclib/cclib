# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for the time attribute."""

import inspect
import os
import unittest

import cclib

import utils


class TimeTest(utils.DataTestCase):

    @utils.for_job_type('BOMD', nsteps=35)
    def test_shape(self, data, nsteps):
        self.assertEqual(data.time.shape, (nsteps, ))


if __name__ == '__main__':
    unittest.main()
