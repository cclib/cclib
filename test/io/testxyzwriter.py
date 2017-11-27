# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for writer xyzwriter module."""

import os
import unittest

import cclib


__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class XYZTest(unittest.TestCase):

    def setUp(self):
        self.XYZ = cclib.io.XYZ

    def test_init(self):
        """Does the class initialize correctly?"""
        fpath = os.path.join(__datadir__, "data/ADF/basicADF2007.01/dvb_gopt.adfout")
        data = cclib.io.ccopen(fpath).parse()
        xyz = cclib.io.xyzwriter.XYZ(data)

        # The object should keep the ccData instance passed to its constructor.
        self.assertEqual(xyz.ccdata, data)


if __name__ == "__main__":
    unittest.main()
