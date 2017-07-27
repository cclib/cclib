# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for writer WFXwriter module."""

import os
import unittest

import cclib

from cclib.io.filewriter import MissingAttributeError

__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class WFXTest(unittest.TestCase):

    def test_missing_attribute_error(self):
        """Check if MissingAttributeError is raised as expected."""
        fpath = os.path.join(__datadir__,
                             "data/GAMESS/basicGAMESS-US2014/C_bigbasis.out")
        required_attrs = ('atomcoords', 'atomnos', 'gbasis', 'charge',
                          'homos', 'mult')
        for attr in required_attrs:
            data = cclib.io.ccopen(fpath).parse()
            delattr(data, attr)

            # WFX files cannot be written if required attrs are missing.
            with self.assertRaises(MissingAttributeError):
                cclib.io.wfxwriter.WFXWriter(data)


if __name__ == "__main__":
    unittest.main()
