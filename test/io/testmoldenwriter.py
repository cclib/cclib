# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for writer moldenwriter module."""

import os
import unittest

import cclib
from cclib.io.filewriter import MissingAttributeError

__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class MOLDENTest(unittest.TestCase):

    def test_missing_attributes(self):
        """Check if MissingAttributeError is raised as expected."""
        fpath = os.path.join(__datadir__,
                             "data/ADF/basicADF2007.01/dvb_gopt.adfout")
        data = cclib.io.ccopen(fpath).parse()
        del data.atomcoords

        # Molden files cannot be wriiten if atomcoords are missing.
        with self.assertRaises(MissingAttributeError):
            cclib.io.moldenwriter.MOLDEN(data)


if __name__ == "__main__":
    unittest.main()
