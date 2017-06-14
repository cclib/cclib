# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for writer moldenwriter module."""

import os
import unittest

import cclib


__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class MOLDENTest(unittest.TestCase):

    def setUp(self):
        self.molden = cclib.io.MOLDEN

    def test_init(self):
        """Does the class initialize correctly?"""
        fpath = os.path.join(__datadir__,
                             "data/ADF/basicADF2007.01/dvb_gopt.adfout")
        data = cclib.io.ccopen(fpath).parse()
        molden = cclib.io.moldenwriter.MOLDEN(data)

        # The object should keep the ccData instance passed to its constructor.
        self.assertEqual(molden.ccdata, data)


if __name__ == "__main__":
    unittest.main()
