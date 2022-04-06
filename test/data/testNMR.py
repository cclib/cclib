# -*- coding: utf-8 -*-
#
# Copyright (c) 2021, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test NMR logfiles in cclib."""

import os
import unittest

from skip import skipForParser

__filedir__ = os.path.realpath(os.path.dirname(__file__))


class GenericNMRTest(unittest.TestCase):
    """Generic NMR unittest"""

    def testsize(self):
        """Check to make sure there are the correct number of tensors parsed"""
        self.assertEqual(len(self.data.nmrtensors), self.data.natom)
        self.assertEqual(len(self.data.nmrtensors[0]), 3)
        self.assertEqual(self.data.nmrtensors[0]["total"].shape, (3, 3))


if __name__ == "__main__":
    import sys

    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite

    suite = DataSuite(["NMR"])
    suite.testall()
