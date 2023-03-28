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
        assert len(self.data.nmrtensors) == self.data.natom
        assert len(self.data.nmrtensors[0]) == 3
        assert self.data.nmrtensors[0]["total"].shape == (3, 3)


class GenericNMRCouplingTest(unittest.TestCase):
    """Generic NMR spin-spin coupling unittest"""

    def testsize(self):
        """Check to make sure there are the correct number of tensors parsed"""
        assert len(self.data.nmrcouplingtensors) == 139
        assert len(list(self.data.nmrcouplingtensors.values())[0]) == 6
        assert list(self.data.nmrcouplingtensors.values())[0]["total"].shape == (3, 3)


if __name__ == "__main__":
    import sys

    sys.path.insert(1, os.path.join(__filedir__, ".."))

    from test_data import DataSuite

    suite = DataSuite(['NMR'])
    suite.testall()
