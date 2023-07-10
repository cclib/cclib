# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import unittest

from cclib.bridge import cclib2rdkit
from cclib.parser.utils import find_package

from ..test_data import getdatafile


class RdkitTest(unittest.TestCase):
    """Test for the cclib2rdkit bridge in cclib."""

    def test_makerdkit(self) -> None:
        if not find_package("rdkit"):
            raise ImportError("Must install rdkit to run this test")

        data, _ = getdatafile("Gaussian", "basicGaussian16", ["water_ccsd.log"])

        mol = cclib2rdkit.makerdkit(data)


if __name__ == "__main__":
    unittest.main()
