# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import os
import unittest

import numpy

from cclib.bridge import cclib2openbabel


class OpenbabelTest(unittest.TestCase):
    """Tests for the cclib2openbabel bridge in cclib."""

    def setUp(self):
        self.path = os.path.abspath(os.path.dirname(__file__))

    def test_makeopenbabel(self):
        import openbabel
        atomnos = numpy.array([1, 8, 1], "i")
        atomcoords = numpy.array([[[-1., 1., 0.], [0., 0., 0.], [1., 1., 0.]]])
        obmol = cclib2openbabel.makeopenbabel(atomcoords, atomnos)
        obconversion = openbabel.OBConversion()
        formatok = obconversion.SetOutFormat("inchi")
        assert obconversion.WriteString(obmol).strip() == "InChI=1S/H2O/h1H2"

    def test_readfile(self):
        """Try to load an XYZ file with uracyl through Openbabel"""
        data = cclib2openbabel.readfile(self.path + "/uracil.xyz", "XYZ")
        assert data.natom == 12


if __name__ == "__main__":
    unittest.main()
