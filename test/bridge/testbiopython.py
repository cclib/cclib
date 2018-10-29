# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import unittest

import numpy

from cclib.bridge import cclib2biopython


class BiopythonTest(unittest.TestCase):
    """Tests for the cclib2biopython bridge in cclib."""

    def test_makebiopython(self):
        from Bio.PDB.Superimposer import Superimposer
        atomnos = numpy.array([1, 8, 1], "i")
        a = numpy.array([[-1, 1, 0], [0, 0, 0], [1, 1, 0]], "f")
        b = numpy.array([[1.1, 2, 0], [1, 1, 0], [2, 1, 0]], "f")
        si = Superimposer()
        si.set_atoms(cclib2biopython.makebiopython(a, atomnos),
                     cclib2biopython.makebiopython(b, atomnos))
        ref = 0.29337859596
        assert abs(si.rms - ref) < 1.0e-6


if __name__ == "__main__":
    unittest.main()
