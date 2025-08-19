# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

from cclib.bridge import cclib2biopython
from cclib.parser.utils import find_package

import numpy


class BiopythonTest:
    """Tests for the cclib2biopython bridge in cclib."""

    def setup_method(self) -> None:
        if not find_package("Bio"):
            raise ImportError("Must install biopython to run this test")

    def test_makebiopython(self) -> None:
        from Bio.PDB.Superimposer import Superimposer

        atomnos = numpy.array([1, 8, 1], "i")
        a = numpy.array([[-1, 1, 0], [0, 0, 0], [1, 1, 0]], "f")
        b = numpy.array([[1.1, 2, 0], [1, 1, 0], [2, 1, 0]], "f")
        si = Superimposer()
        si.set_atoms(
            cclib2biopython.makebiopython(a, atomnos), cclib2biopython.makebiopython(b, atomnos)
        )
        ref = 0.29337859596
        assert abs(si.rms - ref) < 1.0e-6
