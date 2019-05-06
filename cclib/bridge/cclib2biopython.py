# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Bridge for using cclib data in biopython (http://biopython.org)."""

from cclib.parser.utils import PeriodicTable
from cclib.parser.utils import find_package

_found_biopython = find_package("Bio")
if _found_biopython:
    from Bio.PDB.Atom import Atom


def makebiopython(atomcoords, atomnos):
    """Create a list of BioPython Atoms.

    This creates a list of BioPython Atoms suitable for use by
    Bio.PDB.Superimposer, for example.
    """
    if not _found_biopython:
        raise ImportError("You must install `biopython` to use this function")
    pt = PeriodicTable()
    bioatoms = []
    for coords, atomno in zip(atomcoords, atomnos):
        symbol = pt.element[atomno]
        bioatoms.append(Atom(symbol, coords, 0, 0, 0, symbol, 0, symbol.upper()))
    return bioatoms


del find_package
