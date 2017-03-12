# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""Bridge for using cclib data in biopython (http://biopython.org)."""

try:
    from Bio.PDB.Atom import Atom
except ImportError:
    # Fail silently for now.
    pass

from cclib.parser.utils import PeriodicTable


def makebiopython(atomcoords, atomnos):
    """Create a list of BioPython Atoms.

    This creates a list of BioPython Atoms suitable for use
    by Bio.PDB.Superimposer, for example.

    >>> import numpy
    >>> from Bio.PDB.Superimposer import Superimposer
    >>> atomnos = numpy.array([1,8,1],"i")
    >>> a = numpy.array([[-1,1,0],[0,0,0],[1,1,0]],"f")
    >>> b = numpy.array([[1.1,2,0],[1,1,0],[2,1,0]],"f")
    >>> si = Superimposer()
    >>> si.set_atoms(makebiopython(a,atomnos),makebiopython(b,atomnos))
    >>> print si.rms
    0.29337859596
    """
    pt = PeriodicTable()
    bioatoms = []
    for coords, atomno in zip(atomcoords, atomnos):
        symbol = pt.element[atomno]
        bioatoms.append(Atom(symbol, coords, 0, 0, 0, symbol, 0, symbol.upper()))
    return bioatoms


if __name__ == "__main__":
    import doctest
    doctest.testmod()
