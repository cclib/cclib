# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge between cclib data and Open Babel (http://openbabel.org)."""

try:
    import openbabel as ob
    import pybel as pb
except ImportError:
    # Fail silently for now.
    pass

from cclib.parser.data import ccData


def makeopenbabel(atomcoords, atomnos, charge=0, mult=1):
    """Create an Open Babel molecule.

    >>> import numpy, openbabel
    >>> atomnos = numpy.array([1, 8, 1], "i")
    >>> coords = numpy.array([[-1., 1., 0.], [0., 0., 0.], [1., 1., 0.]])
    >>> obmol = makeopenbabel(coords, atomnos)
    >>> obconversion = openbabel.OBConversion()
    >>> formatok = obconversion.SetOutFormat("inchi")
    >>> print obconversion.WriteString(obmol).strip()
    InChI=1/H2O/h1H2
    """
    obmol = ob.OBMol()
    for i in range(len(atomnos)):
        # Note that list(atomcoords[i]) is not equivalent!!!
        # For now, only take the last geometry.
        # TODO: option to export last geometry or all geometries?
        coords = atomcoords[-1][i].tolist()
        atomno = int(atomnos[i])
        obatom = ob.OBAtom()
        obatom.SetAtomicNum(atomno)
        obatom.SetVector(*coords)
        obmol.AddAtom(obatom)
    obmol.ConnectTheDots()
    obmol.PerceiveBondOrders()
    obmol.SetTotalSpinMultiplicity(mult)
    obmol.SetTotalCharge(int(charge))
    return obmol


def makecclib(mols):
    """Create cclib attributes and return a single ccData from a list of
    Pybel (Open Babel wrapper) molecules. Assume that all molecules
    are identical except for their coordinates.

    Beyond the numbers, masses and coordinates, we could also set the
    total charge and multiplicity, but often these are calculated from
    atomic formal charges so it is better to assume that would not be
    correct.
    """

    if mols:
        attributes = {
            'atomcoords': [[atom.coords for atom in pbmol.atoms]
                           for pbmol in mols],
            'atommasses': [atom.atomicmass for atom in mols[0].atoms],
            'atomnos':    [atom.atomicnum for atom in mols[0].atoms],
            'natom':      len(mols[0].atoms),
        }
        return ccData(attributes)


def readfile(fname, format):
    """Read a file with Open Babel and extract cclib attributes."""

    mols = list(pb.readfile(format, fname))
    if mols:
        return makecclib(mols)
    else:
        print("Unable to load the %s reader from Open Babel." % format)
        return {}


if __name__ == "__main__":
    import doctest
    doctest.testmod()
