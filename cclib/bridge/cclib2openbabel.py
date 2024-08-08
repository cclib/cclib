# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge between cclib data and openbabel (http://openbabel.org)."""

from typing import Optional, Set, Union

from cclib.parser.data import ccData
from cclib.parser.utils import find_package

import numpy as np

_found_openbabel = find_package("openbabel")
if _found_openbabel:
    # The exception idiom is needed because OB < 3.0 doesn't set __version__.
    # The `try` block is for OB >= 3.0, and `except` is for 2.4.x and older.
    try:
        from openbabel import openbabel as ob
    except:  # noqa: E722
        import openbabel as ob


def _check_openbabel(found_openbabel: bool) -> None:
    if not found_openbabel:
        raise ImportError("You must install `openbabel` to use this function")


def makecclib(mol: "ob.OBMol") -> ccData:
    """Create cclib attributes and return a ccData from an OpenBabel molecule.

    Beyond the numbers, masses and coordinates, we could also set the total charge
    and multiplicity, but often these are calculated from atomic formal charges
    so it is better to assume that would not be correct.
    """
    _check_openbabel(_found_openbabel)
    attributes = {"atomcoords": [], "atommasses": [], "atomnos": [], "natom": mol.NumAtoms()}
    for atom in ob.OBMolAtomIter(mol):
        attributes["atomcoords"].append([atom.GetX(), atom.GetY(), atom.GetZ()])
        attributes["atommasses"].append(atom.GetAtomicMass())
        attributes["atomnos"].append(atom.GetAtomicNum())
    attributes["atomcoords"] = [attributes["atomcoords"]]
    return ccData(attributes)


def makeopenbabel(
    *,
    atomcoords: Optional[np.ndarray] = None,
    atomnos: Optional[np.ndarray] = None,
    charge: int = 0,
    mult: int = 1,
) -> "ob.OBMol":
    """Create an Open Babel molecule."""
    _check_openbabel(_found_openbabel)
    if atomcoords is None and atomnos is None:
        raise RuntimeError("Must pass at least one of atomcoords or atomnos")
    elif atomcoords is not None and atomnos is not None:
        assert atomcoords.shape[1] == atomnos.shape[0]
        natom = atomnos.shape[0]
    elif atomcoords is None:
        natom = atomnos.shape[0]
    else:
        natom = atomcoords.shape[1]
    obmol = ob.OBMol()
    for i in range(natom):
        obatom = ob.OBAtom()
        if atomcoords is not None:
            # Note that list(atomcoords[i]) is not equivalent!!!
            # For now, only take the last geometry.
            # TODO: option to export last geometry or all geometries?
            coords = atomcoords[-1][i].tolist()
            obatom.SetVector(*coords)
        if atomnos is not None:
            atomno = int(atomnos[i])
            obatom.SetAtomicNum(atomno)
        obmol.AddAtom(obatom)
    obmol.ConnectTheDots()
    obmol.PerceiveBondOrders()
    obmol.SetTotalSpinMultiplicity(mult)
    obmol.SetTotalCharge(int(charge))
    return obmol


def readfile(fname: str, fmt: str) -> Union[ccData, Set]:
    """Read a file with OpenBabel and extract cclib attributes."""
    _check_openbabel(_found_openbabel)
    obc = ob.OBConversion()
    if obc.SetInFormat(fmt):
        mol = ob.OBMol()
        obc.ReadFile(mol, fname)
        return makecclib(mol)
    else:
        print(f"Unable to load the {fmt} reader from OpenBabel.")
        return {}


del find_package
