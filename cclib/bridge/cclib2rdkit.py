# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in rdkit (https://rdkit.org)."""

from cclib.parser.data import ccData
from cclib.parser.utils import find_package

_found_rdkit = find_package("rdkit")

if _found_rdkit:
    from rdkit import Chem

def _check_rdkit():
    if not _found_rdkit:
        raise ImportError("You must install `rdkit >= 2022.09` to use this function.")
    else:
        try:
            from rdkit.Chem import rdDetermineBonds
        except ImportError:
            raise ImportError('This script requires RDKit version 2022.09 or later')


def makerdkit(data, explicit=True):
    """ Create rdkit Mol object from ccData object """

    _check_rdkit()

    mol = Chem.RWMol()
    for atom in data.atomnos:
        mol.AddAtom(Chem.Atom(int(atom)))
    # we'll need space for the 3D coordinates
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for i in range(len(data.atomcoords[-1])):
        conformer.SetAtomPosition(i, data.atomcoords[-1][i])
    mol.AddConformer(conformer)
    # now the new bond connectivity and bond orders
    bonded_mol = Chem.Mol(mol)
    rdDetermineBonds.DetermineConnectivity(bonded_mol)
    rdDetermineBonds.DetermineBondOrders(bonded_mol, data.charge)

    if not explicit:
        return Chem.RemoveHs(bonded_mol, implicitOnly=False)
    else:
        return bonded_mol

del find_package
