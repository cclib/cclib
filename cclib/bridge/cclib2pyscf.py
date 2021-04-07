# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in PySCF (https://github.com/pyscf/pyscf)."""

from cclib.parser.utils import find_package, PeriodicTable
import numpy as np

l_sym2num = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}


_found_pyscf = find_package("pyscf")
if _found_pyscf:
    from pyscf import gto


def _check_pyscf(found_pyscf):
    if not found_pyscf:
        raise ImportError("You must install `pyscf` to use this function")


def makepyscf(data, charge=0, mult=1):
    """Create a Pyscf Molecule."""
    _check_pyscf(_found_pyscf)
    mol = gto.Mole(
        atom=[
            ["{}".format(data.atomnos[i]), data.atomcoords[-1][i]]
            for i in range(data.natom)
        ],
        unit="Angstrom",
        charge=charge,
        multiplicity=mult,
    )
    inputattr = data.__dict__
    pt = PeriodicTable()
    if "gbasis" in inputattr:
        basis = {}  # object for internal PySCF format
        uatoms, uatoms_idx = np.unique(
            data.atomnos, return_index=True
        )  # find unique atoms
        for idx, i in enumerate(uatoms_idx):
            curr_atom_basis = data.gbasis[i]
            for jdx, j in enumerate(curr_atom_basis):
                curr_l = j[0]
                curr_e_prim = j[1]
                new_list = [l_sym2num["{}".format(curr_l)]]
                new_list += curr_e_prim
                if not "{}".format(pt.element[uatoms[idx]]) in basis:
                    basis["{}".format(pt.element[uatoms[idx]])] = [new_list]
                else:
                    basis["{}".format(pt.element[uatoms[idx]])].append(new_list)
        mol.basis = basis
    return mol


del find_package
