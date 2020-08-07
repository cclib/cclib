# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in PySCF (https://github.com/pyscf/pyscf)."""

from cclib.parser.utils import find_package

_found_pyscf = find_package("pyscf")
if _found_pyscf:
    from pyscf import gto


def _check_pyscf(found_pyscf):
    if not found_pyscf:
        raise ImportError("You must install `pyscf` to use this function")


def makepyscf(atomcoords, atomnos, charge=0, mult=1):
    """Create a Pyscf Molecule."""
    _check_pyscf(_found_pyscf)
    mol = gto.Mole(
        atom = [['{}'.format(atomnos[i]),atomcoords[i]] for i in range(len(atomcoords))],
        unit="Angstrom",
        charge=charge,
        multiplicity=mult
    )
    return  mol
  
del find_package
