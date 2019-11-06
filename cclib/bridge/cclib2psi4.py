# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in Psi4 (http://www.psicode.org/)."""

from cclib.parser.utils import find_package

_found_psi4 = find_package("psi4")
if _found_psi4:
    from psi4.core import Molecule


def _check_psi4(found_psi4):
    if not found_psi4:
        raise ImportError("You must install `psi4` to use this function")


def makepsi4(atomcoords, atomnos, charge=0, mult=1):
    """Create a Psi4 Molecule."""
    _check_psi4(_found_psi4)
    return Molecule.from_arrays(
        name="notitle",
        elez=atomnos,
        geom=atomcoords,
        units="Angstrom",
        molecular_charge=charge,
        molecular_multiplicity=mult,
    )


del find_package
