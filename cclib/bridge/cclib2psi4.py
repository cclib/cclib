# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in Psi4 (http://www.psicode.org/)."""

from cclib.parser.utils import find_package

import numpy as np

_found_psi4 = find_package("psi4")
if _found_psi4:
    from psi4.core import Molecule


def _check_psi4(found_psi4: bool) -> None:
    if not found_psi4:
        raise ImportError("You must install `psi4` to use this function")


def makepsi4(
    atomcoords: np.ndarray, atomnos: np.ndarray, charge: int = 0, mult: int = 1
) -> "Molecule":
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
