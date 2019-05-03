# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in PyQuante (http://pyquante.sourceforge.net)."""

from cclib.parser.utils import find_package

_found_pyquante = find_package("PyQuante")
if _found_pyquante:
    from PyQuante.Molecule import Molecule


def _check_pyquante(found_pyquante):
    if not found_pyquante:
        raise ImportError("You must install `PyQuante` to use this function")


def makepyquante(atomcoords, atomnos, charge=0, mult=1):
    """Create a PyQuante Molecule."""
    _check_pyquante(_found_pyquante)
    return Molecule("notitle",
                    list(zip(atomnos, atomcoords)),
                    units="Angstrom",
                    charge=charge, multiplicity=mult)


del find_package
