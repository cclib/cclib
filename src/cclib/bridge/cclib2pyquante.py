# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in PyQuante (http://pyquante.sourceforge.net)."""

from __future__ import print_function

import sys

try:
    from PyQuante.Molecule import Molecule
except ImportError:
    # Fail silently for now.
    pass


def makepyquante(atomcoords, atomnos, charge=0, mult=1):
    """Create a PyQuante Molecule."""
    return Molecule("notitle", list(zip(atomnos, atomcoords)), units="Angstrom",
                    charge=charge, multiplicity=mult)
