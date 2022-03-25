# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in PyQuante (http://pyquante.sourceforge.net)."""

import numpy
from cclib.parser.utils import find_package


class MissingAttributeError(Exception):
    pass


_found_pyquante2 = find_package("pyquante2")
if _found_pyquante2:
    from pyquante2 import molecule


def _check_pyquante():
    if not _found_pyquante2:
        raise ImportError("You must install `pyquante2` to use this function")


def makepyquante(data):
    """Create a PyQuante Molecule from ccData object."""
    _check_pyquante()

    # Check required attributes.
    required_attrs = {"atomcoords", "atomnos"}
    missing = [x for x in required_attrs if not hasattr(data, x)]

    if missing:
        missing = " ".join(missing)
        raise MissingAttributeError(
            f"Could not create pyquante molecule due to missing attribute: {missing}"
        )

    # In pyquante2, molecular geometry is specified in a format of:
    # [(3,.0000000000, .0000000000, .0000000000), (1, .0000000000, .0000000000,1.629912)]
    moldesc = numpy.insert(data.atomcoords[-1], 0, data.atomnos, 1).tolist()

    return molecule(
        [tuple(x) for x in moldesc],
        units="Angstroms",
        charge=data.charge,
        multiplicity=data.mult,
    )


del find_package
