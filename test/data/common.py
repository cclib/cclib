# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Functions used across multiple data tests."""

import itertools

from cclib.parser.data import ccData

import numpy


def get_minimum_carbon_separation(data: "ccData") -> float:
    """Returns minimum carbon distance for any coordinates.

    Note that atomcoords is 3D, and we will take the minimum
    over all coordinates and combinations of carbon atoms.
    """
    # mypy: disable-error-code="attr-defined"
    icarbons = numpy.arange(data.natom)[data.atomnos == 6]
    mindist = 999
    for i, j in itertools.combinations(icarbons, 2):
        vectors = data.atomcoords[:, i] - data.atomcoords[:, j]
        distances = numpy.linalg.norm(vectors, axis=1)
        mindist = min(mindist, min(distances))
    return mindist


def is_optnew(optstatus_value: int) -> bool:
    """Does the given optstatus value represent a new optimization?

    This does not require the optstatus value to only be this status
    exclusively.
    """
    return optstatus_value & ccData.OPT_NEW == ccData.OPT_NEW


def is_optunknown(optstatus_value: int) -> bool:
    """Does the given optstatus value represent an optimization in progress?

    This does not require the optstatus value to only be this status
    exclusively.  'Unknown' is not mutually exclusive with any other status.
    """
    return optstatus_value & ccData.OPT_UNKNOWN == ccData.OPT_UNKNOWN


def is_optdone(optstatus_value: int) -> bool:
    """Does the given optstatus value represent a converged optimization?

    This does not require the optstatus value to only be this status
    exclusively.  However, it should be exclusive with being unconverged.
    """
    return optstatus_value & ccData.OPT_DONE == ccData.OPT_DONE


def is_optunconverged(optstatus_value: int) -> bool:
    """Does the given optstatus value represent an unconverged optimization?

    This does not require the optstatus value to only be this status
    exclusively.  However, it should be exclusive with being converged.
    """
    return optstatus_value & ccData.OPT_UNCONVERGED == ccData.OPT_UNCONVERGED
