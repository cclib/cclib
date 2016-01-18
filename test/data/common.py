# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2016, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Functions used across multiple data tests."""

import itertools

import numpy


def get_minimum_carbon_separation(data):
    """Returns minimum carbon distance for any coordinates.

    Note that atomcoords is 3D, and we will take the minimum
    over all coordinates and combinations of carbon atoms.
    """

    icarbons = numpy.arange(data.natom)[data.atomnos == 6]
    mindist = 999
    for i, j in itertools.combinations(icarbons, 2):
        vectors = data.atomcoords[:, i] - data.atomcoords[:, j]
        distances = numpy.linalg.norm(vectors, axis=1)
        mindist = min(mindist, min(distances))
    return mindist
