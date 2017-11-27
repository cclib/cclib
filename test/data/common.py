# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

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
