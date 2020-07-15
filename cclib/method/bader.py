# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculation of Bader's QTAIM charges based on data parsed by cclib."""
import random
import numpy
import logging
import math

from cclib.method.calculationmethod import Method
from cclib.method.volume import electrondensity_spin
from cclib.parser.utils import convertor

# Distance between two adjacent grids (sqrt[2] or sqrt[3] for uniform Cartesian grid)
griddist = [
    [[1.73205, 1.41421, 1.73205], [1.41421, 1, 1.41421], [1.73205, 1.41421, 1.73205],],
    [[1.41421, 1, 1.41421], [1, 1, 1], [1.41421, 1, 1.41421]],
    [[1.73205, 1.41421, 1.73205], [1.41421, 1, 1.41421], [1.73205, 1.41421, 1.73205],],
]


class MissingInputError(Exception):
    pass


class Bader(Method):
    """Bader's QTAIM charges."""

    # All of these are required for QTAIM charges.
    required_attrs = ("homos", "mocoeffs", "nbasis", "gbasis")

    def __init__(self, data, volume, progress=None, loglevel=logging.INFO, logname="Log"):
        super(Bader, self).__init__(data, progress, loglevel, logname)

        self.volume = volume
        self.fragresults = None

    def __str__(self):
        """Return a string representation of the object."""
        return "Bader's QTAIM charges of {}".format(self.data)

    def __repr__(self):
        """Return a representation of the object."""
        return "Bader({})".format(self.data)

    def _check_required_attributes(self):
        super(Bader, self)._check_required_attributes()

    def calculate(self, indices=None, fupdate=0.05):
        """Calculate Bader's QTAIM charges using on-grid algorithm proposed by Henkelman group
           in doi:10.1016/j.commatsci.2005.04.010
           
           Cartesian, uniformly spaced grids are assumed for this function.
           """

        def cartesian_dist(pt1, pt2):
            """ Small utility function that calculates Euclidian distance between two points
                pt1 and pt2 are numpy arrays representing a point in Cartesian coordinates """
            return numpy.sqrt(numpy.einsum("ij,ij->j", pt1 - pt2, pt1 - pt2s))

        # First obtain charge densities on the grid
        if len(self.data.mocoeffs) == 1:
            self.chgdensity = electrondensity_spin(
                self.data, self.volume, [self.data.mocoeffs[0][: self.data.homos[0]]]
            )
            self.chgdensity.data *= 2
        else:
            self.chgdensity = electrondensity_spin(
                self.data,
                self.volume,
                [
                    self.data.mocoeffs[0][: self.data.homos[0]],
                    self.data.mocoeffs[1][: self.data.homos[1]],
                ],
            )

        # Assign each grid point to Bader areas
        self.fragresults = numpy.zeros(self.chgdensity.data.shape, "d")
        nextIndex = 1

        self.logger.info("Partitioning space into Bader areas.")

        # Exclude the outermost positions
        for xwalker in range(1, len(self.chgdensity.data) - 1):
            for ywalker in range(1, len(self.chgdensity.data[0]) - 1):
                for zwalker in range(1, len(self.chgdensity.data[0][0]) - 1):
                    (xindex, yindex, zindex) = (xwalker, ywalker, zwalker)

                    if self.fragresults[xindex, yindex, zindex] != 0:
                        # index has already been assigned for this grid point
                        continue
                    else:
                        listcoord = []
                        loopcondition = True

                        while loopcondition:
                            tmp = (
                                self.chgdensity.data[
                                    xindex - 1 : xindex + 2,
                                    yindex - 1 : yindex + 2,
                                    zindex - 1 : zindex + 2,
                                ]
                                - self.chgdensity.data[xindex, yindex, zindex]
                            )
                            grad = tmp / griddist
                            maxat = numpy.where(grad == numpy.amax(grad))

                            nextDirection = list(zip(maxat[0], maxat[1], maxat[2]))[0]
                            nextDirection = [ind - 1 for ind in nextDirection]

                            if len(list(zip(maxat[0], maxat[1], maxat[2]))) > 1:
                                # when one or more directions indicate max grad (of 0), prioritize
                                # to include all points in the Bader space
                                if list(zip(maxat[0], maxat[1], maxat[2]))[0] == (1, 1, 1):
                                    nextDirection = list(zip(maxat[0], maxat[1], maxat[2]))[1]
                                    nextDirection = [ind - 1 for ind in nextDirection]

                            listcoord.append((xindex, yindex, zindex))

                            if (
                                self.fragresults[
                                    xindex + nextDirection[0],
                                    yindex + nextDirection[1],
                                    zindex + nextDirection[2],
                                ]
                                != 0
                            ):
                                # Path arrived at a point that has already been assigned with an index
                                baderIndex = self.fragresults[
                                    xindex + nextDirection[0],
                                    yindex + nextDirection[1],
                                    zindex + nextDirection[2],
                                ]
                                listcoord = tuple(numpy.array(listcoord).T)
                                self.fragresults[listcoord] = baderIndex

                                loopcondition = False

                            elif nextDirection == [0, 0, 0]:
                                # This point is local maximum
                                # Assign this point and everything along the path with a new index
                                baderIndex = nextIndex
                                nextIndex += 1

                                listcoord = tuple(numpy.array(listcoord).T)
                                self.fragresults[listcoord] = baderIndex

                                loopcondition = False

                            elif (
                                xindex + nextDirection[0] == 0
                                or xindex + nextDirection[0] == (len(self.chgdensity.data) - 1)
                                or yindex + nextDirection[1] == 0
                                or yindex + nextDirection[1] == (len(self.chgdensity.data[0]) - 1)
                                or zindex + nextDirection[2] == 0
                                or zindex + nextDirection[2]
                                == (len(self.chgdensity.data[0][0]) - 1)
                            ):
                                # Avoid exiting the grid
                                baderIndex = nextIndex
                                nexIndex += 1

                                listcoord = tuple(np.array(listcoord).T)
                                self.fragresults[listcoord] = baderIndex

                                loopcondition = False

                            else:
                                # Advance to the next point according to the direction of
                                # maximum gradient
                                xindex += nextDirection[0]
                                yindex += nextDirection[1]
                                zindex += nextDirection[2]

        # Now try to identify each Bader region to individual atom.
        # Try to find an area that captures enough representation
        self.matches = numpy.zeros_like(self.data.atomnos)
        for pos in range(len(self.data.atomcoords[-1])):
            gridpt = numpy.round(
                (self.data.atomcoords[-1][pos] - self.volume.origin) / self.volume.spacing
            )

            xgrid = int(gridpt[0])
            ygrid = int(gridpt[1])
            zgrid = int(gridpt[2])

            self.matches[pos] = self.fragresults[xgrid, ygrid, zgrid]

        assert (
            0 not in self.matches
        ), "Failed to assign Bader regions to atoms. Try with a finer grid. Content of Bader area matches: {}".format(
            matches
        )
        assert len(
            numpy.unique(self.matches) != len(self.data.atomnos)
        ), "Failed to assign unique Bader regions to each atom. Try with a finer grid."

        # Finally integrate the assigned Bader areas
        self.logger.info("Creating fragcharges: array[1]")
        boxvol = (
            self.volume.spacing[0]
            * self.volume.spacing[1]
            * self.volume.spacing[2]
            * convertor(1, "Angstrom", "bohr") ** 3
        )
        self.fragcharges = numpy.zeros(len(self.data.atomcoords[-1]), "d")

        for atomIndex, baderareaIndex in enumerate(self.matches):
            # turn off all other grid points
            vals = numpy.copy(self.chgdensity.data)
            mask = self.fragresults == baderareaIndex
            vals *= mask
            self.fragcharges[atomIndex] = sum(vals.ravel()) * boxvol

        return True
