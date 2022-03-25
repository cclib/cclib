# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculation of Bader's QTAIM charges based on data parsed by cclib."""
import copy
import random
import numpy
import logging
import math

from cclib.method.calculationmethod import Method
from cclib.method.volume import electrondensity_spin
from cclib.parser.utils import convertor

# Distance between two adjacent grids (sqrt[2] or sqrt[3] for uniform Cartesian grid).
_griddist = numpy.array(
    [
        [[1.73205, 1.41421, 1.73205], [1.41421, 1, 1.41421], [1.73205, 1.41421, 1.73205],],
        [[1.41421, 1, 1.41421], [1, 1, 1], [1.41421, 1, 1.41421]],
        [[1.73205, 1.41421, 1.73205], [1.41421, 1, 1.41421], [1.73205, 1.41421, 1.73205],],
    ],
    dtype=float,
)


class MissingInputError(Exception):
    pass


def __cartesian_dist(pt1, pt2):
    """ Small utility function that calculates Euclidian distance between two points
        pt1 and pt2 are numpy arrays representing a point in Cartesian coordinates. """
    return numpy.sqrt(numpy.einsum("ij,ij->j", pt1 - pt2, pt1 - pt2))


class Bader(Method):
    """Bader's QTAIM charges."""

    # All of these are required for QTAIM charges.
    required_attrs = ("homos", "mocoeffs", "nbasis", "gbasis")

    def __init__(self, data, volume, progress=None, loglevel=logging.INFO, logname="Log"):
        super().__init__(data, progress, loglevel, logname)

        self.volume = volume
        self.fragresults = None

        if numpy.sum(self.data.coreelectrons) != 0:
            # Pseudopotentials can cause Bader spaces to be inaccurate, as suggested by the
            # original paper.
            self.logger.info(
                "It looks like pseudopotentials were used to generate this output. Please note that the Bader charges may not be accurate and may report unexpected results. Consult the original paper (doi:10.1016/j.commatsci.2005.04.010) for more information."
            )

    def __str__(self):
        """Return a string representation of the object."""
        return f"Bader's QTAIM charges of {self.data}"

    def __repr__(self):
        """Return a representation of the object."""
        return f"Bader({self.data})"

    def _check_required_attributes(self):
        super()._check_required_attributes()

    def calculate(self, indices=None, fupdate=0.05):
        """Calculate Bader's QTAIM charges using on-grid algorithm proposed by Henkelman group
           in doi:10.1016/j.commatsci.2005.04.010
           
           Cartesian, uniformly spaced grids are assumed for this function.
           """

        # Obtain charge densities on the grid if it does not contain one.
        if not numpy.any(self.volume.data):
            self.logger.info("Calculating charge densities on the provided empty grid.")
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

        # If charge densities are provided beforehand, log this information
        # `Volume` object does not contain (nor rely on) information about the constituent atoms.
        else:
            self.logger.info("Using charge densities from the provided Volume object.")
            self.chgdensity = self.volume

        # Assign each grid point to Bader areas
        self.fragresults = numpy.zeros(self.chgdensity.data.shape, "d")
        next_index = 1

        self.logger.info("Partitioning space into Bader areas.")

        # Generator to iterate over the elements excluding the outermost positions
        xshape, yshape, zshape = self.chgdensity.data.shape
        indices = (
            (x, y, z)
            for x in range(1, xshape - 1)
            for y in range(1, yshape - 1)
            for z in range(1, zshape - 1)
        )

        for xindex, yindex, zindex in indices:
            if self.fragresults[xindex, yindex, zindex] != 0:
                # index has already been assigned for this grid point
                continue
            else:
                listcoord = []
                local_max_reached = False

                while not local_max_reached:
                    # Here, `delta_rho` corresponds to equation 2,
                    # and `grad_rho_dot_r` corresponds to equation 1 in the aforementioned
                    # paper (doi:10.1016/j.commatsci.2005.04.010)
                    delta_rho = (
                        self.chgdensity.data[
                            xindex - 1 : xindex + 2,
                            yindex - 1 : yindex + 2,
                            zindex - 1 : zindex + 2,
                        ]
                        - self.chgdensity.data[xindex, yindex, zindex]
                    )
                    grad_rho_dot_r = delta_rho / _griddist
                    maxat = numpy.where(grad_rho_dot_r == numpy.amax(grad_rho_dot_r))

                    directions = list(zip(maxat[0], maxat[1], maxat[2]))
                    next_direction = [ind - 1 for ind in directions[0]]

                    if len(directions) > 1:
                        # when one or more directions indicate max grad (of 0), prioritize
                        # to include all points in the Bader space
                        if directions[0] == [1, 1, 1]:
                            next_direction = [ind - 1 for ind in directions[1]]

                    listcoord.append((xindex, yindex, zindex))
                    bader_candidate_index = self.fragresults[
                        xindex + next_direction[0],
                        yindex + next_direction[1],
                        zindex + next_direction[2],
                    ]

                    if bader_candidate_index != 0:
                        # Path arrived at a point that has already been assigned with an index
                        bader_index = bader_candidate_index
                        listcoord = tuple(numpy.array(listcoord).T)
                        self.fragresults[listcoord] = bader_index

                        local_max_reached = True

                    elif (
                        next_direction == [0, 0, 0]
                        or xindex + next_direction[0] == 0
                        or xindex + next_direction[0] == (len(self.chgdensity.data) - 1)
                        or yindex + next_direction[1] == 0
                        or yindex + next_direction[1] == (len(self.chgdensity.data[0]) - 1)
                        or zindex + next_direction[2] == 0
                        or zindex + next_direction[2] == (len(self.chgdensity.data[0][0]) - 1)
                    ):
                        # When next_direction is [0, 0, 0] -- local maximum
                        # Other conditions indicate that the path is heading out to edge of
                        # the grid. Here, assign new Bader space to avoid exiting the grid.
                        bader_index = next_index
                        next_index += 1

                        listcoord = tuple(numpy.array(listcoord).T)
                        self.fragresults[listcoord] = bader_index

                        local_max_reached = True

                    else:
                        # Advance to the next point according to the direction of
                        # maximum gradient
                        xindex += next_direction[0]
                        yindex += next_direction[1]
                        zindex += next_direction[2]

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
        ), f"Failed to assign Bader regions to atoms. Try with a finer grid. Content of Bader area matches: {self.matches}"
        assert len(
            numpy.unique(self.matches) != len(self.data.atomnos)
        ), "Failed to assign unique Bader regions to each atom. Try with a finer grid."

        # Finally integrate the assigned Bader areas
        self.logger.info("Creating fragcharges: array[1]")
        self.fragcharges = numpy.zeros(len(self.data.atomcoords[-1]), "d")

        for atom_index, baderarea_index in enumerate(self.matches):
            self.fragcharges[atom_index] = self.chgdensity.integrate(
                weights=(self.fragresults == baderarea_index)
            )

        return True
