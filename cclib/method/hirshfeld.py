# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculation of Hirshfeld charges based on data parsed by cclib."""
import copy
import random
import numpy
import logging
import math
import os
import sys

from cclib.method.calculationmethod import Method
from cclib.method.stockholder import Stockholder
from cclib.method.volume import electrondensity_spin
from cclib.parser.utils import find_package

from typing import List


class MissingInputError(Exception):
    pass


class ConvergenceError(Exception):
    pass


class Hirshfeld(Stockholder):
    """Hirshfeld charges."""

    # All of these are required for DDEC6 charges.
    required_attrs = ("homos", "mocoeffs", "nbasis", "gbasis")

    def __init__(
        self,
        data,
        volume,
        proatom_path=None,
        progress=None,
        loglevel=logging.INFO,
        logname="Log",
    ):
        """ Initialize Hirshfeld object.
            Inputs are:
                data -- ccData object that describe target molecule.
                volume -- Volume object that describe target Cartesian grid.
                proatom_path -- path to proatom densities
                (directory containing atoms.h5 in horton or c2_001_001_000_400_075.txt in chargemol)
        """
        super().__init__(data, volume, proatom_path, progress, loglevel, logname)

    def __str__(self):
        """Return a string representation of the object."""
        return f"Hirshfeld charges of {self.data}"

    def __repr__(self):
        """Return a representation of the object."""
        return f"Hirshfeld({self.data})"

    def _check_required_attributes(self):
        super()._check_required_attributes()

    def _cartesian_dist(self, pt1, pt2):
        """Small utility function that calculates Euclidian distance between two points.
        
        Arguments pt1 and pt2 are NumPy arrays representing points in Cartesian coordinates.
        """
        return numpy.sqrt(numpy.dot(pt1 - pt2, pt1 - pt2))

    def _read_proatom(
        self, directory, atom_num, charge  # type = str  # type = int  # type = float
    ):
        return super()._read_proatom(directory, atom_num, charge)

    def calculate(self):
        """Calculate Hirshfeld charges."""
        super().calculate()

        # Generator object to iterate over the grid.
        ngridx, ngridy, ngridz = self.charge_density.data.shape
        indices = (
            (i, x, y, z)
            for i in range(self.data.natom)
            for x in range(ngridx)
            for y in range(ngridy)
            for z in range(ngridz)
        )
        grid_shape = (self.data.natom, ngridx, ngridy, ngridz)

        stockholder_w = numpy.zeros(grid_shape)
        self.closest_r_index = numpy.zeros(grid_shape, dtype=int)

        for atomi, xindex, yindex, zindex in indices:
            # Distance of the grid from atom grid.
            dist_r = self._cartesian_dist(
                self.data.atomcoords[-1][atomi],
                self.charge_density.coordinates([xindex, yindex, zindex]),
            )
            self.closest_r_index[atomi][xindex][yindex][zindex] = numpy.abs(
                self.radial_grid_r[atomi] - dist_r
            ).argmin()

            # stockholder_w is radial proatom density projected on Cartesian grid
            stockholder_w[atomi][xindex][yindex][zindex] = self.proatom_density[atomi][
                self.closest_r_index[atomi][xindex][yindex][zindex]
            ]

        # Equation 3 in doi:10.1007/BF01113058
        # stockholder_bigW represents the denominator in this equation
        # \sum{eta_chi(r_chi)}_chi
        # where eta is proatom density, chi is index of atoms in the molecule, and r_chi is
        # radial distance from center of atom chi.
        stockholder_bigW = numpy.sum(stockholder_w, axis=0)

        self.fragcharges = numpy.zeros((self.data.natom))
        self.logger.info("Creating fragcharges: array[1]")

        for atomi in range(self.data.natom):
            # Equation 1 in doi:10.1007/BF01113058
            # The weights supplied for integrate function correspond to W_A in the equation.
            # Q_A = Z_A - \integral(W_A(r) * q_A(r) dr)
            # where Q_A is Hirshfeld charges, W_A is weights determined on each grid and on each
            # atom, and q_A is total charge densities on the Cartesian grid.
            self.fragcharges[atomi] = self.data.atomnos[atomi] - self.charge_density.integrate(
                weights=(stockholder_w[atomi] / stockholder_bigW)
            )

        return True
