# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Iterative Stockholder partitioning based on cclib data."""
import random
import numpy

from cclib.method.stockholder import Stockholder

class IterativeStockholder(Stockholder):
    """Iterative Stockholder Method"""

    def __init__(self, *args):

        # Call the __init__ method of the superclass.
        super(IterativeStockholder, self).__init__(logname="Iterative Stockholder", *args)

    def calculate(self, indices=None, fupdate=0.05, convergence_level = 1e-5):
        """
        Calculate Hirshfeld-I charges based on doi: 10.1063/1.2715563 paper.
        Cartesian, uniformly spaced grids are assumed for this function.
        """

        super(IterativeStockholder, self).calculate()
        
        # Create initial reference charges
        self.reference_charges = []
        self.reference_charges.append(self.calculate_reference_charges())
        
        # Load new proatom densities and create second reference charges.
        self.proatom_density = []
        self.radial_grid_r = []
        for i, atom_number in enumerate(self.data.atomnos):
            density, r = self._read_proatom(
                self.proatom_path, atom_number, float(self.reference_charges[0][i])
            )
            self.proatom_density.append(density)
            self.radial_grid_r.append(r)
        self.reference_charges.append(self.calculate_reference_charges())

        # Load new proatom densities based on the reference charge,
        # Determine new reference charges until convergence.
        while (numpy.dot(self.referencecharges[-1], self.referencecharges[-2]) > convergence_level):
            self.proatom_density = []
            self.radial_grid_r = []
            for i, atom_number in enumerate(self.data.atomnos):
                density, r = self._read_proatom(
                    self.proatom_path, atom_number, float(self.reference_charges[0][i])
                )
                self.proatom_density.append(density)
                self.radial_grid_r.append(r)
            self.reference_charges.append(self.calculate_reference_charges())

        # Store calculated charges in fragcharges
        self.logger.info("Creating fragcharges: array[1]")
        self.fragcharges = self.reference_charges[-1]

    def calculate_reference_charges(self):
        """ Calculate stockholder charges from proatom density and molecular density
        """
        # Generator object to iterate over the grid
        ngridx, ngridy, ngridz = self.charge_density.data.shape
        grid_shape = (self.data.natom, ngridx, ngridy, ngridz)

        stockholder_w = numpy.zeros(grid_shape)
        self.closest_r_index = numpy.zeros(grid_shape, dtype=int)

        indices = numpy.asanyarray(
            tuple(
                (x, y, z)
                for x in range(ngridx)
                for y in range(ngridy)
                for z in range(ngridz)
            )
        )
        coordinates = self.charge_density.coordinates(indices)

        for atomi in range(self.data.natom):
            # Distance of the grid from atom grid
            self.closest_r_index[atomi] = numpy.argmin(
                numpy.abs(
                    self.radial_grid_r[atomi][..., numpy.newaxis] - numpy.linalg.norm(self.data.atomcoords[-1][atomi] - coordinates, axis=1)
                ),
                axis=0
            ).reshape((ngridx, ngridy, ngridz))

            # Equation 54 in doi: 10.1039/c6ra04656h
            stockholder_w[atomi] = self.proatom_density[atomi][self.closest_r_index[atomi]]

        # Equation 53 in doi: 10.1039/c6ra04656h
        stockholder_bigW = numpy.sum(stockholder_w, axis=0)
        stockholdercharges = numpy.zeros((self.data.natom))

        for atomi in range(self.data.natom):
            stockholdercharges[atomi] = self.data.atomnos[atomi] - self.charge_density.integrate(
                weights=(stockholder_w[atomi] / stockholder_bigW)
            )

        return stockholdercharges


