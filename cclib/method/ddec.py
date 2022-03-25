# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculation of DDEC charges based on data parsed by cclib."""
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
from cclib.parser.utils import convertor
from cclib.parser.utils import find_package

from typing import List


class MissingInputError(Exception):
    pass


class ConvergenceError(Exception):
    pass


class DDEC6(Stockholder):
    """DDEC6 charges."""

    # All of these are required for DDEC6 charges.
    required_attrs = ("homos", "mocoeffs", "nbasis", "gbasis")

    def __init__(
        self,
        data,
        volume,
        proatom_path=None,
        progress=None,
        convergence_level=1e-10,
        max_iteration=50,
        loglevel=logging.INFO,
        logname="Log",
    ):
        """ Initialize DDEC6 object.
            Inputs are:
                data -- ccData object that describe target molecule.
                volume -- Volume object that describe target Cartesian grid.
                proatom_path -- path to proatom densities
                (directory containing atoms.h5 in horton or c2_001_001_000_400_075.txt in chargemol)
                convergence_level -- convergence level to use for conditioning densities in step 3
                max_iteration -- maximum iteration to optimize phi in step 3-6
            
            Note:
                Proatom densities are used in DDEC6 algorithm in a similar way to other stockholder
                partitioning methods. They are used as references to appropriately partition the
                total densities (which is the density stored in cube files). Proatom densities are
                densities obtained for single atom or ion in a radial grid that originates from
                the atom or ion.
                In DDEC6 algorithm, stockholder partitioning is heavily modified to ensure that the
                total densities that are partitioned resemble the proatom densities and to prevent
                the numerical algorithm from failing to converge.
        """
        super().__init__(data, volume, proatom_path, progress, loglevel, logname)

        self.convergence_level = convergence_level
        self.max_iteration = max_iteration

        if numpy.sum(self.data.coreelectrons) != 0:
            # TODO: Pseudopotentials should be added back
            pass

    def __str__(self):
        """Return a string representation of the object."""
        return f"DDEC6 charges of {self.data}"

    def __repr__(self):
        """Return a representation of the object."""
        return f"DDEC6({self.data})"

    def _check_required_attributes(self):
        super()._check_required_attributes()

    def _cartesian_dist(self, pt1, pt2):
        """ Small utility function that calculates Euclidian distance between two points
            pt1 and pt2 are numpy arrays representing a point in Cartesian coordinates. """
        return numpy.sqrt(numpy.dot(pt1 - pt2, pt1 - pt2))

    def _read_proatom(
        self, directory, atom_num, charge  # type = str  # type = int  # type = float
    ):
        return super()._read_proatom(directory, atom_num, charge)

    def calculate(self, indices=None, fupdate=0.05):
        """
        Calculate DDEC6 charges based on doi: 10.1039/c6ra04656h paper.
        Cartesian, uniformly spaced grids are assumed for this function.
        """
        super().calculate()

        # Notify user about the total charge in the density grid
        integrated_density = self.charge_density.integrate()
        self.logger.info(
            f"Total charge density in the grid is {integrated_density}. If this does not match what is expected, using a finer grid may help."
        )


        # * STEP 1 *
        # Carry out step 1 of DDEC6 algorithm [Determining reference charge value]
        # Refer to equations 49-57 in doi: 10.1039/c6ra04656h
        self.logger.info("Creating first reference charges. (Step 1/7)")
        (
            reference_charges,
            localized_charges,
            stockholder_charges,
        ) = self.calculate_reference_charges()
        self.reference_charges = [reference_charges]
        self._localized_charges = [localized_charges]
        self._stockholder_charges = [stockholder_charges]

        # * STEP 2 *
        # Load new proatom densities based on the reference charges determined in step 1.
        self.logger.info("Creating second reference charges. (Step 2/7)")
        self.proatom_density = []
        self.radial_grid_r = []
        for i, atom_number in enumerate(self.data.atomnos):
            density, r = self._read_proatom(
                self.proatom_path, atom_number, float(self.reference_charges[0][i])
            )
            self.proatom_density.append(density)
            self.radial_grid_r.append(r)

        # Carry out step 2 of DDEC6 algorithm [Determining ion charge value again]
        ref, loc, stock = self.calculate_reference_charges()
        self.reference_charges.append(ref)
        self._localized_charges.append(loc)
        self._stockholder_charges.append(stock)

        # * STEP 3 *
        # Load new proatom densities based on the reference charges determined in step 2.
        self.proatom_density = []
        self.radial_grid_r = []
        self._cond_density = []
        for i, atom_number in enumerate(self.data.atomnos):
            density, r = self._read_proatom(
                self.proatom_path, atom_number, float(self.reference_charges[1][i])
            )
            self.proatom_density.append(density)
            self.radial_grid_r.append(r)

        # Carry out step 3 of DDEC6 algorithm [Determine conditioned charge density and tau]
        self.logger.info("Conditioning charge densities. (Step 3/7)")
        self.condition_densities()

        # Steps 4 through 7 contain similar routines. Comments the precede each step explain the
        # differences among them.

        # * STEP 4 *
        # In step 4, calculate w, u, g, and h but always skip kappa updates.
        self.logger.info("Optimizing grid weights. (Step 4/7)")
        self._kappa = [0.0] * self.data.atomnos
        # N_A is assigned number of electrons determined using equation 72.
        self.N_A = []
        # u_A is determined using equation 77.
        self.u_A = []

        self.N_A.append(self._calculate_w_and_u())

        # Calculate G_A and H_A based on S4.3 in doi: 10.1039/c6ra04656h
        self.reshape_G()
        self.calculate_H()

        # Update weights (w_A) using equation 96 in doi: 10.1039/c6ra04656h
        # self._cond_density is first created in step 3 by conditioning on the total densities
        # as described in Figure S1. Then, this quantity is updated in every step that follows
        # until the last step in the algorithm when the weights placed on the grid is iteratively
        # updated.
        for atomi in range(self.data.natom):
            self._cond_density[atomi] = math.exp(self._kappa[atomi]) * self._h[atomi]

        # Update rho_cond for next iteration
        self._update_rho_cond()

        # * STEPS 5 and 6 *
        # In step 5 and 6, calculate w and u. Then if update_kappa is found to be true, do not
        # continue to next step. Otherwise, calculate g and h. In both cases, calculate new
        # rho_cond.
        steps = 5
        self._update_kappa = False
        while steps < 7:
            self.logger.info(f"Optimizing grid weights. (Step {steps}/7)")
            self.N_A.append(self._calculate_w_and_u())

            # Determine whether kappa needs to be updated or not based on Figure S4.2
            # of doi: 10.1039/c6ra04656h
            self._update_kappa = self._check_kappa()

            if not self._update_kappa:
                # Increment steps
                steps = steps + 1
                # Calculate G_A and H_A based on S4.3 in doi: 10.1039/c6ra04656h
                self.reshape_G()
                self.calculate_H()

            else:
                # `steps` not incremented in this case
                # First update kappa based on equation 93
                kappa_new = self._kappa - self.N_A[-1] / self.u_A[-1]
                self._kappa = [x if x > 0 else 0.0 for x in kappa_new]

            # Update weights (w_A) using equation 96 in doi: 10.1039/c6ra04656h
            # self._cond_density is first created in step 3 by conditioning on the total densities
            # as described in Figure S1. Then, this quantity is updated in every step that follows
            # until the last step in the algorithm when the weights placed on the grid is
            # iteratively updated.
            for atomi in range(self.data.natom):
                self._cond_density[atomi] = math.exp(self._kappa[atomi]) * self._h[atomi]

            # Update rho_cond for next iteration
            self._update_rho_cond()

        # * STEP 7 *
        # In step 7, calculate w and u. Then, stop to calculate reference_charges.
        self.logger.info("Optimizing grid weights. (Step 7/7)")
        self.N_A.append(self._calculate_w_and_u())

        # Finally, store calculated DDEC6 charges in fragcharges
        self.logger.info("Creating fragcharges: array[1]")
        self.fragcharges = numpy.array(self.data.atomnos - self.N_A[-1], dtype=float)

    def _check_kappa(self):
        """ Return whether kappa needs to be updated or not based on Figure S4.2
            of doi: 10.1039/c6ra04656h
        """
        if numpy.any([x if x < -1e-5 else 0 for x in self.N_A[-1]]):
            return True
        elif (
            self._update_kappa
            and numpy.any(numpy.diff(self.N_A)[-1] < 1e-5)  # change in N_A during last cycle
            and numpy.any(numpy.diff(self.N_A)[-2] < 1e-5)  # change in N_A during last to
            # second cycle
        ):
            self._kappa = [0.0 for x in self.data.atomnos]
            return False
        else:
            return self._update_kappa

    def calculate_reference_charges(self):
        """ Calculate reference charges from proatom density and molecular density
            [STEP 1 and 2]
            
            Function returns calculated reference charges, localized charges, and stockholder
            charges.
        """
        # Generator object to iterate over the grid
        ngridx, ngridy, ngridz = self.charge_density.data.shape
        grid_shape = (self.data.natom, ngridx, ngridy, ngridz)

        stockholder_w = numpy.zeros(grid_shape)
        localized_w = numpy.zeros(grid_shape)
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

        # Equation 55 in doi: 10.1039/c6ra04656h
        localized_w = numpy.power(stockholder_w, 4)

        # Equation 53 in doi: 10.1039/c6ra04656h
        stockholder_bigW = numpy.sum(stockholder_w, axis=0)
        localized_bigW = numpy.sum(localized_w, axis=0)

        reference_charges = numpy.zeros((self.data.natom))
        localizedcharges = numpy.zeros((self.data.natom))
        stockholdercharges = numpy.zeros((self.data.natom))

        for atomi in range(self.data.natom):
            # Equation 52 and 51 in doi: 10.1039/c6ra04656h
            localizedcharges[atomi] = self.data.atomnos[atomi] - self.charge_density.integrate(
                weights=(localized_w[atomi] / localized_bigW)
            )
            stockholdercharges[atomi] = self.data.atomnos[atomi] - self.charge_density.integrate(
                weights=(stockholder_w[atomi] / stockholder_bigW)
            )

            # In DDEC6, weights of 1/3 and 2/3 are assigned for stockholder and localized charges.
            # (Equation 50 and 58 in doi: 10.1039/c6ra04656h)
            reference_charges[atomi] = (stockholdercharges[atomi] / 3.0) + (
                localizedcharges[atomi] * 2.0 / 3.0
            )

        return reference_charges, localizedcharges, stockholdercharges

    def condition_densities(self):
        """ Calculate conditioned densities
            [STEP 3]
        """
        # Generator object to iterate over the grid
        ngridx, ngridy, ngridz = self.charge_density.data.shape

        self._rho_ref = numpy.zeros((ngridx, ngridy, ngridz))

        for atomi in range(self.data.natom):
            # rho_ref -- Equation 41 in doi: 10.1039/c6ra04656h
            self._rho_ref += self.proatom_density[atomi][
                self.closest_r_index[atomi]
            ]

        self._candidates_bigPhi = []
        self._candidates_phi = []

        # Initial conditions are detailed in Figure S1 in doi: 10.1039/c6ra04656h
        phiAI = numpy.zeros_like(self.data.atomnos, dtype=float)
        bigphiAI = numpy.zeros_like(self.data.atomnos, dtype=float)
        self._y_a = []

        for atomi in range(self.data.natom):
            # y_a -- equation 40 in doi: 10.1039/c6ra04656h
            self._y_a.append(self._ya(self.proatom_density, atomi))
            # rho_A^cond -- equation S97 in doi: 10.1039/c6ra04656h
            self._cond_density.append(
                self._y_a[atomi] + bigphiAI[atomi] * numpy.sqrt(self._y_a[atomi])
            )
            # Monotonic Decrease Condition (as expressed in equation S99)
            self._cond_density[atomi] = numpy.minimum.accumulate(self._cond_density[atomi])
            # phi_A^I -- Equation S100 in doi: 10.1039/c6ra04656h
            phiAI[atomi] = (
                self._integrate_from_radial([self._cond_density[atomi]], [atomi])
                - self.data.atomnos[atomi]
                + self.reference_charges[-1][atomi]
            )

            self._candidates_bigPhi.append([bigphiAI[atomi]])
            self._candidates_phi.append([phiAI[atomi]])

            # Attempt to find the point where phiAI is zero iteratively
            # Refer to S101 in doi: 10.1039/c6ra04656h
            self._candidates_phi[atomi], self._candidates_bigPhi[atomi] = self._converge_phi(
                phiAI[atomi], 1, atomi
            )

            # Perform parabolic fit to find optimized phiAI
            # Refer to Figure S1 in doi: 10.1039/c6ra04656h
            bigphiAI[atomi] = self._parabolic_fit(self._y_a[atomi], 1, atomi)

            # Set final conditioned density using chosen Phi
            self._cond_density[atomi] = self._update_phiai(
                self._y_a[atomi], bigphiAI[atomi], atomi
            )[1]

        self.logger.info("Calculating tau and combined conditioned densities.")

        # Calculate tau(r) and rho^cond(r)
        # Refer to equation 65 and 66 in doi: 10.1039/c6ra04656h
        # Assign rho^cond on grid using generator object
        self.rho_cond = copy.deepcopy(self.charge_density)
        self.rho_cond.data = numpy.zeros_like(self.rho_cond.data, dtype=float)
        rho_cond_sqrt = numpy.zeros_like(self.rho_cond.data, dtype=float)

        # Generator object to iterate over the grid
        ngridx, ngridy, ngridz = self.charge_density.data.shape

        self._leftterm = numpy.zeros((self.data.natom, ngridx, ngridy, ngridz), dtype=float)
        # rho_cond_cartesian is rho^cond projected on Cartesian grid
        # (used for Step 4 calculation)
        self._rho_cond_cartesian = numpy.zeros(
            (self.data.natom, ngridx, ngridy, ngridz), dtype=float
        )
        self.tau = []

        # rho_cond -- equation 65 in doi: 10.1039/c6ra04656h
        for atomi in range(self.data.natom):
            self.rho_cond.data += self._cond_density[atomi][
                self.closest_r_index[atomi]
            ]

        rho_cond_sqrt = numpy.sqrt(self.rho_cond.data)

        for atomi in range(self.data.natom):
            self.tau.append(numpy.zeros_like(self.proatom_density[atomi], dtype=float))
            # leftterm is the first spherical average term in equation 66.
            # <rho^cond_A(r_A) / sqrt(rho^cond(r))>
            self._rho_cond_cartesian[atomi] = self._cond_density[atomi][
                self.closest_r_index[atomi]
            ]
            self._leftterm[atomi] = self._rho_cond_cartesian[atomi] / rho_cond_sqrt
            for radiusi in range(len(self.tau[atomi])):
                grid_filter = self.closest_r_index[atomi] == radiusi
                num_grid_filter = numpy.count_nonzero(grid_filter)
                if num_grid_filter < 1:
                    self.tau[atomi][radiusi] = 0.0
                else:
                    leftaverage = numpy.sum(grid_filter * self._leftterm[atomi]) / num_grid_filter
                    rightaverage = numpy.sum(grid_filter * rho_cond_sqrt) / num_grid_filter
                    if leftaverage < 1e-20:
                        self.tau[atomi][radiusi] = 0.0
                    else:
                        self.tau[atomi][radiusi] = numpy.divide(
                            leftaverage,
                            rightaverage,
                            out=numpy.zeros_like(leftaverage),
                            where=rightaverage != 0.0,
                        )
            # Make tau monotonic decreasing
            self.tau[atomi] = numpy.maximum.accumulate(self.tau[atomi][::-1])[::-1]

    def _ya(self, proatom_density, atomi):
        # Function that calculates Y_a^avg
        # See Eq. 40-41 in doi: 10.1039/c6ra04656h
        rho_ref = self._rho_ref
        # Y_a^avg -- Equation 40 in doi: 10.1039/c6ra04656h
        ya = numpy.zeros_like(proatom_density[atomi], dtype=float)
        weights = self.charge_density.data / rho_ref

        for radiusi in range(len(ya)):
            grid_filter = self.closest_r_index[atomi] == radiusi
            num_grid_filter = numpy.count_nonzero(grid_filter)
            if num_grid_filter < 1:
                ya[radiusi] = 0.0
            else:
                spherical_avg = numpy.sum(grid_filter * weights) / num_grid_filter
                ya[radiusi] = proatom_density[atomi][radiusi] * spherical_avg

        # Make y_a monotonic decreasing
        # Refer to module_reshaping_functions.f08::77-79
        ya = numpy.maximum.accumulate(ya[::-1])[::-1]
        ya = numpy.minimum.accumulate(ya)

        # Normalize y_a (see module_DDEC6_valence_iterator.f08::284)
        nelec = self._integrate_from_radial([ya], [atomi])
        ya *= (self.data.atomnos[atomi] - self.reference_charges[-1][atomi]) / nelec

        return ya

    def _calculate_w_and_u(self):
        """ Calculate weights placed on each integration grid point
            [STEP 4-7]
        """
        # From equation 67, w_A(r_A) = self._cond_density
        # From equation 8, W(r) = self.rho_cond for the rest of this function
        ngridx, ngridy, ngridz = self.charge_density.data.shape

        # Evaluate rho_A(r_A) from equation 68, rho_A^avg(r_A) from equation 69,
        # theta(r_A) from equation 70, _wavg from equation 71,
        # N_A from equation 72, rho_wavg from equation 73, and u_A from equation 77.
        self._rho_A = []
        self._rho_A_avg = []
        self._theta = []
        self._wavg = []
        N_A = []
        u_A = []
        self.rho_wavg = []

        for atomi in range(self.data.natom):
            self._rho_A.append(copy.deepcopy(self.charge_density))
            # Equation 68
            self._rho_A[atomi].data = numpy.divide(
                self.charge_density.data * self._rho_cond_cartesian[atomi],
                self.rho_cond.data,
                out=numpy.zeros_like(self.charge_density.data, dtype=float),
                where=self.rho_cond.data != 0,
            )
            self._rho_A_avg.append(numpy.zeros_like(self.proatom_density[atomi], dtype=float))
            self._theta.append(numpy.zeros_like(self.proatom_density[atomi], dtype=float))
            self._wavg.append(numpy.zeros_like(self.proatom_density[atomi], dtype=float))
            # Equation 69, 70, and 71
            self._rho_A_avg[atomi] = self._spherical_average_from_cartesian(
                self._rho_A[atomi].data, atomi, self.radial_grid_r[atomi]
            )
            self._rho_A_avg[atomi] = numpy.maximum.accumulate(self._rho_A_avg[atomi][::-1])[::-1]
            self._theta[atomi] = self._spherical_average_from_cartesian(
                (
                    1
                    - numpy.divide(
                        self._rho_cond_cartesian[atomi],
                        self.rho_cond.data,
                        out=numpy.zeros_like(self.rho_cond.data, dtype=float),
                        where=self.rho_cond.data != 0,
                    )
                )
                * self._rho_A[atomi].data,
                atomi,
                self.radial_grid_r[atomi],
            )
            self._theta[atomi] = numpy.maximum.accumulate(self._theta[atomi][::-1])[::-1]
            self._wavg[atomi] = self._spherical_average_from_cartesian(
                numpy.divide(
                    self._rho_cond_cartesian[atomi],
                    self.rho_cond.data,
                    out=numpy.zeros_like(self.rho_cond.data, dtype=float),
                    where=self.rho_cond.data != 0,
                ),
                atomi,
                self.radial_grid_r[atomi],
            )
            self._wavg[atomi] = numpy.maximum.accumulate(self._wavg[atomi][::-1])[::-1]
            # Equation 72, 73, and 77
            N_A.append(self._rho_A[atomi].integrate())
            self.rho_wavg.append(
                (self._theta[atomi] + self._rho_A_avg[atomi] * self._wavg[atomi] / 5.0)
                / (1.0 - (4.0 / 5.0) * self._wavg[atomi])
            )
            u_A.append(self._integrate_from_radial([self._theta[atomi]], [atomi]))
        self.u_A.append(u_A)

        return N_A

    def reshape_G(self):
        """ Calculate G_A(r_A) and reshape densities
        
            This is a quantity introduced in DDEC6 as a constraint preventing the tails from being
            too diffuse.
            [STEP 4-7]
        """
        self._candidates_bigPhi = []
        self._candidates_phi = []

        # Initial conditions are detailed in Figure S3 in doi: 10.1039/c6ra04656h
        phiAII = numpy.zeros_like(self.data.atomnos, dtype=float)
        bigphiAII = numpy.zeros_like(self.data.atomnos, dtype=float)
        self._g = []
        self._eta = []

        for atomi in range(self.data.natom):
            # G_A -- equation S102 in doi: 10.1039/c6ra04656h
            self._g.append(numpy.zeros_like(self.proatom_density[atomi], dtype=float))
            self._g[atomi] = self.rho_wavg[atomi] + bigphiAII[atomi] * numpy.sqrt(
                self.rho_wavg[atomi]
            )
            # Exponential constraint (as expressed in equation S105)
            self._eta.append((1 - (self.tau[atomi]) ** 2) * 1.75 * convertor(1, "Angstrom", "bohr"))
            exp_applied = self._g[atomi][:-1] * numpy.exp(
                -1 * self._eta[atomi][1:] * numpy.diff(self.radial_grid_r[atomi])
            )
            for radiusi in range(1, len(self._g[atomi])):
                self._g[atomi][radiusi] = min(self._g[atomi][radiusi], exp_applied[radiusi - 1])
            # phi_A^II -- Equation S106 in doi: 10.1039/c6ra04656h
            phiAII[atomi] = self._integrate_from_radial(
                [self._g[atomi] - self.rho_wavg[atomi]], [atomi]
            )

            self._candidates_bigPhi.append([bigphiAII[atomi]])
            self._candidates_phi.append([phiAII[atomi]])

            # Attempt to find the point where phiAI is zero iteratively
            # Refer to S101 in doi: 10.1039/c6ra04656h
            self._candidates_phi[atomi], self._candidates_bigPhi[atomi] = self._converge_phi(
                phiAII[atomi], 1, atomi
            )

            # Perform parabolic fit to find optimized phiAI
            # Refer to Figure S1 in doi: 10.1039/c6ra04656h
            bigphiAII[atomi] = self._parabolic_fit(self.rho_wavg[atomi], 2, atomi)

            # Set final G_A value using chosen Phi
            self._g[atomi] = self._update_phiaii(self.rho_wavg[atomi], bigphiAII[atomi], atomi)[1]

    def calculate_H(self):
        """ Calculate H_A(r_A)
            
            This is a quantity introduced in DDEC6 as a constraint preventing the tails from being
            too contracted.
            [STEP 4-7]
        """
        self._h = []
        for atomi in range(len(self._g)):
            # First set H_est as G_A
            self._h.append(self._g[atomi])

            # Determine eta_upper using equation 86 in doi: 10.1039/c6ra04656h
            # and apply upper limit using equation 91.
            temp = (
                1 - (self.tau[atomi]) ** 2 + self.convergence_level
            )  # convergence_level is added to avoid divide-by-zero in next line for highly polar molecules.
            eta = 2.5 * convertor(1, "Angstrom", "bohr") / temp
            exp_applied = self._h[atomi][:-1] * numpy.exp(
                -1 * eta[1:] * numpy.diff(self.radial_grid_r[atomi])
            )
            for radiusi in range(1, len(self._h[atomi])):
                self._h[atomi][radiusi] = max(self._h[atomi][radiusi], exp_applied[radiusi - 1])

            # Normalize using equation 92 in doi: 10.1039/c6ra04656h.
            self._h[atomi] = (
                self._h[atomi]
                * self._integrate_from_radial([self._g[atomi]], [atomi])
                / self._integrate_from_radial([self._h[atomi]], [atomi])
            )

    def _update_phiai(self, ya, bigphiAI, atomi):
        # Update phi^a_i and quantity that directly follows (cond_density) in each step of
        # iterative optimization. (Refer to Figure S1)
        # Re-evaluate cond_density
        if isinstance(bigphiAI, float) and not numpy.isinf(bigphiAI):
            cond_density = ya + bigphiAI * numpy.sqrt(ya)
        else:
            cond_density = ya

        # Monotonic Decrease Condition
        cond_density = numpy.minimum.accumulate(cond_density)

        # Re-evaluate phi_AI
        phiAI = (
            self._integrate_from_radial([cond_density], [atomi])
            - self.data.atomnos[atomi]
            + self.reference_charges[-1][atomi]
        )

        return phiAI, cond_density

    def _update_phiaii(self, rhowavg, bigphiAI, atomi):
        # Update phi^a_ii and quantity that directly follows (G_A) in each step of
        # iterative optimization. (Refer to Figure S3)
        # Re-evaluate g_a
        # Equations can be found in Figure S3.
        ga = rhowavg + bigphiAI * numpy.sqrt(rhowavg)

        # Exponential Decrease Condition
        exp_applied = ga[:-1] * numpy.exp(
            -1 * self._eta[atomi][1:] * numpy.diff(self.radial_grid_r[atomi])
        )
        for radiusi in range(1, len(ga)):
            ga[radiusi] = min(ga[radiusi], exp_applied[radiusi - 1])

        # Re-evaluate phi_AII
        phiAII = self._integrate_from_radial([ga - rhowavg], [atomi])

        return phiAII, ga

    def _integrate_from_radial(self, radial_density_list, atom_list):
        # Function that reads in list of radial densities, projects it on Cartesian grid,
        # and returns integrated value
        grid = copy.deepcopy(self.charge_density)
        grid.data = numpy.zeros_like(grid.data, dtype=float)

        for density, atomi in zip(radial_density_list, atom_list):
            grid.data += density[self.closest_r_index[atomi]]

        return grid.integrate()

    def _spherical_average_from_cartesian(self, cartesian_grid, atom_index, radius_list):
        spherical_average = numpy.zeros(len(radius_list))
        for radiusi in range(len(radius_list)):
            grid_filter = self.closest_r_index[atom_index] == radiusi
            num_grid_filter = numpy.count_nonzero(grid_filter)
            if num_grid_filter < 1:
                average = 0.0
            else:
                average = numpy.sum(grid_filter * cartesian_grid) / num_grid_filter
                if average < self.convergence_level:
                    average = 0.0
            spherical_average[radiusi] = average
        return spherical_average

    def _update_rho_cond(self):
        # Update total weights on Cartesian grid using equation 65 in doi: 10.1039/c6ra04656h
        ngridx, ngridy, ngridz = self.charge_density.data.shape

        self.rho_cond.data = numpy.zeros_like(self.rho_cond.data, dtype=float)
        self._rho_cond_cartesian = numpy.zeros(
            (self.data.natom, ngridx, ngridy, ngridz), dtype=float
        )

        for atomi in range(self.data.natom):
            self.rho_cond.data += self._cond_density[atomi][
                self.closest_r_index[atomi]
            ]
            self._rho_cond_cartesian[atomi] = self._cond_density[atomi][
                self.closest_r_index[atomi]
            ]

    def _converge_phi(self, phiA, superscript, atomi):
        """ Update phi until it is positive.
            This is used in step 3 (for phi_A^I) and in steps 4-6 (for phi_A^II).
            
            --- Inputs ---
            phiA            Either phi_A^I or phi_A^II
            superscript     1 when calculating phi_I (STEP 3)
                            2 when calculating phi_II (STEPS 4-6)
            atomi           Index of target atom as in ccData object
            
            Refer to Equation S101, Figure S1 and S3 for an overview.
        """
        # Initial value of bigphi is zero (Equation S100)
        bigphiA = 0.0

        # List to store candidate values for parabolic fitting
        candidates_phi = [phiA]
        candidates_bigphi = [bigphiA]

        while phiA <= 0:
            # Iterative algorithm until convergence
            # Refer to S101 in doi: 10.1039/c6ra04656h
            if superscript == 1:
                temp = self._integrate_from_radial([numpy.sqrt(self._y_a[atomi])], [atomi])
            elif superscript == 2:
                temp = self._integrate_from_radial([numpy.sqrt(self.rho_wavg[atomi])], [atomi])

            bigphiA = 2 * bigphiA - phiA / temp

            # When Phi is updated, related quantities are updated as well
            # Refer to S100 in doi: 10.1039/c6ra04656h [Step 3]
            #       or S102 in doi: 10.1039/c6ra04656h [Steps 4-6]

            if superscript == 1:
                phiA, self._cond_density[atomi] = self._update_phiai(
                    self._y_a[atomi], bigphiA, atomi
                )
            elif superscript == 2:
                phiA, self._g[atomi] = self._update_phiaii(self.rho_wavg[atomi], bigphiA, atomi)

            candidates_phi.append(phiA)
            candidates_bigphi.append(bigphiA)

        return candidates_phi, candidates_bigphi

    def _parabolic_fit(self, pseudodensity, superscript, atomi):
        """ Optimize phi using parabolic fitting.
            This is used in step 3 (for phi_A^I) and in steps 4-6 (for phi_A^II).
            
            --- Inputs ---
            phiA            Either phi_A^I or phi_A^II
            superscript     1 when calculating phi_I (STEP 3)
                            2 when calculating phi_II (STEPS 4-6)
            atomi           Index of target atom as in ccData object
            
            Refer to Figure S1 and S3 for an overview.
        """
        # Set update methods for phi_A^I or phi_A^II
        if superscript == 1:

            def update(pdens, bigPhi, atomi):
                return self._update_phiai(pdens, bigPhi, atomi)

        elif superscript == 2:

            def update(pdens, bigPhi, atomi):
                return self._update_phiaii(pseudodensity, bigPhi, atomi)

        # lowerbigPhi is bigPhi that yields biggest negative phi.
        # upperbigPhi is bigPhi that yields smallest positive phi.
        # The point here is to find two phi values that are closest to zero (from positive side
        # and negative side respectively).
        self._candidates_phi[atomi] = numpy.array(self._candidates_phi[atomi], dtype=float)
        self._candidates_bigPhi[atomi] = numpy.array(self._candidates_bigPhi[atomi], dtype=float)
        if numpy.count_nonzero(self._candidates_phi[atomi] < 0) > 0:
            # If there is at least one candidate phi that is negative
            lower_ind = numpy.where(
                self._candidates_phi[atomi]
                == self._candidates_phi[atomi][self._candidates_phi[atomi] < 0].max()
            )[0][0]
            lowerbigPhi = self._candidates_bigPhi[atomi][lower_ind]
            lowerphi = self._candidates_phi[atomi][lower_ind]
        else:  # assign some large negative number otherwise
            lowerbigPhi = numpy.NINF
            lowerphi = numpy.NINF
        if numpy.count_nonzero(self._candidates_phi[atomi] > 0) > 0:
            # If there is at least one candidate phi that is positive
            upper_ind = numpy.where(
                self._candidates_phi[atomi]
                == self._candidates_phi[atomi][self._candidates_phi[atomi] > 0].min()
            )[0][0]
            upperbigPhi = self._candidates_bigPhi[atomi][upper_ind]
            upperphi = self._candidates_phi[atomi][upper_ind]
        else:  # assign some large positive number otherwise
            upperbigPhi = numpy.PINF
            upperphi = numpy.PINF

        for iteration in range(self.max_iteration):
            # Flow diagram on Figure S1 in doi: 10.1039/c6ra04656h details the procedure.
            # Find midpoint between positive bigPhi that yields phi closest to zero and negative
            # bigPhi closest to zero. Then, evaluate phi.
            # This can be thought as linear fitting compared to parabolic fitting below.
            midbigPhi = (lowerbigPhi + upperbigPhi) / 2.0
            midphi = update(pseudodensity, midbigPhi, atomi)[0]
            # Exit conditions -- if any of three phi values are within the convergence level.
            if abs(lowerphi) < self.convergence_level:
                return lowerbigPhi
            elif abs(upperphi) < self.convergence_level:
                return upperbigPhi
            elif abs(midphi) < self.convergence_level:
                return midbigPhi

            # Parabolic fitting as described on Figure S1 in doi: 10.1039/c6ra04656h
            # Type casting here converts from size 1 numpy.ndarray to float
            xpts = numpy.array(
                [float(lowerbigPhi), float(midbigPhi), float(upperbigPhi)], dtype=float
            )
            ypts = numpy.array([float(lowerphi), float(midphi), float(upperphi)], dtype=float)
            fit = numpy.polyfit(xpts, ypts, 2)
            roots = numpy.roots(fit)  # max two roots (bigPhi) from parabolic fitting

            # Find phi for two bigPhis that were obtained from parabolic fitting.
            belowphi = update(pseudodensity, roots.min(), atomi)[0]
            abovephi = update(pseudodensity, roots.max(), atomi)[0]

            # If phi values from parabolically fitted bigPhis lie within the convergence level,
            # exit the iterative algorithm.
            if abs(abovephi) < self.convergence_level:
                return roots.min()
            elif abs(belowphi) < self.convergence_level:
                return roots.max()
            else:
                # Otherwise, corrected phi value is obtained in a way that cuts the numerical
                # search domain in half in each iteration.
                if 3 * abs(abovephi) < abs(belowphi):
                    corbigPhi = roots.max() - 2.0 * abovephi * (roots.max() - roots.min()) / (
                        abovephi - belowphi
                    )
                elif 3 * abs(belowphi) < abs(abovephi):
                    corbigPhi = roots.min() - 2.0 * belowphi * (roots.max() - roots.min()) / (
                        abovephi - belowphi
                    )
                else:
                    corbigPhi = (roots.max() + roots.min()) / 2.0
                # New candidates of phi and bigPhi are determined as bigPhi yielding largest
                # negative phi and bigPhi yielding smallest positve phi. This is analogous to how
                # the first candidiate phi values are evaluated.
                corphi = update(pseudodensity, corbigPhi, atomi)[0]
                self._candidates_bigPhi[atomi] = numpy.array(
                    [lowerbigPhi, midbigPhi, upperbigPhi, roots.max(), roots.min(), corbigPhi,],
                    dtype=float,
                )
                self._candidates_phi[atomi] = numpy.array(
                    [lowerphi, midphi, upperphi, abovephi, belowphi, corphi], dtype=float
                )

                # Set new upperphi and lowerphi
                lower_ind = numpy.where(
                    self._candidates_phi[atomi]
                    == self._candidates_phi[atomi][self._candidates_phi[atomi] < 0].max()
                )[0][0]
                upper_ind = numpy.where(
                    self._candidates_phi[atomi]
                    == self._candidates_phi[atomi][self._candidates_phi[atomi] > 0].min()
                )[0][0]

                lowerphi = self._candidates_phi[atomi][lower_ind]
                upperphi = self._candidates_phi[atomi][upper_ind]

                # If new lowerphi or upperphi values are within convergence level, exit the
                # iterative algorithm. Otherwise, start new linear/parabolic fitting.
                if abs(lowerphi) < self.convergence_level:
                    return self._candidates_bigPhi[atomi][lower_ind]
                elif abs(upperphi) < self.convergence_level:
                    return self._candidates_bigPhi[atomi][upper_ind]
                else:
                    # Fitting needs to continue in this case.
                    lowerbigPhi = self._candidates_bigPhi[atomi][lower_ind]
                    lowerphi = self._candidates_phi[atomi][lower_ind]
                    upperbigPhi = self._candidates_bigPhi[atomi][upper_ind]
                    upperphi = self._candidates_phi[atomi][upper_ind]

        # Raise Exception if convergence is not achieved within max_iteration.
        raise ConvergenceError("Iterative conditioning failed to converge.")
