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
from cclib.method.volume import electrondensity_spin
from cclib.parser.utils import convertor
from cclib.parser.utils import find_package

from typing import List


class MissingInputError(Exception):
    pass

class ConvergenceError(Exception):
    pass

class DDEC6(Method):
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
        super(DDEC6, self).__init__(data, progress, loglevel, logname)

        self.volume = volume
        self.fragresults = None
        self.proatom_path = proatom_path
        self.convergence_level = convergence_level
        self.max_iteration = max_iteration

        if numpy.sum(self.data.coreelectrons) != 0:
            # TODO: Pseudopotentials should be added back
            pass

        # Check whether proatom_path is a valid directory or not.
        assert os.path.isdir(
            proatom_path
        ), "Directory that contains proatom densities should be added as an input."

        # Read in reference charges.
        self.proatom_density = []
        self.radial_grid_r = []
        for atom_number in self.data.atomnos:
            density, r = self._read_proatom(proatom_path, atom_number, 0)
            self.proatom_density.append(density)
            self.radial_grid_r.append(r)

    def __str__(self):
        """Return a string representation of the object."""
        return "DDEC6 charges of {}".format(self.data)

    def __repr__(self):
        """Return a representation of the object."""
        return "DDEC6({})".format(self.data)

    def _check_required_attributes(self):
        super(DDEC6, self)._check_required_attributes()

    def _cartesian_dist(self, pt1, pt2):
        """ Small utility function that calculates Euclidian distance between two points
            pt1 and pt2 are numpy arrays representing a point in Cartesian coordinates. """
        return numpy.sqrt(numpy.dot(pt1 - pt2, pt1 - pt2))

    def _read_proatom(
        self, directory, atom_num, charge  # type = str  # type = int  # type = float
    ):
        # type: (...) -> numpy.ndarray, numpy.ndarray
        """Return a list containing proatom reference densities."""
        # TODO: Treat calculations with psuedopotentials
        # TODO: Modify so that proatom densities are read only once for horton
        #       [https://github.com/cclib/cclib/pull/914#discussion_r464039991]
        # File name format:
        #   ** Chargemol **
        #       c2_[atom number]_[nuclear charge]_[electron count]_[cutoff radius]_[# shells]
        #   ** Horton **
        #       atoms.h5
        # File format:
        #   Starting from line 13, each line contains the charge densities for each shell
        # If `charge` is not an integer, proatom densities have to be linearly interpolated between
        # the densities of the ion/atom with floor(charge) and ceiling(charge)
        charge_floor = int(math.floor(charge))
        charge_ceil = int(math.ceil(charge))

        chargemol_path_floor = os.path.join(
            directory,
            "c2_{:03d}_{:03d}_{:03d}_500_100.txt".format(
                atom_num, atom_num, atom_num - charge_floor
            ),
        )
        chargemol_path_ceil = os.path.join(
            directory,
            "c2_{:03d}_{:03d}_{:03d}_500_100.txt".format(
                atom_num, atom_num, atom_num - charge_ceil
            ),
        )
        horton_path = os.path.join(directory, "atoms.h5")

        if os.path.isfile(chargemol_path_floor) or os.path.isfile(chargemol_path_ceil):
            # Use chargemol proatom densities
            # Each shell is .05 angstroms apart (uniform).
            # *scalefactor* = 10.58354497764173 bohrs in module_global_parameter.f08
            if atom_num <= charge_floor:
                density_floor = numpy.array([0])
            else:
                density_floor = numpy.loadtxt(chargemol_path_floor, skiprows=12, dtype=float)
            if atom_num >= charge_ceil:
                density_ceil = numpy.array([0])
            else:
                density_ceil = numpy.loadtxt(chargemol_path_ceil, skiprows=12, dtype=float)

            density = (charge_ceil - charge) * density_floor + (
                charge - charge_floor
            ) * density_ceil
            radiusgrid = numpy.arange(1, len(density) + 1) * 0.05

        elif os.path.isfile(horton_path):
            # Use horton proatom densities
            assert find_package("h5py"), "h5py is needed to read in proatom densities from horton."

            import h5py

            with h5py.File(horton_path, "r") as proatomdb:
                if atom_num <= charge_floor:
                    density_floor = numpy.array([0])
                    radiusgrid = numpy.array([0])
                else:
                    keystring_floor = "Z={}_Q={:+d}".format(atom_num, charge_floor)
                    density_floor = numpy.asanyarray(list(proatomdb[keystring_floor]["rho"]))

                    # gridspec is specification of integration grid for proatom densities in horton.
                    # Example -- ['PowerRTransform', '1.1774580743206259e-07', '20.140888089596444', '41']
                    #   is constructed using PowerRTransform grid
                    #   with rmin = 1.1774580743206259e-07
                    #        rmax = 20.140888089596444
                    #   and  ngrid = 41
                    # PowerRTransform is default in horton-atomdb.py.
                    gridtype, gridmin, gridmax, gridn = (
                        proatomdb[keystring_floor].attrs["rtransform"].split()
                    )
                    gridmin = convertor(float(gridmin), "bohr", "Angstrom")
                    gridmax = convertor(float(gridmax), "bohr", "Angstrom")
                    gridn = int(gridn)
                    # Convert byte to string in Python3
                    if sys.version[0] == "3":
                        gridtype = gridtype.decode("UTF-8")

                    # First verify that it is one of recognized grids
                    assert gridtype in [
                        "LinearRTransform",
                        "ExpRTransform",
                        "PowerRTransform",
                    ], "Grid type not recognized."

                    if gridtype == "LinearRTransform":
                        # Linear transformation. r(t) = rmin + t*(rmax - rmin)/(npoint - 1)
                        gridcoeff = (gridmax - gridmin) / (gridn - 1)
                        radiusgrid = gridmin + numpy.arange(1, gridn + 1) * gridcoeff
                    elif gridtype == "ExpRTransform":
                        # Exponential transformation. r(t) = rmin*exp(t*log(rmax/rmin)/(npoint - 1))
                        gridcoeff = math.log(gridmax / gridmin) / (gridn - 1)
                        radiusgrid = gridmin * numpy.exp(numpy.arange(1, gridn + 1) * gridcoeff)
                    elif gridtype == "PowerRTransform":
                        # Power transformation. r(t) = rmin*t^power
                        # with  power = log(rmax/rmin)/log(npoint)
                        gridcoeff = math.log(gridmax / gridmin) / math.log(gridn)
                        radiusgrid = gridmin * numpy.power(numpy.arange(1, gridn + 1), gridcoeff)

                if atom_num <= charge_ceil:
                    density_ceil = numpy.array([0])
                else:
                    keystring_ceil = "Z={}_Q={:+d}".format(atom_num, charge_ceil)
                    density_ceil = numpy.asanyarray(list(proatomdb[keystring_ceil]["rho"]))

                density = (charge_ceil - charge) * density_floor + (
                    charge - charge_floor
                ) * density_ceil

                del h5py

        else:
            raise MissingInputError("Pro-atom densities were not found in the specified path.")

        if charge == charge_floor:
            density = density_floor

        return density, radiusgrid

    def calculate(self, indices=None, fupdate=0.05):
        """
        Calculate DDEC6 charges based on doi: 10.1039/c6ra04656h paper.
        Cartesian, uniformly spaced grids are assumed for this function.
        """

        # Obtain charge densities on the grid if it does not contain one.
        if not numpy.any(self.volume.data):
            self.logger.info("Calculating charge densities on the provided empty grid.")
            if len(self.data.mocoeffs) == 1:
                self.charge_density = electrondensity_spin(
                    self.data, self.volume, [self.data.mocoeffs[0][: self.data.homos[0]]]
                )
                self.charge_density.data *= 2
            else:
                self.charge_density = electrondensity_spin(
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
            self.charge_density = self.volume

        # * STEP 1 *
        # Carry out step 1 of DDEC6 algorithm [Determining reference charge value]
        # Refer to equations 49-57 in doi: 10.1039/c6ra04656h
        self.logger.info("Creating first reference charges. (Step 1/7)")
        reference_charges, localized_charges, stockholder_charges = self.calculate_reference_charges()
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
        update_kappa = False
        while steps < 7:
            self.logger.info("Optimizing grid weights. (Step {}/7)".format(steps))
            self.N_A.append(self._calculate_w_and_u())

            # Update kappa as described in Figure S4.2 of doi: 10.1039/c6ra04656h
            if numpy.any([x if x < -1e-5 else 0 for x in self.N_A[-1]]):
                update_kappa = True
            if (
                update_kappa
                and numpy.any(numpy.diff(self.N_A)[-1] < 1e-5)  # change in N_A during last cycle
                and numpy.any(numpy.diff(self.N_A)[-2] < 1e-5)  # change in N_A during last to
                # second cycle
            ):
                update_kappa = False
                self._kappa = [0.0 for x in self.data.atomnos]

            if not update_kappa:
                # Increment steps
                steps = steps + 1
                # Calculate G_A and H_A based on S4.3 in doi: 10.1039/c6ra04656h
                self.reshape_G()
                self.calculate_H()

            else:
                # `steps` not incremented in this case
                # First update kappa based on equation 93
                kappa_new = self._kappa - self.N_A[cycle] / self.u_A[cycle]
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

    def calculate_reference_charges(self):
        """ Calculate reference charges from proatom density and molecular density
            [STEP 1 and 2]
            
            Function returns calculated reference charges, localized charges, and stockholder
            charges.
        """
        # Generator object to iterate over the grid
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
        localized_w = numpy.zeros(grid_shape)
        self.closest_r_index = numpy.zeros(grid_shape, dtype=int)

        for atomi, xindex, yindex, zindex in indices:
            # Distance of the grid from atom grid
            dist_r = self._cartesian_dist(
                self.data.atomcoords[-1][atomi],
                self.charge_density.coordinates([xindex, yindex, zindex]),
            )
            self.closest_r_index[atomi][xindex][yindex][zindex] = numpy.abs(
                self.radial_grid_r[atomi] - dist_r
            ).argmin()

            # Equation 54 in doi: 10.1039/c6ra04656h
            stockholder_w[atomi][xindex][yindex][zindex] = self.proatom_density[atomi][
                self.closest_r_index[atomi][xindex][yindex][zindex]
            ]

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
        indices = (
            (i, x, y, z)
            for i in range(self.data.natom)
            for x in range(ngridx)
            for y in range(ngridy)
            for z in range(ngridz)
        )

        self._rho_ref = numpy.zeros((ngridx, ngridy, ngridz))

        for atomi, xindex, yindex, zindex in indices:
            # rho_ref -- Equation 41 in doi: 10.1039/c6ra04656h
            self._rho_ref[xindex][yindex][zindex] += self.proatom_density[atomi][
                self.closest_r_index[atomi][xindex][yindex][zindex]
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

            while phiAI[atomi] <= 0:
                # Iterative algorithm until convergence
                # Refer to S101 in doi: 10.1039/c6ra04656h
                bigphiAI[atomi] = 2 * bigphiAI[atomi] - phiAI[atomi] / self._integrate_from_radial(
                    [numpy.sqrt(self._y_a[atomi])], [atomi]
                )

                # When Phi is updated, related quantities are updated as well
                # Refer to S100 in doi: 10.1039/c6ra04656h
                phiAI[atomi], self._cond_density[atomi] = self._update_phiai(
                    self._y_a[atomi], bigphiAI[atomi], atomi
                )

                self._candidates_phi[atomi].append(phiAI[atomi])
                self._candidates_bigPhi[atomi].append(bigphiAI[atomi])

            # lowerbigPhi is largest negative Phi.
            # upperbigPhi is smallest positive Phi.
            self._candidates_phi[atomi] = numpy.array(self._candidates_phi[atomi], dtype=float)
            self._candidates_bigPhi[atomi] = numpy.array(
                self._candidates_bigPhi[atomi], dtype=float
            )
            if numpy.count_nonzero(self._candidates_phi[atomi] < 0) > 0:
                lower_ind = numpy.where(
                    self._candidates_phi[atomi]
                    == self._candidates_phi[atomi][self._candidates_phi[atomi] < 0].max()
                )[0][0]
                lowerbigPhi = self._candidates_bigPhi[atomi][lower_ind]
                lowerphi = self._candidates_phi[atomi][lower_ind]
            else:  # assign some large negative number
                lowerbigPhi = numpy.NINF
                lowerphi = numpy.NINF
            if numpy.count_nonzero(self._candidates_phi[atomi] > 0) > 0:
                upper_ind = numpy.where(
                    self._candidates_phi[atomi]
                    == self._candidates_phi[atomi][self._candidates_phi[atomi] > 0].min()
                )[0][0]
                upperbigPhi = self._candidates_bigPhi[atomi][upper_ind]
                upperphi = self._candidates_phi[atomi][upper_ind]
            else:  # assign some large positive number
                upperbigPhi = numpy.PINF
                upperphi = numpy.PINF

            for iteration in range(50):
                # Flow diagram on Figure S1 in doi: 10.1039/c6ra04656h details the procedure.
                midbigPhi = (lowerbigPhi + upperbigPhi) / 2.0
                midphi = self._update_phiai(self._y_a[atomi], midbigPhi, atomi)[0]
                # Exit conditions
                if abs(lowerphi) < self.convergence_level:
                    bigphiAI[atomi] = lowerbigPhi
                    break
                elif abs(upperphi) < self.convergence_level:
                    bigphiAI[atomi] = upperbigPhi
                    break
                elif abs(midphi) < self.convergence_level:
                    bigphiAI[atomi] = midbigPhi
                    break

                # Parabolic fitting as described on Figure S1 in doi: 10.1039/c6ra04656h
                # Type casting here converts from size 1 numpy.ndarray to float
                xpts = numpy.array(
                    [float(lowerbigPhi), float(midbigPhi), float(upperbigPhi)], dtype=float
                )
                ypts = numpy.array(
                    [float(lowerphi), float(midphi), float(upperphi)], dtype=float
                )
                fit = numpy.polyfit(xpts, ypts, 2)
                roots = numpy.roots(fit)  # max two roots (bigPhi)

                belowphi = self._update_phiai(self._y_a[atomi], roots.min(), atomi)[0]
                abovephi = self._update_phiai(self._y_a[atomi], roots.max(), atomi)[0]

                if abs(abovephi) < self.convergence_level:
                    bigphiAI[atomi] = roots.min()
                    break
                elif abs(belowphi) < self.convergence_level:
                    bigphiAI[atomi] = roots.max()
                    break
                else:
                    if 3 * abs(abovephi) < abs(belowphi):
                        corbigPhi = roots.max() - 2.0 * abovephi * (
                            roots.max() - roots.min()
                        ) / (abovephi - belowphi)
                    elif 3 * abs(belowphi) < abs(abovephi):
                        corbigPhi = roots.min() - 2.0 * belowphi * (
                            roots.max() - roots.min()
                        ) / (abovephi - belowphi)
                    else:
                        corbigPhi = (roots.max() + roots.min()) / 2.0
                    # New candidates
                    corphi = self._update_phiai(self._y_a[atomi], corbigPhi, atomi)[0]
                    self._candidates_bigPhi[atomi] = numpy.array(
                        [
                            lowerbigPhi,
                            midbigPhi,
                            upperbigPhi,
                            roots.max(),
                            roots.min(),
                            corbigPhi,
                        ],
                        dtype=float,
                    )
                    self._candidates_phi[atomi] = numpy.array(
                        [lowerphi, midphi, upperphi, abovephi, belowphi, corphi], dtype=float
                    )

                    # Update upperphi and lowerphi
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

                    if abs(lowerphi) < self.convergence_level:
                        bigphiAI[atomi] = self._candidates_bigPhi[atomi][lower_ind]
                        break
                    elif abs(upperphi) < self.convergence_level:
                        bigphiAI[atomi] = self._candidates_bigPhi[atomi][upper_ind]
                        break
                    else:
                        # Fitting needs to continue in this case.
                        lowerbigPhi = self._candidates_bigPhi[atomi][lower_ind]
                        lowerphi = self._candidates_phi[atomi][lower_ind]
                        upperbigPhi = self._candidates_bigPhi[atomi][upper_ind]
                        upperphi = self._candidates_phi[atomi][upper_ind]

                if iteration == self.max_iteration - 1:
                    raise ConvergenceError("Iterative conditioning failed to converge.")

            # Set final conditioned density using chosen Phi
            self._cond_density[atomi] = self._update_phiai(self._y_a[atomi], bigphiAI[atomi], atomi)[1]

        self.logger.info("Calculating tau and combined conditioned densities.")

        # Calculate tau(r) and rho^cond(r)
        # Refer to equation 65 and 66 in doi: 10.1039/c6ra04656h
        # Assign rho^cond on grid using generator object
        self.rho_cond = copy.deepcopy(self.charge_density)
        self.rho_cond.data = numpy.zeros_like(self.rho_cond.data, dtype=float)
        rho_cond_sqrt = numpy.zeros_like(self.rho_cond.data, dtype=float)

        # Generator object to iterate over the grid
        ngridx, ngridy, ngridz = self.charge_density.data.shape
        indices = (
            (i, x, y, z)
            for i in range(self.data.natom)
            for x in range(ngridx)
            for y in range(ngridy)
            for z in range(ngridz)
        )

        self._leftterm = numpy.zeros((self.data.natom, ngridx, ngridy, ngridz), dtype=float)
        # rho_cond_cartesian is rho^cond projected on Cartesian grid
        # (used for Step 4 calculation)
        self._rho_cond_cartesian = numpy.zeros(
            (self.data.natom, ngridx, ngridy, ngridz), dtype=float
        )
        self.tau = []

        # rho_cond -- equation 65 in doi: 10.1039/c6ra04656h
        for atomi, xindex, yindex, zindex in indices:
            self.rho_cond.data[xindex][yindex][zindex] += self._cond_density[atomi][
                self.closest_r_index[atomi][xindex][yindex][zindex]
            ]

        rho_cond_sqrt = numpy.sqrt(self.rho_cond.data)

        for atomi in range(self.data.natom):
            self.tau.append(numpy.zeros_like(self.proatom_density[atomi], dtype=float))
            for xindex, yindex, zindex in numpy.ndindex(ngridx, ngridy, ngridz):
                # leftterm is the first spherical average term in equation 66.
                # <rho^cond_A(r_A) / sqrt(rho^cond(r))>
                self._rho_cond_cartesian[atomi][xindex][yindex][zindex] = self._cond_density[atomi][
                    self.closest_r_index[atomi][xindex][yindex][zindex]
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
            self._rho_A[atomi].data = (
                self.charge_density.data * self._rho_cond_cartesian[atomi] / self.rho_cond.data
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
                (1 - self._rho_cond_cartesian[atomi] / self.rho_cond.data)
                * self._rho_A[atomi].data,
                atomi,
                self.radial_grid_r[atomi],
            )
            self._theta[atomi] = numpy.maximum.accumulate(self._theta[atomi][::-1])[::-1]
            self._wavg[atomi] = self._spherical_average_from_cartesian(
                self._rho_cond_cartesian[atomi] / self.rho_cond.data,
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
            TODO: some similarities exist with step 3 routine. code refactoring may improve readability.
        """
        self._candidates_bigPhi = []
        self._candidates_phi = []

        # Initial conditions are detailed in Figure S3 in doi: 10.1039/c6ra04656h
        phiAII = numpy.zeros_like(self.data.atomnos, dtype=float)
        bigphiAII = numpy.zeros_like(self.data.atomnos, dtype=float)
        self._g = []

        for atomi in range(self.data.natom):
            # G_A -- equation S102 in doi: 10.1039/c6ra04656h
            self._g.append(numpy.zeros_like(self.proatom_density[atomi], dtype=float))
            self._g[atomi] = self.rho_wavg[atomi] + bigphiAII[atomi] * numpy.sqrt(
                self.rho_wavg[atomi]
            )
            # Exponential constraint (as expressed in equation S105)
            eta = (1 - (self.tau[atomi]) ** 2) * 1.75 * convertor(1, "Angstrom", "bohr")
            exp_applied = self._g[atomi][:-1] * numpy.exp(
                -1 * eta[1:] * numpy.diff(self.radial_grid_r[atomi])
            )
            for radiusi in range(1, len(self._g[atomi])):
                self._g[atomi][radiusi] = min(self._g[atomi][radiusi], exp_applied[radiusi - 1])
            # phi_A^II -- Equation S106 in doi: 10.1039/c6ra04656h
            phiAII[atomi] = self._integrate_from_radial(
                [self._g[atomi] - self.rho_wavg[atomi]], [atomi]
            )

            self._candidates_bigPhi.append([bigphiAII[atomi]])
            self._candidates_phi.append([phiAII[atomi]])

            while phiAII[atomi] <= 0:
                # Iterative algorithm until convergence
                # Refer to Figure S3 in doi: 10.1039/c6ra04656h
                temp = self._integrate_from_radial([numpy.sqrt(self.rho_wavg[atomi])], [atomi])
                bigphiAII[atomi] = 2 * bigphiAII[atomi] - phiAII[atomi] / temp

                # When Phi is updated, related quantities are updated as well
                # Refer to S102 in doi: 10.1039/c6ra04656h
                phiAII[atomi], self._g[atomi] = self._update_phiaii(
                    self.rho_wavg[atomi], eta, bigphiAII[atomi], atomi
                )

                self._candidates_phi[atomi].append(phiAII[atomi])
                self._candidates_bigPhi[atomi].append(bigphiAII[atomi])

            # lowerbigPhi is negative Phi with the smallest magnitude.
            # upperbigPhi is positive Phi with the smallest magnitude.
            self._candidates_phi[atomi] = numpy.array(self._candidates_phi[atomi], dtype=float)
            self._candidates_bigPhi[atomi] = numpy.array(
                self._candidates_bigPhi[atomi], dtype=float
            )
            if numpy.count_nonzero(self._candidates_phi[atomi] < 0) > 0:
                lower_ind = numpy.where(
                    self._candidates_phi[atomi]
                    == self._candidates_phi[atomi][self._candidates_phi[atomi] < 0].max()
                )[0][0]
                lowerbigPhi = self._candidates_bigPhi[atomi][lower_ind]
                lowerphi = self._candidates_phi[atomi][lower_ind]
            else:  # assign some large negative number
                lowerbigPhi = numpy.NINF
                lowerphi = numpy.NINF
            if numpy.count_nonzero(self._candidates_phi[atomi] > 0) > 0:
                upper_ind = numpy.where(
                    self._candidates_phi[atomi]
                    == self._candidates_phi[atomi][self._candidates_phi[atomi] > 0].min()
                )[0][0]
                upperbigPhi = self._candidates_bigPhi[atomi][upper_ind]
                upperphi = self._candidates_phi[atomi][upper_ind]
            else:  # assign some large positive number
                upperbigPhi = numpy.PINF
                upperphi = numpy.PINF

            for iteration in range(50):
                # Flow diagram on Figure S1 in doi: 10.1039/c6ra04656h details the procedure.
                midbigPhi = (lowerbigPhi + upperbigPhi) / 2.0
                midphi = self._update_phiai(self._y_a[atomi], midbigPhi, atomi)[0]
                # Exit conditions
                if abs(lowerphi) < self.convergence_level:
                    bigphiAII[atomi] = lowerbigPhi
                    break
                elif abs(upperphi) < self.convergence_level:
                    bigphiAII[atomi] = upperbigPhi
                    break
                elif abs(midphi) < self.convergence_level:
                    bigphiAII[atomi] = midbigPhi
                    break

                # Parabolic fitting as described on Figure S1 in doi: 10.1039/c6ra04656h
                # Type casting here converts from size 1 numpy.ndarray to float
                xpts = numpy.array(
                    [float(lowerbigPhi), float(midbigPhi), float(upperbigPhi)], dtype=float
                )
                ypts = numpy.array(
                    [float(lowerphi), float(midphi), float(upperphi)], dtype=float
                )
                fit = numpy.polyfit(xpts, ypts, 2)
                roots = numpy.roots(fit)  # max two roots (bigPhi)

                belowphi = self._update_phiaii(self.rho_wavg[atomi], eta, roots.min(), atomi)[0]
                abovephi = self._update_phiaii(self.rho_wavg[atomi], eta, roots.max(), atomi)[0]

                if abs(abovephi) < self.convergence_level:
                    bigphiAII[atomi] = roots.min()
                    break
                elif abs(belowphi) < self.convergence_level:
                    bigphiAII[atomi] = roots.max()
                    break
                else:
                    if 3 * abs(abovephi) < abs(belowphi):
                        corbigPhi = roots.max() - 2.0 * abovephi * (
                            roots.max() - roots.min()
                        ) / (abovephi - belowphi)
                    elif 3 * abs(belowphi) < abs(abovephi):
                        corbigPhi = roots.min() - 2.0 * belowphi * (
                            roots.max() - roots.min()
                        ) / (abovephi - belowphi)
                    else:
                        corbigPhi = (roots.max() + roots.min()) / 2.0
                    # New candidates
                    corphi = self._update_phiaii(self.rho_wavg[atomi], eta, corbigPhi, atomi)[0]
                    self._candidates_bigPhi[atomi] = numpy.array(
                        [
                            lowerbigPhi,
                            midbigPhi,
                            upperbigPhi,
                            roots.max(),
                            roots.min(),
                            corbigPhi,
                        ],
                        dtype=float,
                    )
                    self._candidates_phi[atomi] = numpy.array(
                        [lowerphi, midphi, upperphi, abovephi, belowphi, corphi], dtype=float
                    )

                    # Update upperphi and lowerphi
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

                    if abs(lowerphi) < self.convergence_level:
                        bigphiAII[atomi] = self._candidates_bigPhi[atomi][lower_ind]
                        break
                    elif abs(upperphi) < self.convergence_level:
                        bigphiAII[atomi] = self._candidates_bigPhi[atomi][upper_ind]
                        break
                    else:
                        # Fitting needs to continue in this case.
                        lowerbigPhi = self._candidates_bigPhi[atomi][lower_ind]
                        lowerphi = self._candidates_phi[atomi][lower_ind]
                        upperbigPhi = self._candidates_bigPhi[atomi][upper_ind]
                        upperphi = self._candidates_phi[atomi][upper_ind]

                if iteration == self.max_iteration - 1:
                    raise ConvergenceError("Iterative conditioning failed to converge.")

            # Set final G_A value using chosen Phi
            self._g[atomi] = self._update_phiaii(self.rho_wavg[atomi], eta, bigphiAII[atomi], atomi)[1]

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
            eta = 2.5 * convertor(1, "Angstrom", "bohr") / (1 - (self.tau[atomi]) ** 2)
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

    def _update_phiaii(self, rhowavg, eta, bigphiAI, atomi):
        # Update phi^a_ii and quantity that directly follows (G_A) in each step of
        # iterative optimization. (Refer to Figure S3)
        # Re-evaluate g_a
        # Equations can be found in Figure S3.
        ga = rhowavg + bigphiAI * numpy.sqrt(rhowavg)

        # Exponential Decrease Condition
        exp_applied = ga[:-1] * numpy.exp(-1 * eta[1:] * numpy.diff(self.radial_grid_r[atomi]))
        for radiusi in range(1, len(ga)):
            ga[radiusi] = min(ga[radiusi], exp_applied[radiusi - 1])

        # Re-evaluate phi_AII
        phiAII = self._integrate_from_radial([ga - rhowavg], [atomi])

        return phiAII, ga

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

    def _integrate_from_radial(self, radial_density_list, atom_list):
        # Function that reads in list of radial densities, projects it on Cartesian grid,
        # and returns integrated value
        grid = copy.deepcopy(self.charge_density)
        grid.data = numpy.zeros_like(grid.data, dtype=float)

        ngridx, ngridy, ngridz = self.charge_density.data.shape
        indices = ((x, y, z) for x in range(ngridx) for y in range(ngridy) for z in range(ngridz))

        for density, atomi in zip(radial_density_list, atom_list):
            for x, y, z in indices:
                grid.data[x][y][z] = (
                    grid.data[x][y][z] + density[self.closest_r_index[atomi][x][y][z]]
                )

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
        # Generator object to iterate over the grid
        ngridx, ngridy, ngridz = self.charge_density.data.shape

        self.rho_cond.data = numpy.zeros_like(self.rho_cond.data, dtype=float)
        self._rho_cond_cartesian = numpy.zeros(
            (self.data.natom, ngridx, ngridy, ngridz), dtype=float
        )

        for atomi in range(self.data.natom):
            grid = ((x, y, z) for x in range(ngridx) for y in range(ngridy) for z in range(ngridz))
            for xindex, yindex, zindex in grid:
                self.rho_cond.data[xindex][yindex][zindex] += self._cond_density[atomi][
                    self.closest_r_index[atomi][xindex][yindex][zindex]
                ]
                self._rho_cond_cartesian[atomi][xindex][yindex][zindex] = self._cond_density[atomi][
                    self.closest_r_index[atomi][xindex][yindex][zindex]
                ]
