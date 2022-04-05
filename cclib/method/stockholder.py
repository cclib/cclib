# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Stockholder partitioning based on cclib data."""
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


class Stockholder(Method):
    """An abstract base class for stockholder-type methods."""

    # All of these are required for stockholder-type charges.
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
        """ Initialize Stockholder-type method object.
            Inputs are:
                data -- ccData object that describe target molecule.
                volume -- Volume object that describe target Cartesian grid.
                proatom_path -- path to proatom densities
                (directory containing atoms.h5 in horton or c2_001_001_000_400_075.txt in chargemol)
        """
        super().__init__(data, progress, loglevel, logname)

        self.volume = volume
        self.proatom_path = proatom_path

        # Check whether proatom_path is a valid directory or not.
        assert os.path.isdir(
            proatom_path
        ), "Directory that contains proatom densities should be added as an input."

        # Read in reference charges.
        self.proatom_density = []
        self.radial_grid_r = []
        for atom_number in self.data.atomnos:
            density, r = self._read_proatom(self.proatom_path, atom_number, 0)
            self.proatom_density.append(density)
            self.radial_grid_r.append(r)

    def __str__(self):
        """Return a string representation of the object."""
        return "Stockholder"

    def __repr__(self):
        """Return a representation of the object."""
        return "Stockholder"

    def _check_required_attributes(self):
        super()._check_required_attributes()

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
            f"c2_{atom_num:03d}_{atom_num:03d}_{atom_num - charge_floor:03d}_500_100.txt",
        )
        chargemol_path_ceil = os.path.join(
            directory,
            f"c2_{atom_num:03d}_{atom_num:03d}_{atom_num - charge_ceil:03d}_500_100.txt",
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
                    keystring_floor = f"Z={atom_num}_Q={charge_floor:+d}"
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
                    if isinstance(gridtype, bytes):
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
                    keystring_ceil = f"Z={atom_num}_Q={charge_ceil:+d}"
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
        """ Charge density on a Cartesian grid is a common routine required for Stockholder-type
            and related methods. This abstract class prepares the grid if input Volume object
            is empty.
        """
        # Obtain charge densities on the grid if it does not contain one.
        if not numpy.any(self.volume.data):
            self.logger.info("Calculating charge densities on the provided empty grid.")
            if len(self.data.mocoeffs) == 1:
                self.charge_density = electrondensity_spin(
                    self.data, self.volume, [self.data.mocoeffs[0][: self.data.homos[0] + 1]]
                )
                self.charge_density.data *= 2
            else:
                self.charge_density = electrondensity_spin(
                    self.data,
                    self.volume,
                    [
                        self.data.mocoeffs[0][: self.data.homos[0] + 1],
                        self.data.mocoeffs[1][: self.data.homos[1] + 1],
                    ],
                )
        # If charge densities are provided beforehand, log this information
        # `Volume` object does not contain (nor rely on) information about the constituent atoms.
        else:
            self.logger.info("Using charge densities from the provided Volume object.")
            self.charge_density = self.volume
