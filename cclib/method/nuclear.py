# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculate properties of nuclei based on data parsed by cclib."""

import logging
from enum import Enum, auto
from typing import Tuple

from cclib.method.calculationmethod import Method
from cclib.parser.utils import PeriodicTable, convertor, find_package

import numpy as np
import qcelemental as qcel
from qcelemental import constants

_CENTI = 0.01
_GIGA = 1000000000.0
_KILO = 1000.0
_ANGSTROM = 1e-10


def get_isotopic_masses(charges):
    """Return the masses for the given nuclei, respresented by their
    nuclear charges.
    """
    return np.array([qcel.periodictable.to_mass(charge) for charge in charges])


class Nuclear(Method):
    """A container for methods pertaining to atomic nuclei."""

    def __init__(self, data, progress=None, loglevel=logging.INFO, logname="Log") -> None:
        super().__init__(data, progress, loglevel, logname)
        self.required_attrs = ("natom", "atomcoords", "atomnos", "charge")

    def __str__(self) -> str:
        """Return a string representation of the object."""
        return "Nuclear"

    def __repr__(self) -> str:
        """Return a representation of the object."""
        return "Nuclear"

    def stoichiometry(self) -> str:
        """Return the stoichemistry of the object according to the Hill system"""
        cclib_pt = PeriodicTable()
        elements = [cclib_pt.element[ano] for ano in self.data.atomnos]
        counts = {el: elements.count(el) for el in set(elements)}

        formula = ""
        elcount = lambda el, c: f"{el}{int(c)}" if c > 1 else el
        if "C" in elements:
            formula += elcount("C", counts["C"])
            counts.pop("C")
            if "H" in elements:
                formula += elcount("H", counts["H"])
                counts.pop("H")
        for el, c in sorted(counts.items()):
            formula += elcount(el, c)

        if getattr(self.data, "charge", 0):
            magnitude = abs(self.data.charge)
            sign = "+" if self.data.charge > 0 else "-"
            formula += f"({sign}{int(magnitude)})"
        return formula

    def repulsion_energy(self, atomcoords_index: int = -1) -> float:
        """Return the nuclear repulsion energy."""
        nre = 0.0
        atomcoords = convertor(self.data.atomcoords[atomcoords_index], "Angstrom", "bohr")
        for i in range(self.data.natom):
            ri = atomcoords[i]
            zi = self.data.atomnos[i]
            for j in range(i + 1, self.data.natom):
                rj = atomcoords[j]
                zj = self.data.atomnos[j]
                d = np.linalg.norm(ri - rj)
                nre += zi * zj / d
        return convertor(nre, "hartree", "eV")

    def center_of_mass(self, atomcoords_index: int = -1) -> float:
        """Return the center of mass."""
        charges = self.data.atomnos
        coords = self.data.atomcoords[atomcoords_index]
        masses = get_isotopic_masses(charges)

        mwc = coords * masses[:, np.newaxis]
        numerator = np.sum(mwc, axis=0)
        denominator = np.sum(masses)

        return numerator / denominator

    class MOITensorUnits(Enum):
        amu_bohr_squared = auto()
        amu_angstrom_squared = auto()

    def moment_of_inertia_tensor(
        self,
        units: MOITensorUnits = MOITensorUnits.amu_angstrom_squared,
        atomcoords_index: int = -1,
    ) -> np.ndarray:
        """Return the moment of inertia tensor."""

        charges = self.data.atomnos
        coords = self.data.atomcoords[atomcoords_index]
        masses = get_isotopic_masses(charges)

        allowed_units = Nuclear.MOITensorUnits
        # Put input components in the right units before assembling the tensor.
        if units == allowed_units.amu_bohr_squared:
            bohr2ang = constants.atomic_unit_of_length / _ANGSTROM
            coords /= bohr2ang
        elif units == allowed_units.amu_angstrom_squared:
            # No need to do anything, already in the correct units
            pass
        else:
            raise ValueError(f"Invalid units, pick one of Nuclear.MOITensorUnits")

        moi_tensor = np.empty((3, 3))

        moi_tensor[0][0] = np.sum(masses * (coords[:, 1] ** 2 + coords[:, 2] ** 2))
        moi_tensor[1][1] = np.sum(masses * (coords[:, 0] ** 2 + coords[:, 2] ** 2))
        moi_tensor[2][2] = np.sum(masses * (coords[:, 0] ** 2 + coords[:, 1] ** 2))

        moi_tensor[0][1] = np.sum(masses * coords[:, 0] * coords[:, 1])
        moi_tensor[0][2] = np.sum(masses * coords[:, 0] * coords[:, 2])
        moi_tensor[1][2] = np.sum(masses * coords[:, 1] * coords[:, 2])

        moi_tensor[1][0] = moi_tensor[0][1]
        moi_tensor[2][0] = moi_tensor[0][2]
        moi_tensor[2][1] = moi_tensor[1][2]

        return moi_tensor

    def principal_moments_of_inertia(
        self, units: str = "amu_bohr_2", atomcoords_index: int = -1
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return the principal moments of inertia in 3 kinds of units:
        1. [amu][bohr]^2
        2. [amu][angstrom]^2
        3. [g][cm]^2
        and the principal axes.
        """
        choices = ("amu_bohr_2", "amu_angstrom_2", "g_cm_2")
        units = units.lower()
        if units not in choices:
            raise ValueError(f"Invalid units, pick one of {choices}")
        moi_tensor = self.moment_of_inertia_tensor(atomcoords_index=atomcoords_index)
        principal_moments, principal_axes = np.linalg.eigh(moi_tensor)
        if units == "amu_bohr_2":
            bohr2ang = constants.atomic_unit_of_length / _ANGSTROM
            conv = 1 / bohr2ang**2
        elif units == "amu_angstrom_2":
            conv = 1
        elif units == "g_cm_2":
            amu2g = constants.unified_atomic_mass_unit * _KILO
            conv = amu2g * (constants.atomic_unit_of_length * _CENTI) ** 2
        return conv * principal_moments, principal_axes

    def rotational_constants(self, units: str = "ghz", atomcoords_index: int = -1) -> np.ndarray:
        """Compute the rotational constants in 1/cm or GHz."""
        choices = ("invcm", "ghz")
        units = units.lower()
        if units not in choices:
            raise ValueError(f"Invalid units, pick one of {choices}")
        principal_moments = self.principal_moments_of_inertia(
            "amu_angstrom_2", atomcoords_index=atomcoords_index
        )[0]
        bohr2ang = constants.atomic_unit_of_length / _ANGSTROM
        xfamu = 1 / constants.electron_mass_in_u
        rotghz = constants.hartree_hertz_relationship * (bohr2ang**2) / (2 * xfamu * _GIGA)
        if units == "ghz":
            conv = rotghz
        elif units == "invcm":
            ghz2invcm = _GIGA * _CENTI / constants.c
            conv = rotghz * ghz2invcm
        else:
            raise ValueError("Invalid units, pick one of {}".format(choices))
        return conv / principal_moments
