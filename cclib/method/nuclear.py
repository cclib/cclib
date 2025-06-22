# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculate properties of nuclei based on data parsed by cclib."""

import logging
from typing import Tuple

from cclib.method.calculationmethod import Method
from cclib.parser.utils import PeriodicTable, convertor, find_package

import numpy as np
import periodictable as pt

_found_scipy = find_package("scipy")
if _found_scipy:
    import scipy.constants


def _check_periodictable(found_periodictable: bool) -> None:
    if not _found_periodictable:
        raise ImportError("You must install `periodictable` to use this function")


def _check_scipy(found_scipy: bool) -> None:
    if not _found_scipy:
        raise ImportError("You must install `scipy` to use this function")


def get_most_abundant_isotope(element: "pt.Element") -> "pt.Isotope":
    """Given a `periodictable` element, return the most abundant
    isotope.
    """
    most_abundant_isotope = element.isotopes[0]
    abundance = 0
    for iso in element:
        if iso.abundance > abundance:
            most_abundant_isotope = iso
            abundance = iso.abundance
    return most_abundant_isotope


def get_isotopic_masses(charges):
    """Return the masses for the given nuclei, respresented by their
    nuclear charges.
    """
    return np.array([get_most_abundant_isotope(pt.elements[charge]).mass for charge in charges])


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

        def elcount(el: str, c: int) -> str:
            return f"{el}{int(c)}" if c > 1 else el

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

    def moment_of_inertia_tensor(self, atomcoords_index: int = -1) -> np.ndarray:
        """Return the moment of inertia tensor."""
        charges = self.data.atomnos
        coords = self.data.atomcoords[atomcoords_index]
        masses = get_isotopic_masses(charges)

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
        self, units: str = "amu_bohr_2"
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
        moi_tensor = self.moment_of_inertia_tensor()
        principal_moments, principal_axes = np.linalg.eigh(moi_tensor)
        if units == "amu_bohr_2":
            _check_scipy(_found_scipy)
            bohr2ang = scipy.constants.value("atomic unit of length") / scipy.constants.angstrom
            conv = 1 / bohr2ang**2
        elif units == "amu_angstrom_2":
            conv = 1
        elif units == "g_cm_2":
            _check_scipy(_found_scipy)
            amu2g = scipy.constants.value("unified atomic mass unit") * scipy.constants.kilo
            conv = amu2g * (scipy.constants.angstrom / scipy.constants.centi) ** 2
        return conv * principal_moments, principal_axes

    def rotational_constants(self, units: str = "ghz") -> np.ndarray:
        """Compute the rotational constants in 1/cm or GHz."""
        choices = ("invcm", "ghz")
        units = units.lower()
        if units not in choices:
            raise ValueError(f"Invalid units, pick one of {choices}")
        principal_moments = self.principal_moments_of_inertia("amu_angstrom_2")[0]
        _check_scipy(_found_scipy)
        bohr2ang = scipy.constants.value("atomic unit of length") / scipy.constants.angstrom
        xfamu = 1 / scipy.constants.value("electron mass in u")
        xthz = scipy.constants.value("hartree-hertz relationship")
        rotghz = xthz * (bohr2ang**2) / (2 * xfamu * scipy.constants.giga)
        if units == "ghz":
            conv = rotghz
        if units == "invcm":
            ghz2invcm = scipy.constants.giga * scipy.constants.centi / scipy.constants.c
            conv = rotghz * ghz2invcm
        return conv / principal_moments


del find_package
