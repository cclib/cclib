# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculate properties of nuclei based on data parsed by cclib."""

import logging

import numpy as np

from cclib.method.calculationmethod import Method
from cclib.parser.utils import PeriodicTable


class Nuclear(Method):
    """A container for methods pertaining to atomic nuclei."""

    def __init__(self, data, progress=None, loglevel=logging.INFO, logname="Log"):

        self.required_attrs = ('natom','atomcoords','atomnos','charge')

        super(Nuclear, self).__init__(data, progress, loglevel, logname)

    def __str__(self):
        """Return a string representation of the object."""
        return "Nuclear"

    def __repr__(self):
        """Return a representation of the object."""
        return "Nuclear"

    def stoichiometry(self):
        """Return the stoichemistry of the object according to the Hill system"""
        pt = PeriodicTable()
        elements = [pt.element[ano] for ano in self.data.atomnos]
        counts = {el: elements.count(el) for el in set(elements)}

        formula = ""
        elcount = lambda el, c: "%s%i" % (el, c) if c > 1 else el
        if 'C' in elements:
            formula += elcount('C', counts['C'])
            counts.pop('C')
            if 'H' in elements:
              formula += elcount('H', counts['H'])
              counts.pop('H')
        for el, c in sorted(counts.items()):
            formula += elcount(el, c)

        if getattr(self.data, 'charge', 0):
            magnitude = abs(self.data.charge)
            sign = "+" if self.data.charge > 0 else "-"
            formula += "(%s%i)" % (sign, magnitude)
        return formula

    def repulsion_energy(self):
        """Return the nuclear repulsion energy."""

        nre = 0.0
        for i in range(self.data.natom):
            ri = self.data.atomcoords[0][i]
            zi = self.data.atomnos[i]
            for j in range(i+1, self.data.natom):
                rj = self.data.atomcoords[0][j]
                zj = self.data.atomnos[j]
                d = np.linalg.norm(ri-rj)
                nre += zi*zj/d
        return nre

    @staticmethod
    def _get_most_abundant_isotope(element):
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

    @staticmethod
    def get_isotopic_masses(charges):
        """Return the masses for the given nuclei, respresented by their
        nuclear charges.
        """
        import periodictable
        masses = []
        for charge in charges:
            el = periodictable.elements[charge]
            isotope = Nuclear._get_most_abundant_isotope(el)
            mass = isotope.mass
            masses.append(mass)
        return np.array(masses)

    def center_of_mass(self):
        """Return the center of mass."""
        charges = self.data.atomnos
        masses = Nuclear.get_isotopic_masses(charges)
        coords = self.data.atomcoords[-1]

        mwc = coords * masses[:, np.newaxis]
        numerator = np.sum(mwc, axis=0)
        denominator = np.sum(masses)

        return numerator / denominator

    def moment_of_inertia_tensor(self):
        """Return the moment of inertia tensor."""
        charges = self.data.atomnos
        masses = Nuclear.get_isotopic_masses(charges)
        coords = self.data.atomcoords[-1]

        moi_tensor = np.empty((3, 3))

        moi_tensor[0][0] = np.sum(masses * (coords[:, 1]**2 + coords[:, 2]**2))
        moi_tensor[1][1] = np.sum(masses * (coords[:, 0]**2 + coords[:, 2]**2))
        moi_tensor[2][2] = np.sum(masses * (coords[:, 0]**2 + coords[:, 1]**2))

        moi_tensor[0][1] = np.sum(masses * coords[:, 0] * coords[:, 1])
        moi_tensor[0][2] = np.sum(masses * coords[:, 0] * coords[:, 2])
        moi_tensor[1][2] = np.sum(masses * coords[:, 1] * coords[:, 2])

        moi_tensor[1][0] = moi_tensor[0][1]
        moi_tensor[2][0] = moi_tensor[0][2]
        moi_tensor[2][1] = moi_tensor[1][2]

        return moi_tensor

    def principal_moments_of_inertia(self):
        """Return the principal moments of inertia in 3 kinds of units:
        1. [amu][bohr]^2
        2. [amu][angstrom]^2
        3. [g][cm]^2
        and the principal axes.
        """
        import scipy.constants as spc
        moi_tensor = self.moment_of_inertia_tensor()
        principal_moments, principal_axes = np.linalg.eigh(moi_tensor)
        amu2g = spc.value('unified atomic mass unit') * spc.kilo
        bohr2ang = spc.value('atomic unit of length') / spc.angstrom
        conv1 = bohr2ang ** 2
        conv2 = amu2g * (spc.value('atomic unit of length') * spc.centi) ** 2
        return (principal_moments,
                principal_moments * conv1,
                principal_moments * conv2,
                principal_axes)

    def rotational_constants(self, units='ghz'):
        """Compute the rotational constants in 1/cm or GHz."""
        choices = ('invcm', 'ghz')
        if units.lower() not in choices:
            raise ValueError("Invalid units, pick one of {}".format(choices))
        units = units.lower()
        import scipy.constants as spc
        principal_moments = self.principal_moments_of_inertia()[0]
        bohr2ang = spc.value('atomic unit of length') / spc.angstrom
        xfamu = 1 / spc.value('electron mass in u')
        xthz = spc.value('hartree-hertz relationship')
        rotghz = xthz * (bohr2ang ** 2) / (2 * xfamu * spc.giga)
        ghz2invcm = spc.giga * spc.centi / spc.c
        if units == 'ghz':
            return rotghz / principal_moments
        if units == 'invcm':
            return rotghz * ghz2invcm / principal_moments
