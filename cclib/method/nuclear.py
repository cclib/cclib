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
from cclib.parser.utils import find_package

_found_periodictable = find_package("periodictable")
if _found_periodictable:
    import periodictable as pt

_found_scipy = find_package("scipy")
if _found_scipy:
    import scipy.constants as spc


def _check_periodictable(found_periodictable):
    if not _found_periodictable:
        raise ImportError("You must install `periodictable` to use this function")


def _check_scipy(found_scipy):
    if not _found_scipy:
        raise ImportError("You must install `scipy` to use this function")


def get_most_abundant_isotope(element):
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
    _check_periodictable(_found_periodictable)
    masses = []
    for charge in charges:
        el = pt.elements[charge]
        isotope = get_most_abundant_isotope(el)
        mass = isotope.mass
        masses.append(mass)
    return np.array(masses)


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
        cclib_pt = PeriodicTable()
        elements = [cclib_pt.element[ano] for ano in self.data.atomnos]
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
    
    def distances(self, atomcoords_index = -1):
        '''Returns a matrix with distances between i_th and j_th atoms at index [i][j]'''

        coords = self.data.atomcoords[atomcoords_index]
        number_of_atoms = self.data.natom
        distance_matrix = np.zeros([number_of_atoms, number_of_atoms])

        for i in range(0, number_of_atoms - 1):
            for j in range(i + 1, number_of_atoms):
                distance_matrix[i][j] = np.linalg.norm(np.subtract(coords[j], coords[i]))
                distance_matrix[j][i] = distance_matrix[i][j]

        return distance_matrix

    def angles(self, atomcoords_index = -1):
        """Returns a matrix with Angles(in radians) between atoms i, j, k at index [i][j][k]"""

        coords = self.data.atomcoords[atomcoords_index]
        number_of_atoms = self.data.natom
        angle_matrix = np.zeros([number_of_atoms, number_of_atoms, number_of_atoms])
    
        for i in range(0, number_of_atoms - 1):
            for j in range(i + 1, number_of_atoms - 1):
                for k in range(j + 1, number_of_atoms):
                    angle_matrix[i][j][k] = np.arccos(((coords[j, 0]**2 + coords[j, 1]**2 + coords[j, 2]**2) - (((coords[i, 0] + coords[k, 0]) * coords[j, 0]) + ((coords[i, 1] + coords[k, 1]) * coords[j, 1]) + ((coords[i, 2] + coords[k, 2]) * coords[j, 2]))) / (np.linalg.norm(np.subtract(coords[i], coords[k])) * np.linalg.norm(np.subtract(coords[j], coords[k]))))
                    angle_matrix[i][j][k] = angle_matrix[k][j][i]
        return angle_matrix

    def dihedral_angles(self, atomcoords_index = -1):
        """Returns a matrix with Dihedral angles(in radians) between the planes containing the atoms i, j and k, l at index [i][j][k][l]"""
        """Based on this formula -> [https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python]
        1 sqrt, 1 cross product"""

        coords = self.data.atomcoords[atomcoords_index]
        number_of_atoms = self.data.natom
        dihedral_matrix = np.zeros([number_of_atoms, number_of_atoms, number_of_atoms, number_of_atoms])
        
        for i in range(0, number_of_atoms - 1):
            for j in range(i + 1, number_of_atoms - 1):
                for k in range(j + 1, number_of_atoms - 1):
                    for l in range(k + 1, number_of_atoms):
                        b0 = -1.0*(coords[j] - coords[i])
                        b1 = coords[k] - coords[j]
                        b2 = coords[l] - coords[k]

                        # normalize b1 so that it does not influence magnitude of vector
                        # rejections that come next
                        b1 /= np.linalg.norm(b1)

                        # vector rejections
                        # v = projection of b0 onto plane perpendicular to b1
                        # w = projection of b2 onto plane perpendicular to b1
                        v = b0 - np.dot(b0, b1)*b1
                        w = b2 - np.dot(b2, b1)*b1

                        # angle between v and w in a plane is the torsion angle
                        x = np.dot(v, w)
                        y = np.dot(np.cross(b1, v), w)

                        dihedral_matrix[i][j][k][l] = np.arctan2(y, x)
                        
                        # All permutations
                        dihedral_matrix[i][j][l][k] = dihedral_matrix[i][j][k][l]
                        dihedral_matrix[j][i][k][l] = dihedral_matrix[i][j][k][l]
                        dihedral_matrix[j][i][l][k] = dihedral_matrix[i][j][k][l]
                        dihedral_matrix[k][l][i][j] = dihedral_matrix[i][j][k][l]
                        dihedral_matrix[k][l][j][i] = dihedral_matrix[i][j][k][l]
                        dihedral_matrix[l][k][i][j] = dihedral_matrix[i][j][k][l]
                        dihedral_matrix[l][k][j][i] = dihedral_matrix[i][j][k][l]

        return dihedral_matrix

    def repulsion_energy(self, atomcoords_index=-1):
        """Return the nuclear repulsion energy."""
        nre = 0.0
        for i in range(self.data.natom):
            ri = self.data.atomcoords[atomcoords_index][i]
            zi = self.data.atomnos[i]
            for j in range(i+1, self.data.natom):
                rj = self.data.atomcoords[0][j]
                zj = self.data.atomnos[j]
                d = np.linalg.norm(ri-rj)
                nre += zi*zj/d
        return nre

    def center_of_mass(self, atomcoords_index=-1):
        """Return the center of mass."""
        charges = self.data.atomnos
        coords = self.data.atomcoords[atomcoords_index]
        masses = get_isotopic_masses(charges)

        mwc = coords * masses[:, np.newaxis]
        numerator = np.sum(mwc, axis=0)
        denominator = np.sum(masses)

        return numerator / denominator

    def moment_of_inertia_tensor(self, atomcoords_index=-1):
        """Return the moment of inertia tensor."""
        charges = self.data.atomnos
        coords = self.data.atomcoords[atomcoords_index]
        masses = get_isotopic_masses(charges)

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

    def principal_moments_of_inertia(self, units='amu_bohr_2'):
        """Return the principal moments of inertia in 3 kinds of units:
        1. [amu][bohr]^2
        2. [amu][angstrom]^2
        3. [g][cm]^2
        and the principal axes.
        """
        choices = ('amu_bohr_2', 'amu_angstrom_2', 'g_cm_2')
        units = units.lower()
        if units not in choices:
            raise ValueError("Invalid units, pick one of {}".format(choices))
        moi_tensor = self.moment_of_inertia_tensor()
        principal_moments, principal_axes = np.linalg.eigh(moi_tensor)
        if units == "amu_bohr_2":
            _check_scipy(_found_scipy)
            bohr2ang = spc.value("atomic unit of length") / spc.angstrom
            conv = 1 / bohr2ang ** 2
        elif units == "amu_angstrom_2":
            conv = 1
        elif units == "g_cm_2":
            _check_scipy(_found_scipy)
            amu2g = spc.value("unified atomic mass unit") * spc.kilo
            conv = amu2g * (spc.angstrom / spc.centi) ** 2
        return conv * principal_moments, principal_axes

    def rotational_constants(self, units='ghz'):
        """Compute the rotational constants in 1/cm or GHz."""
        choices = ('invcm', 'ghz')
        units = units.lower()
        if units not in choices:
            raise ValueError("Invalid units, pick one of {}".format(choices))
        principal_moments = self.principal_moments_of_inertia("amu_angstrom_2")[0]
        _check_scipy(_found_scipy)
        bohr2ang = spc.value('atomic unit of length') / spc.angstrom
        xfamu = 1 / spc.value('electron mass in u')
        xthz = spc.value('hartree-hertz relationship')
        rotghz = xthz * (bohr2ang ** 2) / (2 * xfamu * spc.giga)
        if units == 'ghz':
            conv = rotghz
        if units == 'invcm':
            ghz2invcm = spc.giga * spc.centi / spc.c
            conv = rotghz * ghz2invcm
        return conv / principal_moments


del find_package
