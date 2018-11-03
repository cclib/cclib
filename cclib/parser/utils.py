# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Utilities often used by cclib parsers and scripts"""

import sys
import numpy


# Define for any Python version <= 3.3,
# See https://github.com/kachayev/fn.py/commit/391824c43fb388e0eca94e568ff62cc35b543ecb
if sys.version_info.major == 2 or sys.version_info.minor <= 3:
    import operator
    def accumulate(iterable, func=operator.add):
        """Return running totals"""
        # accumulate([1,2,3,4,5]) --> 1 3 6 10 15
        # accumulate([1,2,3,4,5], operator.mul) --> 1 2 6 24 120
        it = iter(iterable)
        try:
            total = next(it)
        except StopIteration:
            return
        yield total
        for element in it:
            total = func(total, element)
            yield total
else:
    from itertools import accumulate


def symmetrize(m, use_triangle='lower'):
    """Symmetrize a square NumPy array by reflecting one triangular
    section across the diagonal to the other.
    """

    if use_triangle not in ('lower', 'upper'):
        raise ValueError
    if not len(m.shape) == 2:
        raise ValueError
    if not (m.shape[0] == m.shape[1]):
        raise ValueError

    dim = m.shape[0]

    lower_indices = numpy.tril_indices(dim, k=-1)
    upper_indices = numpy.triu_indices(dim, k=1)

    ms = m.copy()

    if use_triangle == 'lower':
        ms[upper_indices] = ms[lower_indices]
    if use_triangle == 'upper':
        ms[lower_indices] = ms[upper_indices]

    return ms


def convertor(value, fromunits, tounits):
    """Convert from one set of units to another.
    
    Sources:
        NIST 2010 CODATA (http://physics.nist.gov/cuu/Constants/index.html)
        Documentation of GAMESS-US or other programs as noted
    """

    _convertor = {

        "time_au_to_fs": lambda x: x * 0.02418884,
        "fs_to_time_au": lambda x: x / 0.02418884,

        "Angstrom_to_bohr": lambda x: x * 1.8897261245,
        "bohr_to_Angstrom": lambda x: x * 0.5291772109,

        "wavenumber_to_eV":       lambda x: x / 8065.54429,
        "wavenumber_to_hartree":  lambda x: x / 219474.6313708,
        "wavenumber_to_kcal/mol": lambda x: x / 349.7550112,
        "wavenumber_to_kJ/mol":   lambda x: x / 83.5934722814,
        "wavenumber_to_nm":       lambda x: 1e7 / x,
        "wavenumber_to_Hz":       lambda x: x * 29.9792458,

        "eV_to_wavenumber": lambda x: x * 8065.54429,
        "eV_to_hartree":    lambda x: x / 27.21138505,
        "eV_to_kcal/mol":   lambda x: x * 23.060548867,
        "eV_to_kJ/mol":     lambda x: x * 96.4853364596,

        "hartree_to_wavenumber": lambda x: x * 219474.6313708,
        "hartree_to_eV":         lambda x: x * 27.21138505,
        "hartree_to_kcal/mol":   lambda x: x * 627.50947414,
        "hartree_to_kJ/mol":     lambda x: x * 2625.4996398,

        "kcal/mol_to_wavenumber": lambda x: x * 349.7550112,
        "kcal/mol_to_eV":         lambda x: x / 23.060548867,
        "kcal/mol_to_hartree":    lambda x: x / 627.50947414,
        "kcal/mol_to_kJ/mol":     lambda x: x * 4.184,

        "kJ/mol_to_wavenumber": lambda x: x * 83.5934722814,
        "kJ/mol_to_eV":         lambda x: x / 96.4853364596,
        "kJ/mol_to_hartree":    lambda x: x / 2625.49963978,
        "kJ/mol_to_kcal/mol":   lambda x: x / 4.184,
        "nm_to_wavenumber":     lambda x: 1e7 / x,

        # Taken from GAMESS docs, "Further information",
        # "Molecular Properties and Conversion Factors"
        "Debye^2/amu-Angstrom^2_to_km/mol": lambda x: x * 42.255,

        # Conversion for charges and multipole moments.
        "e_to_coulomb":         lambda x: x * 1.602176565  * 1e-19,
        "e_to_statcoulomb":     lambda x: x * 4.80320425   * 1e-10,
        "coulomb_to_e":         lambda x: x * 0.6241509343 * 1e19,
        "statcoulomb_to_e":     lambda x: x * 0.2081943527 * 1e10,
        "ebohr_to_Debye":       lambda x: x * 2.5417462300,
        "ebohr2_to_Buckingham": lambda x: x * 1.3450341749,
        "ebohr2_to_Debye.ang":  lambda x: x * 1.3450341749,
        "ebohr3_to_Debye.ang2": lambda x: x * 0.7117614302,
        "ebohr4_to_Debye.ang3": lambda x: x * 0.3766479268,
        "ebohr5_to_Debye.ang4": lambda x: x * 0.1993134985,

    }

    return _convertor["%s_to_%s" % (fromunits, tounits)](value)


class PeriodicTable(object):
    """Allows conversion between element name and atomic no."""

    def __init__(self):
        self.element = [
            None,
            'H', 'He',
            'Li', 'Be',
            'B', 'C', 'N', 'O', 'F', 'Ne',
            'Na', 'Mg',
            'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
            'K', 'Ca',
            'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
            'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
            'Rb', 'Sr',
            'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
            'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
            'Cs', 'Ba',
            'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
            'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
            'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
            'Fr', 'Ra',
            'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
            'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
            'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo']
        self.number = {}
        for i in range(1, len(self.element)):
            self.number[self.element[i]] = i


class WidthSplitter:
    """Split a line based not on a character, but a given number of field
    widths.
    """

    def __init__(self, widths):
        self.start_indices = [0] + list(accumulate(widths))[:-1]
        self.end_indices = list(accumulate(widths))

    def split(self, line, truncate=True):
        """Split the given line using the field widths passed in on class
        initialization.
        """
        elements = [line[start:end].strip()
                    for (start, end) in zip(self.start_indices, self.end_indices)]
        # Handle lines that contain fewer fields than specified in the
        # widths; they are added as empty strings, so remove them.
        if truncate:
            while len(elements) and elements[-1] == '':
                elements.pop()
        return elements
