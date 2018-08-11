# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculate properties of nuclei based on data parsed by cclib."""

import logging

import numpy

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
                d = numpy.linalg.norm(ri-rj)
                nre += zi*zj/d
        return nre
