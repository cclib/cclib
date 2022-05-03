# -*- coding: utf-8 -*-
#
# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2017, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Analyses related to orbitals."""

import logging

import numpy

from cclib.method.calculationmethod import Method


class Orbitals(Method):
    """A class for orbital related methods."""

    def __init__(self, data, progress=None, \
                 loglevel=logging.INFO, logname="Log"):
        super().__init__(data, progress, loglevel, logname)
        self.required_attrs = ('mocoeffs','moenergies','homos')
        self.fragresults = None

    def __str__(self):
        """Return a string representation of the object."""
        return "Orbitals"

    def __repr__(self):
        """Return a representation of the object."""
        return "Orbitals"

    def closed_shell(self):
        """Return Boolean indicating if system is closed shell."""

        # If there are beta orbitals, we can assume the system is closed
        # shell if the orbital energies are identical within numerical accuracy.
        if len(self.data.mocoeffs) == 2:
            precision = 10e-6
            return numpy.allclose(*self.data.moenergies, atol=precision)

        # Restricted open shell will have one set of MOs but two HOMO indices,
        # and the indices should be different (otherwise it's still closed shell).
        if len(self.data.homos) == 2 and self.data.homos[0] != self.data.homos[1]:
            return False

        return True
