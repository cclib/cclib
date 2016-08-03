# -*- coding: utf-8 -*-
#
# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2016, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Analyses related to orbitals."""

import logging

import numpy

from .calculationmethod import Method


class Orbitals(Method):
    """A class for orbital related methods."""
    
    def __init__(self, data, progress=None, \
                 loglevel=logging.INFO, logname="Log"):

        # Call the __init__ method of the superclass.
        super(Orbitals, self).__init__(data, progress, loglevel, logname)
        self.fragresults = None
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Orbitals"

    def __repr__(self):
        """Return a representation of the object."""
        return "Orbitals"

    def closed_shell(self):
        """Return Boolean indicating if system is closed shell."""

        # Restricted calculation will have one set of orbitals and are
        # therefore closed shell by definition.
        if len(self.data.mocoeffs) == 1:
            assert len(self.data.moenergies) == 1
            return True

        # If there are beta orbitals, we can assume the system is
        # closed shell if the orbital energies are identical within
        # numerical accuracy.
        precision = 10e-6
        return numpy.allclose(*self.data.moenergies, atol=precision)


if __name__ == "__main__":
    import doctest, orbitals
    doctest.testmod(orbitals, verbose=False)
