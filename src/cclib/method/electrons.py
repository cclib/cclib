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

"""Calculate properties for electrons."""

import logging

import numpy

from .calculationmethod import Method


class Electrons(Method):
    """A container for methods pertaining to electrons."""
    
    def __init__(self, data, progress=None, loglevel=logging.INFO, logname="Log"):

        super(Electrons, self).__init__(data, progress, loglevel, logname)
        
    def __str__(self):
        """Returns a string representation of the object."""
        return "Electrons"

    def __repr__(self):
        """Returns a representation of the object."""
        return "Electrons"
    
    def count(self, core=False):
        """Returns the electron count in system.
        
        Normally returns electrons used in calculation, but will include
        core electrons in pseudopotentials if core is True.
        """
        nelectrons = sum(self.data.atomnos) - self.data.charge
        if core:
            nelectrons += sum(self.data.coreelectrons)
        return nelectrons


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
