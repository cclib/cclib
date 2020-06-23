# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculation of Becke charges using bridge to horton."""

import random

import numpy

from cclib.parser.utils import find_package
from cclib.method.population import Population
from cclib.method.density import Density
from cclib.bridge.cclib2horton import makehorton

class Becke(Method):
    """Becke Population Analysis"""

    def __init__(self, *args):

        # Call the __init__ method of the superclass.
        super(Becke, self).__init__(logname="Becke Population Analysis", *args)
        
        # horton 3 does not have denspart released yet
        try:
            from horton import __version__
        except:
            raise ImportError(
                "You must have horton 2 installed to use this function."
            )
        else:
            from horton import BeckeMolGrid
            from horton import getgobasis
            from horton.matrix.dense import DenseTwoIndex

    def __str__(self):
        """Return a string representation of the object."""
        return "Becke charges of {}".format(self.data)

    def __repr__(self):
        """Return a representation of the object."""
        return "Becke({})".format(self.data)

    def calculate(self, indices=None, fupdate=0.05):
        """Perform a Becke population analysis."""
        # Based on horton documentation
        # (https://theochem.github.io/horton/2.1.0/user_postproc_aim.html#horton-wpart-py-aim-analysis-based-on-a-wavefunction-file)
        # Use bridging function to create horton IOData object
        iodat = makehorton(self.data)
        
        # Calculate density matrix and convert to object used in horton
        d = Density(self.data)
        d.calculate()
        den = DenseTwoIndex(len(iodat.orb_alpha))
        
        if len(d.density) == 1:
            # restricted case
            den._array = d.density[0]
        elif len(d.density) == 2:
            # unrestricted case
            den._array = d.density[0] + d.density[1]
        
        # Create integration grid
        grid = BeckeMolGrid(iodat.coordinates, iodat.numbers, iodat.pseudo_numbers, mode='only')
        
        # Define Gaussian basis set
        # -- This needs to be modified so that it first generates GOBasisFamily instance from cclib's gbasis attribute and passes it to 'default' argument
        gob = get_gobasis(iodat.coordinates, iodat.numbers, default = 'STO-3G')

        # Calculate grid density
        moldens = gob.compute_grid_density_dm(den, grid.points)
        
        # Partition charges
        self.logger.info("Creating fragcharges: array[1]")
        wpart = BeckeWPart(ht.coordinates, ht.numbers, ht.pseudo_numbers, grid, moldens, local=True)
        wpart.do_charges()
        
        # Save charges into applicable format
        self.fragcharges = wpart['charges']

        return True
