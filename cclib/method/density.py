# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Building the density matrix from data parsed by cclib."""

import logging
import random

import numpy

from cclib.method.calculationmethod import Method


class Density(Method):
    """Calculate the density matrix"""
    def __init__(self, data, progress=None, loglevel=logging.INFO,
                 logname="Density"):
        super().__init__(data, progress, loglevel, logname)

    def __str__(self):
        """Return a string representation of the object."""
        return f"Density matrix of {self.data}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'Density matrix("{self.data}")'

    def calculate(self, fupdate=0.05):
        """Calculate the density matrix."""

        # Do we have the needed info in the data object?
        if not hasattr(self.data, "mocoeffs"):
            self.logger.error("Missing mocoeffs")
            return False
        if not hasattr(self.data,"nbasis"):
            self.logger.error("Missing nbasis")
            return False
        if not hasattr(self.data,"homos"):
            self.logger.error("Missing homos")
            return False

        self.logger.info("Creating attribute density: array[3]")
        size = self.data.nbasis
        unrestricted = (len(self.data.mocoeffs) == 2)

        #determine number of steps, and whether process involves beta orbitals
        nstep = self.data.homos[0] + 1
        if unrestricted:
            self.density = numpy.zeros([2, size, size], "d")
            nstep += self.data.homos[1] + 1
        else:
            self.density = numpy.zeros([1, size, size], "d")

        #intialize progress if available
        if self.progress:
            self.progress.initialize(nstep)

        step = 0
        for spin in range(len(self.data.mocoeffs)):

            for i in range(self.data.homos[spin] + 1):

                if self.progress and random.random() < fupdate:
                    self.progress.update(step, "Density Matrix")

                col = numpy.reshape(self.data.mocoeffs[spin][i], (size, 1))
                colt = numpy.reshape(col, (1, size))

                tempdensity = numpy.dot(col, colt)
                self.density[spin] = numpy.add(self.density[spin],
                                                 tempdensity)

                step += 1

        if not unrestricted: #multiply by two to account for second electron
            self.density[0] = numpy.add(self.density[0], self.density[0])

        if self.progress:
            self.progress.update(nstep, "Done")

        return True #let caller know we finished density
