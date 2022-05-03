# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""C-squared population analysis."""

import random

import numpy

from cclib.method.population import Population


class CSPA(Population):
    """The C-squared population analysis."""

    # Overlaps are not required for CSPA.
    overlap_attributes = ()

    def __init__(self, *args):
        super().__init__(logname="CSPA", *args)

    def __str__(self):
        """Return a string representation of the object."""
        return f"CSPA of {self.data}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'CSPA("{self.data}")'

    def calculate(self, indices=None, fupdate=0.05):
        """Perform the C squared population analysis.

        Inputs:
           indices - list of lists containing atomic orbital indices of fragments
        """
        self.logger.info("Creating attribute aoresults: array[3]")

        # Determine number of steps, and whether process involves beta orbitals.
        unrestricted = (len(self.data.mocoeffs)==2)
        nbasis = self.data.nbasis
        self.aoresults = []
        alpha = len(self.data.mocoeffs[0])
        self.aoresults.append(numpy.zeros([alpha, nbasis], "d"))
        nstep = alpha
        if unrestricted:
            beta = len(self.data.mocoeffs[1])
            self.aoresults.append(numpy.zeros([beta, nbasis], "d"))
            nstep += beta

        # Intialize progress if available.
        if self.progress:
            self.progress.initialize(nstep)

        step = 0
        for spin in range(len(self.data.mocoeffs)):

            for i in range(len(self.data.mocoeffs[spin])):

                if self.progress and random.random() < fupdate:
                    self.progress.update(step, "C^2 Population Analysis")

                submocoeffs = self.data.mocoeffs[spin][i]
                scale = numpy.inner(submocoeffs, submocoeffs)
                tempcoeffs = numpy.multiply(submocoeffs, submocoeffs)
                tempvec = tempcoeffs/scale
                self.aoresults[spin][i] = numpy.divide(tempcoeffs, scale).astype("d")

                step += 1

        if self.progress:
            self.progress.update(nstep, "Done")

        retval = super().partition(indices)

        if not retval:
            self.logger.error("Error in partitioning results")
            return False

        self.logger.info("Creating fragcharges: array[1]")
        size = len(self.fragresults[0][0])
        self.fragcharges = numpy.zeros([size], "d")
        alpha = numpy.zeros([size], "d")
        if unrestricted:
            beta = numpy.zeros([size], "d")

        for spin in range(len(self.fragresults)):

            for i in range(self.data.homos[spin] + 1):

                temp = numpy.reshape(self.fragresults[spin][i], (size,))
                self.fragcharges = numpy.add(self.fragcharges, temp)
                if spin == 0:
                    alpha = numpy.add(alpha, temp)
                elif spin == 1:
                    beta = numpy.add(beta, temp)

        if not unrestricted:
            self.fragcharges = numpy.multiply(self.fragcharges, 2)
        else:
            self.logger.info("Creating fragspins: array[1]")
            self.fragspins = numpy.subtract(alpha, beta)

        return True
