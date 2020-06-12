# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculation of Bickelhaupt population analysis based on data parsed by cclib."""

import random

import numpy

from cclib.method.population import Population


class Bickelhaupt(Population):
    """Bickelhaupt population analysis."""

    def __init__(self, *args):

        # Call the __init__ method of the superclass.
        super(Bickelhaupt, self).__init__(logname="Bickelhaupt Population Analysis", *args)

    def __str__(self):
        """Return a string representation of the object."""
        return "Bickelhaupt charges of {}".format(self.data)

    def __repr__(self):
        """Return a representation of the object."""
        return "Bickelhaupt({})".format(self.data)

    def calculate(self, indices=None, fupdate=0.05):
        """Perform a Bickelhaupt population analysis."""

        # Determine number of steps, and whether process involves beta orbitals.
        self.logger.info("Creating attribute aoresults: [array[2]]")
        nbasis = self.data.nbasis
        alpha = len(self.data.mocoeffs[0])
        self.aoresults = [ numpy.zeros([alpha, nbasis], "d") ]
        w = [ numpy.zeros([self.data.nbasis, self.data.nbasis], "d") ]
        nstep = alpha + nbasis
        unrestricted = (len(self.data.mocoeffs) == 2)
        if unrestricted:
            beta = len(self.data.mocoeffs[1])
            self.aoresults.append(numpy.zeros([beta, nbasis], "d"))
            w.append(numpy.zeros([self.data.nbasis, self.data.nbasis], "d"))
            nstep += (beta + nbasis)

        # Intialize progress if available.
        if self.progress:
            self.progress.initialize(nstep)

        step = 0

        # First determine the weight matrix
        for spin in range(len(self.data.mocoeffs)):
            for i in range(self.data.nbasis):
                if self.progress and random.random() < fupdate:
                    self.progress.update(step, "Bickelhaupt Population Analysis")
                for j in range(self.data.nbasis):
                    # W_{ab} = 2 * \sum_k c_{ak}^2 / (\sum_i c_{ai}^2 + \sum_j c_{bj}^2)
                    # The dimension of W is (nbasis)x(nbasis)
                    aTerm = numpy.dot(self.data.mocoeffs[spin][0:(self.data.homos[spin]+1),i], self.data.mocoeffs[spin][0:(self.data.homos[spin]+1),i])
                    bTerm = numpy.dot(self.data.mocoeffs[spin][0:(self.data.homos[spin]+1),j], self.data.mocoeffs[spin][0:(self.data.homos[spin]+1),j])
                    w[spin][i][j] = 2 * aTerm / (aTerm + bTerm)
                step += 1

        for spin in range(len(self.data.mocoeffs)):
            for i in range(len(self.data.mocoeffs[spin])): # for each mo
                for a in range(self.data.nbasis): # for each basis
                    # X_{ai} = \sum_b w_{ab} c_{ai} c_{bi} S_{ab}
                    if hasattr(self.data, "aooverlaps"):
                        overlaps = self.data.aooverlaps[a]
                    # handle spin-unrestricted beta case
                    elif hasattr(self.data, "fooverlaps2") and spin == 1:
                        overlaps = self.data.fooverlaps2[a]

                    elif hasattr(self.data, "fooverlaps"):
                        overlaps = self.data.fooverlaps[a]
                    
                    temp = numpy.multiply(w[spin][a], self.data.mocoeffs[spin][i][a])
                    temp = numpy.multiply(temp, self.data.mocoeffs[spin][i])
                    self.aoresults[spin][i][a] = numpy.sum(numpy.multiply(temp, overlaps))

                step += 1

        if self.progress:
            self.progress.update(nstep, "Done")

        retval = super(Bickelhaupt, self).partition(indices)

        if not retval:
            self.logger.error("Error in partitioning results")
            return False

        # Create array for Bickelhaupt charges.
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
