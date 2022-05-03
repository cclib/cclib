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
        super().__init__(logname="Bickelhaupt Population Analysis", *args)

    def __str__(self):
        """Return a string representation of the object."""
        return f"Bickelhaupt charges of {self.data}"

    def __repr__(self):
        """Return a representation of the object."""
        return f"Bickelhaupt({self.data})"

    def calculate(self, indices=None, fupdate=0.05):
        """Perform a Bickelhaupt population analysis."""
        # Bickelhaupt population analysis uses the relative magnitude of the diagonal terms to
        # partition off-diagonal terms. The weights are therefore calculated using the following
        # formula:
        # W_{ab} = 2 * \sum_k c_{ak}^2 / (\sum_i c_{ai}^2 + \sum_j c_{bj}^2)
        # The dimension of W is (nbasis)x(nbasis).
        # Then, each term is calculated similarly to Mulliken charges, in a form of
        # X_{ai} = \sum_b W_{ab} c_{ai} c_{bi} S_{ab}
        # This agrees with quation 11 on 10.1021/om950966x, with the charge (q) terms substituted
        # with MO coefficients.

        # Determine number of steps, and whether process involves beta orbitals.
        self.logger.info("Creating attribute aoresults: [array[2]]")
        nbasis = self.data.nbasis
        alpha = len(self.data.mocoeffs[0])
        self.aoresults = [numpy.zeros([alpha, nbasis], "d")]
        w = [numpy.zeros([self.data.nbasis, self.data.nbasis], "d")]
        nstep = alpha + nbasis
        unrestricted = len(self.data.mocoeffs) == 2
        if unrestricted:
            beta = len(self.data.mocoeffs[1])
            self.aoresults.append(numpy.zeros([beta, nbasis], "d"))
            w.append(numpy.zeros([self.data.nbasis, self.data.nbasis], "d"))
            nstep += beta + nbasis

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
                    aTerm = numpy.dot(
                        self.data.mocoeffs[spin][0 : (self.data.homos[spin] + 1), i],
                        self.data.mocoeffs[spin][0 : (self.data.homos[spin] + 1), i],
                    )
                    bTerm = numpy.dot(
                        self.data.mocoeffs[spin][0 : (self.data.homos[spin] + 1), j],
                        self.data.mocoeffs[spin][0 : (self.data.homos[spin] + 1), j],
                    )
                    w[spin][i][j] = 2 * aTerm / (aTerm + bTerm)
                step += 1

        for spin in range(len(self.data.mocoeffs)):
            for i in range(len(self.data.mocoeffs[spin])):  # for each mo
                for a in range(self.data.nbasis):  # for each basis
                    if hasattr(self.data, "aooverlaps"):
                        overlaps = self.data.aooverlaps[a]
                    # handle spin-unrestricted beta case
                    elif hasattr(self.data, "fooverlaps2") and spin == 1:
                        overlaps = self.data.fooverlaps2[a]

                    elif hasattr(self.data, "fooverlaps"):
                        overlaps = self.data.fooverlaps[a]

                    self.aoresults[spin][i][a] = numpy.sum(
                        (
                            (w[spin][a] * self.data.mocoeffs[spin][i][a])
                            * self.data.mocoeffs[spin][i]
                        )
                        * overlaps
                    )

                step += 1

        if self.progress:
            self.progress.update(nstep, "Done")

        retval = super().partition(indices)

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
                self.fragcharges = self.fragcharges + temp
                if spin == 0:
                    alpha += temp
                elif spin == 1:
                    beta += temp

        if not unrestricted:
            self.fragcharges = self.fragcharges * 2
        else:
            self.logger.info("Creating fragspins: array[1]")
            self.fragspins = alpha - beta

        return True
