# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Population analyses based on cclib data."""

import logging

import numpy

from cclib.method.calculationmethod import Method, MissingAttributeError


class Population(Method):
    """An abstract base class for population-type methods."""

    # All of these are typically required for population analyses.
    required_attrs = ('homos', 'mocoeffs', 'nbasis')

    # At least one of these are typically required.
    overlap_attributes = ('aooverlaps', 'fooverlaps')

    def __init__(self, data, progress=None, \
                 loglevel=logging.INFO, logname="Log"):
        super().__init__(data, progress, loglevel, logname)

        self.fragresults = None

    def __str__(self):
        """Return a string representation of the object."""
        return "Population"

    def __repr__(self):
        """Return a representation of the object."""
        return "Population"

    def _check_required_attributes(self):
        super()._check_required_attributes()
        
        if self.overlap_attributes and not any(hasattr(self.data, a) for a in self.overlap_attributes):
            raise MissingAttributeError(
                    'Need overlap matrix (aooverlaps or fooverlaps attribute) for Population methods')

    def partition(self, indices=None):

        if not hasattr(self, "aoresults"):
            self.calculate()

        if not indices:

            # Build list of groups of orbitals in each atom for atomresults.
            if hasattr(self.data, "aonames"):
                names = self.data.aonames
            elif hasattr(self.data, "fonames"):
                names = self.data.fonames

            atoms = []
            indices = []

            name = names[0].split('_')[0]
            atoms.append(name)
            indices.append([0])

            for i in range(1, len(names)):
                name = names[i].split('_')[0]
                try:
                    index = atoms.index(name)
                except ValueError: #not found in atom list
                    atoms.append(name)
                    indices.append([i])
                else:
                    indices[index].append(i)

        natoms = len(indices)
        nmocoeffs = len(self.aoresults[0])

        # Build results numpy array[3].
        alpha = len(self.aoresults[0])
        results = []
        results.append(numpy.zeros([alpha, natoms], "d"))

        if len(self.aoresults) == 2:
            beta = len(self.aoresults[1])
            results.append(numpy.zeros([beta, natoms], "d"))

        # For each spin, splice numpy array at ao index,
        #   and add to correct result row.
        for spin in range(len(results)):

            for i in range(natoms): # Number of groups.

                for j in range(len(indices[i])): # For each group.

                    temp = self.aoresults[spin][:, indices[i][j]]
                    results[spin][:, i] = numpy.add(results[spin][:, i], temp)

        self.logger.info("Saving partitioned results in fragresults: [array[2]]")
        self.fragresults = results

        return True
