# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Formatted Checkpoint files"""


from __future__ import print_function

import re

import numpy

from cclib.parser import data
from cclib.parser import logfileparser
from cclib.parser import utils

class FChk(logfileparser.Logfile):
    """A Formatted checkpoint file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(FChk, self).__init__(logname="FCHK", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "Formatted checkpoint file %s" % self.filename

    def __repr__(self):
        """Return a representation of the object."""
        return 'FCHK("%s")' % self.filename

    def normalisesym(self, symlabel):
        """Just return label"""
        return symlabel

    def extract(self, inputfile, line):

        if line[0:14] == 'Atomic numbers':
            self.updateprogress(inputfile, "Basic Information", self.fupdate)

            self.natom = int(line.split()[-1])
            atomnos = self._parse_block(inputfile, self.natom, int)
            self.set_attribute('atomnos', atomnos)

        if line[0:19] == 'Number of electrons':
            alpha = next(inputfile)
            alpha_homo = int(alpha.split()[-1]) - 1

            beta = next(inputfile)
            beta_homo = int(beta.split()[-1]) - 1

            self.homos = [alpha_homo, beta_homo]

        if line[0:29] == 'Current cartesian coordinates':
            count = int(line.split()[-1])
            assert count % 3 == 0

            coords = numpy.array(self._parse_block(inputfile, count, float))
            coords.shape = (1, int(count / 3), 3)
            self.set_attribute('atomcoords', coords)

        if line[0:25] == 'Number of basis functions':
            self.nbasis = int(line.split()[-1])

        if line[0:14] == 'Overlap Matrix':
            self.updateprogress(inputfile, "Overlap Matrix", self.fupdate)

            count = int(line.split()[-1])

            # triangle matrix, with number of elements in a row:
            # 1 + 2 + 3 + .... + self.nbasis
            assert count == (self.nbasis + 1) * self.nbasis / 2
            raw_overlaps = self._parse_block(inputfile, count, float)

            # now turn into matrix
            overlaps = numpy.zeros((self.nbasis, self.nbasis))
            raw_index = 0
            for row in range(self.nbasis):
                for col in range(row, self.nbasis):
                    overlaps[row, col] = raw_overlaps[raw_index]
                    overlaps[col, row] = raw_overlaps[raw_index]
                    raw_index += 1

            self.set_attribute('aooverlaps', overlaps)

    def _parse_block(self, inputfile, count, type):
        atomnos = []
        while len(atomnos) < count :
            line = next(inputfile)
            atomnos.extend([ type(x) for x in line.split()])
        return atomnos