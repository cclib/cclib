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

SHELL_ORBITALS = {
    0: ['S'],
    1: ['PX', 'PY', 'PZ'],
    -1: ['S', 'PX', 'PY', 'PZ'],
    2: ['D1', 'D2', 'D3', 'D4', 'D5', 'D6'],
    -2: ['D1', 'D2', 'D3', 'D4', 'D5'],
    3:  ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10'],
    -3: ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7'],
    4: ['G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12','G13'],
    -4: ['G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9']

}

SHELL_START = {
    0: 1,
    1: 2,
    -1: 1,
    2: 3,
    -2: 3,
    3: 4,
    -3: 4
}


def shell_to_orbitals(type, offset):
    return ['{}{}'.format(SHELL_START[type] + offset, x) for x in SHELL_ORBITALS[type]]


class FChk(logfileparser.Logfile):
    """A Formatted checkpoint file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(FChk, self).__init__(logname="FChk", *args, **kwargs)

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

        if line[0:6] == 'Charge':
            self.set_attribute('charge', int(line.split()[-1]))

        if line[0:12] == 'Multiplicity':
            self.set_attribute('mult', int(line.split()[-1]))

        if line[0:14] == 'Atomic numbers':
            self.natom = int(line.split()[-1])
            atomnos = self._parse_block(inputfile, self.natom, int, 'Basic Information')
            self.set_attribute('atomnos', atomnos)

        if line[0:19] == 'Number of electrons':
            alpha = next(inputfile)
            alpha_homo = int(alpha.split()[-1]) - 1

            beta = next(inputfile)
            beta_homo = int(beta.split()[-1]) - 1

            self.set_attribute('homos', [alpha_homo, beta_homo])

        if line[0:29] == 'Current cartesian coordinates':
            count = int(line.split()[-1])
            assert count % 3 == 0

            coords = numpy.array(self._parse_block(inputfile, count, float, 'Coordinates'))
            coords.shape = (1, int(count / 3), 3)
            self.set_attribute('atomcoords', coords)

        if line[0:25] == 'Number of basis functions':
            self.set_attribute('nbasis', int(line.split()[-1]))

        if line[0:14] == 'Overlap Matrix':
            count = int(line.split()[-1])

            # triangle matrix, with number of elements in a row:
            # 1 + 2 + 3 + .... + self.nbasis
            assert count == (self.nbasis + 1) * self.nbasis / 2
            raw_overlaps = self._parse_block(inputfile, count, float, 'Overlap Matrix')

            # now turn into matrix
            overlaps = numpy.zeros((self.nbasis, self.nbasis))
            raw_index = 0
            for row in range(self.nbasis):
                for col in range(row, self.nbasis):
                    overlaps[row, col] = raw_overlaps[raw_index]
                    overlaps[col, row] = raw_overlaps[raw_index]
                    raw_index += 1

            self.set_attribute('aooverlaps', overlaps)

        if line[0:31] == 'Number of independent functions':
            self.set_attribute('nmo', int(line.split()[-1]))

        if line[0:21] == 'Alpha MO coefficients':
            count = int(line.split()[-1])
            assert count == self.nbasis * self.nmo

            coeffs = numpy.array(self._parse_block(inputfile, count, float, 'Alpha Coefficients'))
            coeffs.shape = (self.nmo, self.nbasis)
            self.set_attribute('mocoeffs', [coeffs])

        if line[0:22] == 'Alpha Orbital Energies':
            count = int(line.split()[-1])
            assert count == self.nmo

            energies = numpy.array(self._parse_block(inputfile, count, float, 'Alpha MO Energies'))
            self.set_attribute('moenergies', [energies])

        if line[0:20] == 'Beta MO coefficients':
            count = int(line.split()[-1])
            assert count == self.nbasis * self.nmo

            coeffs = numpy.array(self._parse_block(inputfile, count, float, 'Beta Coefficients'))
            coeffs.shape = (self.nmo, self.nbasis)
            self.append_attribute('mocoeffs', coeffs)

        if line[0:21] == 'Beta Orbital Energies':
            count = int(line.split()[-1])
            assert count == self.nmo

            energies = numpy.array(self._parse_block(inputfile, count, float, 'Alpha MO Energies'))
            self.append_attribute('moenergies', energies)

        if line[0:11] == 'Shell types':
            self.parse_aonames(line, inputfile)

    def parse_aonames(self, line, inputfile):
        # e.g.: Shell types                                I   N=          28
        count = int(line.split()[-1])
        shell_types = self._parse_block(inputfile, count, int, 'Atomic Orbital Names')

        # e.g.: Number of primitives per shell             I   N=          28
        next(inputfile)
        self._parse_block(inputfile, count, int, 'Atomic Orbital Names')

        # e.g. Shell to atom map                          I   N=          28
        next(inputfile)
        shell_map = self._parse_block(inputfile, count, int, 'Atomic Orbital Names')
        table = utils.PeriodicTable()

        elements = (table.element[x] for x in self.atomnos)
        atom_labels = ["{}{}".format(y, x) for x, y in enumerate(elements, 1)]

        # get orbitals for first atom
        atom = shell_map[0] - 1
        offset = 0
        aonames = ["{}_{}".format(atom_labels[atom], x) for x in shell_to_orbitals(shell_types[0], offset)]

        # get rest
        for i in range(1, len(shell_types)):
            _type = shell_types[i]
            atom = shell_map[i] - 1
            offset += 1

            # determine if we move to next atom
            if atom != shell_map[i - 1] - 1:
                offset = 0

            # determine if we've changed shell type
            if _type != shell_types[i - 1]:
                offset = 0

            orbitals = ["{}_{}".format(atom_labels[atom], x) for x in shell_to_orbitals(_type, offset)]
            aonames.extend(orbitals)

        assert len(aonames) == self.nbasis, 'Length of aonames != nbasis: {} != {}'.format(len(aonames), self.nbasis)
        self.set_attribute('aonames', aonames)

    def _parse_block(self, inputfile, count, type, msg):
        atomnos = []
        while len(atomnos) < count :
            self.updateprogress(inputfile, msg, self.fupdate)
            line = next(inputfile)
            atomnos.extend([ type(x) for x in line.split()])
        return atomnos
