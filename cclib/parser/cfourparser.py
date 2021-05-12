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
    -1: 2,
    2: 3,
    -2: 3,
    3: 4,
    -3: 4
}


def _shell_to_orbitals(type, offset):
    """Convert a Fchk shell type and offset to a list of string representations.

    For example, shell type = -2 corresponds to d orbitals (spherical basis) with
    an offset = 1 would correspond to the 4d orbitals, so this function returns
    `['4D1', '4D2', '4D3', '4D4', '4D5']`.
    """

    return ['{}{}'.format(SHELL_START[type] + offset, x) for x in SHELL_ORBITALS[type]]


class CFOUR(logfileparser.Logfile):
    """ A CFOUR output file
    """

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(CFOUR, self).__init__(logname="CFOUR", *args, **kwargs)
        self.start = True

    def __str__(self):
        """Return a string representation of the object."""
        return "CFOUR output file {}".format(self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'CFOUR("{}")'.format(self.filename)

    def normalisesym(self, symlabel):
        """Just return label"""
        return symlabel

    def extract(self, inputfile, line):

        # just opened file, skip first line to get basis
        if 'Version' in line:
            package_version = line.split()[1]
            self.metadata["package_version"] = package_version

        if '@GETXYZ-I' in line:
            self.natom = int(line.split()[1])
            print(self.natom)

        if 'Coordinates (in bohr)' in line:
            # skip labeling lines
            #         Coordinates used in calculation (QCOMP)
            # ----------------------------------------------------------------
            # Z-matrix   Atomic            Coordinates (in bohr)
            #  Symbol    Number           X              Y              Z
            # ----------------------------------------------------------------
            #     O         8         0.00000000     0.00000000     0.00000000
            for i in range(2):
                 line = next(inputfile)
            coord_block = numpy.array(self._parse_block(inputfile, self.natom,str, 'parsing coordinate blocks'))
            self.atomnos = numpy.array(coord_block[:,1],dtype=int)
            self.atomcoords = numpy.array(coord_block[:,2:],dtype=float)
            
    def _parse_block(self, inputfile, count, type, msg):
        """
        pases through a block
        PARAMETERS
        ---
        inputfile: str
          input file currently being parsed
        count: int
          the number of lines to parse
        type: object type definition
          type of data returned
        msg: str
          message to pass to the update log.
        """
        property = []
        while len(property) < count :
            self.updateprogress(inputfile, msg, self.fupdate)
            line = next(inputfile)
            property.append([ type(x) for x in line.split()])
            print(property)
        print(numpy.shape(property))
        return property
