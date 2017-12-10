# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A reader for XYZ (Cartesian coordinate) files."""

from ..parser.data import ccData
from ..parser.utils import PeriodicTable

from . import filereader


class XYZ(filereader.Reader):
    """A reader for XYZ (Cartesian coordinate) files."""

    def __init__(self, source, *args, **kwargs):
        super(XYZ, self).__init__(source, *args, **kwargs)
        self.pt = PeriodicTable()

    def generate_repr(self):
        """Convert the raw contents of the source into the internal representation."""

        assert hasattr(self, 'filecontents')

        it = iter(self.filecontents.splitlines())

        # Ordering of lines:
        # 1. number of atoms
        # 2. comment line
        # 3. line of at least 4 columns: 1 is atomic symbol (str), 2-4 are atomic coordinates (float)
        #    repeat for numver of atoms
        # (4. optional blank line)
        # repeat for multiple sets of coordinates (TODO)

        line = next(it)
        if line.strip() == '':
            line = next(it)
        tokens = line.split()
        assert len(tokens) >= 1
        natom = int(tokens[0])

        comment = next(it)

        lines = []
        for line in it:
            tokens = line.split()
            assert len(tokens) >= 4
            lines.append(tokens)
        assert len(lines) == natom

        atomsyms = [line[0] for line in lines]
        atomnos = [self.pt.number[atomsym] for atomsym in atomsyms]
        atomcoords = [line[1:4] for line in lines]
        # Everything beyond the fourth column is ignored.

        attributes = {
            'natom': natom,
            'atomnos': atomnos,
            'atomcoords': [atomcoords],
        }

        data = ccData(attributes)

        return data
