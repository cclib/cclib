# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A reader for XYZ (Cartesian coordinate) files."""

from cclib.io import filereader
from cclib.parser.data import ccData
from cclib.parser.utils import PeriodicTable


class XYZ(filereader.Reader):
    """A reader for XYZ (Cartesian coordinate) files."""

    def __init__(self, source, *args, **kwargs):
        super().__init__(source, *args, **kwargs)

        self.pt = PeriodicTable()

    def parse(self):
        super().parse()

        self.generate_repr()

        return self.data

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
        # repeat for multiple sets of coordinates

        all_atomcoords = []
        comments = []

        while True:

            try:
                line = next(it)
                if line.strip() == '':
                    line = next(it)
                tokens = line.split()
                assert len(tokens) >= 1
                natom = int(tokens[0])

                comments.append(next(it))

                lines = []
                for _ in range(natom):
                    line = next(it)
                    tokens = line.split()
                    assert len(tokens) >= 4
                    lines.append(tokens)
                assert len(lines) == natom

                atomsyms = [line[0] for line in lines]
                atomnos = [self.pt.number[atomsym] for atomsym in atomsyms]
                atomcoords = [line[1:4] for line in lines]
                # Everything beyond the fourth column is ignored.
                all_atomcoords.append(atomcoords)

            except StopIteration:
                break

        attributes = {
            'natom': natom,
            'atomnos': atomnos,
            'atomcoords': all_atomcoords,
            'metadata': {"comments": comments},
        }

        self.data = ccData(attributes)
