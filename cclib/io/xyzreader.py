# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A reader for XYZ (Cartesian coordinate) files."""

from cclib.io import filereader
from cclib.parser.data import ccData
from cclib.parser.utils import PeriodicTable


class XYZ(filereader.Reader):
    """A reader for XYZ (Cartesian coordinate) files."""

    def __init__(self, source, *args, **kwargs) -> None:
        super().__init__(source, *args, **kwargs)

        self.pt = PeriodicTable()

    def parse(self) -> ccData:
        super().parse()

        self.generate_repr()

        return self.data

    def generate_repr(self) -> None:
        """Convert the raw contents of the source into the internal representation."""

        assert hasattr(self, "filecontents")

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
        atomnos = None
        natom = None

        while True:
            line = next((line for line in it if line.strip() != ""), None)
            if line is None:
                break

            tokens = line.split()
            if len(tokens) < 1:
                raise ValueError("XYZ atom-count line is empty")
            block_natom = int(tokens[0])
            if block_natom < 0:
                raise ValueError("XYZ atom count cannot be negative")

            try:
                comment = next(it)
            except StopIteration as exc:
                raise ValueError("XYZ block is missing its comment line") from exc

            lines = []
            for _ in range(block_natom):
                try:
                    line = next(it)
                except StopIteration as exc:
                    raise ValueError("XYZ block is truncated before all atoms are read") from exc
                tokens = line.split()
                if len(tokens) < 4:
                    raise ValueError("XYZ atom line must contain an element and three coordinates")
                lines.append(tokens)

            atomsyms = [line[0] for line in lines]
            block_atomnos = [self.pt.number[atomsym] for atomsym in atomsyms]
            if atomnos is None:
                atomnos = block_atomnos
                natom = block_natom
            elif block_natom != natom:
                raise ValueError("XYZ trajectory blocks must contain the same number of atoms")
            atomcoords = [line[1:4] for line in lines]
            comments.append(comment)
            # Everything beyond the fourth column is ignored.
            all_atomcoords.append(atomcoords)

        if natom is None:
            raise ValueError("XYZ file does not contain a geometry block")

        attributes = {
            "natom": natom,
            "atomnos": atomnos,
            "atomcoords": all_atomcoords,
            "metadata": {"comments": comments},
        }

        self.data = ccData(attributes)
