# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for MOLDEN format files."""

import os.path
import ntpath
import numpy as np

from . import filewriter
from cclib.parser.data import ccData

class MOLDEN(filewriter.Writer):
    """A writer for chemical JSON (MOLDEN) files."""
    def __init__(self, ccdata, *args, **kwargs):
        """Initialize the chemical JSON writer object.

        Inputs:
          ccdata - An instance of ccData, parsed from a logfile.
        """
        # Call the __init__ method of the superclass
        super(MOLDEN, self).__init__(ccdata, *args, **kwargs)

    def pathname(self, path):
        """
        This function is OS independent and returns the file name irrespective of
        the file path containing forward slash or backward slash - which is valid
        in Windows.
        """
        head, tail = ntpath.split(path)
        return tail or ntpath.basename(head)

    def _coords_from_ccdata(self, index):
        """Create an [Atoms] section using geometry at the given index."""

        element_list = [self.pt.element[Z] for Z in self.ccdata.atomnos]
        atomcoords = self.ccdata.atomcoords[index]
        atomnos = self.ccdata.atomnos
        nos = range(1, self.ccdata.natom)
        
        # element_name number atomic_number x y z
        atom_template = '{:3s} {:3d} {:3d} {:15.10f} {:15.10f} {:15.10f}'
        block = []
        for element, no, atomno, (x, y, z) in zip(element_list, nos, atomnos,\
            atomcoords):
            block.append(atom_template.format(element, no, atomno, x, y, z))

        return block

    def generate_repr(self):
        """Generate the MOLDEN representation of the logfile data."""

        molden_block = ["[MOLDEN FORMAT]"]

        # Title of file
        molden_block.append("[Title]")
        molden_block.append(self.pathname(os.path.splitext(self.jobfilename)[0]))

        # Coordinates for the Electron Density/Molecular orbitals
        # [Atoms] (Angs|AU)
        unit = "Angs"
        molden_block.append("[Atoms] %s" % unit)
        # Last set of coordinates for geometry optimization runs.
        index = -1
        molden_block.extend(self._coords_from_ccdata(index))

        return '\n'.join(molden_block)


if __name__ == "__main__":
    pass
