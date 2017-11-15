# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for XYZ (Cartesian coordinate) files."""

from collections import Iterable

from . import filewriter


class XYZ(filewriter.Writer):
    """A writer for XYZ (Cartesian coordinate) files."""

    def __init__(self, ccdata, splitfiles=False,
                 firstgeom=False, lastgeom=False, allgeom=False,
                 *args, **kwargs):
        """Initialize the XYZ writer object.

        Inputs:
          ccdata - An instance of ccData, parse from a logfile.
          splitfiles - Boolean to write multiple files if multiple files are requested. [TODO]
          firstgeom - Boolean to write the first available geometry from the logfile.
          lastgeom - Boolean to write the last available geometry from the logfile.
          allgeom - Boolean to write all available geometries from the logfile.
        """

        # Call the __init__ method of the superclass
        super(XYZ, self).__init__(ccdata, *args, **kwargs)

        self.do_firstgeom = firstgeom
        self.do_lastgeom = lastgeom
        self.do_allgeom = allgeom

        self.indices = kwargs.get('indices')

    def generate_repr(self):
        """Generate the XYZ representation of the logfile data."""

        # Options for output (to a single file):
        # 1. Write all geometries from an optimization, which programs like VMD
        #   can read in like a trajectory.
        # 2. Write the final converged geometry, which for any job other than
        #   a geometry optimization would be the single/only geometry.
        # 3. Write the very first geometry, which for any job other than a
        #   geometry optimization would be the single/only geometry.
        # 4. Write the first and last geometries from a geometry optimization.
        # 5. Write arbitrary structures via zero-based indexing.
        # TODO: Options for output (to multiple files)

        xyzblock = []

        if hasattr(self.ccdata, 'atomcoords'):

            lencoords = len(self.ccdata.atomcoords)

            indices = set()

            # Collect the indices.
            if lencoords == 1 or self.do_firstgeom:
                indices.add(0)
            if self.do_lastgeom:
                indices.add(lencoords - 1)
            if self.do_allgeom:
                for index in range(lencoords):
                    indices.add(index)

            if self.indices:
                if isinstance(self.indices, Iterable):
                    for i in self.indices:
                        if i < 0:
                            i += lencoords
                        indices.add(i)
                else:
                    assert isinstance(self.indices, int)
                    if self.indices < 0:
                        self.indices += lencoords
                    indices.add(self.indices)

            # Generate the XYZ string for each index.
            indices = sorted(indices)
            for i in indices:
                xyzblock.append(self._xyz_from_ccdata(i))

        return '\n'.join(xyzblock)

    def _xyz_from_ccdata(self, index):
        """Create an XYZ file of the geometry at the given index."""

        natom = str(self.ccdata.natom)
        element_list = [self.pt.element[Z] for Z in self.ccdata.atomnos]
        atomcoords = self.ccdata.atomcoords[index]

        # Create a comment derived from the filename and the index.
        if index == -1:
            geometry_num = len(self.ccdata.atomcoords)
        else:
            geometry_num = index + 1
        if self.jobfilename is not None:
            comment = "{}: Geometry {}".format(self.jobfilename, geometry_num)
        else:
            comment = "Geometry {}".format(geometry_num)

        atom_template = '{:3s} {:15.10f} {:15.10f} {:15.10f}'
        block = []
        block.append(natom)
        block.append(comment)
        for element, (x, y, z) in zip(element_list, atomcoords):
            block.append(atom_template.format(element, x, y, z))
        return '\n'.join(block)


if __name__ == "__main__":
    pass
