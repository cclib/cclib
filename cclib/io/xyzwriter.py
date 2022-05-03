# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for XYZ (Cartesian coordinate) files."""

from cclib.io import filewriter


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
        super().__init__(ccdata, *args, **kwargs)

        self.required_attrs = ('natom', 'atomcoords', 'atomnos')

        self.do_firstgeom = firstgeom
        self.do_lastgeom = lastgeom
        self.do_allgeom = allgeom

        self.natom = str(self.ccdata.natom)
        self.element_list = [self.pt.element[Z] for Z in self.ccdata.atomnos]

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

        lencoords = len(self.ccdata.atomcoords)

        # Collect the indices.
        if lencoords == 1 or self.do_firstgeom:
            self.indices.add(0)
        if self.do_lastgeom:
            self.indices.add(lencoords - 1)
        if self.do_allgeom:
            for i in range(lencoords):
                self.indices.add(i)

        # Generate the XYZ string for each index.
        indices = sorted(self.indices)
        if not indices:
            indices = [-1]
        for i in indices:
            xyzblock.append(self._xyz_from_ccdata(i))

        # Ensure an extra newline at the very end.
        xyzblock.append('')

        return '\n'.join(xyzblock)

    def _xyz_from_ccdata(self, index):
        """Create an XYZ file of the geometry at the given index."""

        atomcoords = self.ccdata.atomcoords[index]
        existing_comment = "" if "comments" not in self.ccdata.metadata \
            else self.ccdata.metadata["comments"][index]

        # Create a comment derived from the filename and the index.
        if index == -1:
            geometry_num = len(self.ccdata.atomcoords)
        else:
            geometry_num = index + 1
        if self.jobfilename is not None:
            comment = f"{self.jobfilename}: Geometry {geometry_num}"
        else:
            comment = f"Geometry {geometry_num}"
        # Wrap the geometry number part of the comment in square brackets,
        # prefixing it with one previously parsed if it existed.
        if existing_comment:
            comment = f"{existing_comment} [{comment}]"
        else:
            comment = f"[{comment}]"

        atom_template = '{:3s} {:15.10f} {:15.10f} {:15.10f}'
        block = []
        block.append(self.natom)
        block.append(comment)
        for element, (x, y, z) in zip(self.element_list, atomcoords):
            block.append(atom_template.format(element, x, y, z))
        return '\n'.join(block)
