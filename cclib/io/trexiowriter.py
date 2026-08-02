# Copyright (c) 2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
"""A writer for trexio format files (hdf5)."""

from cclib.io import filewriter


class Trexio(filewriter.Writer):
    """A writer for Trexio files."""

    def __init__(self, ccdata, *args, **kwargs):
        """Initialize the Trexio writer object.

        Inputs:
          ccdata - An instance of ccData, parse from a logfile.
        """

        self.required_attrs = ("natom", "atomcoords", "atomnos")

        # Call the __init__ method of the superclass
        super().__init__(ccdata, *args, **kwargs)

        self.natom = str(self.ccdata.natom)
        self.element_list = [self.pt.element[Z] for Z in self.ccdata.atomnos]

    def generate_repr(self):
        """Generate the Trexio representation of the logfile data."""

        # Options for output (to a single file):
        # 1. Write the single geometry at a time

        xyzblock = []

        lencoords = len(self.ccdata.atomcoords)

        # Collect the indices.
        if lencoords == 1:
            self.indices.add(0)

        # Generate the XYZ string for each index.
        indices = sorted(self.indices)
        if not indices:
            indices = [-1]
        for i in indices:
            xyzblock.append(self._xyz_from_ccdata(i))

        # Ensure an extra newline at the very end.
        xyzblock.append("")

        return "\n".join(xyzblock)

    def _xyz_from_ccdata(self, index):
        """Write to trexio hdf5 format file of the geometry"""

        atomcoords = self.ccdata.atomcoords[index]
        existing_comment = (
            ""
            if "comments" not in self.ccdata.metadata
            else self.ccdata.metadata["comments"][index]
        )

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

        atom_template = "{:3s} {:15.10f} {:15.10f} {:15.10f}"
        block = []
        block.append(self.natom)
        block.append(comment)
        for element, (x, y, z) in zip(self.element_list, atomcoords):
            block.append(atom_template.format(element, x, y, z))
        return "\n".join(block)
