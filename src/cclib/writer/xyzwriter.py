# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""A writer for XYZ (Cartesian coordinate) files."""

from . import filewriter
from cclib.parser.utils import PeriodicTable


class XYZ(filewriter.Writer):
    """A writer for XYZ (Cartesian coordinate) files."""

    def __init__(self, ccdata, jobfilename=None, splitfiles=False,
                 firstgeom=False, lastgeom=True, allgeom=False,
                 *args, **kwargs):

        # Call the __init__ method of the superclass
        super(XYZ, self).__init__(ccdata, *args, **kwargs)

        self.ccdata = ccdata
        self.jobfilename = jobfilename
        self.pt = PeriodicTable()

        self.do_firstgeom = firstgeom
        self.do_lastgeom = lastgeom
        self.do_allgeom = allgeom

        self.generate_repr()

    def generate_repr(self):
        """Generate the XYZ representation of the logfile data."""

        # Options for output (to a single file):
        # 1. Write the final converged geometry, which for any job other than
        #   a geometry optimization would be the single/only geometry.
        # 2. Write the very first geometry, which for any job other than a
        #   geometry optimization would be the single/only geometry.
        # 3. Write all geometries from an optimization, which programs like VMD
        #   can read in like a trajectory.
        # Options for ouput (to multiple files):
        # 1. Write all geometries from an optimization, to suitably named files. [TODO]

        xyzblock = []

        lencoords = len(self.ccdata.atomcoords)

        if lencoords == 1:
            xyzblock.append(self._xyz_from_ccdata(self.ccdata, -1))
        elif self.do_allgeom:
            for index in range(lencoords):
                xyzblock.append(self._xyz_from_ccdata(self.ccdata, index))
        elif self.do_firstgeom:
            xyzblock.append(self._xyz_from_ccdata(self.ccdata, 0))
        elif self.do_lastgeom:
            xyzblock.append(self._xyz_from_ccdata(self.ccdata, -1))
        # If none of the options are set, return the empty string.
        else:
            xyzblock.append("")

        return '\n'.join(xyzblock)

    def _xyz_from_ccdata(self, ccdata, index):
        """Create an XYZ file of the geometry at the given index."""
        natom = str(ccdata.natom)
        element_list = [self.pt.element[Z] for Z in ccdata.atomnos]
        atomcoords = ccdata.atomcoords[index]
        # Create a comment derived from the filename and the index.
        comment = "{}: Geometry {}".format(self.jobfilename, index + 1)

        atom_template = '{:3s} {:15.10f} {:15.10f} {:15.10f}'
        block = []
        block.append(natom)
        block.append(comment)
        for element, (x, y, z) in zip(element_list, atomcoords):
            block.append(atom_template.format(element, x, y, z))
        return '\n'.join(block)


if __name__ == "__main__":
    pass
