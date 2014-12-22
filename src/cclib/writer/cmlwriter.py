# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""A writer for chemical markup language (CML) files."""

from . import filewriter


class CML(filewriter.Writer):
    """A writer for chemical markup language (CML) files."""

    def __init__(self, ccdata, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(CML, self).__init__(ccdata, *args, **kwargs)

    def generate_repr(self):
        """Generate the CML representation of the logfile data."""
        pass


if __name__ == "__main__":
    pass
