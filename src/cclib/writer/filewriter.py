# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Generic file writer and related tools"""


class Writer(object):
    """Abstract class for writer objects.

    Subclasses defined by cclib:
        CJSON, CML, XYZ
    """

    def __init__(self, ccdata, *args, **kwargs):
        """Initialize the Writer object.

        This should be called by a subclass in its own __init__ method.

        Inputs:
            ccdata - ...
        """

    def generate_repr(self):
        """Generate the written representation of the logfile data.
        """
        pass

if __name__ == "__main__":
    pass
