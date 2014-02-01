"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 796 $"


import re

import numpy

from .import logfileparser
from . import utils


class NWChem(logfileparser.Logfile):
    """An NWChem log file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(NWChem, self).__init__(logname="NWChem", *args, **kwargs)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "NWChem log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'NWChem("%s")' % (self.filename)
    
    def normalisesym(self, label):
        """Use standard symmetry labels instead of NWChem labels.

        To normalise:
        (1) If label is one of [SG, PI, PHI, DLTA], replace by [sigma, pi, phi, delta]
        (2) replace any G or U by their lowercase equivalent

        >>> sym = NWChem("dummyfile").normalisesym
        >>> labels = ['A1', 'AG', 'A1G', "SG", "PI", "PHI", "DLTA", 'DLTU', 'SGG']
        >>> map(sym, labels)
        ['A1', 'Ag', 'A1g', 'sigma', 'pi', 'phi', 'delta', 'delta.u', 'sigma.g']
        """
        # FIXME if necessary
        return label

    def before_parsing(self):

        # Set any global variables for the parser here
        pass

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""
        
        # Number of atoms (THIS IS FROM THE GAUSSIAN PARSER, PLEASE REPLACE)
        if line[1:8] == "NAtoms=":

            self.updateprogress(inputfile, "Attributes", self.fupdate)
                    
            natom = int(line.split()[1])
            if hasattr(self, "natom"):
                assert self.natom == natom
            else:
                # I wonder whether this code will ever be executed.
                self.natom = natom

if __name__ == "__main__":
    import doctest, nwchemparser
    doctest.testmod(nwchemparser, verbose=False)
