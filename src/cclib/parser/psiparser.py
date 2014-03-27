# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2008-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.


import re

import numpy

from . import logfileparser
from . import utils


class Psi(logfileparser.Logfile):
    """A Psi log file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(Psi, self).__init__(logname="Psi", *args, **kwargs)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Psi log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Psi("%s")' % (self.filename)
    
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

        if line[2:16] == "charge       =":
            self.charge = int(line.split()[-1])

        if line[2:16] == "multiplicity =":
            self.mult = int(line.split()[-1])

        if "Number of atoms" in line:
            self.natom = int(line.split()[-1])

        if "Number of atomic orbitals" in line:
            self.nbasis = int(line.split()[-1])

        if "SCF total energy" in line:
            self.scfenergies = [float(line.split()[-1])]

        if line.strip() == "Orbital energies (a.u.):":
            self.moenergies = []
            blank = next(inputfile)
            doubly = next(inputfile)
            line = next(inputfile)
            while line.strip():
                for i in range(len(line.split())//2):
                    self.moenergies.append(line.split()[i*2+1])
                line = next(inputfile)
            blank = next(inputfile)
            blank = next(inputfile)
            unoccupied = next(inputfile)
            while line.strip():
                for i in range(len(line.split())//2):
                    self.moenergies.append(line.split()[i*2+1])
                line = next(inputfile)


if __name__ == "__main__":
    import doctest, psiparser
    doctest.testmod(psiparser, verbose=False)
