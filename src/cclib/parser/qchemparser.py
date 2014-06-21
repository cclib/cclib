# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

from __future__ import print_function

from . import logfileparser
from . import utils


class QChem(logfileparser.Logfile):
    """A Q-Chem 4 log file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(QChem, self).__init__(logname="QChem", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "QChem log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'QChem("%s")' % (self.filename)

    def normalisesym(self, label):
        pass

    def before_parsing(self):
        pass

    def after_parsing(self):
        pass

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # Extract the atomic numbers and coordinates of the atoms.
        if line.find("Standard Nuclear Orientation (Angstroms)") > -1:
            self.skip_lines(inputfile, ['cols', 'd'])
            atomelements = []
            atomcoords = []
            line = next(inputfile)
            while list(set(line.strip())) != ["-"]:
                entry = line.split()
                atomelements.append(entry[1])
                atomcoords.append(list(map(float, entry[2:])))
                line = next(inputfile)

            atomnos = list(utils.PeriodicTable().number[i] for i in atomelements)
            self.set_attribute('natom', len(atomnos))
            self.set_attribute('atomnos', atomnos)
            self.set_attribute('atomcoords', atomcoords)

        # TODO:
        # 'aonames'
        # 'aooverlaps'
        # 'atombasis'
        # 'atomcharges'
        # 'atommasses'
        # 'atomspins'
        # 'ccenergies'
        # 'charge'
        # 'coreelectrons'
        # 'enthalpy'
        # 'entropy'
        # 'etenergies'
        # 'etoscs'
        # 'etrotats'
        # 'etsecs'
        # 'etsyms'
        # 'freeenergy'
        # 'fonames'
        # 'fooverlaps'
        # 'fragnames'
        # 'frags'
        # 'gbasis'
        # 'geotargets'
        # 'geovalues'
        # 'grads'
        # 'hessian'
        # 'homos'
        # 'mocoeffs'
        # 'moenergies'
        # 'mosyms'
        # 'mpenergies'
        # 'mult'
        # 'nbasis'
        # 'nmo'
        # 'nocoeffs'
        # 'optdone'
        # 'scancoords'
        # 'scanenergies'
        # 'scannames'
        # 'scanparm'
        # 'scfenergies'
        # 'scftargets'
        # 'scfvalues'
        # 'temperature'
        # 'vibanharms'
        # 'vibdisps'
        # 'vibfreqs'
        # 'vibirs'
        # 'vibramans'
        # 'vibsyms'


if __name__ == "__main__":
    import sys
    import doctest, qchemparser

    if len(sys.argv) == 1:
        doctest.testmod(qchemparser, verbose=False)

    if len(sys.argv) == 2:
        parser = qchemparser.QChem(sys.argv[1])
        data = parser.parse()

    if len(sys.argv) > 2:
        for i in range(len(sys.argv[2:])):
            if hasattr(data, sys.argv[2 + i]):
                print(getattr(data, sys.argv[2 + i]))

