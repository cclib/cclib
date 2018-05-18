# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Molcas output files"""

from __future__ import print_function

import re

from cclib.parser import logfileparser
from cclib.parser import utils


class Molcas(logfileparser.Logfile):
    """A Molcas log file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(Molcas, self).__init__(logname="Molcas", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "Molcas log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Molcas("%s")' % (self.filename)

    #These are yet to be implemented.
    def normalisesym(self, label):
        """Does Molcas require symmetry label normalization?"""

    def before_parsing(self):
        # Compile the regex for extracting the element symbol from the
        # atom label in the "Molecular structure info" block.
        self.re_atomelement = re.compile('([a-zA-Z]+)\d+')   

        # Compile the dashes-and-or-spaces-only regex.
        self.re_dashes_and_spaces = re.compile('^[\s-]+$')

    def after_parsing(self):
        pass

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        ## This section is present when executing &GATEWAY.
        # ++    Molecular structure info:
        #       -------------------------

        #                     ************************************************ 
        #                     **** Cartesian Coordinates / Bohr, Angstrom **** 
        #                     ************************************************ 

        #      Center  Label                x              y              z                     x              y              z
        #         1      C1               0.526628      -2.582937       0.000000              0.278679      -1.366832       0.000000
        #         2      C2               2.500165      -0.834760       0.000000              1.323030      -0.441736       0.000000
        if line[25:63] == 'Cartesian Coordinates / Bohr, Angstrom':
            if not hasattr(self, 'atomnos'):
                self.atomnos = []
            if not hasattr(self, 'atomcoords'):
                self.atomcoords = []

            self.skip_lines(inputfile, ['stars', 'blank', 'header'])

            line = next(inputfile)

            atomelements = []
            atomcoords = []

            while not self.re_dashes_and_spaces.search(line):
                sline = line.split()
                atomelements.append(self.re_atomelement.search(sline[1]).groups()[0])
                atomcoords.append(list(map(float, sline[5:])))
                line = next(inputfile)

            self.atomcoords.append(atomcoords)

            if self.atomnos == []:
                self.atomnos = [utils.PeriodicTable().number[atomelement] for atomelement in atomelements]

            if not hasattr(self, 'natom'):
                self.set_attribute('natom', len(self.atomnos))

        ## This section is present when executing &SCF.
        if line[:29] == '++    Orbital specifications:':

            self.skip_lines(inputfile, ['dashes', 'blank'])

            line = next(inputfile)

            while line[:2] != '--':

                if line[6:30] == 'Total number of orbitals':
                    self.set_attribute('nmo', int(line.split()[-1]))
                if line[6:31] == 'Number of basis functions':
                    self.set_attribute('nbasis', int(line.split()[-1]))

                line = next(inputfile)

        #Parsing the molecular charge
        if line[6:23] == 'Molecular charge ':
            self.set_attribute('charge', int(float(line.split()[-1])))


        if line[:24] == '++    Molecular charges:':

            atomcharges = []

            while line[6:29] != 'Total electronic charge':
                line = next(inputfile)
                if line[6:9] == 'N-E':
                    atomcharges.extend(map(float, line.split()[1:]))

            # Molcas only performs Mulliken population analysis.
            self.set_attribute('atomcharges', {'mullikoen': atomcharges})

            # Ensure the charge printed here is identical to the
            # charge printed before entering the SCF.
            self.skip_line(inputfile, 'blank')
            line = next(inputfile)
            assert line[6:30] == 'Total            charge='
            if hasattr(self, 'charge'):
                assert int(float(line.split()[2])) == self.charge

        # This section is present when executing &SCF
        # This section parses the total SCF Energy.
        if line[:22] == '::    Total SCF energy':
            if not hasattr(self, 'scfenergies'):
                self.scfenergies = []
            scfenergy = float(line.split()[-1])
            self.scfenergies.append(utils.convertor(scfenergy, 'hartree', 'eV'))


if __name__ == '__main__':
    import sys
    import doctest, molcasparser

    if len(sys.argv) == 1:
        doctest.testmod(molcasparser, verbose=False)

    if len(sys.argv) == 2:
        parser = molcasparser.Molcas(sys.argv[1])
        data = parser.parse()

    if len(sys.argv) > 2:
        for i in range(len(sys.argv[2:])):
            if hasattr(data, sys.argv[2 + i]):
                print(getattr(data, sys.argv[2 + i]))