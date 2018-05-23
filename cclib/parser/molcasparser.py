# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
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
        #  ++    Orbital specifications:
        #  -----------------------

        #  Symmetry species               1

        #  Frozen orbitals                0
        #  Occupied orbitals              3
        #  Secondary orbitals            77
        #  Deleted orbitals               0
        #  Total number of orbitals      80
        #  Number of basis functions     80
        #  --
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

        #  ++    Molecular charges:
        #  ------------------

        #  Mulliken charges per centre and basis function type
        #  ---------------------------------------------------

        #         C1    
        #  1s     2.0005
        #  2s     2.0207
        #  2px    0.0253
        #  2pz    0.1147
        #  2py    1.8198
        #  *s    -0.0215
        #  *px    0.0005
        #  *pz    0.0023
        #  *py    0.0368
        #  *d2+   0.0002
        #  *d1+   0.0000
        #  *d0    0.0000
        #  *d1-   0.0000
        #  *d2-   0.0000
        #  *f3+   0.0000
        #  *f2+   0.0001
        #  *f1+   0.0000
        #  *f0    0.0001
        #  *f1-   0.0001
        #  *f2-   0.0000
        #  *f3-   0.0003
        #  *g4+   0.0000
        #  *g3+   0.0000
        #  *g2+   0.0000
        #  *g1+   0.0000
        #  *g0    0.0000
        #  *g1-   0.0000
        #  *g2-   0.0000
        #  *g3-   0.0000
        #  *g4-   0.0000
        #  Total  6.0000

        #  N-E    0.0000

        #  Total electronic charge=    6.000000

        #  Total            charge=    0.000000
        #--
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
        # *****************************************************************************************************************************
        # *                                                                                                                           *
        # *                                             SCF/KS-DFT Program, Final results                                             *
        # *                                                                                                                           *
        # *                                                                                                                           *
        # *                                                                                                                           *
        # *                                                       Final Results                                                       *
        # *                                                                                                                           *
        # *****************************************************************************************************************************

        # ::    Total SCF energy                                -37.6045426484
        if line[:22] == '::    Total SCF energy' or line[:25] == '::    Total KS-DFT energy':
            if not hasattr(self, 'scfenergies'):
                self.scfenergies = []
            scfenergy = float(line.split()[-1])
            self.scfenergies.append(utils.convertor(scfenergy, 'hartree', 'eV'))

        ## Parsing the scftargets in this section
        #  ++    Optimization specifications:
        #  ----------------------------

        #  SCF Algorithm: Conventional
        #  Minimized density differences are used
        #  Number of density matrices in core                9
        #  Maximum number of NDDO SCF iterations           400
        #  Maximum number of HF  SCF iterations            400
        #  Threshold for SCF energy change            0.10E-08
        #  Threshold for density matrix               0.10E-03
        #  Threshold for Fock matrix                  0.15E-03
        if line[:34] == '++    Optimization specifications:':

            while not line[6:37] == 'Threshold for SCF energy change':
                line = next(inputfile)

            if line[6:37] == 'Threshold for SCF energy change':
                if not hasattr(self, 'scftargets'):
                    self.scftargets = [[]]
                target = float(line.split()[-1])
                self.scftargets[0].append(target)

            line = next(inputfile)

            if line[6:34] == 'Threshold for density matrix':
                if not hasattr(self, 'scftargets'):
                    self.scftargets = [[]]
                target = float(line.split()[-1])
                self.scftargets[0].append(target)

            line = next(inputfile)

            if line[6:31] == 'Threshold for Fock matrix':
                if not hasattr(self, 'scftargets'):
                    self.scftargets = [[]]
                target = float(line.split()[-1])
                self.scftargets[0].append(target)


        #  ++ Convergence information
        #                                     SCF        iterations: Energy and convergence statistics
  
        #  Iter     Tot. SCF       One-electron     Two-electron   Energy   Max Dij or  Max Fij    DNorm      TNorm     AccCon    Time
        #             Energy          Energy          Energy       Change   Delta Norm                                          in Sec.
        #     1    -36.83817703    -50.43096166     13.59278464  0.00E+00   0.16E+00*  0.27E+01*   0.30E+01   0.33E+02   NoneDa    0.
        #     2    -36.03405202    -45.74525152      9.71119950  0.80E+00*  0.14E+00*  0.93E-02*   0.26E+01   0.43E+01   Damp      0.
        #     3    -37.08936118    -48.41536598     11.32600480 -0.11E+01*  0.12E+00*  0.91E-01*   0.97E+00   0.16E+01   Damp      0.
        #     4    -37.31610460    -50.54103969     13.22493509 -0.23E+00*  0.11E+00*  0.96E-01*   0.72E+00   0.27E+01   Damp      0.
        #     5    -37.33596239    -49.47021484     12.13425245 -0.20E-01*  0.59E-01*  0.59E-01*   0.37E+00   0.16E+01   Damp      0.
        if line[46:91] == 'iterations: Energy and convergence statistics':

            self.skip_line(inputfile,'blank')

            while line.split() != ['Energy', 'Energy', 'Energy', 'Change', 'Delta', 'Norm', 'in', 'Sec.']: 
                line = next(inputfile)

            if not hasattr(self, "scfvalues"):
                self.scfvalues = []

            line = next(inputfile)
            while line.split()[0] != 'Convergence':
                if line.split()[0].isdigit():
                    energy = float(line.split()[4].replace('*',''))
                    density = float(line.split()[5].replace('*',''))
                    fock = float(line.split()[6].replace('*',''))
                    self.scfvalues.append([energy, density, fock])

                try:
                    line = next(inputfile)
                    if line.split() == []:
                        line = next(inputfile)
                except StopIteration:
                    self.logger.warning('File terminated before end of last SCF!')
                    break



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