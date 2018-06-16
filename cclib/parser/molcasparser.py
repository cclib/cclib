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
            self.set_attribute('atomcharges', {'mulliken': atomcharges})

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
        #  Threshold for linear dependence            0.10E-08
        #  Threshold at which DIIS is turned on       0.15E+00
        #  Threshold at which QNR/C2DIIS is turned on 0.75E-01
        #  Threshold for Norm(delta) (QNR/C2DIIS)     0.20E-04
        if line[:34] == '++    Optimization specifications:':

            scftargets = []
            while not line[6:37] == 'Threshold for SCF energy change':
                line = next(inputfile)

            if line[6:37] == 'Threshold for SCF energy change':
                target = float(line.split()[-1])
                scftargets.append(target)

            line = next(inputfile)

            if line[6:34] == 'Threshold for density matrix':
                target = float(line.split()[-1])
                scftargets.append(target)

            line = next(inputfile)

            if line[6:31] == 'Threshold for Fock matrix':
                target = float(line.split()[-1])
                scftargets.append(target)

            self.skip_lines(inputfile,['Threshold for linear dependence', 'Threshold at which DIIS is turned on', \
                                         'Threshold at which QNR/C2DIIS is turned on'])
            line = next(inputfile)

            if line[6:31] == 'Threshold for Norm(delta)':
                target = float(line.split()[-1])
                scftargets.append(target)

            if not hasattr(self, 'scftargets'):
                self.scftargets = []
            self.scftargets.append(scftargets)

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

            scfvalues = []
            line = next(inputfile)
            while line.split()[0] != 'Convergence':
                if line.split()[0].isdigit():
                    energy = float(line.split()[4].replace('*',''))
                    density = float(line.split()[5].replace('*',''))
                    fock = float(line.split()[6].replace('*',''))
                    dnorm = float(line.split()[7].replace('*',''))
                    scfvalues.append([energy, density, fock, dnorm])

                try:
                    line = next(inputfile)
                    if line.split() == []:
                        line = next(inputfile)
                except StopIteration:
                    self.logger.warning('File terminated before end of last SCF!')
                    break

            if not hasattr(self, "scfvalues"):
                self.scfvalues = []
            self.scfvalues.append(scfvalues)

        #  Harmonic frequencies in cm-1
        #
        #  IR Intensities in km/mol
        #
        #                         1         2         3         4         5         6
        #
        #      Frequency:       i60.14    i57.39    128.18    210.06    298.24    309.65
        #
        #      Intensity:    3.177E-03 2.129E-06 4.767E-01 2.056E-01 6.983E-07 1.753E-07
        #      Red. mass:      2.42030   2.34024   2.68044   3.66414   2.61721   3.34904
        #
        #      C1         x   -0.00000   0.00000   0.00000  -0.05921   0.00000  -0.06807
        #      C1         y    0.00001  -0.00001  -0.00001   0.00889   0.00001  -0.02479
        #      C1         z   -0.03190   0.04096  -0.03872   0.00001  -0.12398  -0.00002
        #      C2         x   -0.00000   0.00001   0.00000  -0.06504   0.00000  -0.03487
        #      C2         y    0.00000  -0.00000  -0.00000   0.01045   0.00001  -0.05659
        #      C2         z   -0.03703  -0.03449  -0.07269   0.00000  -0.07416  -0.00001
        #      C3         x   -0.00000   0.00001   0.00000  -0.06409  -0.00001   0.05110
        #      C3         y   -0.00000   0.00001   0.00000   0.00152   0.00000  -0.03263
        #      C3         z   -0.03808  -0.08037  -0.07267  -0.00001   0.07305   0.00000
        # ...
        #      H20        y    0.00245  -0.00394   0.03215   0.03444  -0.10424  -0.10517
        #      H20        z    0.00002  -0.00001   0.00000  -0.00000  -0.00000   0.00000
        #
        #
        #
        # ++ Thermochemistry
        if line[1:29] == 'Harmonic frequencies in cm-1':

            self.skip_line(inputfile,'blank')
            line = next(inputfile)

            while 'Thermochemistry' not in line:

                if 'Frequency:' in line:
                    if not hasattr(self, 'vibfreqs'):
                        self.vibfreqs = []
                    vibfreqs = [float(i.replace('i','-')) for i in line.split()[1:]]
                    self.vibfreqs.extend(vibfreqs)

                if 'Intensity:' in line:
                    if not hasattr(self, 'vibirs'):
                        self.vibirs = []
                    vibirs = map(float, line.split()[1:])
                    self.vibirs.extend(vibirs)

                if 'Red.' in line: 
                    self.skip_line(inputfile,'blank')
                    line = next(inputfile)
                    if not hasattr(self, 'vibdisps'):
                        self.vibdisps = []
                    disps = []
                    for n in range(3*self.natom):
                        numbers = [float(s) for s in line[17:].split()]
                        # The atomindex should start at 0 instead of 1.
                        atomindex = int(re.search(r'\d+$', line.split()[0]).group()) - 1
                        numbermodes = len(numbers)
                        if len(disps) == 0:
                            # Appends empty array of the following 
                            # dimensions (numbermodes, natom, 0) to disps.
                            for mode in range(numbermodes):
                                disps.append([[] for x in range(0, self.natom)])
                        for mode in range(numbermodes):
                            disps[mode][atomindex].append(numbers[mode])
                        line = next(inputfile)
                    self.vibdisps.extend(disps)

                line = next(inputfile)

        ## Parsing thermochemistry attributes here
        #  ++ Thermochemistry
        #
        #   *********************
        #   *                   *
        #   *  THERMOCHEMISTRY  *
        #   *                   *
        #   *********************
        #
        #   Mass-centered Coordinates (Angstrom):
        #   ***********************************************************
        # ...
        #   *****************************************************
        #   Temperature =     0.00 Kelvin, Pressure =   1.00 atm
        #   -----------------------------------------------------
        #   Molecular Partition Function and Molar Entropy:
        #                          q/V (M**-3)    S(kcal/mol*K)
        #   Electronic            0.100000D+01        0.000
        #   Translational         0.100000D+01        0.000
        #   Rotational            0.100000D+01        2.981
        #   Vibrational           0.100000D+01        0.000
        #   TOTAL                 0.100000D+01        2.981
        #
        #   Thermal contributions to INTERNAL ENERGY:
        #   Electronic           0.000 kcal/mol      0.000000 au.
        #   Translational        0.000 kcal/mol      0.000000 au.
        #   Rotational           0.000 kcal/mol      0.000000 au.
        #   Vibrational        111.885 kcal/mol      0.178300 au.
        #   TOTAL              111.885 kcal/mol      0.178300 au.
        #
        #   Thermal contributions to
        #   ENTHALPY           111.885 kcal/mol      0.178300 au.
        #   GIBBS FREE ENERGY  111.885 kcal/mol      0.178300 au.
        #
        #   Sum of energy and thermal contributions
        #   INTERNAL ENERGY                       -382.121931 au.
        #   ENTHALPY                              -382.121931 au.
        #   GIBBS FREE ENERGY                     -382.121931 au.
        #   -----------------------------------------------------
        # ...
        #   ENTHALPY                              -382.102619 au.
        #   GIBBS FREE ENERGY                     -382.179819 au.
        #   -----------------------------------------------------
        #  --
        #
        #  ++    Isotopic shifts:
        if line[4:19] == 'THERMOCHEMISTRY':

            temperature_values = []
            pressure_values = []
            entropy_values = []
            internal_energy_values = []
            enthalpy_values = []
            free_energy_values = []

            while 'Isotopic' not in line:

                if line[1:12] == 'Temperature':
                    temperature_values.append(float(line.split()[2]))
                    pressure_values.append(float(line.split()[6]))

                if line[1:48] == 'Molecular Partition Function and Molar Entropy:':
                    while 'TOTAL' not in line:
                        line = next(inputfile)
                    entropy_values.append(utils.convertor(float(line.split()[2]), 'kcal', 'hartree'))

                if line[1:40] == 'Sum of energy and thermal contributions':
                    internal_energy_values.append(float(next(inputfile).split()[2]))
                    enthalpy_values.append(float(next(inputfile).split()[1]))
                    free_energy_values.append(float(next(inputfile).split()[3]))

                line=next(inputfile)
            # When calculations for more than one temperature value are
            # performed, the values corresponding to room temperature (298.15 K)
            # are returned and if no calculations are performed for 298.15 K, then
            # the values corresponding last temperature value are returned.
            index  = -1
            if 298.15 in temperature_values:
            	index = temperature_values.index(298.15)

            self.set_attribute('temperature', temperature_values[index])
            if len(temperature_values) > 1:
                self.logger.warning('More than 1 values of temperature found')

            self.set_attribute('pressure', pressure_values[index])
            if len(pressure_values) > 1:
                self.logger.warning('More than 1 values of pressure found')

            self.set_attribute('entropy', entropy_values[index])
            if len(entropy_values) > 1:
                self.logger.warning('More than 1 values of entropy found')

            self.set_attribute('enthalpy', enthalpy_values[index])
            if len(enthalpy_values) > 1:
                self.logger.warning('More than 1 values of enthalpy found')

            self.set_attribute('freeenergy', free_energy_values[index])
            if len(free_energy_values) > 1:
                self.logger.warning('More than 1 values of freeenergy found')

        ## Parsing Geometrical Optimization attributes in this section.
        #  **********************************************************************************************************************
        #  *                                    Energy Statistics for Geometry Optimization                                     *
        #  **********************************************************************************************************************
        #                          Energy     Grad      Grad              Step                 Estimated   Geom       Hessian
        #  Iter      Energy       Change     Norm      Max    Element    Max     Element     Final Energy Update Update   Index
        #    1   -382.30023222  0.00000000 0.107221  0.039531 nrc047   0.085726  nrc047     -382.30533799 RS-RFO  None      0
        #    2   -382.30702964 -0.00679742 0.043573  0.014908 nrc001   0.068195  nrc001     -382.30871333 RS-RFO  BFGS      0
        #    3   -382.30805348 -0.00102384 0.014883  0.005458 nrc010  -0.020973  nrc001     -382.30822089 RS-RFO  BFGS      0
        # ...
        # ...
        #   18   -382.30823419 -0.00000136 0.001032  0.000100 nrc053   0.012319  nrc053     -382.30823452 RS-RFO  BFGS      0
        #   19   -382.30823198  0.00000221 0.001051 -0.000092 nrc054   0.066565  nrc053     -382.30823822 RS-RFO  BFGS      0
        #   20   -382.30820252  0.00002946 0.001132 -0.000167 nrc021  -0.064003  nrc053     -382.30823244 RS-RFO  BFGS      0
        #
        #         +----------------------------------+----------------------------------+
        #         +    Cartesian Displacements       +    Gradient in internals         +
        #         +  Value      Threshold Converged? +  Value      Threshold Converged? +
        #   +-----+----------------------------------+----------------------------------+
        #   + RMS + 5.7330E-02  1.2000E-03     No    + 1.6508E-04  3.0000E-04     Yes   +
        #   +-----+----------------------------------+----------------------------------+
        #   + Max + 1.2039E-01  1.8000E-03     No    + 1.6711E-04  4.5000E-04     Yes   +
        #   +-----+----------------------------------+----------------------------------+
        if 'Convergence criterion on energy change' in line:
            self.energy_threshold =  float(line.split()[6])

        if 'Energy Statistics for Geometry Optimization' in line:
            if not hasattr(self, 'geovalues'):
                self.geovalues = []

            self.skip_lines(inputfile, ['stars', 'header'])
            line = next(inputfile)
            assert 'Iter      Energy       Change     Norm' in line
            # A variable keeping track of ongoing iteration.
            iter_number = len(self.geovalues) + 1
            # Iterate till blank line.
            while line.split() != []:
                for i in range(iter_number):
                    line = next(inputfile)
                self.geovalues.append([float(line.split()[2])])
                line = next(inputfile)
            # Along with energy change, RMS and Max values of change in 
            # Cartesian Diaplacement and Gradients are used as optimization
            # criteria. 
            self.skip_lines(inputfile, ['border', 'header', 'header', 'border'])
            line = next(inputfile)
            assert '+ RMS +' in line
            line_rms = line.split()
            line = next(inputfile)
            line_max = next(inputfile).split()
            if not hasattr(self, 'geotargets'):
                # The attribute geotargets is an array consisting of the following
                # values: [Energy threshold, RMS Displacements threshold, RMS Gradient threshold, \
                #          Max Displacements threshold, Max Gradient threshold].
                self.geotargets = [self.energy_threshold, float(line_rms[4]), float(line_rms[8]), float(line_max[4]), float(line_max[8])]

            self.geovalues[iter_number - 1].extend([float(line_rms[3]), float(line_rms[7]), float(line_max[3]), float(line_max[7])])

        #   *********************************************************
        #   * Nuclear coordinates for the next iteration / Angstrom *
        #   *********************************************************
        #    ATOM              X               Y               Z
        #    C1               0.235560       -1.415847        0.012012
        #    C2               1.313797       -0.488199        0.015149
        #    C3               1.087050        0.895510        0.014200
        # ...
        # ...
        #    H19             -0.021327       -4.934915       -0.029355
        #    H20             -1.432030       -3.721047       -0.039835
        #
        #  --
        if 'Nuclear coordinates for the next iteration / Angstrom' in line:
            self.skip_lines(inputfile, ['s','header'])
            line = next(inputfile)

            atomcoords = []
            while line.split() != []:
                atomcoords.append([float(c) for c in line.split()[1:]])
                line = next(inputfile)

            self.atomcoords.append(atomcoords)

        #  **********************************************************************************************************************
        #  *                                    Energy Statistics for Geometry Optimization                                     *
        #  **********************************************************************************************************************
        #                         Energy     Grad      Grad              Step                 Estimated   Geom       Hessian
        #  Iter      Energy       Change     Norm      Max    Element    Max     Element     Final Energy Update Update   Index
        #    1   -382.30023222  0.00000000 0.107221  0.039531 nrc047   0.085726  nrc047     -382.30533799 RS-RFO  None      0
        # ...
        # ...
        #   23   -382.30823115 -0.00000089 0.001030  0.000088 nrc053   0.000955  nrc053     -382.30823118 RS-RFO  BFGS      0
        #
        #         +----------------------------------+----------------------------------+
        #         +    Cartesian Displacements       +    Gradient in internals         +
        #         +  Value      Threshold Converged? +  Value      Threshold Converged? +
        #   +-----+----------------------------------+----------------------------------+
        #   + RMS + 7.2395E-04  1.2000E-03     Yes   + 2.7516E-04  3.0000E-04     Yes   +
        #   +-----+----------------------------------+----------------------------------+
        #   + Max + 1.6918E-03  1.8000E-03     Yes   + 8.7768E-05  4.5000E-04     Yes   +
        #   +-----+----------------------------------+----------------------------------+
        #
        #   Geometry is converged in  23 iterations to a Minimum Structure
        if 'Geometry is converged' in line:
            if not hasattr(self, 'optdone'):
                self.optdone = []
            self.optdone.append(len(self.atomcoords))

        #   *********************************************************
        #   * Nuclear coordinates of the final structure / Angstrom *
        #   *********************************************************
        #    ATOM              X               Y               Z
        #    C1               0.235547       -1.415838        0.012193
        #    C2               1.313784       -0.488201        0.015297
        #    C3               1.087036        0.895508        0.014333
        # ...
        # ...
        #    H19             -0.021315       -4.934913       -0.029666
        #    H20             -1.431994       -3.721026       -0.041078
        if 'Nuclear coordinates of the final structure / Angstrom' in line:
            self.skip_lines(inputfile, ['s','header'])
            line = next(inputfile)

            atomcoords = []

            while line.split() != []:
                atomcoords.append([float(c) for c in line.split()[1:]])
                line = next(inputfile)

            self.atomcoords.append(atomcoords)



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