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

import numpy

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

        # Charge and multiplicity.
        # Only present in the input file.

        if line.find("$molecule") > -1:
            line = next(inputfile)
            charge, mult = map(int, line.split())
            self.set_attribute('charge', charge)
            self.set_attribute('mult', mult)

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

        # Number of basis functions.
        # Because Q-Chem's integral recursion scheme is defined using
        # Cartesian basis functions, there is often a distinction between the
        # two in the output. We only parse for *pure* functions.
        # Examples:
        #  Only one type:
        #   There are 30 shells and 60 basis functions
        #  Both Cartesian and pure:
        #   ...

        if line.find("basis functions") > -1:
            self.nbasis = int(line.split()[-3])

        # Section with SCF iterations.

        if line.find("SCF converges when DIIS error is below") > -1:
            if not hasattr(self, 'scftargets'):
                self.scftargets = []
            diis_target = float(line.split()[-1])
            self.scftargets.append([diis_target])

        if line.find("Cycle       Energy         DIIS Error") > -1:
            self.skip_lines(inputfile, ['d'])
            line = next(inputfile)
            scfvalues = []
            values = []
            while list(set(line.strip())) != ['-']:
                if not hasattr(self, 'scfvalues'):
                    self.scfvalues = []
                # Q-Chem only outputs the DIIS error, even if one chooses
                # (G)DM, RCA, ...
                diis_error = float(line.split()[2])
                values.append(diis_error)
                line = next(inputfile)
            # Only so we remain a list of arrays of rank 2.
            scfvalues.append(values)
            self.scfvalues.append(numpy.array(scfvalues))

        if line.find("Total energy in the final basis set") > -1:
            if not hasattr(self, 'scfenergies'):
                self.scfenergies = []
            scfenergy = float(line.split()[-1])
            self.scfenergies.append(utils.convertor(scfenergy, "hartree", "eV"))

        # For IR-related jobs, the Hessian is printed (dim: 3*natom, 3*natom).
        # if line.find("Hessian of the SCF Energy") > -1:
        #     # A maximum of 6 columns/block.
        #     if not hasattr(self, 'hessian'):
        #         self.hessian = []

        # Start of the IR/Raman frequency section.

        if line.find("VIBRATIONAL ANALYSIS") > -1:

            while "STANDARD THERMODYNAMIC QUANTITIES" not in line:

                ## IR, optional Raman:

                # **********************************************************************
                # **                                                                  **
                # **                       VIBRATIONAL ANALYSIS                       **
                # **                       --------------------                       **
                # **                                                                  **
                # **        VIBRATIONAL FREQUENCIES (CM**-1) AND NORMAL MODES         **
                # **     FORCE CONSTANTS (mDYN/ANGSTROM) AND REDUCED MASSES (AMU)     **
                # **                  INFRARED INTENSITIES (KM/MOL)                   **
                ##** RAMAN SCATTERING ACTIVITIES (A**4/AMU) AND DEPOLARIZATION RATIOS **
                # **                                                                  **
                # **********************************************************************


                # Mode:                 1                      2                      3
                # Frequency:      -106.88                -102.91                 161.77
                # Force Cnst:      0.0185                 0.0178                 0.0380
                # Red. Mass:       2.7502                 2.8542                 2.4660
                # IR Active:          NO                     YES                    YES
                # IR Intens:        0.000                  0.000                  0.419
                # Raman Active:       YES                    NO                     NO
                ##Raman Intens:     2.048                  0.000                  0.000
                ##Depolar:          0.750                  0.000                  0.000
                #               X      Y      Z        X      Y      Z        X      Y      Z
                # C          0.000  0.000 -0.100   -0.000  0.000 -0.070   -0.000 -0.000 -0.027
                # C          0.000  0.000  0.045   -0.000  0.000 -0.074    0.000 -0.000 -0.109
                # C          0.000  0.000  0.148   -0.000 -0.000 -0.074    0.000  0.000 -0.121
                # C          0.000  0.000  0.100   -0.000 -0.000 -0.070    0.000  0.000 -0.027
                # C          0.000  0.000 -0.045    0.000 -0.000 -0.074   -0.000 -0.000 -0.109
                # C          0.000  0.000 -0.148    0.000  0.000 -0.074   -0.000 -0.000 -0.121
                # H         -0.000  0.000  0.086   -0.000  0.000 -0.082    0.000 -0.000 -0.102
                # H          0.000  0.000  0.269   -0.000 -0.000 -0.091    0.000  0.000 -0.118
                # H          0.000  0.000 -0.086    0.000 -0.000 -0.082   -0.000  0.000 -0.102
                # H         -0.000  0.000 -0.269    0.000  0.000 -0.091   -0.000 -0.000 -0.118
                # C          0.000 -0.000  0.141   -0.000 -0.000 -0.062   -0.000  0.000  0.193
                # C         -0.000 -0.000 -0.160    0.000  0.000  0.254   -0.000  0.000  0.043
                # H          0.000 -0.000  0.378   -0.000  0.000 -0.289    0.000  0.000  0.519
                # H         -0.000 -0.000 -0.140    0.000  0.000  0.261   -0.000 -0.000  0.241
                # H         -0.000 -0.000 -0.422    0.000  0.000  0.499   -0.000  0.000 -0.285
                # C          0.000 -0.000 -0.141    0.000  0.000 -0.062   -0.000 -0.000  0.193
                # C         -0.000 -0.000  0.160   -0.000 -0.000  0.254    0.000  0.000  0.043
                # H          0.000 -0.000 -0.378    0.000 -0.000 -0.289   -0.000  0.000  0.519
                # H         -0.000 -0.000  0.140   -0.000 -0.000  0.261    0.000  0.000  0.241
                # H         -0.000 -0.000  0.422   -0.000 -0.000  0.499    0.000  0.000 -0.285
                # TransDip   0.000 -0.000 -0.000    0.000 -0.000 -0.000   -0.000  0.000  0.021

                # Mode:                 4                      5                      6
                # ...

                # There isn't any symmetry information for normal modes present
                # in Q-Chem.
                # if not hasattr(self, 'vibsyms'):
                #     self.vibsyms = []

                if "Frequency:" in line:
                    if not hasattr(self, 'vibfreqs'):
                        self.vibfreqs = []
                    vibfreqs = map(float, line.split()[1:])
                    self.vibfreqs.extend(vibfreqs)

                if "IR Intens:" in line:
                    if not hasattr(self, 'vibirs'):
                        self.vibirs = []
                    vibirs = map(float, line.split()[2:])
                    self.vibirs.extend(vibirs)

                if "Raman Intens:" in line:
                    if not hasattr(self, 'vibramans'):
                        self.vibramans = []
                    vibramans = map(float, line.split()[2:])
                    self.vibramans.extend(vibramans)

                # This is the start of the displacement block.
                if line.split()[0:3] == ["X", "Y", "Z"]:
                    if not hasattr(self, 'vibdisps'):
                        self.vibdisps = []
                    disps = []
                    for k in range(self.natom):
                        line = next(inputfile)
                        numbers = list(map(float, line.split()[1:]))
                        N = len(numbers) // 3
                        if not disps:
                            for n in range(N):
                                disps.append([])
                        for n in range(N):
                            disps[n].append(numbers[3*n:(3*n)+3])
                    self.vibdisps.extend(disps)

                line = next(inputfile)

                # Anharmonic vibrational analysis.
                # Q-Chem includes 3 theories: VPT2, TOSH, and VCI.
                # For now, just take the VPT2 results.

                if line.find("VIBRATIONAL ANHARMONIC ANALYSIS") > -1:

                    while list(set(line.strip())) != ["="]:
                        if "VPT2" in line:
                            if not hasattr(self, 'vibanharms'):
                                self.vibanharms = []
                            self.vibanharms.append(float(line.split()[-1]))
                        line = next(inputfile)

        # TODO:
        # 'aonames'
        # 'atombasis'
        # 'atomcharges'
        # 'atommasses'
        # 'atomspins'
        # 'ccenergies'
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
        # 'nocoeffs'
        # 'optdone'
        # 'scancoords'
        # 'scanenergies'
        # 'scannames'
        # 'scanparm'
        # 'temperature'


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

