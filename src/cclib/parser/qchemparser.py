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

        if '$molecule' in line:
            line = next(inputfile)
            charge, mult = map(int, line.split())
            self.set_attribute('charge', charge)
            self.set_attribute('mult', mult)

        # Extract the atomic numbers and coordinates of the atoms.

        if 'Standard Nuclear Orientation (Angstroms)' in line:
            if not hasattr(self, 'atomcoords'):
                self.atomcoords = []
            self.skip_lines(inputfile, ['cols', 'dashes'])
            atomelements = []
            atomcoords = []
            line = next(inputfile)
            while list(set(line.strip())) != ['-']:
                entry = line.split()
                atomelements.append(entry[1])
                atomcoords.append(list(map(float, entry[2:])))
                line = next(inputfile)

            self.atomcoords.append(atomcoords)

            atomnos = list(utils.PeriodicTable().number[i]
                           for i in atomelements)

            if not hasattr(self, 'natom'):
                self.natom = len(atomnos)
            if not hasattr(self, 'atomnos'):
                self.atomnos = atomnos

        # Number of basis functions.
        # Because Q-Chem's integral recursion scheme is defined using
        # Cartesian basis functions, there is often a distinction between the
        # two in the output. We only parse for *pure* functions.
        # Examples:
        #  Only one type:
        #   There are 30 shells and 60 basis functions
        #  Both Cartesian and pure:
        #   ...

        if 'basis functions' in line:
            self.nbasis = int(line.split()[-3])

        # Section with SCF iterations.

        if 'SCF converges when DIIS error is below' in line:
            if not hasattr(self, 'scftargets'):
                self.scftargets = []
            diis_target = float(line.split()[-1])
            self.scftargets.append([diis_target])

        if 'Cycle       Energy         DIIS Error' in line:
            self.skip_lines(inputfile, ['d'])
            line = next(inputfile)
            values = []
            while list(set(line.strip())) != ['-']:
                if not hasattr(self, 'scfvalues'):
                    self.scfvalues = []
                # Q-Chem only outputs the DIIS error, even if one chooses
                # (G)DM, RCA, ...
                diis_error = float(line.split()[2])
                values.append([diis_error])
                line = next(inputfile)
            # Only so we remain a list of arrays of rank 2.
            self.scfvalues.append(numpy.array(values))

        if 'Total energy in the final basis set' in line:
            if not hasattr(self, 'scfenergies'):
                self.scfenergies = []
            scfenergy = float(line.split()[-1])
            self.scfenergies.append(utils.convertor(scfenergy, 'hartree', 'eV'))

        # Geometry optimization.

        if 'Maximum     Tolerance    Cnvgd?' in line:
            line_g = list(map(float, next(inputfile).split()[1:3]))
            line_d = list(map(float, next(inputfile).split()[1:3]))
            line_e = next(inputfile).split()[2:4]

            if not hasattr(self, 'geotargets'):
                self.geotargets = [line_g[1], line_d[1], self.float(line_e[1])]
            if not hasattr(self, 'geovalues'):
                self.geovalues = []
            try:
                ediff = abs(self.float(line_e[0]))
            except ValueError:
                ediff = numpy.nan
            geovalues = [line_g[0], line_d[0], ediff]
            self.geovalues.append(geovalues)

        if '**  OPTIMIZATION CONVERGED  **' in line:
            if not hasattr(self, 'optdone'):
                self.optdone = []
            self.optdone.append(len(self.atomcoords))

        if '**  MAXIMUM OPTIMIZATION CYCLES REACHED  **' in line:
            if not hasattr(self, 'optdone'):
                self.optdone = []

        # Moller-Plesset corrections.

        # There are multiple modules in Q-Chem for calculating MPn energies:
        # cdman, ccman, and ccman2, all with different output.
        #
        # MP2, RI-MP2, and local MP2 all default to cdman, which has a simple
        # block of output after the regular SCF iterations.
        #
        # MP3 is handled by ccman2.
        #
        # MP4 and variants are handled by ccman.

        if 'MP2         total energy' in line:
            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []
            mp2energy = float(line.split()[4])
            self.mpenergies.append([mp2energy])

        # This is the MP3 case.
        if 'MP2 energy' in line:
            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []
            mp2energy = float(line.split()[3])
            line = next(inputfile)
            line = next(inputfile)
            # Just a safe check.
            if 'MP3 energy' in line:
                mp3energy = float(line.split()[3])
            self.mpenergies.append([mp2energy, mp3energy])

        # Molecular orbital energies and symmetries.

        if 'Orbital Energies (a.u.) and Symmetries' in line:

            #  --------------------------------------------------------------
            #              Orbital Energies (a.u.) and Symmetries
            #  --------------------------------------------------------------
            #
            #  Alpha MOs, Restricted
            #  -- Occupied --
            # -10.018 -10.018 -10.008 -10.008 -10.007 -10.007 -10.006 -10.005
            #   1 Bu    1 Ag    2 Bu    2 Ag    3 Bu    3 Ag    4 Bu    4 Ag
            #  -9.992  -9.992  -0.818  -0.755  -0.721  -0.704  -0.670  -0.585
            #   5 Ag    5 Bu    6 Ag    6 Bu    7 Ag    7 Bu    8 Bu    8 Ag
            #  -0.561  -0.532  -0.512  -0.462  -0.439  -0.410  -0.400  -0.397
            #   9 Ag    9 Bu   10 Ag   11 Ag   10 Bu   11 Bu   12 Bu   12 Ag
            #  -0.376  -0.358  -0.349  -0.330  -0.305  -0.295  -0.281  -0.263
            #  13 Bu   14 Bu   13 Ag    1 Au   15 Bu   14 Ag   15 Ag    1 Bg
            #  -0.216  -0.198  -0.160
            #   2 Au    2 Bg    3 Bg
            #  -- Virtual --
            #   0.050   0.091   0.116   0.181   0.280   0.319   0.330   0.365
            #   3 Au    4 Au    4 Bg    5 Au    5 Bg   16 Ag   16 Bu   17 Bu
            #   0.370   0.413   0.416   0.422   0.446   0.469   0.496   0.539
            #  17 Ag   18 Bu   18 Ag   19 Bu   19 Ag   20 Bu   20 Ag   21 Ag
            #   0.571   0.587   0.610   0.627   0.646   0.693   0.743   0.806
            #  21 Bu   22 Ag   22 Bu   23 Bu   23 Ag   24 Ag   24 Bu   25 Ag
            #   0.816
            #  25 Bu
            #
            #  Beta MOs, Restricted
            #  -- Occupied --
            # -10.018 -10.018 -10.008 -10.008 -10.007 -10.007 -10.006 -10.005
            #   1 Bu    1 Ag    2 Bu    2 Ag    3 Bu    3 Ag    4 Bu    4 Ag
            #  -9.992  -9.992  -0.818  -0.755  -0.721  -0.704  -0.670  -0.585
            #   5 Ag    5 Bu    6 Ag    6 Bu    7 Ag    7 Bu    8 Bu    8 Ag
            #  -0.561  -0.532  -0.512  -0.462  -0.439  -0.410  -0.400  -0.397
            #   9 Ag    9 Bu   10 Ag   11 Ag   10 Bu   11 Bu   12 Bu   12 Ag
            #  -0.376  -0.358  -0.349  -0.330  -0.305  -0.295  -0.281  -0.263
            #  13 Bu   14 Bu   13 Ag    1 Au   15 Bu   14 Ag   15 Ag    1 Bg
            #  -0.216  -0.198  -0.160
            #   2 Au    2 Bg    3 Bg
            #  -- Virtual --
            #   0.050   0.091   0.116   0.181   0.280   0.319   0.330   0.365
            #   3 Au    4 Au    4 Bg    5 Au    5 Bg   16 Ag   16 Bu   17 Bu
            #   0.370   0.413   0.416   0.422   0.446   0.469   0.496   0.539
            #  17 Ag   18 Bu   18 Ag   19 Bu   19 Ag   20 Bu   20 Ag   21 Ag
            #   0.571   0.587   0.610   0.627   0.646   0.693   0.743   0.806
            #  21 Bu   22 Ag   22 Bu   23 Bu   23 Ag   24 Ag   24 Bu   25 Ag
            #   0.816
            #  25 Bu
            #  --------------------------------------------------------------

            self.skip_line(inputfile, 'dashes')
            line = next(inputfile)
            # Sometimes Q-Chem gets a little confused...
            while 'Warning : Irrep of orbital' in line:
                line = next(inputfile)

            line = next(inputfile)
            unres = False
            energies_alpha = []
            symbols_alpha = []
            if line.split()[2] == 'Unrestricted':
                unres = True
                energies_beta = []
                symbols_beta = []
            line = next(inputfile)

            while len(energies_alpha) < self.nbasis:
                if 'Occupied' in line or 'Virtual' in line:
                    # A nice trick to find where the HOMO is.
                    if 'Virtual' in line:
                        if not hasattr(self, 'homos'):
                            self.homos = [len(energies_alpha)-1]
                    line = next(inputfile)
                # Parse the energies and symmetries in pairs of lines.
                # energies = [utils.convertor(energy, 'hartree', 'eV')
                #             for energy in map(float, line.split())]
                # This convoluted bit handles '*******' when present.
                energies = []
                energy_line = line.split()
                for e in energy_line:
                    try:
                        energy = utils.convertor(self.float(e), 'hartree', 'eV')
                    except ValueError:
                        energy = numpy.nan
                    energies.append(energy)
                energies_alpha.extend(energies)
                line = next(inputfile)
                symbols = line.split()[1::2]
                symbols_alpha.extend(symbols)
                line = next(inputfile)

            # Only look at the second block if doing an unrestricted calculation.
            # This might be a problem for ROHF/ROKS.
            if unres:
                self.skip_line(inputfile, 'header')
                line = next(inputfile)
                while len(energies_beta) < self.nbasis:
                    if 'Occupied' in line or 'Virtual' in line:
                        # This will definitely exist, thanks to the above block.
                        if 'Virtual' in line:
                            if len(self.homos) == 1:
                                self.homos.append(len(energies_beta)-1)
                        line = next(inputfile)
                    energies = []
                    energy_line = line.split()
                    for e in energy_line:
                        try:
                            energy = utils.convertor(self.float(e), 'hartree', 'eV')
                        except ValueError:
                            energy = numpy.nan
                        energies.append(energy)
                    energies_beta.extend(energies)
                    line = next(inputfile)
                    symbols = line.split()[1::2]
                    symbols_beta.extend(symbols)
                    line = next(inputfile)

            if not hasattr(self, 'moenergies'):
                self.moenergies = []
            if not hasattr(self, 'mosyms'):
                self.mosyms = []
            self.moenergies.append(numpy.array(energies_alpha))
            self.mosyms.append(symbols_alpha)
            if unres:
                self.moenergies.append(numpy.array(energies_beta))
                self.mosyms.append(symbols_beta)

        # Molecular orbital energies, no symmetries.

        if line.strip() == 'Orbital Energies (a.u.)':

            # In the case of no orbital symmetries, the beta spin block is not
            # present for restricted calculations.

            #  --------------------------------------------------------------
            #                     Orbital Energies (a.u.)
            #  --------------------------------------------------------------
            #
            #  Alpha MOs
            #  -- Occupied --
            # ******* -38.595 -34.580 -34.579 -34.578 -19.372 -19.372 -19.364
            # -19.363 -19.362 -19.362  -4.738  -3.252  -3.250  -3.250  -1.379
            #  -1.371  -1.369  -1.365  -1.364  -1.362  -0.859  -0.855  -0.849
            #  -0.846  -0.840  -0.836  -0.810  -0.759  -0.732  -0.729  -0.704
            #  -0.701  -0.621  -0.610  -0.595  -0.587  -0.584  -0.578  -0.411
            #  -0.403  -0.355  -0.354  -0.352
            #  -- Virtual --
            #  -0.201  -0.117  -0.099  -0.086   0.020   0.031   0.055   0.067
            #   0.075   0.082   0.086   0.092   0.096   0.105   0.114   0.148
            #
            #  Beta MOs
            #  -- Occupied --
            # ******* -38.561 -34.550 -34.549 -34.549 -19.375 -19.375 -19.367
            # -19.367 -19.365 -19.365  -4.605  -3.105  -3.103  -3.102  -1.385
            #  -1.376  -1.376  -1.371  -1.370  -1.368  -0.863  -0.858  -0.853
            #  -0.849  -0.843  -0.839  -0.818  -0.765  -0.738  -0.737  -0.706
            #  -0.702  -0.624  -0.613  -0.600  -0.591  -0.588  -0.585  -0.291
            #  -0.291  -0.288  -0.275
            #  -- Virtual --
            #  -0.139  -0.122  -0.103   0.003   0.014   0.049   0.049   0.059
            #   0.061   0.070   0.076   0.081   0.086   0.090   0.098   0.106
            #   0.138
            #  --------------------------------------------------------------

            self.skip_lines(inputfile, ['dashes', 'blank'])
            line = next(inputfile)
            unres = False
            energies_alpha = []
            line = next(inputfile)

            while len(energies_alpha) < self.nbasis:
                if 'Occupied' in line or 'Virtual' in line:
                    # A nice trick to find where the HOMO is.
                    if 'Virtual' in line:
                        if not hasattr(self, 'homos'):
                            self.homos = [len(energies_alpha)-1]
                    line = next(inputfile)
                # Parse the energies and symmetries in pairs of lines.
                # energies = [utils.convertor(energy, 'hartree', 'eV')
                #             for energy in map(float, line.split())]
                # This convoluted bit handles '*******' when present.
                energies = []
                energy_line = line.split()
                for e in energy_line:
                    try:
                        energy = utils.convertor(self.float(e), 'hartree', 'eV')
                        energy = self.float(e)
                    except ValueError:
                        energy = numpy.nan
                    energies.append(energy)
                energies_alpha.extend(energies)
                line = next(inputfile)

            line = next(inputfile)
            # Only look at the second block if doing an unrestricted calculation.
            # This might be a problem for ROHF/ROKS.
            if line.strip() == 'Beta MOs':
                unres = True
                energies_beta = []

                self.skip_lines(inputfile, ['blank'])
                line = next(inputfile)
                while len(energies_beta) < self.nbasis:
                    if 'Occupied' in line or 'Virtual' in line:
                        # This will definitely exist, thanks to the above block.
                        if 'Virtual' in line:
                            if len(self.homos) == 1:
                                self.homos.append(len(energies_beta)-1)
                        line = next(inputfile)
                    energies = []
                    energy_line = line.split()
                    for e in energy_line:
                        try:
                            energy = utils.convertor(self.float(e), 'hartree', 'eV')
                            energy = self.float(e)
                        except ValueError:
                            energy = numpy.nan
                        energies.append(energy)
                    energies_beta.extend(energies)
                    line = next(inputfile)

            if not hasattr(self, 'moenergies'):
                self.moenergies = []
            self.moenergies.append(numpy.array(energies_alpha))
            if unres:
                self.moenergies.append(numpy.array(energies_beta))

        # Population analysis.

        if 'Ground-State Mulliken Net Atomic Charges' in line:

            self.skip_line(inputfile, 'blank')
            line = next(inputfile)
            has_spins = False
            if 'Spin' in line:
                if not hasattr(self, 'atomspins'):
                    self.atomspins = dict()
                has_spins = True
                spins = []
            self.skip_line(inputfile, 'dashes')
            if not hasattr(self, 'atomcharges'):
                self.atomcharges = dict()
            charges = []
            line = next(inputfile)

            while list(set(line.strip())) != ['-']:
                elements = line.split()
                charge = self.float(elements[2])
                charges.append(charge)
                if has_spins:
                    spin = self.float(elements[3])
                    spins.append(spin)
                line = next(inputfile)

            self.atomcharges['mulliken'] = numpy.array(charges)
            if has_spins:
                self.atomspins['mulliken'] = numpy.array(spins)

        # For IR-related jobs, the Hessian is printed (dim: 3*natom, 3*natom).
        # if 'Hessian of the SCF Energy' in line:
        #     # A maximum of 6 columns/block.
        #     if not hasattr(self, 'hessian'):
        #         self.hessian = []

        # Start of the IR/Raman frequency section.

        if 'VIBRATIONAL ANALYSIS' in line:

            while 'STANDARD THERMODYNAMIC QUANTITIES' not in line:

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

                if 'Frequency:' in line:
                    if not hasattr(self, 'vibfreqs'):
                        self.vibfreqs = []
                    vibfreqs = map(float, line.split()[1:])
                    self.vibfreqs.extend(vibfreqs)

                if 'IR Intens:' in line:
                    if not hasattr(self, 'vibirs'):
                        self.vibirs = []
                    vibirs = map(float, line.split()[2:])
                    self.vibirs.extend(vibirs)

                if 'Raman Intens:' in line:
                    if not hasattr(self, 'vibramans'):
                        self.vibramans = []
                    vibramans = map(float, line.split()[2:])
                    self.vibramans.extend(vibramans)

                # This is the start of the displacement block.
                if line.split()[0:3] == ['X', 'Y', 'Z']:
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

                if 'VIBRATIONAL ANHARMONIC ANALYSIS' in line:

                    while list(set(line.strip())) != ['=']:
                        if 'VPT2' in line:
                            if not hasattr(self, 'vibanharms'):
                                self.vibanharms = []
                            self.vibanharms.append(float(line.split()[-1]))
                        line = next(inputfile)

        # TODO:
        # 'aonames'
        # 'atombasis'
        # 'atommasses'
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
        # 'grads'
        # 'hessian'
        # 'mocoeffs'
        # 'mpenergies'
        # 'nocoeffs'
        # 'scancoords'
        # 'scanenergies'
        # 'scannames'
        # 'scanparm'
        # 'temperature'


if __name__ == '__main__':
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

