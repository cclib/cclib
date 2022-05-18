# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for ORCA output files"""


import re
from itertools import zip_longest

import numpy
from packaging.version import parse as parse_version

from cclib.parser import logfileparser
from cclib.parser import utils


class ORCA(logfileparser.Logfile):
    """An ORCA log file."""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="ORCA", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return f"ORCA log file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'ORCA("{self.filename}")'

    def normalisesym(self, label):
        """ORCA does not require normalizing symmetry labels."""
        return label

    def before_parsing(self):

        self.uses_symmetry = False

        # A geometry optimization is started only when
        # we parse a cycle (so it will be larger than zero().
        self.gopt_cycle = 0

        # Keep track of whether this is a relaxed scan calculation
        self.is_relaxed_scan = False

    def after_parsing(self):
        # ORCA doesn't add the dispersion energy to the "Total energy" (which
        # we parse), only to the "FINAL SINGLE POINT ENERGY" (which we don't
        # parse).
        if hasattr(self, "scfenergies") and hasattr(self, "dispersionenergies"):
            for i, (scfenergy, dispersionenergy) in enumerate(
                zip_longest(self.scfenergies, self.dispersionenergies)
            ):
                # It isn't as problematic if there are more dispersion than
                # SCF energies, since all dispersion energies can still be
                # added to the SCF energies, hence the difference in log level.
                if dispersionenergy is None:
                    self.logger.error(
                        "The number of SCF and dispersion energies are not equal: %d vs. %d, "
                        "can't add dispersion energy to all SCF energies",
                        len(self.scfenergies),
                        len(self.dispersionenergies)
                    )
                    break
                if scfenergy is None:
                    self.logger.warning(
                        "The number of SCF and dispersion energies are not equal: %d vs. %d, "
                        "can't add dispersion energy to all SCF energies",
                        len(self.scfenergies),
                        len(self.dispersionenergies)
                    )
                    break
                self.scfenergies[i] += dispersionenergy

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # Extract the version number.
        if "Program Version" == line.strip()[:15]:
            # Handle development versions.
            self.metadata["legacy_package_version"] = line.split()[2]
            self.metadata["package_version"] = self.metadata["legacy_package_version"].replace(".x", "dev")
            possible_revision_line = next(inputfile)
            if "SVN: $Rev" in possible_revision_line:
                version = re.search(r'\d+', possible_revision_line).group()
                self.metadata["package_version"] += f"+{version}"

        # ================================================================================
        #                                         WARNINGS
        #                        Please study these warnings very carefully!
        # ================================================================================
        #
        # Warning: TCutStore was < 0. Adjusted to Thresh (uncritical)
        #
        # WARNING: your system is open-shell and RHF/RKS was chosen
        #   ===> : WILL SWITCH to UHF/UKS
        #
        #
        # INFO   : the flag for use of LIBINT has been found!
        #
        # ================================================================================
        if "WARNINGS" == line.strip():
            self.skip_lines(inputfile, ['text', '=', 'blank'])
            if 'warnings' not in self.metadata:
                self.metadata['warnings'] = []
            if 'info' not in self.metadata:
                self.metadata['info'] = []

            line = next(inputfile)
            while line[0] != '=':
                if line.lower()[:7] == 'warning':
                    self.metadata['warnings'].append('')
                    while len(line) > 1 and set(line.strip()) != {'='}:
                        self.metadata['warnings'][-1] += line[9:].strip()
                        line = next(inputfile)
                elif line.lower()[:4] == 'info':
                    self.metadata['info'].append('')
                    while len(line) > 1 and set(line.strip()) != {'='}:
                        self.metadata['info'][-1] += line[9:].strip()
                        line = next(inputfile)
                else:
                    line = next(inputfile)

        # ================================================================================
        #                                        INPUT FILE
        # ================================================================================
        # NAME = input.dat
        # |  1> %pal nprocs 4 end
        # |  2> ! B3LYP def2-svp
        # |  3> ! Grid4
        # |  4>
        # |  5> *xyz 0 3
        # |  6>     O   0   0   0
        # |  7>     O   0   0   1.5
        # |  8> *
        # |  9>
        # | 10>                          ****END OF INPUT****
        # ================================================================================
        if "INPUT FILE" == line.strip():
            self.skip_line(inputfile, '=')
            self.metadata['input_file_name'] = next(inputfile).split()[-1]

            # First, collect all the lines...
            lines = []
            for line in inputfile:
                if line[0] != '|':
                    break
                lines.append(line[line.find('> ')+2:])

            self.metadata['input_file_contents'] = ''.join(lines[:-1])
            lines_iter = iter(lines[:-1])

            keywords = []
            coords = []
            # ...then parse them separately.
            for line in lines_iter:
                line = line.strip()
                if not line:
                    continue

                # Keywords block
                if line[0] == '!':
                    keywords += line[1:].split()

                # Impossible to parse without knowing whether a keyword opens a new block
                elif line[0] == '%':
                    pass
                # Geometry block
                elif line[0] == '*':
                    coord_type, charge, multiplicity = line[1:].split()[:3]
                    self.set_attribute('charge', int(charge))
                    self.set_attribute('multiplicity', int(multiplicity))
                    coord_type = coord_type.lower()
                    self.metadata['coord_type'] = coord_type
                    if coord_type == 'xyz':
                        def splitter(line):
                            atom, x, y, z = line.split()[:4]
                            return [atom, float(x), float(y), float(z)]
                    elif coord_type in ['int', 'internal']:
                        def splitter(line):
                            atom, a1, a2, a3, bond, angle, dihedral = line.split()[:7]
                            # This could be some combination of floats and variables
                            # C                  0    0    0       0.0                  0.0                   0.0
                            # C                  3    2    1       {B3}                 {A2}                  {D1}
                            try:
                                return [atom, int(a1), int(a2), int(a3), float(bond), float(angle), float(dihedral)]
                            except:
                                return [atom, int(a1), int(a2), int(a3), str(bond), str(angle), str(dihedral)]
                    elif coord_type == 'gzmt':
                        def splitter(line):
                            vals = line.split()[:7]
                            if len(vals) == 7:
                                atom, a1, bond, a2, angle, a3, dihedral = vals
                                return [atom, int(a1), float(bond), int(a2), float(angle), int(a3), float(dihedral)]
                            elif len(vals) == 5:
                                return [vals[0], int(vals[1]), float(vals[2]), int(vals[3]), float(vals[4])]
                            elif len(vals) == 3:
                                return [vals[0], int(vals[1]), float(vals[2])]
                            elif len(vals) == 1:
                                return [vals[0]]
                            self.logger.warning('Incorrect number of atoms in input geometry.')
                    elif 'file' in coord_type:
                        pass
                    else:
                        self.logger.warning('Invalid coordinate type.')

                    if 'file' not in coord_type:
                        for line in lines_iter:
                            if not line:
                                continue
                            if line[0] == '#' or line.strip(' ') == '\n':
                                continue
                            if line.strip()[0] == '*' or line.strip() == "end":
                                break
                            # Strip basis specification that can appear after coordinates
                            line = line.split('newGTO')[0].strip()
                            coords.append(splitter(line))
            self.metadata['keywords'] = keywords
            self.metadata['coords'] = coords
        # If the calculations is a unrelaxed parameter scan then immediately following the 
        # input file block is the following section:
                
        # | 12> **                         ****END OF INPUT****
        # ================================================================================
        # 
        #                        ******************************
        #                        * Parameter Scan Calculation *
        #                        ******************************
        # 
        # Trajectory settings:
        #     -> SCF surface will be mapped
        # 
        # There are 1 parameter(s) to be scanned
        #              R: range=   0.58220000 ..   5.08220000  steps=   46
        # There will be   46 energy evaluations
        #
        #
        # following this each calculation has the following block at the start
        #
        #         *************************************************************
        #                                TRAJECTORY STEP   1
        #                  R  :   0.58220000
        #         *************************************************************

        if 'Parameter Scan Calculation' in line:
            self.skip_lines(inputfile,['s', 'b', 'Trajectory settings', 'Surface information', 'b'])
            line = next(inputfile)
            num_params = int(line.strip().split()[2])
            for i in range(num_params):
                line = next(inputfile).strip()
                self.append_attribute('scannames', line.split(':')[0])
        if 'TRAJECTORY STEP' in line:
            current_params = []
            for i in range(len(self.scannames)):
                line = next(inputfile)
                current_params.append(float(line.split(':')[-1].strip()))
            self.append_attribute('scanparm', tuple(current_params))

        # If the calculations is a relaxed parameter scan then immediately following the 
        # input file block is the following section:

        #                        ******************************
        #                        *    Relaxed Surface Scan    *
        #                        ******************************
        # 
        #         Dihedral (  9,   8,   3,   2):   range=   0.00000000 .. 360.00000000  steps =   12
        # 
        # There is 1 parameter to be scanned.
        # There will be   12 constrained geometry optimizations.
        # 
        # 
        #          *************************************************************
        #          *               RELAXED SURFACE SCAN STEP   1               *
        #          *                                                           *
        #          *   Dihedral (  9,   8,   3,   2)  :   0.00000000           *
        #          *************************************************************

        if 'Relaxed Surface Scan' in line:
            self.skip_lines(inputfile,['s', 'b'])
            line = next(inputfile)
            while not line.isspace():
                line = line.strip()
                self.append_attribute('scannames', line.split(':')[0])
                line = next(inputfile)
            line = next(inputfile)
            num_params = int(line.strip().split()[2])
        
        if line[0:15] == "Number of atoms":

            natom = int(line.split()[-1])
            self.set_attribute('natom', natom)

        if line[1:13] == "Total Charge":

            charge = int(line.split()[-1])
            self.set_attribute('charge', charge)

            line = next(inputfile)

            mult = int(line.split()[-1])
            self.set_attribute('mult', mult)

        if line[1:18] == "Symmetry handling":
            self.uses_symmetry = True

            line = next(inputfile)
            assert "Point group" in line
            point_group_full = line.split()[3].lower()
            line = next(inputfile)
            assert "Used point group" in line
            point_group_abelian = line.split()[4].lower()
            line = next(inputfile)
            assert "Number of irreps" in line
            nirrep = int(line.split()[4])
            for n in range(nirrep):
                line = next(inputfile)
                assert "symmetry adapted basis functions" in line
                irrep = line[8:13]
                if not hasattr(self, 'symlabels'):
                    self.symlabels = []
                self.symlabels.append(self.normalisesym(irrep))

            self.metadata['symmetry_detected'] = point_group_full
            self.metadata['symmetry_used'] = point_group_abelian

        # SCF convergence output begins with:
        #
        # --------------
        # SCF ITERATIONS
        # --------------
        #
        # However, there are two common formats which need to be handled, implemented as separate functions.
        if line.strip() == "SCF ITERATIONS":

            self.skip_line(inputfile, 'dashes')

            line = next(inputfile)
            columns = line.split()
            # "Starting incremental Fock matrix formation" doesn't
            # necessarily appear before the extended format.
            if not columns:
                self.parse_scf_expanded_format(inputfile, columns)
            # A header with distinct columns indicates the condensed
            # format.
            elif columns[1] == "Energy":
                self.parse_scf_condensed_format(inputfile, columns)
            # Assume the extended format.
            else:
                self.parse_scf_expanded_format(inputfile, columns)

        # Information about the final iteration, which also includes the convergence
        # targets and the convergence values, is printed separately, in a section like this:
        #
        #       *****************************************************
        #       *                     SUCCESS                       *
        #       *           SCF CONVERGED AFTER   9 CYCLES          *
        #       *****************************************************
        #
        # ...
        #
        # Total Energy       :         -382.04963064 Eh          -10396.09898 eV
        #
        # ...
        #
        # -------------------------   ----------------
        # FINAL SINGLE POINT ENERGY     -382.049630637
        # -------------------------   ----------------
        #
        # We cannot use this last message as a stop condition in general, because
        # often there is vibrational output before it. So we use the 'Total Energy'
        # line. However, what comes after that is different for single point calculations
        # and in the inner steps of geometry optimizations.
        if "SCF CONVERGED AFTER" in line:

            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            if not hasattr(self, "scfvalues"):
                self.scfvalues = []
            if not hasattr(self, "scftargets"):
                self.scftargets = []

            while not "Total Energy       :" in line:
                line = next(inputfile)
            energy = utils.convertor(float(line.split()[3]), "hartree", "eV")
            self.scfenergies.append(energy)

            self._append_scfvalues_scftargets(inputfile, line)

        # Sometimes the SCF does not converge, but does not halt the
        # the run (like in bug 3184890). In this this case, we should
        # remain consistent and use the energy from the last reported
        # SCF cycle. In this case, ORCA print a banner like this:
        #
        #       *****************************************************
        #       *                     ERROR                         *
        #       *           SCF NOT CONVERGED AFTER   8 CYCLES      *
        #       *****************************************************
        if "SCF NOT CONVERGED AFTER" in line:

            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            if not hasattr(self, "scfvalues"):
                self.scfvalues = []
            if not hasattr(self, "scftargets"):
                self.scftargets = []

            energy = utils.convertor(self.scfvalues[-1][-1][0], "hartree", "eV")
            self.scfenergies.append(energy)

            self._append_scfvalues_scftargets(inputfile, line)

        """
-------------------------------------------------------------------------------
                          DFT DISPERSION CORRECTION

                              DFTD3 V3.1  Rev 1
                              USING zero damping
-------------------------------------------------------------------------------
The omegaB97X-D3 functional is recognized. Fit by Chai et al.
Active option DFTDOPT                   ...         3

molecular C6(AA) [au] = 9563.878941


            DFT-D V3
 parameters
 s6 scaling factor         :     1.0000
 rs6 scaling factor        :     1.2810
 s8 scaling factor         :     1.0000
 rs8 scaling factor        :     1.0940
 Damping factor alpha6     :    14.0000
 Damping factor alpha8     :    16.0000
 ad hoc parameters k1-k3   :    16.0000     1.3333    -4.0000

 Edisp/kcal,au: -10.165629059768  -0.016199959356
 E6   /kcal   :  -4.994512983
 E8   /kcal   :  -5.171116077
 % E8         :  50.868628459

-------------------------   ----------------
Dispersion correction           -0.016199959
-------------------------   ----------------
"""
        if "DFT DISPERSION CORRECTION" in line:
            # A bunch of parameters are printed the first time dispersion is called
            # However, they vary wildly in form and number, making parsing problematic
            line = next(inputfile)
            while 'Dispersion correction' not in line:
                line = next(inputfile)
            dispersion = utils.convertor(float(line.split()[-1]), "hartree", "eV")
            self.append_attribute("dispersionenergies", dispersion)

        # The convergence targets for geometry optimizations are printed at the
        # beginning of the output, although the order and their description is
        # different than later on. So, try to standardize the names of the criteria
        # and save them for later so that we can get the order right.
        #
        #                        *****************************
        #                        * Geometry Optimization Run *
        #                        *****************************
        #
        # Geometry optimization settings:
        # Update method            Update   .... BFGS
        # Choice of coordinates    CoordSys .... Redundant Internals
        # Initial Hessian          InHess   .... Almoef's Model
        #
        # Convergence Tolerances:
        # Energy Change            TolE     ....  5.0000e-06 Eh
        # Max. Gradient            TolMAXG  ....  3.0000e-04 Eh/bohr
        # RMS Gradient             TolRMSG  ....  1.0000e-04 Eh/bohr
        # Max. Displacement        TolMAXD  ....  4.0000e-03 bohr
        # RMS Displacement         TolRMSD  ....  2.0000e-03 bohr
        #
        if line[25:50] == "Geometry Optimization Run":

            stars = next(inputfile)
            blank = next(inputfile)

            line = next(inputfile)
            while line[0:23] != "Convergence Tolerances:":
                line = next(inputfile)

            if hasattr(self, 'geotargets'):
                self.logger.warning('The geotargets attribute should not exist yet. There is a problem in the parser.')
            self.geotargets = []
            self.geotargets_names = []

            # There should always be five tolerance values printed here.
            for i in range(5):
                line = next(inputfile)
                name = line[:25].strip().lower().replace('.', '').replace('displacement', 'step')
                target = float(line.split()[-2])
                self.geotargets_names.append(name)
                self.geotargets.append(target)

        # The convergence targets for relaxed surface scan steps are printed at the
        # beginning of the output, although the order and their description is
        # different than later on. So, try to standardize the names of the criteria
        # and save them for later so that we can get the order right.
        #
        #         *************************************************************
        #         *               RELAXED SURFACE SCAN STEP  12               *
        #         *                                                           *
        #         *   Dihedral ( 11,  10,   3,   4)  : 180.00000000           *
        #         *************************************************************
        #
        # Geometry optimization settings:
        # Update method            Update   .... BFGS
        # Choice of coordinates    CoordSys .... Redundant Internals
        # Initial Hessian          InHess   .... Almoef's Model
        #
        # Convergence Tolerances:
        # Energy Change            TolE     ....  5.0000e-06 Eh
        # Max. Gradient            TolMAXG  ....  3.0000e-04 Eh/bohr
        # RMS Gradient             TolRMSG  ....  1.0000e-04 Eh/bohr
        # Max. Displacement        TolMAXD  ....  4.0000e-03 bohr
        # RMS Displacement         TolRMSD  ....  2.0000e-03 bohr
        if 'RELAXED SURFACE SCAN STEP' in line:
            self.skip_lines(inputfile,['b'])
            current_params = []
            for i in range(len(self.scannames)):
                line = next(inputfile)
                line = line.replace('*','')
                current_params.append(float(line.split(':')[-1].strip()))
            self.append_attribute('scanparm', tuple(current_params))

            self.is_relaxed_scan = True
            while "Convergence Tolerances:" not in line:
                line = next(inputfile)

            self.geotargets = []
            self.geotargets_names = []

            # There should always be five tolerance values printed here.
            for i in range(5):
                line = next(inputfile)
                name = line[:25].strip().lower().replace('.', '').replace('displacement', 'step')
                target = float(line.split()[-2])
                self.geotargets_names.append(name)
                self.geotargets.append(target)
        

        # Moller-Plesset energies.
        #
        # ---------------------------------------
        # MP2 TOTAL ENERGY:      -76.112119693 Eh
        # ---------------------------------------
        if 'MP2 TOTAL ENERGY' in line[:16]:

            if not hasattr(self, 'mpenergies'):
                self.metadata['methods'].append('MP2')
                self.mpenergies = []

            self.mpenergies.append([])
            mp2energy = utils.float(line.split()[-2])
            self.mpenergies[-1].append(utils.convertor(mp2energy, 'hartree', 'eV'))

        # MP2 energy output line is different for MP3, since it uses the MDCI
        # code, which is also in charge of coupled cluster.
        #
        # MP3 calculation:
        # E(MP2)  =    -76.112119775   EC(MP2)=    -0.128216417
        # E(MP3)  =    -76.113783480   EC(MP3)=    -0.129880122  E3=    -0.001663705
        #
        # CCSD calculation:
        # E(MP2)                                     ...     -0.393722942
        # Initial E(tot)                             ...  -1639.631576169
        # <T|T>                                      ...      0.087231847
        # Number of pairs included                   ... 55
        # Total number of pairs                      ... 55
        if 'E(MP2)' in line:

            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []

            self.mpenergies.append([])
            mp2energy = utils.float(line.split()[-1])
            self.mpenergies[-1].append(utils.convertor(mp2energy, 'hartree', 'eV'))

            line = next(inputfile)
            if line[:6] == 'E(MP3)':
                self.metadata['methods'].append('MP3')
                mp3energy = utils.float(line.split()[2])
                self.mpenergies[-1].append(utils.convertor(mp3energy, 'hartree', 'eV'))
            else:
                assert line[:14] == 'Initial E(tot)'

        # ----------------------
        # COUPLED CLUSTER ENERGY
        # ----------------------
        #
        # E(0)                                       ...  -1639.237853227
        # E(CORR)                                    ...     -0.360153516
        # E(TOT)                                     ...  -1639.598006742
        # Singles Norm <S|S>**1/2                    ...      0.176406354  
        # T1 diagnostic                              ...      0.039445660  
        if line[:22] == 'COUPLED CLUSTER ENERGY':
            self.skip_lines(inputfile, ['d', 'b'])
            line = next(inputfile)
            assert line[:4] == 'E(0)'
            scfenergy = utils.convertor(utils.float(line.split()[-1]), 'hartree', 'eV')
            line = next(inputfile)
            assert line[:7] == 'E(CORR)'
            while 'E(TOT)' not in line:
                line = next(inputfile)
            self.append_attribute(
                'ccenergies',
                utils.convertor(utils.float(line.split()[-1]), 'hartree', 'eV')
            )
            line = next(inputfile)
            assert line[:23] == 'Singles Norm <S|S>**1/2'
            line = next(inputfile)
            self.metadata["t1_diagnostic"] = utils.float(line.split()[-1])

        # ------------------
        # CARTESIAN GRADIENT
        # ------------------
        #
        # 1   H   :    0.000000004    0.019501450   -0.021537091
        # 2   O   :    0.000000054   -0.042431648    0.042431420
        # 3   H   :    0.000000004    0.021537179   -0.019501388
        #
        # ORCA MP2 module has different signal than 'CARTESIAN GRADIENT'.
        #
        # The final MP2 gradient
        # 0:   0.01527469  -0.00292883   0.01125000
        # 1:   0.00098782  -0.00040549   0.00196825
        # 2:  -0.01626251   0.00333431  -0.01321825
        if line[:18] == 'CARTESIAN GRADIENT' or line[:22] == 'The final MP2 gradient':

            grads = []
            if line[:18] == 'CARTESIAN GRADIENT':
                self.skip_lines(inputfile, ['dashes', 'blank'])

            line = next(inputfile).strip()
            if 'CONSTRAINED CARTESIAN COORDINATES' in line:
                self.skip_line(
                    inputfile, 'constrained Cartesian coordinate warning'
                )
                line = next(inputfile).strip()

            while line:
                tokens = line.split()
                x, y, z = float(tokens[-3]), float(tokens[-2]), float(tokens[-1])
                grads.append((x, y, z))
                line = next(inputfile).strip()

            if not hasattr(self, 'grads'):
                self.grads = []
            self.grads.append(grads)


        # After each geometry optimization step, ORCA prints the current convergence
        # parameters and the targets (again), so it is a good idea to check that they
        # have not changed. Note that the order of these criteria here are different
        # than at the beginning of the output, so make use of the geotargets_names created
        # before and save the new geovalues in correct order.
        #
        #          ----------------------|Geometry convergence|---------------------
        #          Item                value                 Tolerance   Converged
        #          -----------------------------------------------------------------
        #          Energy change       0.00006021            0.00000500      NO
        #          RMS gradient        0.00031313            0.00010000      NO
        #          RMS step            0.01596159            0.00200000      NO
        #          MAX step            0.04324586            0.00400000      NO
        #          ....................................................
        #          Max(Bonds)      0.0218      Max(Angles)    2.48
        #          Max(Dihed)        0.00      Max(Improp)    0.00
        #          -----------------------------------------------------------------
        #
        if line[33:53] == "Geometry convergence":

            headers = next(inputfile)
            dashes = next(inputfile)

            names = []
            values = []
            targets = []
            line = next(inputfile)
            # Handle both the dots only and dashes only cases
            while len(list(set(line.strip()))) != 1:
                name = line[10:28].strip().lower()
                tokens = line.split()
                value = float(tokens[2])
                target = float(tokens[3])
                names.append(name)
                values.append(value)
                targets.append(target)
                line = next(inputfile)

            # The energy change is normally not printed in the first iteration, because
            # there was no previous energy -- in that case assume zero. There are also some
            # edge cases where the energy change is not printed, for example when internal
            # angles become improper and internal coordinates are rebuilt as in regression
            # CuI-MePY2-CH3CN_optxes, and in such cases use NaN.
            newvalues = []
            for i, n in enumerate(self.geotargets_names):
                if (n == "energy change") and (n not in names):
                    if self.is_relaxed_scan:
                        newvalues.append(0.0)
                    else:
                        newvalues.append(numpy.nan)
                else:
                    newvalues.append(values[names.index(n)])
                    assert targets[names.index(n)] == self.geotargets[i]

            self.append_attribute("geovalues", newvalues)

        """ Grab cartesian coordinates
        ---------------------------------
        CARTESIAN COORDINATES (ANGSTROEM)
        ---------------------------------
        H      0.000000    0.000000    0.000000
        O      0.000000    0.000000    1.000000
        H      0.000000    1.000000    1.000000
        """
        if line[0:33] == "CARTESIAN COORDINATES (ANGSTROEM)":
            next(inputfile)

            atomnos = []
            atomcoords = []
            line = next(inputfile)
            while len(line) > 1:
                atom, x, y, z = line.split()
                if atom[-1] != ">":
                    atomnos.append(self.table.number[atom])
                    atomcoords.append([float(x), float(y), float(z)])
                line = next(inputfile)

            self.set_attribute('natom', len(atomnos))
            self.set_attribute('atomnos', atomnos)
            self.append_attribute("atomcoords", atomcoords)

        """ Grab atom masses
        ----------------------------
        CARTESIAN COORDINATES (A.U.)
        ----------------------------
        NO LB      ZA    FRAG     MASS         X           Y           Z
        0 H     1.0000    0     1.008    0.000000    0.000000    0.000000
        1 O     8.0000    0    15.999    0.000000    0.000000    1.889726
        2 H     1.0000    0     1.008    0.000000    1.889726    1.889726
        """
        if line[0:28] == "CARTESIAN COORDINATES (A.U.)" and not hasattr(self, 'atommasses'):
            next(inputfile)
            next(inputfile)

            line = next(inputfile)
            self.atommasses = []
            while len(line) > 1:
                if line[:32] == '* core charge reduced due to ECP':
                    break
                if line.strip() == "> coreless ECP center with (optional) point charge":
                    break
                no, lb, za, frag, mass, x, y, z = line.split()
                if lb[-1] != ">":
                    self.atommasses.append(float(mass))
                line = next(inputfile)

        if line[21:68] == "FINAL ENERGY EVALUATION AT THE STATIONARY POINT":
            if not hasattr(self, 'optdone'):
                self.optdone = []
            self.optdone.append(len(self.atomcoords))

        if "The optimization did not converge" in line:
            if not hasattr(self, 'optdone'):
                self.optdone = []

        if line[0:16] == "ORBITAL ENERGIES":

            self.skip_lines(inputfile, ['d', 'text', 'text'])

            self.mooccnos = [[]]
            self.moenergies = [[]]
            self.mosyms = [[]]

            line = next(inputfile)
            while len(line) > 20:  # restricted calcs are terminated by ------
                info = line.split()
                mooccno = int(float(info[1]))
                moenergy = float(info[2])
                mosym = 'A'
                if self.uses_symmetry:
                    mosym = self.normalisesym(info[4].split('-')[1])
                self.mooccnos[0].append(mooccno)
                self.moenergies[0].append(utils.convertor(moenergy, "hartree", "eV"))
                self.mosyms[0].append(mosym)
                line = next(inputfile)

            line = next(inputfile)

            # handle beta orbitals for UHF
            if line[17:35] == "SPIN DOWN ORBITALS":
                text = next(inputfile)

                self.mooccnos.append([])
                self.moenergies.append([])
                self.mosyms.append([])

                line = next(inputfile)
                while len(line) > 20:  # actually terminated by ------
                    info = line.split()
                    mooccno = int(float(info[1]))
                    moenergy = float(info[2])
                    mosym = 'A'
                    if self.uses_symmetry:
                        mosym = self.normalisesym(info[4].split('-')[1])
                    self.mooccnos[1].append(mooccno)
                    self.moenergies[1].append(utils.convertor(moenergy, "hartree", "eV"))
                    self.mosyms[1].append(mosym)
                    line = next(inputfile)

            if not hasattr(self, 'homos'):
                doubly_occupied = self.mooccnos[0].count(2)
                singly_occupied = self.mooccnos[0].count(1)
                # Restricted closed-shell.
                if doubly_occupied > 0 and singly_occupied == 0:
                    self.set_attribute('homos', [doubly_occupied - 1])
                # Restricted open-shell.
                elif doubly_occupied > 0 and singly_occupied > 0:
                    self.set_attribute('homos', [doubly_occupied + singly_occupied - 1,
                                                 doubly_occupied - 1])
                # Unrestricted.
                else:
                    assert len(self.moenergies) == 2
                    assert doubly_occupied == 0
                    assert self.mooccnos[1].count(2) == 0
                    nbeta = self.mooccnos[1].count(1)
                    self.set_attribute('homos', [singly_occupied - 1, nbeta - 1])

        # So nbasis was parsed at first with the first pattern, but it turns out that
        # semiempirical methods (at least AM1 as reported by Julien Idé) do not use this.
        # For this reason, also check for the second patterns, and use it as an assert
        # if nbasis was already parsed. Regression PCB_1_122.out covers this test case.
        if line[1:32] == "# of contracted basis functions":
            self.set_attribute('nbasis', int(line.split()[-1]))
        if line[1:27] == "Basis Dimension        Dim":
            self.set_attribute('nbasis', int(line.split()[-1]))

        if line[0:14] == "OVERLAP MATRIX":

            self.skip_line(inputfile, 'dashes')

            self.aooverlaps = numpy.zeros((self.nbasis, self.nbasis), "d")
            for i in range(0, self.nbasis, 6):
                self.updateprogress(inputfile, "Overlap")

                header = next(inputfile)
                size = len(header.split())

                for j in range(self.nbasis):
                    line = next(inputfile)
                    broken = line.split()
                    self.aooverlaps[j, i:i+size] = list(map(float, broken[1:size+1]))

        # Molecular orbital coefficients are parsed here, but also related things
        #like atombasis and aonames if possible.
        #
        # Normally the output is easy to parse like this:
        # ------------------
        # MOLECULAR ORBITALS
        # ------------------
        #                       0         1         2         3         4         5   
        #                  -19.28527 -19.26828 -19.26356 -19.25801 -19.25765 -19.21471
        #                    2.00000   2.00000   2.00000   2.00000   2.00000   2.00000
        #                   --------  --------  --------  --------  --------  --------
        #   0C   1s         0.000002 -0.000001  0.000000  0.000000 -0.000000  0.000001
        #   0C   2s        -0.000007  0.000006 -0.000002 -0.000000  0.000001 -0.000003
        #   0C   3s        -0.000086 -0.000061  0.000058 -0.000033 -0.000027 -0.000058
        # ...
        #
        # But when the numbers get big, things get yucky since ORCA does not use
        # fixed width formatting for the floats, and does not insert extra spaces
        # when the numbers get wider. So things get stuck together overflowing columns,
        # like this:
        #   12C   6s       -11.608845-53.775398161.302640-76.633779 29.914985 22.083999
        #
        # One assumption that seems to hold is that there are always six significant
        # digits in the coefficients, so we can try to use that to delineate numbers
        # when the parsing gets rough. This is what we do below with a regex, and a case
        # like this is tested in regression ORCA/ORCA4.0/invalid-literal-for-float.out
        # which was reported in https://github.com/cclib/cclib/issues/629
        if line[0:18] == "MOLECULAR ORBITALS":

            self.skip_line(inputfile, 'dashes')

            aonames = []
            atombasis = [[] for i in range(self.natom)]
            mocoeffs = [numpy.zeros((self.nbasis, self.nbasis), "d")]

            for spin in range(len(self.moenergies)):

                if spin == 1:
                    self.skip_line(inputfile, 'blank')
                    mocoeffs.append(numpy.zeros((self.nbasis, self.nbasis), "d"))

                for i in range(0, self.nbasis, 6):

                    self.updateprogress(inputfile, "Coefficients")

                    self.skip_lines(inputfile, ['numbers', 'energies', 'occs'])
                    dashes = next(inputfile)

                    for j in range(self.nbasis):
                        line = next(inputfile)

                        # Only need this in the first iteration.
                        if spin == 0 and i == 0:
                            atomname = line[3:5].split()[0]
                            num = int(line[0:3])
                            orbital = line.split()[1].upper()

                            aonames.append(f"{atomname}{int(num + 1)}_{orbital}")
                            atombasis[num].append(j)

                        # This regex will tease out all number with exactly
                        # six digits after the decimal point.
                        coeffs = re.findall(r'-?\d+\.\d{6}', line)

                        # Something is very wrong if this does not hold.
                        assert len(coeffs) <= 6

                        mocoeffs[spin][i:i+len(coeffs), j] = [float(c) for c in coeffs]

            self.set_attribute('aonames', aonames)
            self.set_attribute('atombasis', atombasis)
            self.set_attribute("mocoeffs", mocoeffs)

        # Basis set information
        # ORCA prints this out in a somewhat indirect fashion.
        # Therefore, parsing occurs in several steps:
        # 1. read which atom belongs to which basis set group
        if line[0:21] == "BASIS SET INFORMATION":
            line = next(inputfile)
            line = next(inputfile)

            self.tmp_atnames = [] # temporary attribute, needed later
            while(not line[0:5] == '-----'):
                if line[0:4] == "Atom":
                    self.tmp_atnames.append(line[8:12].strip())
                line = next(inputfile)

        # 2. Read information for the basis set groups
        if line[0:25] == "BASIS SET IN INPUT FORMAT":
            line = next(inputfile)
            line = next(inputfile)

            # loop over basis set groups
            gbasis_tmp = {}
            while(not line[0:5] == '-----'):
                if line[1:7] == 'NewGTO':
                    bas_atname = line.split()[1]
                    gbasis_tmp[bas_atname] = []

                    line = next(inputfile)
                    # loop over contracted GTOs
                    while(not line[0:6] == '  end;'):
                        words = line.split()
                        ang = words[0]
                        nprim = int(words[1])

                        # loop over primitives
                        coeff = []
                        for iprim in range(nprim):
                            words = next(inputfile).split()
                            coeff.append( (float(words[1]), float(words[2])) )
                        gbasis_tmp[bas_atname].append((ang, coeff))

                        line = next(inputfile)
                line = next(inputfile)

            # 3. Assign the basis sets to gbasis
            self.gbasis = []
            for bas_atname in self.tmp_atnames:
                self.gbasis.append(gbasis_tmp[bas_atname])
            del self.tmp_atnames

        """
        --------------------------
        THERMOCHEMISTRY AT 298.15K
        --------------------------

        Temperature         ... 298.15 K
        Pressure            ... 1.00 atm
        Total Mass          ... 130.19 AMU

        Throughout the following assumptions are being made:
          (1) The electronic state is orbitally nondegenerate
          ...

        freq.      45.75  E(vib)   ...       0.53 
        freq.      78.40  E(vib)   ...       0.49
        ...


        ------------
        INNER ENERGY
        ------------

        The inner energy is: U= E(el) + E(ZPE) + E(vib) + E(rot) + E(trans)
             E(el)   - is the total energy from the electronic structure calc
             ...

        Summary of contributions to the inner energy U:
        Electronic energy                ...   -382.05075804 Eh
        ...
        """
        if line.strip().startswith('THERMOCHEMISTRY AT'):

            self.skip_lines(inputfile, ['dashes', 'blank'])
            self.temperature = float(next(inputfile).split()[2])
            self.pressure = float(next(inputfile).split()[2])
            total_mass = float(next(inputfile).split()[3])

            # Vibrations, rotations, and translations
            line = next(inputfile)
            while line[:17] != 'Electronic energy':
                line = next(inputfile)
            self.electronic_energy = float(line.split()[3])
            self.set_attribute("zpve", float(next(inputfile).split()[4]))
            thermal_vibrational_correction = float(next(inputfile).split()[4])
            thermal_rotional_correction = float(next(inputfile).split()[4])
            thermal_translational_correction = float(next(inputfile).split()[4])
            self.skip_lines(inputfile, ['dashes'])
            total_thermal_energy = float(next(inputfile).split()[3])

            # Enthalpy
            while line[:17] != 'Total free energy':
                line = next(inputfile)
            thermal_enthalpy_correction = float(next(inputfile).split()[4])
            next(inputfile)

            # For a single atom, ORCA provides the total free energy or inner energy
            # which includes a spurious vibrational correction (see #817 for details).
            if self.natom > 1:
                self.enthalpy = float(next(inputfile).split()[3])
            else:
                self.enthalpy = self.electronic_energy + thermal_translational_correction

            # Entropy
            while line[:18] != 'Electronic entropy':
                line = next(inputfile)
            electronic_entropy = float(line.split()[3])
            vibrational_entropy = float(next(inputfile).split()[3])
            rotational_entropy = float(next(inputfile).split()[3])
            translational_entropy = float(next(inputfile).split()[3])
            self.skip_lines(inputfile, ['dashes'])

            # ORCA prints -inf for single atom entropy.
            if self.natom > 1:
                self.entropy = float(next(inputfile).split()[4]) / self.temperature
            else:
                self.entropy = (electronic_entropy + translational_entropy) / self.temperature

            while (line[:25] != 'Final Gibbs free enthalpy') and (line[:23] != 'Final Gibbs free energy'):
                line = next(inputfile)
            self.skip_lines(inputfile, ['dashes'])


            # ORCA prints -inf for sinle atom free energy.
            if self.natom > 1:
                self.freeenergy = float(line.split()[5])
            else:
                self.freeenergy = self.enthalpy - self.temperature * self.entropy

        # Read TDDFT information
        if any(x in line for x in ("TD-DFT/TDA EXCITED", "TD-DFT EXCITED")):
            # Could be singlets or triplets
            if line.find("SINGLETS") >= 0:
                sym = "Singlet"
            elif line.find("TRIPLETS") >= 0:
                sym = "Triplet"
            else:
                sym = "Not specified"

            etsecs = []
            etenergies = []
            etsyms = []

            lookup = {'a': 0, 'b': 1}
            line = next(inputfile)
            while line.find("STATE") < 0:
                line = next(inputfile)
            # Contains STATE or is blank
            while line.find("STATE") >= 0:
                broken = line.split()
                etenergies.append(float(broken[7]))
                etsyms.append(sym)
                line = next(inputfile)
                sec = []
                # Contains SEC or is blank
                while line.strip():
                    start = line[0:8].strip()
                    start = (int(start[:-1]), lookup[start[-1]])
                    end = line[10:17].strip()
                    end = (int(end[:-1]), lookup[end[-1]])
                    # Coeffients are not printed for RPA, only
                    # TDA/CIS.
                    contrib = line[35:47].strip()
                    try:
                        contrib = float(contrib)
                    except ValueError:
                        contrib = numpy.nan
                    sec.append([start, end, contrib])
                    line = next(inputfile)
                    # ORCA 5.0 seems to print symmetry at end of block listing transitions
                    if 'Symmetry' in line:
                        line = next(inputfile)
                etsecs.append(sec)
                line = next(inputfile)

            self.extend_attribute('etenergies', etenergies)
            self.extend_attribute('etsecs', etsecs)
            self.extend_attribute('etsyms', etsyms)

        # Parse the various absorption spectra for TDDFT and ROCIS.
        if 'ABSORPTION SPECTRUM' in line or 'ELECTRIC DIPOLE' in line:
            # CASSCF has an anomalous printing of ABSORPTION SPECTRUM.
            if line[:-1] == 'ABSORPTION SPECTRUM':
                return

            line = line.strip()

            # Standard header, occasionally changes
            header = ['d', 'header', 'header', 'd']
            energy_intensity = None

            if line == "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS":
                def energy_intensity(line):
                    """ TDDFT and related methods standard method of output
-----------------------------------------------------------------------------
         ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
-----------------------------------------------------------------------------
State   Energy  Wavelength   fosc         T2         TX        TY        TZ
        (cm-1)    (nm)                  (au**2)     (au)      (au)      (au)
-----------------------------------------------------------------------------
   1 5184116.7      1.9   0.040578220   0.00258  -0.05076  -0.00000  -0.00000
"""
                    try:
                        state, energy, wavelength, intensity, t2, tx, ty, tz = line.split()
                    except ValueError as e:
                        # Must be spin forbidden and thus no intensity
                        energy = line.split()[1]
                        intensity = 0
                    return energy, intensity

            # Check for variations
            elif line == 'COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM' or \
               line == 'COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (origin adjusted)':
                def energy_intensity(line):
                    """ TDDFT with DoQuad == True
------------------------------------------------------------------------------------------------------
                COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM
------------------------------------------------------------------------------------------------------
State   Energy Wavelength    D2        m2        Q2         D2+m2+Q2       D2/TOT    m2/TOT    Q2/TOT
        (cm-1)   (nm)                (*1e6)    (*1e6)
------------------------------------------------------------------------------------------------------
   1 61784150.6      0.2   0.00000   0.00000   3.23572   0.00000323571519   0.00000   0.00000   1.00000
"""
                    state, energy, wavelength, d2, m2, q2, intensity, d2_contrib, m2_contrib, q2_contrib = line.split()
                    return energy, intensity

            elif line == 'COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (Origin Independent, Length Representation)':
                def energy_intensity(line):
                    """ TDDFT with doQuad == True (Origin Independent Length Representation)
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                    COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (Origin Independent, Length Representation)
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
State  Energy   Wavelength       D2            m2              Q2               DM             DO               D2+m2+Q2+DM+DO          D2/TOT          m2/TOT          Q2/TOT         DM/TOT          DO/TOT
       (cm-1)      (nm)                      (*1e6)          (*1e6)           (*1e6)         (*1e6)
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   1 61784150.6      0.2      0.00000         0.00000         3.23572         0.00000         0.00000         0.00000323571519         0.00000         0.00000         1.00000         0.00000          0.00000
   2 61793079.3      0.2      0.00000         0.00000         2.85949         0.00000        -0.00000         0.00000285948800         0.00000         0.00000         1.00000         0.00000         -0.00000
"""
                    vals = line.split()
                    if len(vals) < 14:
                        return vals[1], 0
                    return vals[1], vals[8]

            elif line[:5] == 'X-RAY' and \
                (line[6:23] == 'EMISSION SPECTRUM' or line[6:25] == 'ABSORPTION SPECTRUM'):
                def energy_intensity(line):
                    """ X-Ray from XES (emission or absorption, electric or velocity dipole moments)
-------------------------------------------------------------------------------------
          X-RAY ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
-------------------------------------------------------------------------------------
       Transition          Energy           INT             TX        TY        TZ
                            (eV)        (normalized)       (au)      (au)      (au)
-------------------------------------------------------------------------------------
    1   90a ->    0a      8748.824     0.000002678629     0.00004  -0.00001   0.00003
"""
                    state, start, arrow, end, energy, intensity, tx, ty, tz = line.split()
                    return energy, intensity

            elif line[:70] == 'COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE X-RAY':
                header = ['header', 'd', 'header', 'd', 'header', 'header', 'd']
                def energy_intensity(line):
                    """ XAS with quadrupole (origin adjusted)
-------------------------------------------------------------------------------------------------------------------------------
          COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE X-RAY ABSORPTION SPECTRUM
                                      (origin adjusted)
-------------------------------------------------------------------------------------------------------------------------------
                                                        INT (normalized)
                                     ---------------------------------------------------------
       Transition         Energy        D2             M2             Q2           D2+M2+Q2       D2/TOT     M2/TOT     Q2/TOT
                           (eV)                      (*1e6)         (*1e6)
-------------------------------------------------------------------------------------------------------------------------------
    1   90a ->    0a     8748.824    0.000000       0.000292       0.003615     0.000000027512   0.858012   0.010602   0.131386
"""
                    state, start, arrow, end, energy, d2, m2, q2, intensity, d2_contrib, m2_contrib, q2_contrib = line.split()
                    return energy, intensity

            elif line[:55] == 'SPIN ORBIT CORRECTED ABSORPTION SPECTRUM VIA TRANSITION':
                def energy_intensity(line):
                    """ ROCIS dipole approximation with SOC == True (electric or velocity dipole moments)
-------------------------------------------------------------------------------
SPIN ORBIT CORRECTED ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
-------------------------------------------------------------------------------
States    Energy  Wavelength   fosc         T2         TX        TY        TZ
          (cm-1)    (nm)                  (au**2)     (au)      (au)      (au)
-------------------------------------------------------------------------------
 0  1       0.0      0.0   0.000000000   0.00000   0.00000   0.00000   0.00000
 0  2 5184116.4      1.9   0.020288451   0.00258   0.05076   0.00003   0.00000
"""
                    state, state2, energy, wavelength, intensity, t2, tx, ty, tz = line.split()
                    return energy, intensity

            elif line[:79] == 'ROCIS COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM' \
                 or line[:87] == 'SOC CORRECTED COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM':
                def energy_intensity(line):
                    """ ROCIS with DoQuad = True and SOC = True (also does origin adjusted)
------------------------------------------------------------------------------------------------------
          ROCIS COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM
------------------------------------------------------------------------------------------------------
States  Energy Wavelength    D2        m2        Q2         D2+m2+Q2       D2/TOT    m2/TOT    Q2/TOT
        (cm-1)   (nm)                (*1e6)    (*1e6)     (*population)
------------------------------------------------------------------------------------------------------
 0  1       0.0      0.0   0.00000   0.00000   0.00000   0.00000000000000   0.00000   0.00000   0.00000
 0  2 669388066.6      0.0   0.00000   0.00000   0.00876   0.00000000437784   0.00000   0.00000   1.00000
"""
                    state, state2, energy, wavelength, d2, m2, q2, intensity, d2_contrib, m2_contrib, q2_contrib = line.split()
                    return energy, intensity

            # Clashes with Orca 2.6 (and presumably before) TDDFT absorption spectrum printing
            elif line == 'ABSORPTION SPECTRUM' and \
                 parse_version(self.metadata['package_version']).release > (2, 6):
                def energy_intensity(line):
                    """ CASSCF absorption spectrum
------------------------------------------------------------------------------------------
                                ABSORPTION SPECTRUM
------------------------------------------------------------------------------------------
  States           Energy   Wavelength   fosc          T2        TX         TY        TZ
                   (cm-1)     (nm)                   (D**2)      (D)        (D)       (D)
------------------------------------------------------------------------------------------
  0( 0)-> 1( 0) 1   83163.2    120.2   0.088250385   2.25340   0.00000   0.00000   1.50113
"""
                    reg = r'(\d+)\( ?(\d+)\)-> ?(\d+)\( ?(\d+)\) (\d+)'+ r'\s+(\d+\.\d+)'*4 + r'\s+(-?\d+\.\d+)'*3
                    res = re.search(reg, line)
                    jstate, jblock, istate, iblock, mult, energy, wavelength, intensity, t2, tx, ty, tz = res.groups()
                    return energy, intensity

            name = line
            self.skip_lines(inputfile, header)

            if not hasattr(self, 'transprop'):
                self.transprop = {}

            if energy_intensity is not None:
                etenergies = []
                etoscs = []
                line = next(inputfile)
                # The sections are occasionally ended with dashed lines
                # other times they are blank (other than a new line)
                while len(line.strip('-')) > 2:
                    energy, intensity = energy_intensity(line)
                    etenergies.append(float(energy))
                    etoscs.append(float(intensity))

                    line = next(inputfile)

                self.set_attribute('etenergies', etenergies)
                self.set_attribute('etoscs', etoscs)
                self.transprop[name] = (numpy.asarray(etenergies), numpy.asarray(etoscs))

        if line.strip() == "CD SPECTRUM":
            # -------------------------------------------------------------------
            #                              CD SPECTRUM
            # -------------------------------------------------------------------
            # State   Energy Wavelength       R         MX        MY        MZ   
            #         (cm-1)   (nm)       (1e40*cgs)   (au)      (au)      (au)  
            # -------------------------------------------------------------------
            #    1   43167.6    231.7      0.00000   0.00000  -0.00000   0.00000
            # ...
            #    6   25291.1    395.4 spin forbidden
            #
            # OR (from 4.2.0 onwards)
            #------------------------------------------------------------------------------
            #                             CD SPECTRUM
            #------------------------------------------------------------------------------
            #      States        Energy   Wavelength   R*T        RX        RY        RZ
            #                    (cm-1)      (nm)   (1e40*sgs)   (au)      (au)      (au)
            #------------------------------------------------------------------------------
            #  0( 1)-> 1( 1) 1   37192.8    268.9     0.00000  -0.00000  -0.34085   0.00000
            # ...
            #------------------------------------------------------------------------------
            etenergies = []
            etrotats = []
            self.skip_lines(inputfile, ["d", "State   Energy Wavelength", "(cm-1)   (nm)", "d"])
            line = next(inputfile)
            while line.strip() and not utils.str_contains_only(line.strip(), ['-']):
                tokens = line.split()
                if "spin forbidden" in line:
                    etrotat, mx, my, mz = 0.0, 0.0, 0.0, 0.0
                    etenergies.append(utils.float(tokens[-4]))
                else:
                    etrotat, mx, my, mz = [utils.float(t) for t in tokens[-4:]]
                    etenergies.append(utils.float(tokens[-6]))
                etrotats.append(etrotat)
                line = next(inputfile)
            self.set_attribute("etrotats", etrotats)
            if not hasattr(self, "etenergies"):
                self.logger.warning("etenergies not parsed before ECD section, "
                                    "the output file may be malformed")
                self.set_attribute("etenergies", etenergies)

        # ---------------
        # CHEMICAL SHIFTS
        # ---------------
        #
        # Note: using conversion factor for au to ppm alpha^2/2 =   26.625677252
        # GIAO: Doing para- and diamagnetic shielding integrals analytically     ...done
        #  --------------
        #  Nucleus   0C :
        #  --------------
        #
        # Diamagnetic contribution to the shielding tensor (ppm) :
        #            261.007         -0.294        -0.000
        #             -0.093        266.836         0.000
        #             -0.000         -0.000       244.502
        #
        # Paramagnetic contribution to the shielding tensor (ppm):
        #           -179.969          6.352        -0.000
        #              3.571       -210.368        -0.000
        #              0.000          0.000       -31.944
        #
        # Total shielding tensor (ppm):
        #             81.038          6.058        -0.000
        #              3.478         56.468        -0.000
        #              0.000          0.000       212.558
        #
        #
        #  Diagonalized sT*s matrix:
        #
        #  sDSO           266.692          261.151          244.502  iso=     257.448
        #  sPSO          -211.114         -179.223          -31.944  iso=    -140.760
        #         ---------------  ---------------  ---------------
        #  Total           55.577           81.929          212.558  iso=     116.688
        # ...
        #
        #--------------------------
        #CHEMICAL SHIELDING SUMMARY (ppm)
        #--------------------------
        #
        #
        #  Nucleus  Element    Isotropic     Anisotropy
        #  -------  -------  ------------   ------------
        #      0       C          116.686        143.809
        #      1       C          122.158        130.692
        # ...
        if line[:15] == 'CHEMICAL SHIFTS':
            nmrtensors = dict()
            while line.strip() != 'CHEMICAL SHIELDING SUMMARY (ppm)':
                if line[:8] == ' Nucleus':
                    atom = int(re.search(r'Nucleus\s+(\d+)\w', line).groups()[0])
                    self.skip_lines(inputfile, ['-', ''])
                    atomtensors = dict()
                    for _ in range(3):
                        t_type = next(inputfile).split()[0].lower()
                        tensor = numpy.zeros((3, 3))
                        for j, row in zip(range(3), inputfile):
                            tensor[j, :] = list(map(float, row.split()))
                        atomtensors[t_type] = tensor
                        self.skip_line(inputfile, '')
                    nmrtensors[atom] = atomtensors
                line = next(inputfile)

            self.skip_lines(inputfile, ['-', '', '', 'text', '-'])

            # Not currently used.
            isotropic, anisotropic = [], []
            for line in inputfile:
                if not line.strip():
                    break
                nucleus, element, iso, aniso = line.split()
                isotropic.append(float(iso))
                anisotropic.append(float(aniso))

            self.set_attribute('nmrtensors', nmrtensors)

        if line[:23] == "VIBRATIONAL FREQUENCIES":

            self.skip_lines(inputfile, ['d', 'b'])

            # Starting with 4.1, a scaling factor for frequencies is printed
            if float(self.metadata["package_version"][:3]) > 4.0:
                self.skip_lines(inputfile, ['Scaling factor for frequencies', 'b'])

            if self.natom > 1:
                vibfreqs = numpy.zeros(3 * self.natom)
                for i, line in zip(range(3 * self.natom), inputfile):
                    vibfreqs[i] = float(line.split()[1])

                nonzero = numpy.nonzero(vibfreqs)[0]
                self.first_mode = nonzero[0]
                # Take all modes after first
                # Mode between imaginary and real modes could be 0
                self.num_modes = 3*self.natom - self.first_mode
                if self.num_modes > 3*self.natom - 6:
                    msg = "Modes corresponding to rotations/translations may be non-zero."
                    if self.num_modes == 3*self.natom - 5:
                        msg += '\n You can ignore this if the molecule is linear.'
                self.set_attribute('vibfreqs', vibfreqs[self.first_mode:])
            else:
                # we have a single atom
                self.set_attribute('vibfreqs', numpy.array([]))

        # NORMAL MODES
        # ------------
        #
        # These modes are the cartesian displacements weighted by the diagonal matrix
        # M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom
        # Thus, these vectors are normalized but *not* orthogonal
        #
        #                   0          1          2          3          4          5
        #       0       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
        #       1       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
        #       2       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
        # ...
        if line[:12] == "NORMAL MODES":
            if self.natom > 1:
                all_vibdisps = numpy.zeros((3 * self.natom, self.natom, 3), "d")

                self.skip_lines(inputfile, ['d', 'b', 'text', 'text', 'text', 'b'])

                for mode in range(0, 3 * self.natom, 6):
                    header = next(inputfile)
                    for atom in range(self.natom):
                        all_vibdisps[mode:mode + 6, atom, 0] = next(inputfile).split()[1:]
                        all_vibdisps[mode:mode + 6, atom, 1] = next(inputfile).split()[1:]
                        all_vibdisps[mode:mode + 6, atom, 2] = next(inputfile).split()[1:]

                self.set_attribute('vibdisps', all_vibdisps[self.first_mode:])
            else:
                # we have a single atom
                self.set_attribute('vibdisps', numpy.array([]))

        # ORCA 4 example
        # -----------
        # IR SPECTRUM
        # -----------
        #
        #  Mode    freq (cm**-1)   T**2         TX         TY         TZ
        # -------------------------------------------------------------------
        #    6:      2069.36    1.674341  ( -0.000000   0.914970  -0.914970)
        #    7:      3978.81   76.856228  (  0.000000   6.199041  -6.199042)
        #    8:      4113.34   61.077784  ( -0.000000   5.526201   5.526200)
        # ...

        # ORCA 5 example
        # -----------
        # IR SPECTRUM
        # -----------
        #
        # Mode   freq       eps      Int      T**2         TX        TY        TZ
        #       cm**-1   L/(mol*cm) km/mol    a.u.
        # ----------------------------------------------------------------------------
        #  6:     45.66   0.000006    0.03  0.000039  ( 0.000000  0.000000  0.006256)
        #  7:     78.63   0.000000    0.00  0.000000  ( 0.000000  0.000000  0.000000)
        # ...
        if line[:11] == "IR SPECTRUM":
            package_version = self.metadata.get('package_version', None)
            if package_version is None:
                self.logger.warn('package_version has not been set, assuming 5.x.x')
                package_version = '5.x.x'
            major_version = int(package_version[0])
            if major_version <= 4:
                self.skip_lines(inputfile, ['d', 'b', 'header', 'd'])
                regex = r'\s+(?P<num>\d+):\s+(?P<frequency>\d+\.\d+)\s+(?P<intensity>\d+\.\d+)'
            else:
                self.skip_lines(inputfile, ['d', 'b', 'header', 'units', 'd'])
                regex = r'\s+(?P<num>\d+):\s+(?P<frequency>\d+\.\d+)\s+(?P<eps>\d+\.\d+)\s+(?P<intensity>\d+\.\d+)'

            if self.natom > 1:
                all_vibirs = numpy.zeros((3 * self.natom,), "d")

                line = next(inputfile)
                matches = re.match(regex, line)
                while matches:
                    num = int(matches.group('num'))
                    intensity = float(matches.group('intensity'))
                    all_vibirs[num] = intensity
                    line = next(inputfile)
                    matches = re.match(regex, line)

                self.set_attribute('vibirs', all_vibirs[self.first_mode:])
            else:
                # we have a single atom
                self.set_attribute('vibirs', numpy.array([]))

        # --------------
        # RAMAN SPECTRUM
        # --------------
        #
        #  Mode    freq (cm**-1)   Activity   Depolarization
        # -------------------------------------------------------------------
        #    6:       296.23      5.291229      0.399982
        #    7:       356.70      0.000000      0.749764
        #    8:       368.27      0.000000      0.202068
        if line[:14] == "RAMAN SPECTRUM":
            self.skip_lines(inputfile, ['d', 'b', 'header', 'd'])

            if self.natom > 1:
                all_vibramans = numpy.zeros(3 * self.natom)

                line = next(inputfile)
                while len(line) > 2:
                    num = int(line[0:4])
                    all_vibramans[num] = float(line.split()[2])
                    line = next(inputfile)

                self.set_attribute('vibramans', all_vibramans[self.first_mode:])
            else:
                # we have a single atom
                self.set_attribute('vibramans', numpy.array([]))

        # ORCA will print atomic charges along with the spin populations,
        #   so care must be taken about choosing the proper column.
        # Population analyses are performed usually only at the end
        #   of a geometry optimization or other run, so we want to
        #   leave just the final atom charges.
        # Here is an example for Mulliken charges:
        # --------------------------------------------
        # MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS
        # --------------------------------------------
        #    0 H :    0.126447    0.002622
        #    1 C :   -0.613018   -0.029484
        #    2 H :    0.189146    0.015452
        #    3 H :    0.320041    0.037434
        # ...
        # Sum of atomic charges         :   -0.0000000
        # Sum of atomic spin populations:    1.0000000
        if line[:23] == "MULLIKEN ATOMIC CHARGES":
            self.parse_charge_section(line, inputfile, 'mulliken')
        # Things are the same for Lowdin populations, except that the sums
        #   are not printed (there is a blank line at the end).
        if line[:22] == "LOEWDIN ATOMIC CHARGES":
            self.parse_charge_section(line, inputfile, 'lowdin')
        #CHELPG Charges
        #--------------------------------
        #  0   C   :       0.363939
        #  1   H   :       0.025695
        # ...
        #--------------------------------
        #Total charge:    -0.000000
        #--------------------------------
        if line.startswith('CHELPG Charges'):
            self.parse_charge_section(line, inputfile, 'chelpg')

        # It is not stated explicitely, but the dipole moment components printed by ORCA
        # seem to be in atomic units, so they will need to be converted. Also, they
        # are most probably calculated with respect to the origin .
        #
        # -------------
        # DIPOLE MOMENT
        # -------------
        #                                 X             Y             Z
        # Electronic contribution:      0.00000      -0.00000      -0.00000
        # Nuclear contribution   :      0.00000       0.00000       0.00000
        #                         -----------------------------------------
        # Total Dipole Moment    :      0.00000      -0.00000      -0.00000
        #                         -----------------------------------------
        # Magnitude (a.u.)       :      0.00000
        # Magnitude (Debye)      :      0.00000
        #
        if line.strip() == "DIPOLE MOMENT":

            self.skip_lines(inputfile, ['d', 'XYZ', 'electronic', 'nuclear', 'd'])
            total = next(inputfile)
            assert "Total Dipole Moment" in total

            reference = [0.0, 0.0, 0.0]
            dipole = numpy.array([float(d) for d in total.split()[-3:]])
            dipole = utils.convertor(dipole, "ebohr", "Debye")

            if not hasattr(self, 'moments'):
                self.set_attribute('moments', [reference, dipole])
            else:
                try:
                    assert numpy.all(self.moments[1] == dipole)
                except AssertionError:
                    self.logger.warning('Overwriting previous multipole moments with new values')
                    self.set_attribute('moments', [reference, dipole])

        if "Molecular Dynamics Iteration" in line:
            self.skip_lines(inputfile, ['d', 'ORCA MD', 'd', 'New Coordinates'])
            line = next(inputfile)
            tokens = line.split()
            assert tokens[0] == "time"
            time = utils.convertor(float(tokens[2]), "time_au", "fs")
            self.append_attribute('time', time)

        # Static polarizability.
        if line.strip() == "THE POLARIZABILITY TENSOR":
            if not hasattr(self, 'polarizabilities'):
                self.polarizabilities = []
            self.skip_lines(inputfile, ['d', 'b'])
            line = next(inputfile)
            assert line.strip() == "The raw cartesian tensor (atomic units):"
            polarizability = []
            for _ in range(3):
                line = next(inputfile)
                polarizability.append(line.split())
            self.polarizabilities.append(numpy.array(polarizability))

        if line.strip() == 'ORCA-CASSCF':
            # -------------------------------------------------------------------------------
            #                               ORCA-CASSCF
            # -------------------------------------------------------------------------------
            #
            # Symmetry handling      UseSym     ... ON
            # Point group                       ... C2
            # Used point group                  ... C2
            # Number of irreps                  ... 2
            #    Irrep    A has   10 SALCs (ofs=   0) #(closed)=   0 #(active)=   2
            #    Irrep    B has   10 SALCs (ofs=  10) #(closed)=   0 #(active)=   2
            #  Symmetries of active orbitals:
            #    MO =    0  IRREP= 0 (A)
            #    MO =    1  IRREP= 1 (B)
            self.skip_lines(inputfile, ['d', 'b'])
            vals = next(inputfile).split()
            # Symmetry section is only printed if symmetry is used.
            if vals[0] == 'Symmetry':
                assert vals[-1] == 'ON'
                point_group = next(inputfile).split()[-1]
                used_point_group = next(inputfile).split()[-1]
                num_irreps = int(next(inputfile).split()[-1])
                num_active = 0
                # Parse the irreps.
                for i, line in zip(range(num_irreps), inputfile):
                    reg = r'Irrep\s+(\w+) has\s+(\d+) SALCs \(ofs=\s*(\d+)\) #\(closed\)=\s*(\d+) #\(active\)=\s*(\d+)'
                    groups = re.search(reg, line).groups()
                    irrep = groups[0]
                    salcs, ofs, closed, active = map(int, groups[1:])
                    num_active += active
                self.skip_line(inputfile, 'Symmetries')
                # Parse the symmetries of the active orbitals.
                for i, line in zip(range(num_active), inputfile):
                    reg = r'(\d+)  IRREP= (\d+) \((\w+)\)'
                    groups = re.search(reg, line).groups()
                    mo, irrep_idx, irrep = groups

            # Skip until the system specific settings.
            # This will align the cases of symmetry on and off.
            line = next(inputfile)
            while line[:25] != 'SYSTEM-SPECIFIC SETTINGS:':
                line = next(inputfile)

            # SYSTEM-SPECIFIC SETTINGS:
            # Number of active electrons          ...    4
            # Number of active orbitals           ...    4
            # Total number of electrons           ...    4
            # Total number of orbitals            ...   20
            num_el = int(next(inputfile).split()[-1])
            num_orbs = int(next(inputfile).split()[-1])
            total_el = int(next(inputfile).split()[-1])
            total_orbs = int(next(inputfile).split()[-1])

            line = utils.skip_until_no_match(inputfile, r'^\s*$|^Total number aux.*$|^Determined.*$')

            # Determined orbital ranges:
            #    Internal       0 -   -1 (   0 orbitals)
            #    Active         0 -    3 (   4 orbitals)
            #    External       4 -   19 (  16 orbitals)
            orbital_ranges = []
            # Parse the orbital ranges for: Internal, Active, and External orbitals.
            for i in range(3):
                vals = line.split()
                start, end, num = int(vals[1]), int(vals[3]), int(vals[5])
                # Change from inclusive to exclusive in order to match python.
                end = end + 1
                assert end - start == num
                orbital_ranges.append((start, end, num))

            line = next(inputfile)
            while line[:8] != 'CI-STEP:':
                line = next(inputfile)

            # CI-STEP:
            # CI strategy                         ... General CI
            # Number of symmetry/multplity blocks ...    1
            # BLOCK  1 WEIGHT=   1.0000
            #   Multiplicity                      ...    1
            #   Irrep                             ...    0 (A)
            #   #(Configurations)                 ...   11
            #   #(CSFs)                           ...   12
            #   #(Roots)                          ...    1
            #     ROOT=0 WEIGHT=    1.000000
            self.skip_line(inputfile, 'CI strategy')
            num_blocks = int(next(inputfile).split()[-1])
            for b in range(1, num_blocks + 1):
                line = utils.skip_until_no_match(inputfile, r'^\s*$')
                vals = line.split()
                block = int(vals[1])
                weight = float(vals[3])
                assert b == block
                mult = int(next(inputfile).split()[-1])
                vals = next(inputfile).split()
                # The irrep will only be printed if using symmetry.
                if vals[0] == 'Irrep':
                    irrep_idx = int(vals[-2])
                    irrep = vals[-1].strip('()')
                    vals = next(inputfile).split()
                num_confs = int(vals[-1])
                num_csfs = int(next(inputfile).split()[-1])
                num_roots = int(next(inputfile).split()[-1])
                # Parse the roots.
                for r, line in zip(range(num_roots), inputfile):
                    reg = r'=(\d+) WEIGHT=\s*(\d\.\d+)'
                    groups = re.search(reg, line).groups()
                    root = int(groups[0])
                    weight = float(groups[1])
                    assert r == root

            # Skip additional setup printing and CASSCF iterations.
            line = next(inputfile).strip()
            while line != 'CASSCF RESULTS':
                line = next(inputfile).strip()

            # --------------
            # CASSCF RESULTS
            # --------------
            #
            # Final CASSCF energy       : -14.597120777 Eh    -397.2078 eV
            self.skip_lines(inputfile, ['d', 'b'])
            casscf_energy = float(next(inputfile).split()[4])

            # This is only printed for first and last step of geometry optimization.
            # ----------------
            # ORBITAL ENERGIES
            # ----------------
            #
            #   NO   OCC          E(Eh)            E(eV)    Irrep
            #    0   0.0868       0.257841         7.0162    1-A
            self.skip_lines(inputfile, ['b', 'd'])
            if next(inputfile).strip() == 'ORBITAL ENERGIES':
                self.skip_lines(inputfile, ['d', 'b', 'NO'])
                orbitals = []
                vals = next(inputfile).split()
                while vals:
                    occ, eh, ev = map(float, vals[1:4])
                    # The irrep will only be printed if using symmetry.
                    if len(vals) == 5:
                        idx, irrep = vals[4].split('-')
                        orbitals.append((occ, ev, int(idx), irrep))
                    else:
                        orbitals.append((occ, ev))
                    vals = next(inputfile).split()
                self.skip_lines(inputfile, ['b', 'd'])

            # Orbital Compositions
            # ---------------------------------------------
            # CAS-SCF STATES FOR BLOCK  1 MULT= 1 IRREP= Ag NROOTS= 2
            # ---------------------------------------------
            #
            # ROOT   0:  E=     -14.5950507665 Eh
            #       0.89724 [     0]: 2000
            for b in range(num_blocks):
                line = utils.skip_until_no_match(inputfile, r'^\s*$|^-*$')
                # Parse the block data.
                reg = r'BLOCK\s+(\d+) MULT=\s*(\d+) (IRREP=\s*\w+ )?(NROOTS=\s*(\d+))?'
                groups = re.search(reg, line).groups()
                block = int(groups[0])
                mult = int(groups[1])
                # The irrep will only be printed if using symmetry.
                if groups[2] is not None:
                    irrep = groups[2].split('=')[1].strip()
                nroots = int(groups[3].split('=')[1])

                self.skip_lines(inputfile, ['d', 'b'])

                line = next(inputfile).strip()
                while line:
                    if line[:4] == 'ROOT':
                        # Parse the root section.
                        reg = r'(\d+):\s*E=\s*(-?\d+.\d+) Eh(\s+\d+\.\d+ eV)?(\s+\d+\.\d+)?'
                        groups = re.search(reg, line).groups()
                        root = int(groups[0])
                        energy = float(groups[1])
                        # Excitation energies are only printed for excited state roots.
                        if groups[2] is not None:
                            excitation_energy_ev = float(groups[2].split()[0])
                            excitation_energy_cm = float(groups[3])
                    else:
                        # Parse the occupations section.
                        reg = r'(\d+\.\d+) \[\s*(\d+)\]: (\d+)'
                        groups = re.search(reg, line).groups()
                        coeff = float(groups[0])
                        number = float(groups[1])
                        occupations = list(map(int, groups[2]))

                    line = next(inputfile).strip()

            # Skip any extended wavefunction printing.
            while line != 'DENSITY MATRIX':
                line = next(inputfile).strip()

            self.skip_lines(inputfile, ['d', 'b'])
            # --------------
            # DENSITY MATRIX
            # --------------
            #
            #                    0          1          2          3
            #       0       0.897244   0.000000   0.000000   0.000000
            #       1       0.000000   0.533964   0.000000   0.000000
            density = numpy.zeros((num_orbs, num_orbs))
            for i in range(0, num_orbs, 6):
                next(inputfile)
                for j, line in zip(range(num_orbs), inputfile):
                    density[j][i:i + 6] = list(map(float, line.split()[1:]))

            line = utils.skip_until_no_match(inputfile, r'^\s*$|^-*$|^Trace.*$|^Extracting.*$')
            
            # This is only printed for open-shells.
            # -------------------
            # SPIN-DENSITY MATRIX
            # -------------------
            #
            #                   0          1          2          3          4          5
            #       0      -0.003709   0.001410   0.000074  -0.000564  -0.007978   0.000735
            #       1       0.001410  -0.001750  -0.000544  -0.003815   0.008462  -0.004529
            if line.strip() == 'SPIN-DENSITY MATRIX':
                self.skip_lines(inputfile, ['d', 'b'])
                spin_density = numpy.zeros((num_orbs, num_orbs))
                for i in range(0, num_orbs, 6):
                    next(inputfile)
                    for j, line in zip(range(num_orbs), inputfile):
                        spin_density[j][i:i + 6] = list(map(float, line.split()[1:]))
                self.skip_lines(inputfile, ['Trace', 'b', 'd', 'ENERGY'])
            self.skip_lines(inputfile, ['d', 'b'])

            # -----------------
            # ENERGY COMPONENTS
            # -----------------
            #
            # One electron energy          :    -18.811767801 Eh        -511.8942 eV
            # Two electron energy          :      4.367616074 Eh         118.8489 eV
            # Nuclear repuslion energy     :      0.000000000 Eh           0.0000 eV
            #                                ----------------
            #                                   -14.444151727
            #
            # Kinetic energy               :     14.371970266 Eh         391.0812 eV
            # Potential energy             :    -28.816121993 Eh        -784.1265 eV
            # Virial ratio                 :     -2.005022378
            #                                ----------------
            #                                   -14.444151727
            #
            # Core energy                  :    -13.604678408 Eh     -370.2021 eV
            one_el_energy = float(next(inputfile).split()[4])
            two_el_energy = float(next(inputfile).split()[4])
            nuclear_repulsion_energy = float(next(inputfile).split()[4])
            self.skip_line(inputfile, 'dashes')
            energy = float(next(inputfile).strip())
            self.skip_line(inputfile, 'blank')
            kinetic_energy = float(next(inputfile).split()[3])
            potential_energy = float(next(inputfile).split()[3])
            virial_ratio = float(next(inputfile).split()[3])
            self.skip_line(inputfile, 'dashes')
            energy = float(next(inputfile).strip())
            self.skip_line(inputfile, 'blank')
            core_energy = float(next(inputfile).split()[3])

        if line[:15] == 'TOTAL RUN TIME:':
            self.metadata['success'] = True

    def parse_charge_section(self, line, inputfile, chargestype):
        """Parse a charge section, modifies class in place

        Parameters
        ----------
        line : str
          the line which triggered entry here
        inputfile : file
          handle to file object
        chargestype : str
          what type of charge we're dealing with, must be one of
          'mulliken', 'lowdin' or 'chelpg'
        """
        has_spins = 'AND SPIN POPULATIONS' in line

        if not hasattr(self, "atomcharges"):
            self.atomcharges = {}
        if has_spins and not hasattr(self, "atomspins"):
            self.atomspins = {}

        self.skip_line(inputfile, 'dashes')

        # depending on chargestype, decide when to stop parsing lines
        # start, stop - indices for slicing lines and grabbing values
        if chargestype == 'mulliken':
            should_stop = lambda x: x.startswith('Sum of atomic charges')
            start, stop = 8, 20
        elif chargestype == 'lowdin':
            # stops when blank line encountered
            should_stop = lambda x: not bool(x.strip())
            start, stop = 8, 20
        elif chargestype == 'chelpg':
            should_stop = lambda x: x.startswith('---')
            start, stop = 11, 26

        charges = []
        if has_spins:
            spins = []

        line = next(inputfile)
        while not should_stop(line):
            # Don't add point charges or embedding potentials.
            if "Q :" not in line:
                charges.append(float(line[start:stop]))
                if has_spins:
                    spins.append(float(line[stop:]))
            line = next(inputfile)

        self.atomcharges[chargestype] = charges
        if has_spins:
            self.atomspins[chargestype] = spins

    def parse_scf_condensed_format(self, inputfile, line):
        """ Parse the SCF convergence information in condensed format """

        # This is what it looks like
        # ITER       Energy         Delta-E        Max-DP      RMS-DP      [F,P]     Damp
        #                ***  Starting incremental Fock matrix formation  ***
        #   0   -384.5203638934   0.000000000000 0.03375012  0.00223249  0.1351565 0.7000
        #   1   -384.5792776162  -0.058913722842 0.02841696  0.00175952  0.0734529 0.7000
        #                                ***Turning on DIIS***
        #   2   -384.6074211837  -0.028143567475 0.04968025  0.00326114  0.0310435 0.0000
        #   3   -384.6479682063  -0.040547022616 0.02097477  0.00121132  0.0361982 0.0000
        #   4   -384.6571124353  -0.009144228947 0.00576471  0.00035160  0.0061205 0.0000
        #   5   -384.6574659959  -0.000353560584 0.00191156  0.00010160  0.0025838 0.0000
        #   6   -384.6574990782  -0.000033082375 0.00052492  0.00003800  0.0002061 0.0000
        #   7   -384.6575005762  -0.000001497987 0.00020257  0.00001146  0.0001652 0.0000
        #   8   -384.6575007321  -0.000000155848 0.00008572  0.00000435  0.0000745 0.0000
        #          **** Energy Check signals convergence ****

        assert line[2] == "Delta-E"
        assert line[3] == "Max-DP"

        if not hasattr(self, "scfvalues"):
            self.scfvalues = []

        self.scfvalues.append([])

        # Try to keep track of the converger (NR, DIIS, SOSCF, etc.).
        diis_active = True
        while line:

            maxDP = None
            if 'Newton-Raphson' in line:
                diis_active = False
            elif 'SOSCF' in line:
                diis_active = False
            elif line[0].isdigit():
                shim = 0
                try:
                    energy = float(line[1])
                    deltaE = float(line[2])
                    maxDP = float(line[3 + int(not diis_active)])
                    rmsDP = float(line[4 + int(not diis_active)])
                except ValueError as e:
                    # Someone in Orca forgot to properly add spaces in the scf printing
                    # code looks like:
                    # %3i %17.10f%12.12f%11.8f %11.8f
                    if line[1].count('.') == 2:
                        integer1, decimal1_integer2, decimal2 = line[1].split('.')
                        decimal1, integer2 = decimal1_integer2[:10], decimal1_integer2[10:]
                        energy = float(integer1 + '.' + decimal1)
                        deltaE = float(integer2 + '.' + decimal2)
                        maxDP = float(line[2 + int(not diis_active)])
                        rmsDP = float(line[3 + int(not diis_active)])
                    elif line[1].count('.') == 3:
                        integer1, decimal1_integer2, decimal2_integer3, decimal3 = line[1].split('.')
                        decimal1, integer2 = decimal1_integer2[:10], decimal1_integer2[10:]
                        decimal2, integer3 = decimal2_integer3[:12], decimal2_integer3[12:]
                        energy = float(integer1 + '.' + decimal1)
                        deltaE = float(integer2 + '.' + decimal2)
                        maxDP = float(integer3 + '.' + decimal3)
                        rmsDP = float(line[2 + int(not diis_active)])
                    elif line[2].count('.') == 2:
                        integer1, decimal1_integer2, decimal2 = line[2].split('.')
                        decimal1, integer2 = decimal1_integer2[:12], decimal1_integer2[12:]
                        deltaE = float(integer1 + '.' + decimal1)
                        maxDP = float(integer2 + '.' + decimal2)
                        rmsDP = float(line[3 + int(not diis_active)])
                    else:
                        raise e

                self.scfvalues[-1].append([deltaE, maxDP, rmsDP])

            try:
                line = next(inputfile).split()
            except StopIteration:
                self.logger.warning(
                    f"File terminated before end of last SCF! Last Max-DP: {maxDP}"
                )
                break

    def parse_scf_expanded_format(self, inputfile, line):
        """ Parse SCF convergence when in expanded format. """

        # The following is an example of the format
        # -----------------------------------------
        #
        #               ***  Starting incremental Fock matrix formation  ***
        #
        #                         ----------------------------
        #                         !        ITERATION     0   !
        #                         ----------------------------
        #   Total Energy        :    -377.960836651297 Eh
        #   Energy Change       :    -377.960836651297 Eh
        #   MAX-DP              :       0.100175793695
        #   RMS-DP              :       0.004437973661
        #   Actual Damping      :       0.7000
        #   Actual Level Shift  :       0.2500 Eh
        #   Int. Num. El.       :    43.99982197 (UP=   21.99991099 DN=   21.99991099)
        #   Exchange            :   -34.27550826
        #   Correlation         :    -2.02540957
        #
        #
        #                         ----------------------------
        #                         !        ITERATION     1   !
        #                         ----------------------------
        #   Total Energy        :    -378.118458080109 Eh
        #   Energy Change       :      -0.157621428812 Eh
        #   MAX-DP              :       0.053240648588
        #   RMS-DP              :       0.002375092508
        #   Actual Damping      :       0.7000
        #   Actual Level Shift  :       0.2500 Eh
        #   Int. Num. El.       :    43.99994143 (UP=   21.99997071 DN=   21.99997071)
        #   Exchange            :   -34.00291075
        #   Correlation         :    -2.01607243
        #
        #                               ***Turning on DIIS***
        #
        #                         ----------------------------
        #                         !        ITERATION     2   !
        #                         ----------------------------
        # ....
        #
        if not hasattr(self, "scfvalues"):
            self.scfvalues = []

        self.scfvalues.append([])

        line = "Foo"  # dummy argument to enter loop
        while line.find("******") < 0:
            try:
                line = next(inputfile)
            except StopIteration:
                self.logger.warning('File terminated before end of last SCF!')
                break
            info = line.split()
            if len(info) > 1 and info[1] == "ITERATION":
                dashes = next(inputfile)
                energy_line = next(inputfile).split()
                energy = float(energy_line[3])
                deltaE_line = next(inputfile).split()
                deltaE = float(deltaE_line[3])
                if energy == deltaE:
                    deltaE = 0
                maxDP_line = next(inputfile).split()
                maxDP = float(maxDP_line[2])
                rmsDP_line = next(inputfile).split()
                rmsDP = float(rmsDP_line[2])
                self.scfvalues[-1].append([deltaE, maxDP, rmsDP])

        return

    # end of parse_scf_expanded_format

    def _append_scfvalues_scftargets(self, inputfile, line):
        # The SCF convergence targets are always printed after this, but apparently
        # not all of them always -- for example the RMS Density is missing for geometry
        # optimization steps. So, assume the previous value is still valid if it is
        # not found. For additional certainty, assert that the other targets are unchanged.
        while not "Last Energy change" in line:
            line = next(inputfile)
        deltaE_value = float(line.split()[4])
        deltaE_target = float(line.split()[7])
        line = next(inputfile)
        if "Last MAX-Density change" in line:
            maxDP_value = float(line.split()[4])
            maxDP_target = float(line.split()[7])
            line = next(inputfile)
            if "Last RMS-Density change" in line:
                rmsDP_value = float(line.split()[4])
                rmsDP_target = float(line.split()[7])
            else:
                rmsDP_value = self.scfvalues[-1][-1][2]
                rmsDP_target = self.scftargets[-1][2]
                assert deltaE_target == self.scftargets[-1][0]
                assert maxDP_target == self.scftargets[-1][1]
            self.scfvalues[-1].append([deltaE_value, maxDP_value, rmsDP_value])
            self.scftargets.append([deltaE_target, maxDP_target, rmsDP_target])
