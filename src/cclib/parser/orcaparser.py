# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for ORCA output files"""


from __future__ import print_function

import numpy

from cclib.parser import logfileparser
from cclib.parser import utils


class ORCA(logfileparser.Logfile):
    """An ORCA log file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(ORCA, self).__init__(logname="ORCA", *args, **kwargs)


    def __str__(self):
        """Return a string representation of the object."""
        return "ORCA log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'ORCA("%s")' % (self.filename)

    def normalisesym(self, label):
        """ORCA does not require normalizing symmetry labels."""
        return label

    def before_parsing(self):

        # A geometry optimization is started only when
        # we parse a cycle (so it will be larger than zero().
        self.gopt_cycle = 0

        # Keep track of whether this is a relaxed scan calculation
        self.is_relaxed_scan = False

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        #  extract the version number first
        if "Program Version" in line:
            self.metadata["package_version"] = line.split()[2]

        if line[0:15] == "Number of atoms":

            natom = int(line.split()[-1])
            self.set_attribute('natom', natom)

        if line[1:13] == "Total Charge":

            charge = int(line.split()[-1])
            self.set_attribute('charge', charge)

            line = next(inputfile)

            mult = int(line.split()[-1])
            self.set_attribute('mult', mult)

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
        if line[25:50] == "RELAXED SURFACE SCAN STEP":

            self.is_relaxed_scan = True
            blank = next(inputfile)
            info = next(inputfile)
            stars = next(inputfile)
            blank = next(inputfile)

            line = next(inputfile)
            while line[0:23] != "Convergence Tolerances:":
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

            if not hasattr(self, "geovalues"):
                self.geovalues = []

            headers = next(inputfile)
            dashes = next(inputfile)

            names = []
            values = []
            targets = []
            line = next(inputfile)
            while list(set(line.strip())) != ["."]:
                name = line[10:28].strip().lower()
                value = float(line.split()[2])
                target = float(line.split()[3])
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

            self.geovalues.append(newvalues)

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
                atomnos.append(self.table.number[atom])
                atomcoords.append([float(x), float(y), float(z)])
                line = next(inputfile)

            self.set_attribute('natom', len(atomnos))
            self.set_attribute('atomnos', atomnos)
            if not hasattr(self, 'atomcoords'):
                self.atomcoords = []
            self.atomcoords.append(atomcoords)

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
                if line == '* core charge reduced due to ECP\n':
                    break
                no, lb, za, frag, mass, x, y, z = line.split()
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

            line = next(inputfile)
            while len(line) > 20:  # restricted calcs are terminated by ------
                info = line.split()
                mooccno = int(float(info[1]))
                moenergy = float(info[2])
                self.mooccnos[0].append(mooccno)
                self.moenergies[0].append(utils.convertor(moenergy, "hartree", "eV"))
                line = next(inputfile)

            line = next(inputfile)

            # handle beta orbitals for UHF
            if line[17:35] == "SPIN DOWN ORBITALS":
                text = next(inputfile)

                self.mooccnos.append([])
                self.moenergies.append([])

                line = next(inputfile)
                while len(line) > 20:  # actually terminated by ------
                    info = line.split()
                    mooccno = int(float(info[1]))
                    moenergy = float(info[2])
                    self.mooccnos[1].append(mooccno)
                    self.moenergies[1].append(utils.convertor(moenergy, "hartree", "eV"))
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

        # Molecular orbital coefficients.
        # This is also where atombasis is parsed.
        if line[0:18] == "MOLECULAR ORBITALS":

            self.skip_line(inputfile, 'dashes')

            mocoeffs = [numpy.zeros((self.nbasis, self.nbasis), "d")]
            self.aonames = []
            self.atombasis = []
            for n in range(self.natom):
                self.atombasis.append([])

            for spin in range(len(self.moenergies)):

                if spin == 1:
                    self.skip_line(inputfile, 'blank')
                    mocoeffs.append(numpy.zeros((self.nbasis, self.nbasis), "d"))

                for i in range(0, self.nbasis, 6):

                    self.updateprogress(inputfile, "Coefficients")

                    self.skip_lines(inputfile, ['numbers', 'energies', 'occs'])

                    dashes = next(inputfile)
                    broken = dashes.split()
                    size = len(broken)

                    for j in range(self.nbasis):
                        line = next(inputfile)
                        broken = line.split()

                        #only need this on the first time through
                        if spin == 0 and i == 0:
                            atomname = line[3:5].split()[0]
                            num = int(line[0:3])
                            orbital = broken[1].upper()

                            self.aonames.append("%s%i_%s" % (atomname, num+1, orbital))
                            self.atombasis[num].append(j)

                        temp = []
                        vals = line[16:-1]  # -1 to remove the last blank space
                        for k in range(0, len(vals), 10):
                            temp.append(float(vals[k:k+10]))
                        mocoeffs[spin][i:i+size, j] = temp

            self.mocoeffs = mocoeffs

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

        """ Banner announcing Thermochemistry
        --------------------------
        THERMOCHEMISTRY AT 298.15K
        --------------------------
        """
        if 'THERMOCHEMISTRY AT' == line[:18]:

            next(inputfile)
            next(inputfile)
            self.temperature = float(next(inputfile).split()[2])
            self.pressure = float(next(inputfile).split()[2])
            total_mass = float(next(inputfile).split()[3])

            # Vibrations, rotations, and translations
            line = next(inputfile)
            while line[:17] != 'Electronic energy':
                line = next(inputfile)
            self.zpe = next(inputfile).split()[4]
            thermal_vibrational_correction = float(next(inputfile).split()[4])
            thermal_rotional_correction = float(next(inputfile).split()[4])
            thermal_translational_correction = float(next(inputfile).split()[4])
            next(inputfile)
            total_thermal_energy = float(next(inputfile).split()[3])

            # Enthalpy
            line = next(inputfile)
            while line[:17] != 'Total free energy':
                line = next(inputfile)
            thermal_enthalpy_correction = float(next(inputfile).split()[4])
            next(inputfile)
            self.enthalpy = float(next(inputfile).split()[3])

            # Entropy
            line = next(inputfile)
            while line[:18] != 'Electronic entropy':
                line = next(inputfile)
            electronic_entropy = float(line.split()[3])
            vibrational_entropy = float(next(inputfile).split()[3])
            rotational_entropy = float(next(inputfile).split()[3])
            translational_entropy = float(next(inputfile).split()[3])
            next(inputfile)
            self.entropy = float(next(inputfile).split()[4])

            line = next(inputfile)
            while line[:25] != 'Final Gibbs free enthalpy':
                line = next(inputfile)
            self.freeenergy = float(line.split()[5])

        # Read TDDFT information
        if any(x in line for x in ("TD-DFT/TDA EXCITED", "TD-DFT EXCITED")):
            # Could be singlets or triplets
            if line.find("SINGLETS") >= 0:
                sym = "Singlet"
            elif line.find("TRIPLETS") >= 0:
                sym = "Triplet"
            else:
                sym = "Not specified"

            if not hasattr(self, "etenergies"):
                self.etsecs = []
                self.etenergies = []
                self.etsyms = []

            lookup = {'a': 0, 'b': 1}
            line = next(inputfile)
            while line.find("STATE") < 0:
                line = next(inputfile)
            # Contains STATE or is blank
            while line.find("STATE") >= 0:
                broken = line.split()
                self.etenergies.append(float(broken[-2]))
                self.etsyms.append(sym)
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
                self.etsecs.append(sec)
                line = next(inputfile)


        # Parse the various absorption spectra for TDDFT and ROCIS
        if 'ABSORPTION SPECTRUM' in line or 'ELECTRIC DIPOLE' in line:
            line = line.strip()

            # Standard header, occasionally changes
            header = ['d', 'header', 'header', 'd']

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
            if line == 'COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM' or \
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
                    """ ROCIS dipole approximation with SOC == True (electric or veloctiy dipole moments)
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

            elif line[:79] == 'ROCIS COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM':
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

            name = line
            self.skip_lines(inputfile, header)

            if not hasattr(self, 'transprop'):
                self.transprop = {}

            etenergies = []
            etoscs = []
            line = next(inputfile)
            while len(line) > 1:
                energy, intensity = energy_intensity(line)
                etenergies.append(float(energy))
                etoscs.append(float(intensity))

                line = next(inputfile)

            self.etenergies = numpy.array(etenergies)
            self.etoscs = numpy.array(etoscs)
            self.transprop[name] = (self.etenergies, self.etoscs)


        if line[0:23] == "VIBRATIONAL FREQUENCIES":

            self.skip_lines(inputfile, ['d', 'b'])

            self.vibfreqs = numpy.zeros((3 * self.natom,), "d")

            for i in range(3 * self.natom):
                line = next(inputfile)
                self.vibfreqs[i] = float(line.split()[1])

            if numpy.any(self.vibfreqs[0:6] != 0):
                msg = "Modes corresponding to rotations/translations "
                msg += "may be non-zero."
                self.logger.warning(msg)

            self.vibfreqs = self.vibfreqs[6:]

        if line[0:12] == "NORMAL MODES":
            """ Format:
            NORMAL MODES
            ------------

            These modes are the cartesian displacements weighted by the diagonal matrix
            M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom
            Thus, these vectors are normalized but *not* orthogonal

                              0          1          2          3          4          5
                  0       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
                  1       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
                  2       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
            ...
            """

            self.vibdisps = numpy.zeros((3 * self.natom, self.natom, 3), "d")

            self.skip_lines(inputfile, ['d', 'b', 'text', 'text', 'text', 'b'])

            for mode in range(0, 3 * self.natom, 6):
                header = next(inputfile)
                for atom in range(self.natom):
                    x = next(inputfile).split()[1:]
                    y = next(inputfile).split()[1:]
                    z = next(inputfile).split()[1:]

                    self.vibdisps[mode:mode + 6, atom, 0] = x
                    self.vibdisps[mode:mode + 6, atom, 1] = y
                    self.vibdisps[mode:mode + 6, atom, 2] = z

            self.vibdisps = self.vibdisps[6:]

        if line[0:11] == "IR SPECTRUM":

            self.skip_lines(inputfile, ['d', 'b', 'header', 'd'])

            self.vibirs = numpy.zeros((3 * self.natom,), "d")

            line = next(inputfile)
            while len(line) > 2:
                num = int(line[0:4])
                self.vibirs[num] = float(line.split()[2])
                line = next(inputfile)

            self.vibirs = self.vibirs[6:]

        if line[0:14] == "RAMAN SPECTRUM":

            self.skip_lines(inputfile, ['d', 'b', 'header', 'd'])

            self.vibramans = numpy.zeros((3 * self.natom,), "d")

            line = next(inputfile)
            while len(line) > 2:
                num = int(line[0:4])
                self.vibramans[num] = float(line.split()[2])
                line = next(inputfile)

            self.vibramans = self.vibramans[6:]

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
                self.moments = [reference, dipole]
            else:
                try:
                    assert numpy.all(self.moments[1] == dipole)
                except AssertionError:
                    self.logger.warning('Overwriting previous multipole moments with new values')
                    self.moments = [reference, dipole]

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
                self.logger.warning('File terminated before end of last SCF! Last Max-DP: {}'.format(maxDP))
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
