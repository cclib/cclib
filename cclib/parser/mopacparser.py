# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for MOPAC output files"""

# Based on parser in RMG-Py by Greg Magoon
# https://github.com/ReactionMechanismGenerator/RMG-Py/blob/master/external/cclib/parser/mopacparser.py
# Also parts from Ben Albrecht
# https://github.com/ben-albrecht/cclib/blob/master/cclib/parser/mopacparser.py
# Merged and modernized by Geoff Hutchison

import re
import math

import numpy

from cclib.parser import data
from cclib.parser import logfileparser
from cclib.parser import utils


def symbol2int(symbol):
    t = utils.PeriodicTable()
    return t.number[symbol]

class MOPAC(logfileparser.Logfile):
    """A MOPAC20XX output file."""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="MOPAC", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return f"MOPAC log file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'MOPAC("{self.filename}")'

    def normalisesym(self, label):
        """MOPAC does not require normalizing symmetry labels."""
        return label

    def before_parsing(self):
        #TODO

        # Defaults
        charge = 0
        self.set_attribute('charge', charge)
        mult = 1
        self.set_attribute('mult', mult)

        # Keep track of whether or not we're performing an
        # (un)restricted calculation.
        self.unrestricted = False
        self.is_rohf = False

        # Keep track of 1SCF vs. gopt since gopt is default
        self.onescf = False
        self.geomdone = False

        # Compile the dashes-and-or-spaces-only regex.
        self.re_dashes_and_spaces = re.compile(r'^[\s-]+$')

        self.star = ' * '
        self.stars = ' *******************************************************************************'

        self.spinstate = {'SINGLET': 1,
                          'DOUBLET': 2,
                          'TRIPLET': 3,
                          'QUARTET': 4,
                          'QUINTET': 5,
                          'SEXTET': 6,
                          'HEPTET': 7,
                          'OCTET': 8,
                          'NONET': 9}

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # Extract the package version.
        if "For non-commercial use only" in line:
            # Ignore the platorm information for now (the last character).
            self.metadata["package_version"] = line.split()[8][:-1]
            # Use the year as the legacy (short) package version.
            self.skip_lines(
                inputfile, ["Stewart Computational Chemistry", "s", "s", "s", "s"]
            )
            self.metadata["legacy_package_version"] = next(inputfile).split()[1][5:]

        # Extract the atomic numbers and coordinates from the optimized geometry
        # note that cartesian coordinates section occurs multiple times in the file, and we want to end up using the last instance
        # also, note that the section labeled cartesian coordinates doesn't have as many decimal places as the one used here
        # Example 1 (not used):
        #          CARTESIAN COORDINATES
        #
        #    NO.       ATOM               X         Y         Z
        #
        #     1         O                  4.7928   -0.8461    0.3641
        #     2         O                  5.8977   -0.3171    0.0092
        # ...
        # Example 2 (used):
        #   ATOM   CHEMICAL          X               Y               Z
        #  NUMBER    SYMBOL      (ANGSTROMS)     (ANGSTROMS)     (ANGSTROMS)
        #
        #     1       O          4.79280259  *  -0.84610232  *   0.36409474  *
        #     2       O          5.89768035  *  -0.31706418  *   0.00917035  *
        # ... etc.
        if line.split() == ["NUMBER", "SYMBOL", "(ANGSTROMS)", "(ANGSTROMS)", "(ANGSTROMS)"]:

            self.updateprogress(inputfile, "Attributes", self.cupdate)

            self.inputcoords = []
            self.inputatoms = []

            blankline = inputfile.next()

            atomcoords = []
            line = inputfile.next()
            while len(line.split()) > 6:
                # MOPAC Version 14.019L 64BITS suddenly appends this block with
                # "CARTESIAN COORDINATES" block with no blank line.
                tokens = line.split()
                self.inputatoms.append(symbol2int(tokens[1]))
                xc = float(tokens[2])
                yc = float(tokens[4])
                zc = float(tokens[6])
                atomcoords.append([xc, yc, zc])
                line = inputfile.next()

            self.inputcoords.append(atomcoords)

            if not hasattr(self, "natom"):
                self.atomnos = numpy.array(self.inputatoms, 'i')
                self.natom = len(self.atomnos)

        if 'CHARGE ON SYSTEM =' in line:
            charge = int(line.split()[5])
            self.set_attribute('charge', charge)

        if 'SPIN STATE DEFINED' in line:
            # find the multiplicity from the line token (SINGLET, DOUBLET, TRIPLET, etc)
            mult = self.spinstate[line.split()[1]]
            self.set_attribute('mult', mult)

        # Read energy (in kcal/mol, converted to eV)
        #
        # FINAL HEAT OF FORMATION =       -333.88606 KCAL =   -1396.97927 KJ
        if 'FINAL HEAT OF FORMATION =' in line:
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            self.scfenergies.append(utils.convertor(utils.float(line.split()[5]), "kcal/mol", "eV"))

        # Molecular mass parsing (units will be amu)
        #
        # MOLECULAR WEIGHT        ==        130.1890
        if line[0:35] == '          MOLECULAR WEIGHT        =':
            self.molmass = utils.float(line.split()[3])

        #rotational constants
        #Example:
        #          ROTATIONAL CONSTANTS IN CM(-1)
        #
        #          A =    0.01757641   B =    0.00739763   C =    0.00712013
        # could also read in moment of inertia, but this should just differ by a constant: rot cons= h/(8*Pi^2*I)
        # note that the last occurence of this in the thermochemistry section has reduced precision,
        # so we will want to use the 2nd to last instance
        if line[0:40] == '          ROTATIONAL CONSTANTS IN CM(-1)':
            blankline = inputfile.next()
            rotinfo = inputfile.next()
            if not hasattr(self, "rotcons"):
                self.rotcons = []
            broken = rotinfo.split()
            # leave the rotational constants in Hz
            a = float(broken[2])
            b = float(broken[5])
            c = float(broken[8])
            self.rotcons.append([a, b, c])

        # Start of the IR/Raman frequency section.
        # Example:
        # VIBRATION    1    1A       ATOM PAIR        ENERGY CONTRIBUTION    RADIAL
        # FREQ.        15.08        C 12 --  C 16           +7.9% (999.0%)     0.0%
        # T-DIPOLE    0.2028        C 16 --  H 34           +5.8% (999.0%)    28.0%
        # TRAVEL      0.0240        C 16 --  H 32           +5.6% (999.0%)    35.0%
        # RED. MASS   1.7712        O  1 --  O  4           +5.2% (999.0%)     0.4%
        # EFF. MASS7752.8338
        #
        # VIBRATION    2    2A       ATOM PAIR        ENERGY CONTRIBUTION    RADIAL
        # FREQ.        42.22        C 11 --  C 15           +9.0% (985.8%)     0.0%
        # T-DIPOLE    0.1675        C 15 --  H 31           +6.6% (843.6%)     3.3%
        # TRAVEL      0.0359        C 15 --  H 29           +6.0% (802.8%)    24.5%
        # RED. MASS   1.7417        C 13 --  C 17           +5.8% (792.7%)     0.0%
        # EFF. MASS1242.2114
        if line[1:10] == 'VIBRATION':
            self.updateprogress(inputfile, "Frequency Information", self.fupdate)

            # get the vib symmetry
            if len(line.split()) >= 3:
                sym = line.split()[2]
                if not hasattr(self, 'vibsyms'):
                    self.vibsyms = []
                self.vibsyms.append(sym)

            line = inputfile.next()
            if 'FREQ' in line:
                if not hasattr(self, 'vibfreqs'):
                    self.vibfreqs = []
                freq = float(line.split()[1])
                self.vibfreqs.append(freq)

            line = inputfile.next()
            if 'T-DIPOLE' in line:
                if not hasattr(self, 'vibirs'):
                    self.vibirs = []
                tdipole = float(line.split()[1])
                # transform to km/mol
                self.vibirs.append(math.sqrt(tdipole))

            line = inputfile.next()
            if 'TRAVEL' in line:
                pass

            line = inputfile.next()
            if 'RED. MASS' in line:
                if not hasattr(self, 'vibrmasses'):
                    self.vibrmasses = []
                rmass = float(line.split()[2])
                self.vibrmasses.append(rmass)

        # Orbital eigenvalues, e.g.
        #           ALPHA EIGENVALUES
        #            BETA EIGENVALUES
        # or just "EIGENVALUES" for closed-shell
        if 'EIGENVALUES' in line:
            if not hasattr(self, 'moenergies'):
                self.moenergies = [] # list of arrays

            energies = []
            line = inputfile.next()
            while len(line.split()) > 0:
                energies.extend([float(i) for i in line.split()])
                line = inputfile.next()
            self.moenergies.append(energies)

        # todo:
        # Partial charges and dipole moments
        # Example:
        # NET ATOMIC CHARGES

        if line[:16] == '== MOPAC DONE ==':
            self.metadata['success'] = True
