# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Psi3 output files."""

import numpy

from cclib.parser import logfileparser
from cclib.parser import utils


class Psi3(logfileparser.Logfile):
    """A Psi3 log file."""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="Psi3", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return f"Psi3 log file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'Psi3("{self.filename}")'

    def normalisesym(self, label):
        """Psi3 does not require normalizing symmetry labels."""
        return label

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        if "Version" in line:
            self.metadata["package_version"] = ''.join(line.split()[1:]).lower()

        # Psi3 prints the coordinates in several configurations, and we will parse the
        # the canonical coordinates system in Angstroms as the first coordinate set,
        # although it is actually somewhere later in the input, after basis set, etc.
        # We can also get or verify the number of atoms and atomic numbers from this block.
        if line.strip() == "-Geometry in the canonical coordinate system (Angstrom):":

            self.skip_lines(inputfile, ['header', 'd'])

            coords = []
            numbers = []
            line = next(inputfile)
            while line.strip():

                tokens = line.split()

                element = tokens[0]
                numbers.append(self.table.number[element])

                x = float(tokens[1])
                y = float(tokens[2])
                z = float(tokens[3])
                coords.append([x, y, z])

                line = next(inputfile)

            self.set_attribute('natom', len(coords))
            self.set_attribute('atomnos', numbers)

            if not hasattr(self, 'atomcoords'):
                self.atomcoords = []
            self.atomcoords.append(coords)

        if line.strip() == '-SYMMETRY INFORMATION:':
            line = next(inputfile)
            assert "Computational point group is" in line
            point_group_abelian = line.split()[4].lower()
            # TODO Psi 3 doesn't print the full point group?
            point_group_full = point_group_abelian
            self.metadata['symmetry_detected'] = point_group_full
            self.metadata['symmetry_used'] = point_group_abelian
            while line.strip():
                if "Number of atoms" in line:
                    self.set_attribute('natom', int(line.split()[-1]))
                line = next(inputfile)
        if line.strip() == "-BASIS SET INFORMATION:":
            line = next(inputfile)
            while line.strip():
                if "Number of SO" in line:
                    self.set_attribute('nbasis', int(line.split()[-1]))
                line = next(inputfile)

        # In Psi3, the section with the contraction scheme can be used to infer atombasis.
        if line.strip() == "-Contraction Scheme:":

            self.skip_lines(inputfile, ['header', 'd'])

            indices = []
            line = next(inputfile)
            while line.strip():
                shells = line.split('//')[-1]
                expression = shells.strip().replace(' ', '+')
                expression = expression.replace('s', '*1')
                expression = expression.replace('p', '*3')
                expression = expression.replace('d', '*6')
                nfuncs = eval(expression)
                if len(indices) == 0:
                    indices.append(range(nfuncs))
                else:
                    start = indices[-1][-1] + 1
                    indices.append(range(start, start+nfuncs))
                line = next(inputfile)

            self.set_attribute('atombasis', indices)

        if line.strip() == "CINTS: An integrals program written in C":

            self.skip_lines(inputfile, ['authors', 'd', 'b', 'b'])

            line = next(inputfile)
            assert line.strip() == "-OPTIONS:"
            while line.strip():
                line = next(inputfile)

            line = next(inputfile)
            assert line.strip() == "-CALCULATION CONSTANTS:"
            while line.strip():
                if "Number of atoms" in line:
                    natom = int(line.split()[-1])
                    self.set_attribute('natom', natom)
                if "Number of symmetry orbitals" in line:
                    nbasis = int(line.split()[-1])
                    self.set_attribute('nbasis', nbasis)
                line = next(inputfile)

        if line.strip() == "CSCF3.0: An SCF program written in C":

            self.skip_lines(inputfile, ['b', 'authors', 'b', 'd', 'b',
                                        'mult', 'mult_comment', 'b'])

            line = next(inputfile)
            while line.strip():
                if line.split()[0] == "multiplicity":
                    mult = int(line.split()[-1])
                    self.set_attribute('mult', mult)
                if line.split()[0] == "charge":
                    charge = int(line.split()[-1])
                    self.set_attribute('charge', charge)
                if line.split()[0] == "convergence":
                    conv = float(line.split()[-1])
                if line.split()[0] == "reference":
                    self.reference = line.split()[-1]
                line = next(inputfile)

            if not hasattr(self, 'scftargets'):
                self.scftargets = []
            self.scftargets.append([conv])

        #  ==> Iterations <==

        # Psi3 converges just the density elements, although it reports in the iterations
        # changes in the energy as well as the DIIS error.
        psi3_iterations_header = "iter       total energy        delta E         delta P          diiser"
        if line.strip() == psi3_iterations_header:

            if not hasattr(self, 'scfvalues'):
                self.scfvalues = []
            self.scfvalues.append([])

            line = next(inputfile)
            while line.strip():
                ddensity = float(line.split()[-2])
                self.scfvalues[-1].append([ddensity])
                line = next(inputfile)

        # This section, from which we parse molecular orbital symmetries and
        # orbital energies, is quite similar for both Psi3 and Psi4, and in fact
        # the format for orbtials is the same, although the headers and spacers
        # are a bit different. Let's try to get both parsed with one code block.
        #
        # Here is how the block looks like for Psi4:
        #
        #	Orbital Energies (a.u.)
        #	-----------------------
        #
        #	Doubly Occupied:
        #
        #	   1Bu   -11.040586     1Ag   -11.040524     2Bu   -11.031589
        #	   2Ag   -11.031589     3Bu   -11.028950     3Ag   -11.028820
        # (...)
        #	  15Ag    -0.415620     1Bg    -0.376962     2Au    -0.315126
        #	   2Bg    -0.278361     3Bg    -0.222189
        #
        #	Virtual:
        #
        #	   3Au     0.198995     4Au     0.268517     4Bg     0.308826
        #	   5Au     0.397078     5Bg     0.521759    16Ag     0.565017
        # (...)
        #	  24Ag     0.990287    24Bu     1.027266    25Ag     1.107702
        #	  25Bu     1.124938
        #
        # The case is different in the trigger string.
        if "orbital energies (a.u.)" in line.lower():

            self.moenergies = [[]]
            self.mosyms = [[]]

            self.skip_line(inputfile, 'blank')

            occupied = next(inputfile)
            if self.reference[0:2] == 'RO' or self.reference[0:1] == 'R':
                assert 'doubly occupied' in occupied.lower()
            elif self.reference[0:1] == 'U':
                assert 'alpha occupied' in occupied.lower()

            # Parse the occupied MO symmetries and energies.
            self._parse_mosyms_moenergies(inputfile, 0)

            # The last orbital energy here represents the HOMO.
            self.homos = [len(self.moenergies[0])-1]
            # For a restricted open-shell calculation, this is the
            # beta HOMO, and we assume the singly-occupied orbitals
            # are all alpha, which are handled next.
            if self.reference[0:2] == 'RO':
                self.homos.append(self.homos[0])

            self.skip_line(inputfile, 'blank')

            unoccupied = next(inputfile)
            if self.reference[0:2] == 'RO':
                assert unoccupied.strip() == 'Singly Occupied:'
            elif self.reference[0:1] == 'R':
                assert unoccupied.strip() == 'Unoccupied orbitals'
            elif self.reference[0:1] == 'U':
                assert unoccupied.strip() == 'Alpha Virtual:'

            # Parse the unoccupied MO symmetries and energies.
            self._parse_mosyms_moenergies(inputfile, 0)

            # Here is where we handle the Beta or Singly occupied orbitals.
            if self.reference[0:1] == 'U':
                self.mosyms.append([])
                self.moenergies.append([])
                line = next(inputfile)
                assert line.strip() == 'Beta Occupied:'
                self.skip_line(inputfile, 'blank')
                self._parse_mosyms_moenergies(inputfile, 1)
                self.homos.append(len(self.moenergies[1])-1)
                line = next(inputfile)
                assert line.strip() == 'Beta Virtual:'
                self.skip_line(inputfile, 'blank')
                self._parse_mosyms_moenergies(inputfile, 1)
            elif self.reference[0:2] == 'RO':
                line = next(inputfile)
                assert line.strip() == 'Virtual:'
                self.skip_line(inputfile, 'blank')
                self._parse_mosyms_moenergies(inputfile, 0)

        # Both Psi3 and Psi4 print the final SCF energy right after
        # the orbital energies, but the label is different. Psi4 also
        # does DFT, and the label is also different in that case.
        if "* SCF total energy" in line:
            e = float(line.split()[-1])
            if not hasattr(self, 'scfenergies'):
                self.scfenergies = []
            self.scfenergies.append(utils.convertor(e, 'hartree', 'eV'))

        # We can also get some higher moments in Psi3, although here the dipole is not printed
        # separately and the order is not lexicographical. However, the numbers seem
        # kind of strange -- the quadrupole seems to be traceless, although I'm not sure
        # whether the standard transformation has been used. So, until we know what kind
        # of moment these are and how to make them raw again, we will only parse the dipole.
        #
        # --------------------------------------------------------------
        #                *** Electric multipole moments ***
        # --------------------------------------------------------------
        #
        #  CAUTION : The system has non-vanishing dipole moment, therefore
        #    quadrupole and higher moments depend on the reference point.
        #
        # -Coordinates of the reference point (a.u.) :
        #           x                     y                     z
        #  --------------------  --------------------  --------------------
        #          0.0000000000          0.0000000000          0.0000000000
        #
        # -Electric dipole moment (expectation values) :
        #
        #    mu(X)  =  -0.00000 D  =  -1.26132433e-43 C*m  =  -0.00000000 a.u.
        #    mu(Y)  =   0.00000 D  =   3.97987832e-44 C*m  =   0.00000000 a.u.
        #    mu(Z)  =   0.00000 D  =   0.00000000e+00 C*m  =   0.00000000 a.u.
        #    |mu|   =   0.00000 D  =   1.32262368e-43 C*m  =   0.00000000 a.u.
        #
        # -Components of electric quadrupole moment (expectation values) (a.u.) :
        #
        #     Q(XX) =   10.62340220   Q(YY) =    1.11816843   Q(ZZ) =  -11.74157063
        #     Q(XY) =    3.64633112   Q(XZ) =    0.00000000   Q(YZ) =    0.00000000
        #
        if line.strip() == "*** Electric multipole moments ***":

            self.skip_lines(inputfile, ['d', 'b', 'caution1', 'caution2', 'b'])

            coordinates = next(inputfile)
            assert coordinates.split()[-2] == "(a.u.)"
            self.skip_lines(inputfile, ['xyz', 'd'])
            line = next(inputfile)
            self.origin = numpy.array([float(x) for x in line.split()])
            self.origin = utils.convertor(self.origin, 'bohr', 'Angstrom')

            self.skip_line(inputfile, "blank")
            line = next(inputfile)
            assert "Electric dipole moment" in line
            self.skip_line(inputfile, "blank")

            # Make sure to use the column that has the value in Debyes.
            dipole = []
            for i in range(3):
                line = next(inputfile)
                dipole.append(float(line.split()[2]))

            if not hasattr(self, 'moments'):
                self.moments = [self.origin, dipole]
            else:
                assert self.moments[1] == dipole

    def _parse_mosyms_moenergies(self, inputfile, spinidx):
        """Parse molecular orbital symmetries and energies from the
        'Post-Iterations' section.
        """
        line = next(inputfile)
        while line.strip():
            for i in range(len(line.split()) // 2):
                self.mosyms[spinidx].append(line.split()[i*2][-2:])
                moenergy = utils.convertor(float(line.split()[i*2+1]), "hartree", "eV")
                self.moenergies[spinidx].append(moenergy)
            line = next(inputfile)
        return
