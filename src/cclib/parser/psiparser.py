# -*- coding: utf-8 -*-
#
# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2008-2016, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Parser for Psi3 and Psi4 output files"""


import re

import numpy

from . import data
from . import logfileparser
from . import utils


class Psi(logfileparser.Logfile):
    """A Psi log file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(Psi, self).__init__(logname="Psi", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "Psi log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Psi("%s")' % (self.filename)

    def before_parsing(self):

        # There are some major differences between the output of Psi3 and Psi4,
        # so it will be useful to register which one we are dealing with.
        self.version = None

        # Early beta versions of Psi4 normalize basis function
        # coefficients when printing.
        self.version_4_beta = False

        # This is just used to track which part of the output we are in for Psi4,
        # with changes triggered by ==> things like this <== (Psi3 does not have this)
        self.section = None

        # Keep track of whether or not we're performing an
        # (un)restricted calculation. Can be ('RHF', 'UHF', 'ROHF',
        # 'CUHF').
        self.reference = 'RHF'

    def after_parsing(self):

        # Newer versions of Psi4 don't explicitly print the number of atoms.
        if not hasattr(self, 'natom'):
            if hasattr(self, 'atomnos'):
                self.set_attribute('natom', len(self.atomnos))

    def normalisesym(self, label):
        """Use standard symmetry labels instead of Psi labels."""
        # Psi uses the correct labels.
        return label

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # The version should always be detected.
        if "PSI3: An Open-Source Ab Initio" in line:
            self.version = 3
        if "PSI4: An Open-Source Ab Initio".lower() in line.lower():
            self.version = 4
            # A more detailed version (minor and patch level) appears
            # on the next line.
            line = next(inputfile)
            # Keep track of early versions of Psi4.
            if "beta" in line:
                self.version_4_beta = True

        # This will automatically change the section attribute for Psi4, when encountering
        # a line that <== looks like this ==>, to whatever is in between.
        if (line.strip()[:3] == "==>") and (line.strip()[-3:] == "<=="):
            self.section = line.strip()[4:-4]

        # Determine whether or not the reference wavefunction is
        # restricted, unrestricted, or restricted open-shell.
        if line.strip() == "SCF":
            self.skip_line(inputfile, 'author list')
            line = next(inputfile)
            self.reference = line.split()[0]
            # Work with a complex reference as if it's real.
            if self.reference[0] == 'C':
                self.reference = self.reference[1:]

        # Psi3 prints the coordinates in several configurations, and we will parse the
        # the canonical coordinates system in Angstroms as the first coordinate set,
        # although it is actually somewhere later in the input, after basis set, etc.
        # We can also get or verify the number of atoms and atomic numbers from this block.
        if (self.version == 3) and (line.strip() == "-Geometry in the canonical coordinate system (Angstrom):"):

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

        #  ==> Geometry <==
        #
        #    Molecular point group: c2h
        #    Full point group: C2h
        #
        #    Geometry (in Angstrom), charge = 0, multiplicity = 1:
        #
        #       Center              X                  Y                   Z
        #    ------------   -----------------  -----------------  -----------------
        #           C         -1.415253322400     0.230221785400     0.000000000000
        #           C          1.415253322400    -0.230221785400     0.000000000000
        # ...
        #
        if (self.section == "Geometry") and ("Geometry (in Angstrom), charge" in line):

            assert line.split()[3] == "charge"
            charge = int(line.split()[5].strip(','))
            self.set_attribute('charge', charge)

            assert line.split()[6] == "multiplicity"
            mult = int(line.split()[8].strip(':'))
            self.set_attribute('mult', mult)

            self.skip_line(inputfile, "blank")
            line = next(inputfile)

            # Usually there is the header and dashes, but, for example, the coordinates
            # printed when a geometry optimization finishes do not have it.
            if line.split()[0] == "Center":
                self.skip_line(inputfile, "dashes")
                line = next(inputfile)

            elements = []
            coords = []
            atommasses = []
            while line.strip():
                chomp = line.split()
                el, x, y, z = chomp[:4]
                if len(el) > 1:
                    el = el[0] + el[1:].lower()
                elements.append(el)
                coords.append([float(x), float(y), float(z)])
                # Newer versions of Psi4 print atomic masses.
                if len(chomp) == 5:
                    atommasses.append(float(chomp[4]))
                line = next(inputfile)

            # The 0 is to handle the presence of ghost atoms.
            self.set_attribute('atomnos', [self.table.number.get(el, 0) for el in elements])

            if not hasattr(self, 'atomcoords'):
                self.atomcoords = []

            # This condition discards any repeated coordinates that Psi print. For example,
            # geometry optimizations will print the coordinates at the beginning of and SCF
            # section and also at the start of the gradient calculation.
            if len(self.atomcoords) == 0 or self.atomcoords[-1] != coords:
                self.atomcoords.append(coords)

            if len(atommasses) > 0:
                if not hasattr(self, 'atommasses'):
                    self.atommasses = atommasses

        # In Psi3 there are these two helpful sections.
        if (self.version == 3) and (line.strip() == '-SYMMETRY INFORMATION:'):
            line = next(inputfile)
            while line.strip():
                if "Number of atoms" in line:
                    self.set_attribute('natom', int(line.split()[-1]))
                line = next(inputfile)
        if (self.version == 3) and (line.strip() == "-BASIS SET INFORMATION:"):
            line = next(inputfile)
            while line.strip():
                if "Number of SO" in line:
                    self.set_attribute('nbasis', int(line.split()[-1]))
                line = next(inputfile)

        # Psi4 repeats the charge and multiplicity after the geometry.
        if (self.section == "Geometry") and (line[2:16].lower() == "charge       ="):
            charge = int(line.split()[-1])
            self.set_attribute('charge', charge)
        if (self.section == "Geometry") and (line[2:16].lower() == "multiplicity ="):
            mult = int(line.split()[-1])
            self.set_attribute('mult', mult)

        # In Psi3, the section with the contraction scheme can be used to infer atombasis.
        if (self.version == 3) and line.strip() == "-Contraction Scheme:":

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

        # In Psi3, the integrals program prints useful information when invoked.
        if (self.version == 3) and (line.strip() == "CINTS: An integrals program written in C"):

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

        # In Psi3, this part contains alot of important data pertaining to the SCF, but not only:
        if (self.version == 3) and (line.strip() == "CSCF3.0: An SCF program written in C"):

            self.skip_lines(inputfile, ['b', 'authors', 'b', 'd', 'b', 'mult', 'mult_comment', 'b'])

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
                line = next(inputfile)

            if not hasattr(self, 'scftargets'):
                self.scftargets = []
            self.scftargets.append([conv])

        # The printout for Psi4 has a more obvious trigger for the SCF parameter printout.
        if (self.section == "Algorithm") and (line.strip() == "==> Algorithm <=="):

            self.skip_line(inputfile, 'blank')

            line = next(inputfile)
            while line.strip():
                if "Energy threshold" in line:
                    etarget = float(line.split()[-1])
                if "Density threshold" in line:
                    dtarget = float(line.split()[-1])
                line = next(inputfile)

            if not hasattr(self, "scftargets"):
                self.scftargets = []
            self.scftargets.append([etarget, dtarget])

        # This section prints contraction information before the atomic basis set functions and
        # is a good place to parse atombasis indices as well as atomnos. However, the section this line
        # is in differs between HF and DFT outputs.
        #
        #  -Contraction Scheme:
        #    Atom   Type   All Primitives // Shells:
        #   ------ ------ --------------------------
        #       1     C     6s 3p // 2s 1p
        #       2     C     6s 3p // 2s 1p
        #       3     C     6s 3p // 2s 1p
        # ...
        if (self.section == "Primary Basis" or self.section == "DFT Potential") and line.strip() == "-Contraction Scheme:":

            self.skip_lines(inputfile, ['headers', 'd'])

            atomnos = []
            atombasis = []
            atombasis_pos = 0
            line = next(inputfile)
            while line.strip():

                element = line.split()[1]
                if len(element) > 1:
                    element = element[0] + element[1:].lower()
                atomnos.append(self.table.number[element])

                # To count the number of atomic orbitals for the atom, sum up the orbitals
                # in each type of shell, times the numbers of shells. Currently, we assume
                # the multiplier is a single digit and that there are only s and p shells,
                # which will need to be extended later when considering larger basis sets,
                # with corrections for the cartesian/spherical cases.
                ao_count = 0
                shells = line.split('//')[1].split()
                for s in shells:
                    count, type = s
                    multiplier = 3*(type == 'p') or 1
                    ao_count += multiplier*int(count)

                if len(atombasis) > 0:
                    atombasis_pos = atombasis[-1][-1] + 1
                atombasis.append(list(range(atombasis_pos, atombasis_pos+ao_count)))

                line = next(inputfile)

            self.set_attribute('natom', len(atomnos))
            self.set_attribute('atomnos', atomnos)
            self.set_attribute('atombasis', atombasis)

        # The atomic basis set is straightforward to parse, but there are some complications
        # when symmetry is used, because in that case Psi4 only print the symmetry-unique atoms,
        # and the list of symmetry-equivalent ones is not printed. Therefore, for simplicity here
        # when an atomic is missing (atom indices are printed) assume the atomic orbitals of the
        # last atom of the same element before it. This might not work if a mixture of basis sets
        # is used somehow... but it should cover almost all cases for now.
        #
        # Note that Psi also print normalized coefficients (details below).
        #
        #  ==> AO Basis Functions <==
        #
        #    [ STO-3G ]
        #    spherical
        #    ****
        #    C   1
        #    S   3 1.00
        #                        71.61683700           2.70781445
        #                        13.04509600           2.61888016
        # ...
        if (self.section == "AO Basis Functions") and (line.strip() == "==> AO Basis Functions <=="):

            def get_symmetry_atom_basis(gbasis):
                """Get symmetry atom by replicating the last atom in gbasis of the same element."""

                missing_index = len(gbasis)
                missing_atomno = self.atomnos[missing_index]

                ngbasis = len(gbasis)
                last_same = ngbasis - self.atomnos[:ngbasis][::-1].index(missing_atomno) - 1
                return gbasis[last_same]

            dfact = lambda n: (n <= 0) or n * dfact(n-2)

            # Early beta versions of Psi4 normalize basis function
            # coefficients when printing.
            if self.version_4_beta:
                def get_normalization_factor(exp, lx, ly, lz):
                    norm_s = (2*exp/numpy.pi)**0.75
                    if lx + ly + lz > 0:
                        nom = (4*exp)**((lx+ly+lz)/2.0)
                        den = numpy.sqrt(dfact(2*lx-1) * dfact(2*ly-1) * dfact(2*lz-1))
                        return norm_s * nom / den
                    else:
                        return norm_s
            else:
                get_normalization_factor = lambda exp, lx, ly, lz: 1

            self.skip_lines(inputfile, ['b', 'basisname'])

            line = next(inputfile)
            spherical = line.strip() == "spherical"
            if hasattr(self, 'spherical_basis'):
                assert self.spherical_basis == spherical
            else:
                self.spherical_basis = spherical

            gbasis = []
            self.skip_line(inputfile, 'stars')
            line = next(inputfile)
            while line.strip():

                element, index = line.split()
                if len(element) > 1:
                    element = element[0] + element[1:].lower()
                atomno = self.table.number[element]
                index = int(index)

                # This is the code that adds missing atoms when symmetry atoms are excluded
                # from the basis set printout. Again, this will work only if all atoms of
                # the same element use the same basis set.
                while index > len(gbasis) + 1:
                    gbasis.append(get_symmetry_atom_basis(gbasis))

                gbasis.append([])
                line = next(inputfile)
                while line.find("*") == -1:

                    # The shell type and primitive count is in the first line.
                    shell_type, nprimitives, smthg = line.split()
                    nprimitives = int(nprimitives)

                    # Get the angular momentum for this shell type.
                    momentum = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4, 'H': 5, 'I': 6}[shell_type.upper()]

                    # Read in the primitives.
                    primitives_lines = [next(inputfile) for i in range(nprimitives)]
                    primitives = [list(map(float, pl.split())) for pl in primitives_lines]

                    # Un-normalize the coefficients. Psi prints the normalized coefficient
                    # of the highest polynomial, namely XX for D orbitals, XXX for F, and so on.
                    for iprim, prim in enumerate(primitives):
                        exp, coef = prim
                        coef = coef / get_normalization_factor(exp, momentum, 0, 0)
                        primitives[iprim] = [exp, coef]

                    primitives = [tuple(p) for p in primitives]
                    shell = [shell_type, primitives]
                    gbasis[-1].append(shell)

                    line = next(inputfile)

                line = next(inputfile)

            # We will also need to add symmetry atoms that are missing from the input
            # at the end of this block, if the symmetry atoms are last.
            while len(gbasis) < self.natom:
                gbasis.append(get_symmetry_atom_basis(gbasis))

            self.gbasis = gbasis

        # A block called 'Calculation Information' prints these before starting the SCF.
        if (self.section == "Pre-Iterations") and ("Number of atoms" in line):
            natom = int(line.split()[-1])
            self.set_attribute('natom', natom)
        if (self.section == "Pre-Iterations") and ("Number of atomic orbitals" in line):
            nbasis = int(line.split()[-1])
            self.set_attribute('nbasis', nbasis)
        if (self.section == "Pre-Iterations") and ("Total" in line):
            chomp = line.split()
            nbasis = int(chomp[1])
            self.set_attribute('nbasis', nbasis)

        #  ==> Iterations <==

        # Psi3 converges just the density elements, although it reports in the iterations
        # changes in the energy as well as the DIIS error.
        psi3_iterations_header = "iter       total energy        delta E         delta P          diiser"
        if (self.version == 3) and (line.strip() == psi3_iterations_header):

            if not hasattr(self, 'scfvalues'):
                self.scfvalues = []
            self.scfvalues.append([])

            line = next(inputfile)
            while line.strip():
                ddensity = float(line.split()[-2])
                self.scfvalues[-1].append([ddensity])
                line = next(inputfile)

        # Psi4 converges both the SCF energy and density elements and reports both in the
        # iterations printout. However, the default convergence scheme involves a density-fitted
        # algorithm for efficiency, and this is often followed by a something with exact electron
        # repulsion integrals. In that case, there are actually two convergence cycles performed,
        # one for the density-fitted algorithm and one for the exact one, and the iterations are
        # printed in two blocks separated by some set-up information.
        if (self.section == "Iterations") and (line.strip() == "==> Iterations <=="):

            if not hasattr(self, 'scfvalues'):
                self.scfvalues = []
            scfvals = []

            self.skip_lines(inputfile, ['b', 'header', 'b'])
            line = next(inputfile)
            # Read each SCF iteration.
            while line.strip() != "==> Post-Iterations <==":
                if line.strip() and line.split()[0][0] == '@':
                    denergy = float(line.split()[4])
                    ddensity = float(line.split()[5])
                    scfvals.append([denergy, ddensity])
                try:
                    line = next(inputfile)
                except StopIteration:
                    self.logger.warning('File terminated before end of last SCF! Last density err: {}'.format(ddensity))
                    break
            self.section = "Post-Iterations"
            self.scfvalues.append(scfvals)

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

            # If this is Psi4, we will be in the appropriate section.
            assert (self.version == 3) or (self.section == "Post-Iterations")

            self.moenergies = [[]]
            self.mosyms = [[]]

            # Psi4 has dashes under the trigger line, but Psi3 did not.
            if self.version == 4:
                self.skip_line(inputfile, 'dashes')
            self.skip_line(inputfile, 'blank')

            # Both versions have this case-insensitive substring.
            occupied = next(inputfile)
            if self.reference[0:2] == 'RO' or self.reference[0:1] == 'R':
                assert 'doubly occupied' in occupied.lower()
            elif self.reference[0:1] == 'U':
                assert 'alpha occupied' in occupied.lower()

            # Psi4 now has a blank line, Psi3 does not.
            if self.version == 4:
                self.skip_line(inputfile, 'blank')

            # Parse the occupied MO symmetries and energies.
            self._parse_mosyms_moenergies(inputfile, 0)

            # The last orbital energy here represents the HOMO.
            self.homos = [len(self.moenergies[0])-1]
            # For a restricted open-shell calculation, this is the
            # beta HOMO, and we assume the singly-occupied orbitals
            # are all alpha, which are handled next.
            if self.reference[0:2] == 'RO':
                self.homos.append(self.homos[0])

            # Different numbers of blank lines in Psi3 and Psi4.
            if self.version == 3:
                self.skip_line(inputfile, 'blank')

            unoccupied = next(inputfile)
            if self.reference[0:2] == 'RO':
                assert unoccupied.strip() == 'Singly Occupied:'
            elif self.reference[0:1] == 'R':
                # The header for virtual orbitals is different for
                # Psi3 and Psi4.
                if self.version == 3:
                    assert unoccupied.strip() == 'Unoccupied orbitals'
                else:
                    assert unoccupied.strip() == 'Virtual:'
            elif self.reference[0:1] == 'U':
                assert unoccupied.strip() == 'Alpha Virtual:'

            # Psi4 now has a blank line, Psi3 does not.
            if self.version == 4:
                self.skip_line(inputfile, 'blank')

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

            if self.version == 4:
                line = next(inputfile)
                assert line.strip() == 'Final Occupation by Irrep:'
                line = next(inputfile)
                irreps = line.split()
                line = next(inputfile)
                tokens = line.split()
                assert tokens[0] == 'DOCC'
                docc = sum([int(x.replace(',', '')) for x in tokens[2:-1]])
                line = next(inputfile)
                if line.strip():
                    tokens = line.split()
                    assert tokens[0] in ('SOCC', 'NA')
                    socc = sum([int(x.replace(',', '')) for x in tokens[2:-1]])
                    # Fix up the restricted open-shell alpha HOMO.
                    if self.reference[0:2] == 'RO':
                        self.homos[0] += socc

        # Both Psi3 and Psi4 print the final SCF energy right after the orbital energies,
        # but the label is different. Psi4 also does DFT, and the label is also different in that case.
        if (self.version == 3 and "* SCF total energy" in line) or \
           (self.section == "Post-Iterations" and "Final Energy:" in line):
            if self.version == 3:
                e = float(line.split()[-1])
            else:
                e = float(line.split()[3])
            if not hasattr(self, 'scfenergies'):
                self.scfenergies = []
            self.scfenergies.append(utils.convertor(e, 'hartree', 'eV'))

        #  ==> Molecular Orbitals <==
        #
        #                 1            2            3            4            5
        #
        #    1    0.7014827    0.7015412    0.0096801    0.0100168    0.0016438
        #    2    0.0252630    0.0251793   -0.0037890   -0.0037346    0.0016447
        # ...
        #   59    0.0000133   -0.0000067    0.0000005   -0.0047455   -0.0047455
        #   60    0.0000133    0.0000067    0.0000005    0.0047455   -0.0047455
        #
        # Ene   -11.0288198  -11.0286067  -11.0285837  -11.0174766  -11.0174764
        # Sym            Ag           Bu           Ag           Bu           Ag
        # Occ             2            2            2            2            2
        #
        #
        #                11           12           13           14           15
        #
        #    1    0.1066946    0.1012709    0.0029709    0.0120562    0.1002765
        #    2   -0.2753689   -0.2708037   -0.0102079   -0.0329973   -0.2790813
        # ...
        #
        if (self.section) and ("Molecular Orbitals" in self.section) \
           and ("Molecular Orbitals" in line):

            self.skip_line(inputfile, 'blank')

            mocoeffs = []
            indices = next(inputfile)
            while indices.strip():

                if indices[:3] == '***':
                    break

                indices = [int(i) for i in indices.split()]

                if len(mocoeffs) < indices[-1]:
                    for i in range(len(indices)):
                        mocoeffs.append([])
                else:
                    assert len(mocoeffs) == indices[-1]

                self.skip_line(inputfile, 'blank')

                line = next(inputfile)
                while line.strip():
                    iao = int(line.split()[0])
                    coeffs = [float(c) for c in line.split()[1:]]
                    for i, c in enumerate(coeffs):
                        mocoeffs[indices[i]-1].append(c)
                    line = next(inputfile)

                energies = next(inputfile)
                symmetries = next(inputfile)
                occupancies = next(inputfile)

                self.skip_lines(inputfile, ['b', 'b'])
                indices = next(inputfile)

            if not hasattr(self, 'mocoeffs'):
                self.mocoeffs = []
            self.mocoeffs.append(mocoeffs)

        # The formats for Mulliken and Lowdin atomic charges are the same, just with
        # the name changes, so use the same code for both.
        #
        # Properties computed using the SCF density density matrix
        #   Mulliken Charges: (a.u.)
        #    Center  Symbol    Alpha    Beta     Spin     Total
        #        1     C     2.99909  2.99909  0.00000  0.00182
        #        2     C     2.99909  2.99909  0.00000  0.00182
        # ...
        for pop_type in ["Mulliken", "Lowdin"]:
            if line.strip() == "%s Charges: (a.u.)" % pop_type:
                if not hasattr(self, 'atomcharges'):
                    self.atomcharges = {}
                header = next(inputfile)

                line = next(inputfile)
                while not line.strip():
                    line = next(inputfile)

                charges = []
                while line.strip():
                    ch = float(line.split()[-1])
                    charges.append(ch)
                    line = next(inputfile)
                self.atomcharges[pop_type.lower()] = charges

        # This is for the older conventional MP2 code in 4.0b5.
        mp_trigger = "MP2 Total Energy (a.u.)"
        if line.strip()[:len(mp_trigger)] == mp_trigger:
            mpenergy = utils.convertor(float(line.split()[-1]), 'hartree', 'eV')
            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []
            self.mpenergies.append([mpenergy])
        # This is for the newer DF-MP2 code in 4.0.
        if 'DF-MP2 Energies' in line:
            while 'Total Energy' not in line:
                line = next(inputfile)
            mpenergy = utils.convertor(float(line.split()[3]), 'hartree', 'eV')
            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []
            self.mpenergies.append([mpenergy])

        # Note this is just a start and needs to be modified for CCSD(T), etc.
        ccsd_trigger = "* CCSD total energy"
        if line.strip()[:len(ccsd_trigger)] == ccsd_trigger:
            ccsd_energy = utils.convertor(float(line.split()[-1]), 'hartree', 'eV')
            if not hasattr(self, "ccenergis"):
                self.ccenergies = []
            self.ccenergies.append(ccsd_energy)

        # The geometry convergence targets and values are printed in a table, with the legends
        # describing the convergence annotation. Probably exact slicing of the line needs
        # to be done in order to extract the numbers correctly. If there are no values for
        # a paritcular target it means they are not used (marked also with an 'o'), and in this case
        # we will set a value of numpy.inf so that any value will be smaller.
        #
        #  ==> Convergence Check <==
        #
        #  Measures of convergence in internal coordinates in au.
        #  Criteria marked as inactive (o), active & met (*), and active & unmet ( ).
        #  ---------------------------------------------------------------------------------------------
        #   Step     Total Energy     Delta E     MAX Force     RMS Force      MAX Disp      RMS Disp
        #  ---------------------------------------------------------------------------------------------
        #    Convergence Criteria    1.00e-06 *    3.00e-04 *             o    1.20e-03 *             o
        #  ---------------------------------------------------------------------------------------------
        #      2    -379.77675264   -7.79e-03      1.88e-02      4.37e-03 o    2.29e-02      6.76e-03 o  ~
        #  ---------------------------------------------------------------------------------------------
        #
        if (self.section == "Convergence Check") and line.strip() == "==> Convergence Check <==":

            self.skip_lines(inputfile, ['b', 'units', 'comment', 'dash+tilde', 'header', 'dash+tilde'])

            # These are the position in the line at which numbers should start.
            starts = [27, 41, 55, 69, 83]

            criteria = next(inputfile)
            geotargets = []
            for istart in starts:
                if criteria[istart:istart+9].strip():
                    geotargets.append(float(criteria[istart:istart+9]))
                else:
                    geotargets.append(numpy.inf)

            self.skip_line(inputfile, 'dashes')

            values = next(inputfile)
            geovalues = []
            for istart in starts:
                if values[istart:istart+9].strip():
                    geovalues.append(float(values[istart:istart+9]))

            if not hasattr(self, "optstatus"):
                self.optstatus = []
            self.optstatus.append(data.ccData.OPT_NEW)

            # This assertion may be too restrictive, but we haven't seen the geotargets change.
            # If such an example comes up, update the value since we're interested in the last ones.
            if not hasattr(self, 'geotargets'):
                self.geotargets = geotargets
            else:
                assert self.geotargets == geotargets

            if not hasattr(self, 'geovalues'):
                self.geovalues = []
            self.geovalues.append(geovalues)

        # This message signals a converged optimization, in which case we want
        # to append the index for this step to optdone, which should be equal
        # to the number of geovalues gathered so far.
        if "Optimization is complete!" in line:

            # This is a workaround for Psi4.0/sample_opt-irc-2.out;
            # IRC calculations currently aren't parsed properly for
            # optimization parameters.
            if hasattr(self, 'geovalues'):

                if not hasattr(self, 'optdone'):
                    self.optdone = []
                self.optdone.append(len(self.geovalues))

                assert hasattr(self, "optstatus") and len(self.optstatus) > 0
                self.optstatus[-1] = data.ccData.OPT_DONE

        # This message means that optimization has stopped for some reason, but we
        # still want optdone to exist in this case, although it will be an empty list.
        if line.strip() == "Optimizer: Did not converge!":
            if not hasattr(self, 'optdone'):
                self.optdone = []

            assert hasattr(self, "optstatus") and len(self.optstatus) > 0
            self.optstatus[-1] = data.ccData.OPT_UNCONVERGED

        # The reference point at which properties are evaluated in Psi4 is explicitely stated,
        # so we can save it for later. It is not, however, a part of the Properties section,
        # but it appears before it and also in other places where properies that might depend
        # on it are printed.
        #
        # Properties will be evaluated at   0.000000,   0.000000,   0.000000 Bohr
        #
        if (self.version == 4) and ("Properties will be evaluated at" in line.strip()):
            self.origin = numpy.array([float(x.strip(',')) for x in line.split()[-4:-1]])
            assert line.split()[-1] == "Bohr"
            self.origin = utils.convertor(self.origin, 'bohr', 'Angstrom')

        # The properties section print the molecular dipole moment:
        #
        #  ==> Properties <==
        #
        #
        #Properties computed using the SCF density density matrix
        #  Nuclear Dipole Moment: (a.u.)
        #     X:     0.0000      Y:     0.0000      Z:     0.0000
        #
        #  Electronic Dipole Moment: (a.u.)
        #     X:     0.0000      Y:     0.0000      Z:     0.0000
        #
        #  Dipole Moment: (a.u.)
        #     X:     0.0000      Y:     0.0000      Z:     0.0000     Total:     0.0000
        #
        if (self.section == "Properties") and line.strip() == "Dipole Moment: (a.u.)":

            line = next(inputfile)
            dipole = numpy.array([float(line.split()[1]), float(line.split()[3]), float(line.split()[5])])
            dipole = utils.convertor(dipole, "ebohr", "Debye")

            if not hasattr(self, 'moments'):
                # Old versions of Psi4 don't print the origin; assume
                # it's at zero.
                if not hasattr(self, 'origin'):
                    self.origin = numpy.array([0.0, 0.0, 0.0])
                self.moments = [self.origin, dipole]
            else:
                try:
                    assert numpy.all(self.moments[1] == dipole)
                except AssertionError:
                    self.logger.warning('Overwriting previous multipole moments with new values')
                    self.logger.warning('This could be from post-HF properties or geometry optimization')
                    self.moments = [self.origin, dipole]

        # Higher multipole moments are printed separately, on demand, in lexicographical order.
        #
        # Multipole Moments:
        #
        # ------------------------------------------------------------------------------------
        #     Multipole             Electric (a.u.)       Nuclear  (a.u.)        Total (a.u.)
        # ------------------------------------------------------------------------------------
        #
        # L = 1.  Multiply by 2.5417462300 to convert to Debye
        # Dipole X            :          0.0000000            0.0000000            0.0000000
        # Dipole Y            :          0.0000000            0.0000000            0.0000000
        # Dipole Z            :          0.0000000            0.0000000            0.0000000
        #
        # L = 2.  Multiply by 1.3450341749 to convert to Debye.ang
        # Quadrupole XX       :      -1535.8888701         1496.8839996          -39.0048704
        # Quadrupole XY       :        -11.5262958           11.4580038           -0.0682920
        # ...
        #
        if line.strip() == "Multipole Moments:":

            self.skip_lines(inputfile, ['b', 'd', 'header', 'd', 'b'])

            # The reference used here should have been printed somewhere
            # before the properties and parsed above.
            moments = [self.origin]

            line = next(inputfile)
            while "----------" not in line.strip():

                rank = int(line.split()[2].strip('.'))

                multipole = []
                line = next(inputfile)
                while line.strip():

                    value = float(line.split()[-1])
                    fromunits = "ebohr" + (rank > 1)*("%i" % rank)
                    tounits = "Debye" + (rank > 1)*".ang" + (rank > 2)*("%i" % (rank-1))
                    value = utils.convertor(value, fromunits, tounits)
                    multipole.append(value)

                    line = next(inputfile)

                multipole = numpy.array(multipole)
                moments.append(multipole)
                line = next(inputfile)

            if not hasattr(self, 'moments'):
                self.moments = moments
            else:
                for im, m in enumerate(moments):
                    if len(self.moments) <= im:
                        self.moments.append(m)
                    else:
                        assert numpy.allclose(self.moments[im], m, atol=1.0e4)

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
        if (self.version == 3) and line.strip() == "*** Electric multipole moments ***":

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

        ## Harmonic frequencies.

        # -------------------------------------------------------------

        #   Computing second-derivative from gradients using projected,
        #   symmetry-adapted, cartesian coordinates (fd_freq_1).

        #   74 gradients passed in, including the reference geometry.
        #   Generating complete list of displacements from unique ones.

        #       Operation 2 takes plus displacements of irrep Bg to minus ones.
        #       Operation 3 takes plus displacements of irrep Au to minus ones.
        #       Operation 2 takes plus displacements of irrep Bu to minus ones.

        #         Irrep      Harmonic Frequency
        #                         (cm-1)
        #       -----------------------------------------------
        #            Au          137.2883
        if line.strip() == 'Irrep      Harmonic Frequency':

            vibsyms = []
            vibfreqs = []

            self.skip_lines(inputfile, ['(cm-1)', 'dashes'])

            ## The first section contains the symmetry of each normal
            ## mode and its frequency.
            line = next(inputfile)
            while '---' not in line:
                chomp = line.split()
                vibsym = chomp[0]
                vibfreq = Psi.parse_vibfreq(chomp[1])
                vibsyms.append(vibsym)
                vibfreqs.append(vibfreq)
                line = next(inputfile)

            self.set_attribute('vibsyms', vibsyms)
            self.set_attribute('vibfreqs', vibfreqs)

            line = next(inputfile)
            assert line.strip() == ''
            line = next(inputfile)
            assert 'Normal Modes' in line
            line = next(inputfile)
            assert 'Molecular mass is' in line
            if hasattr(self, 'atommasses'):
                assert abs(float(line.split()[3]) - sum(self.atommasses)) < 1.0e-4
            line = next(inputfile)
            assert line.strip() == 'Frequencies in cm^-1; force constants in au.'
            line = next(inputfile)
            assert line.strip() == ''
            line = next(inputfile)

            ## The second section contains the frequency, force
            ## constant, and displacement for each normal mode, along
            ## with the atomic masses.

            #       Normal Modes (non-mass-weighted).
            #       Molecular mass is  130.07825 amu.
            #       Frequencies in cm^-1; force constants in au.

            #  Frequency:        137.29
            #  Force constant:   0.0007
            #            X       Y       Z           mass
            # C      0.000   0.000   0.050      12.000000
            # C      0.000   0.000   0.050      12.000000
            for vibfreq in self.vibfreqs:
                _vibfreq = Psi.parse_vibfreq(line[13:].strip())
                assert abs(vibfreq - _vibfreq) < 1.0e-2
                line = next(inputfile)
                # Can't do anything with this for now.
                assert 'Force constant:' in line
                line = next(inputfile)
                assert 'X       Y       Z           mass' in line
                line = next(inputfile)
                if not hasattr(self, 'vibdisps'):
                    self.vibdisps = []
                normal_mode_disps = []
                # for k in range(self.natom):
                while line.strip():
                    chomp = line.split()
                    # Do nothing with this for now.
                    atomsym = chomp[0]
                    atomcoords = [float(x) for x in chomp[1:4]]
                    # Do nothing with this for now.
                    atommass = float(chomp[4])
                    normal_mode_disps.append(atomcoords)
                    line = next(inputfile)
                self.vibdisps.append(normal_mode_disps)
                line = next(inputfile)

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

    @staticmethod
    def parse_vibfreq(vibfreq):
        """Imaginary frequencies are printed as '12.34i', rather than
        '-12.34'.
        """
        is_imag = vibfreq[-1] == 'i'
        if is_imag:
            return -float(vibfreq[:-1])
        else:
            return float(vibfreq)


if __name__ == "__main__":
    import doctest, psiparser
    doctest.testmod(psiparser, verbose=False)
