# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Psi4 output files."""

from collections import namedtuple

import numpy

from cclib.parser import data
from cclib.parser import logfileparser
from cclib.parser import utils


class Psi4(logfileparser.Logfile):
    """A Psi4 log file."""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="Psi4", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return f"Psi4 log file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'Psi4("{self.filename}")'

    def before_parsing(self):

        # Early beta versions of Psi4 normalize basis function
        # coefficients when printing.
        self.version_4_beta = False

        # This is just used to track which part of the output we are in, with
        # changes triggered by ==> things like this <==.
        self.section = None

        # There are also sometimes subsections within each section, denoted
        # with =>/<= rather than ==>/<==.
        self.subsection = None

    def after_parsing(self):
        super(Psi4, self).after_parsing()

        # Newer versions of Psi4 don't explicitly print the number of atoms.
        if not hasattr(self, 'natom'):
            if hasattr(self, 'atomnos'):
                self.set_attribute('natom', len(self.atomnos))

    def normalisesym(self, label):
        """Use standard symmetry labels instead of Psi4 labels.

        To normalise:
        (1) `App` -> `A"`
        (2) `Ap` -> `A'`
        """
        return label.replace("pp", '"').replace("p", "'")

    # Match the number of skipped lines required based on the type of
    # gradient present (determined from the header), as otherwise the
    # parsing is identical.
    GradientInfo = namedtuple('GradientInfo', ['gradient_type', 'header', 'skip_lines'])
    GRADIENT_TYPES = {
        'analytic': GradientInfo('analytic',
                                 '-Total Gradient:',
                                 ['header', 'dash header']),
        'numerical': GradientInfo('numerical',
                                  '## F-D gradient (Symmetry 0) ##',
                                  ['Irrep num and total size', 'b', '123', 'b']),
    }
    GRADIENT_HEADERS = set([gradient_type.header
                            for gradient_type in GRADIENT_TYPES.values()])

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # Extract the version number and the version control
        # information, if it exists.
        if "An Open-Source Ab Initio Electronic Structure Package" in line:
            version_line = next(inputfile)
            tokens = version_line.split()
            package_version = tokens[1].split("-")[-1]
            self.metadata["legacy_package_version"] = package_version
            # Keep track of early versions of Psi4.
            if "beta" in package_version:
                self.version_4_beta = True
                # `beta2+` -> `0!0.beta2`
                package_version = f"0!0.{package_version}"
                if package_version[-1] == "+":
                    # There is no good way to keep the bare plus sign around,
                    # but this version is so old...
                    package_version = package_version[:-1]
            else:
                package_version = f"1!{package_version}"
            self.skip_line(inputfile, "blank")
            line = next(inputfile)
            if "Git:" in line:
                tokens = line.split()
                assert tokens[1] == "Rev"
                revision = '-'.join(tokens[2:]).replace("{", "").replace("}", "")
                dev_flag = "" if "dev" in package_version else ".dev"
                package_version = f"{package_version}{dev_flag}+{revision}"
            self.metadata["package_version"] = package_version

        # This will automatically change the section attribute for Psi4, when encountering
        # a line that <== looks like this ==>, to whatever is in between.
        if (line.strip()[:3] == "==>") and (line.strip()[-3:] == "<=="):
            self.section = line.strip()[4:-4]
            if self.section == "DFT Potential":
                self.metadata["methods"].append("DFT")

        # There is also the possibility of subsections.
        if (line.strip()[:2] == "=>") and (line.strip()[-2:] == "<="):
            self.subsection = line.strip()[3:-3]

        # Determine whether or not the reference wavefunction is
        # restricted, unrestricted, or restricted open-shell.
        if line.strip() == "SCF":
            while "Reference" not in line:
                line = next(inputfile)
            self.reference = line.split()[0]
            # Work with a complex reference as if it's real.
            if self.reference[0] == 'C':
                self.reference = self.reference[1:]

        # Parse the XC density functional
        # => Composite Functional: B3LYP <=
        if self.section == "DFT Potential" and "composite functional" in line.lower():
            chomp = line.split()
            functional = chomp[-2]
            self.metadata["functional"] = functional

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
        if (self.section == "Geometry") and ("Molecular point group" in line):

            point_group_abelian = line.split()[3].lower()
            line = next(inputfile)
            if "Full point group" in line:
                point_group_full = line.split()[3].lower()
            else:
                # TODO this isn't right, need to "know" about symmetry.
                point_group_full = point_group_abelian

            self.metadata['symmetry_detected'] = point_group_full
            self.metadata['symmetry_used'] = point_group_abelian

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
            if len(self.atomcoords) == 0 \
                or (self.atomcoords[-1] != coords and not hasattr(self, 'finite_difference')):
                self.atomcoords.append(coords)

            if len(atommasses) > 0:
                if not hasattr(self, 'atommasses'):
                    self.atommasses = atommasses

        # Psi4 repeats the charge and multiplicity after the geometry.
        if (self.section == "Geometry") and (line[2:16].lower() == "charge       ="):
            charge = int(line.split()[-1])
            self.set_attribute('charge', charge)
        if (self.section == "Geometry") and (line[2:16].lower() == "multiplicity ="):
            mult = int(line.split()[-1])
            self.set_attribute('mult', mult)

        # The printout for Psi4 has a more obvious trigger for the SCF parameter printout.
        if (self.section == "Algorithm") and (line.strip() == "==> Algorithm <==") \
            and not hasattr(self, 'finite_difference'):

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
        if self.section == "Primary Basis":
            if line[2:12] == "Basis Set:":
                self.metadata["basis_set"] = line.split()[2]

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
                    shell_type, nprimitives, _ = line.split()
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

        # Psi4 converges both the SCF energy and density elements and reports both in the
        # iterations printout. However, the default convergence scheme involves a density-fitted
        # algorithm for efficiency, and this is often followed by a something with exact electron
        # repulsion integrals. In that case, there are actually two convergence cycles performed,
        # one for the density-fitted algorithm and one for the exact one, and the iterations are
        # printed in two blocks separated by some set-up information.
        if (self.section == "Iterations") and (line.strip() == "==> Iterations <==") \
            and not hasattr(self, 'finite_difference'):

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
                    self.logger.warning(
                        f"File terminated before end of last SCF! Last density err: {ddensity}"
                    )
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
        if ("orbital energies (a.u.)" in line.lower()  or "orbital energies [eh]" in line.lower()) \
            and not hasattr(self, 'finite_difference'):

            # If this is Psi4, we will be in the appropriate section.
            assert self.section == "Post-Iterations"

            self.moenergies = [[]]
            self.mosyms = [[]]

            # Psi4 has dashes under the trigger line, but Psi3 did not.
            self.skip_line(inputfile, 'dashes')
            self.skip_line(inputfile, 'blank')

            # Both versions have this case-insensitive substring.
            occupied = next(inputfile)
            if self.reference[0:2] == 'RO' or self.reference[0:1] == 'R':
                assert 'doubly occupied' in occupied.lower()
            elif self.reference[0:1] == 'U':
                assert 'alpha occupied' in occupied.lower()

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

            unoccupied = next(inputfile)
            if self.reference[0:2] == 'RO':
                assert unoccupied.strip() == 'Singly Occupied:'
            elif self.reference[0:1] == 'R':
                assert unoccupied.strip() == 'Virtual:'
            elif self.reference[0:1] == 'U':
                assert unoccupied.strip() == 'Alpha Virtual:'

            # Psi4 now has a blank line, Psi3 does not.
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
        if self.section == "Post-Iterations" and "Final Energy:" in line \
            and not hasattr(self, 'finite_difference'):
            e = float(line.split()[3])
            if not hasattr(self, 'scfenergies'):
                self.scfenergies = []
            self.scfenergies.append(utils.convertor(e, 'hartree', 'eV'))

        if self.subsection == "Energetics":
            if "Empirical Dispersion Energy" in line:
                dispersion = utils.convertor(float(line.split()[-1]), "hartree", "eV")
                self.append_attribute("dispersionenergies", dispersion)

        #   ==> Molecular Orbitals <==
        #
        #                     1            2            3            4            5
        #
        # 1    H1 s0         0.1610392    0.1040990    0.0453848    0.0978665    1.0863246
        # 2    H1 s0         0.3066996    0.0742959    0.8227318    1.3460922   -1.6429494
        # 3    H1 s0         0.1669296    1.5494169   -0.8885631   -1.8689490    1.0473633
        # 4    H2 s0         0.1610392   -0.1040990    0.0453848   -0.0978665   -1.0863246
        # 5    H2 s0         0.3066996   -0.0742959    0.8227318   -1.3460922    1.6429494
        # 6    H2 s0         0.1669296   -1.5494169   -0.8885631    1.8689490   -1.0473633
        #
        #             Ene    -0.5279195    0.1235556    0.3277474    0.5523654    2.5371710
        #             Sym            Ag          B3u           Ag          B3u          B3u
        #             Occ             2            0            0            0            0
        #
        #
        #                         6
        #
        # 1    H1 s0         1.1331221
        # 2    H1 s0        -1.2163107
        # 3    H1 s0         0.4695317
        # 4    H2 s0         1.1331221
        # 5    H2 s0        -1.2163107
        # 6    H2 s0         0.4695317
        #
        #            Ene     2.6515637
        #            Sym            Ag
        #            Occ             0

        if (self.section) and ("Molecular Orbitals" in self.section) \
           and ("Molecular Orbitals" in line) and not hasattr(self, 'finite_difference'):

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

                n = len(indices)
                line = next(inputfile)
                while line.strip():
                    chomp = line.split()
                    m = len(chomp)
                    iao = int(chomp[0])
                    coeffs = [float(c) for c in chomp[m - n:]]
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
            if line.strip() == f"{pop_type} Charges: (a.u.)":
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
            self.metadata["methods"].append("MP2")
            mpenergy = utils.convertor(float(line.split()[-1]), 'hartree', 'eV')
            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []
            self.mpenergies.append([mpenergy])
        # This is for the newer DF-MP2 code in 4.0.
        if 'DF-MP2 Energies' in line:
            self.metadata["methods"].append("DF-MP2")
            while 'Total Energy' not in line:
                line = next(inputfile)
            mpenergy = utils.convertor(float(line.split()[3]), 'hartree', 'eV')
            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []
            self.mpenergies.append([mpenergy])

        # Note this is just a start and needs to be modified for CCSD(T), etc.
        ccsd_trigger = "* CCSD total energy"
        if line.strip()[:len(ccsd_trigger)] == ccsd_trigger:
            self.metadata["methods"].append("CCSD")
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
        if (self.section == "Convergence Check") and line.strip() == "==> Convergence Check <==" \
            and not hasattr(self, 'finite_difference'):

            if not hasattr(self, "optstatus"):
                self.optstatus = []
            self.optstatus.append(data.ccData.OPT_UNKNOWN)

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
            step = int(values.split()[0])
            geovalues = []
            for istart in starts:
                if values[istart:istart+9].strip():
                    geovalues.append(float(values[istart:istart+9]))

            if step == 1:
                self.optstatus[-1] += data.ccData.OPT_NEW

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
                self.optstatus[-1] += data.ccData.OPT_DONE

        # This message means that optimization has stopped for some reason, but we
        # still want optdone to exist in this case, although it will be an empty list.
        if line.strip() == "Optimizer: Did not converge!":
            if not hasattr(self, 'optdone'):
                self.optdone = []

            assert hasattr(self, "optstatus") and len(self.optstatus) > 0
            self.optstatus[-1] += data.ccData.OPT_UNCONVERGED

        # The reference point at which properties are evaluated in Psi4 is explicitely stated,
        # so we can save it for later. It is not, however, a part of the Properties section,
        # but it appears before it and also in other places where properies that might depend
        # on it are printed.
        #
        # Properties will be evaluated at   0.000000,   0.000000,   0.000000 Bohr
        #
        # OR
        #
        # Properties will be evaluated at   0.000000,   0.000000,   0.000000 [a0]
        #
        if "Properties will be evaluated at" in line.strip():
            self.origin = numpy.array([float(x.strip(',')) for x in line.split()[-4:-1]])
            assert line.split()[-1] in ["Bohr", "[a0]"]
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
                    fromunits = f"ebohr{(rank > 1) * f'{int(rank)}'}"
                    tounits = (
                        f"Debye{(rank > 1) * '.ang'}{(rank > 2) * f'{int(rank - 1)}'}"
                    )
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

        ## Analytic Gradient
        #        -Total Gradient:
        #   Atom            X                  Y                   Z
        #  ------   -----------------  -----------------  -----------------
        #     1       -0.000000000000     0.000000000000    -0.064527252292
        #     2        0.000000000000    -0.028380539652     0.032263626146
        #     3       -0.000000000000     0.028380539652     0.032263626146

        ## Finite Differences Gradient
        # -------------------------------------------------------------
        #   ## F-D gradient (Symmetry 0) ##
        #   Irrep: 1 Size: 3 x 3
        #
        #                  1                   2                   3
        #
        #     1     0.00000000000000     0.00000000000000    -0.02921303282515
        #     2     0.00000000000000    -0.00979709321487     0.01460651641258
        #     3     0.00000000000000     0.00979709321487     0.01460651641258
        if line.strip() in Psi4.GRADIENT_HEADERS:

            # Handle the different header lines between analytic and
            # numerical gradients.
            gradient_skip_lines = [
                info.skip_lines
                for info in Psi4.GRADIENT_TYPES.values()
                if info.header == line.strip()
            ][0]
            gradient = self.parse_gradient(inputfile, gradient_skip_lines)

            if not hasattr(self, 'grads'):
                self.grads = []
            self.grads.append(gradient)

        # OLD Normal mode output parser (PSI4 < 1)

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
                vibfreq = Psi4.parse_vibfreq(chomp[1])
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
                _vibfreq = Psi4.parse_vibfreq(line[13:].strip())
                assert abs(vibfreq - _vibfreq) < 1.0e-2
                line = next(inputfile)
                assert 'Force constant:' in line
                if not hasattr(self, "vibfconsts"):
                    self.vibfconsts = []
                self.vibfconsts.append(
                    utils.convertor(float(line.split()[2]), "hartree/bohr2", "mDyne/angstrom")
                )
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

        # NEW Normal mode output parser (PSI4 >= 1)

        #   ==> Harmonic Vibrational Analysis <==
        #   ...
        #   Vibration                       7                   8                   9
        #   ...
        #
        #   Vibration                       10                  11                  12
        #   ...

        if line.strip() == '==> Harmonic Vibrational Analysis <==':

            vibsyms = []
            vibfreqs = []
            vibdisps = []
            vibrmasses = []
            vibfconsts = []

            # Skip lines till the first Vibration block
            while not line.strip().startswith('Vibration'):
                line = next(inputfile)

            n_modes = 0
            # Parse all the Vibration blocks
            while line.strip().startswith('Vibration'):
                n = len(line.split()) - 1
                n_modes += n
                vibfreqs_, vibsyms_, vibdisps_, vibrmasses_, vibfconsts_ = self.parse_vibration(n, inputfile)
                vibfreqs.extend(vibfreqs_)
                vibsyms.extend(vibsyms_)
                vibdisps.extend(vibdisps_)
                vibrmasses.extend(vibrmasses_)
                vibfconsts.extend(vibfconsts_)
                line = next(inputfile)

            # It looks like the symmetry of the normal mode may be missing 
            # from some / most. Only include them if they are there for all

            if len(vibfreqs) == n_modes:
                self.set_attribute('vibfreqs', vibfreqs)

            if len(vibsyms) == n_modes:
                self.set_attribute('vibsyms', vibsyms)

            if len(vibdisps) == n_modes:
                self.set_attribute('vibdisps', vibdisps)

            if len(vibdisps) == n_modes:
                self.set_attribute('vibrmasses', vibrmasses)

            if len(vibdisps) == n_modes:
                self.set_attribute('vibfconsts', vibfconsts)

        # Second one is 1.0, first one is 1.2 and newer
        if (self.section == "Thermochemistry Energy Analysis" and "Thermochemistry Energy Analysis" in line) \
           or (self.section == "Energy Analysis" and "Energy Analysis" in line):

            self.skip_lines(
                inputfile,
                [
                    "b",
                    "Raw electronic energy",
                    "Total E0",
                    "b",
                    "Zero-point energy, ZPE_vib = Sum_i nu_i / 2",
                    "Electronic ZPE",
                    "Translational ZPE",
                    "Rotational ZPE"
                ]
            )
            line = next(inputfile)
            assert "Vibrational ZPE" in line
            self.set_attribute("zpve", float(line.split()[6]))

        # If finite difference is used to compute forces (i.e. by displacing
        # slightly all the atoms), a series of additional scf calculations is
        # performed. Orbitals, geometries, energies, etc. for these shouln't be
        # included in the parsed data.

        if line.strip().startswith('Using finite-differences of gradients'):
            self.set_attribute('finite_difference', True)

        # This is the result of calling `print_variables()` and contains all
        # current inner variables known to Psi4.
        if line.strip() == "Variable Map:":
            self.skip_line(inputfile, "d")
            line = next(inputfile)
            while line.strip():
                tokens = line.split()
                # Remove double quotation marks
                name = " ".join(tokens[:-2])[1:-1]
                value = float(tokens[-1])
                if name == "CC T1 DIAGNOSTIC":
                    self.metadata["t1_diagnostic"] = value
                line = next(inputfile)

        if line[:54] == '*** Psi4 exiting successfully. Buy a developer a beer!'\
                or line[:54] == '*** PSI4 exiting successfully. Buy a developer a beer!':
            self.metadata['success'] = True

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

    def parse_gradient(self, inputfile, skip_lines):
        """Parse the nuclear gradient section into a list of lists with shape
        [natom, 3].
        """
        self.skip_lines(inputfile, skip_lines)
        line = next(inputfile)
        gradient = []

        while line.strip():
            idx, x, y, z = line.split()
            gradient.append((float(x), float(y), float(z)))
            line = next(inputfile)
        return gradient

    @staticmethod
    def parse_vibration(n, inputfile):

        #   Freq [cm^-1]                1501.9533           1501.9533           1501.9533
        #   Irrep
        #   Reduced mass [u]              1.1820              1.1820              1.1820
        #   Force const [mDyne/A]         1.5710              1.5710              1.5710
        #   Turning point v=0 [a0]        0.2604              0.2604              0.2604
        #   RMS dev v=0 [a0 u^1/2]        0.2002              0.2002              0.2002
        #   Char temp [K]               2160.9731           2160.9731           2160.9731
        #   ----------------------------------------------------------------------------------
        #       1   C               -0.00  0.01  0.13   -0.00 -0.13  0.01   -0.13  0.00 -0.00
        #       2   H                0.33 -0.03 -0.38    0.02  0.60 -0.02    0.14 -0.01 -0.32
        #       3   H               -0.32 -0.03 -0.37   -0.01  0.60 -0.01    0.15 -0.01  0.33
        #       4   H                0.02  0.32 -0.36    0.01  0.16 -0.34    0.60 -0.01  0.01
        #       5   H                0.02 -0.33 -0.39    0.01  0.13  0.31    0.60  0.01  0.01

        line = next(inputfile)
        assert 'Freq' in line
        chomp = line.split()
        vibfreqs = [Psi4.parse_vibfreq(x) for x in chomp[-n:]]

        line = next(inputfile)
        assert 'Irrep' in line
        chomp = line.split()
        vibsyms = [irrep for irrep in chomp[1:]]

        line = next(inputfile)
        assert 'Reduced mass' in line
        chomp = line.split()
        vibrmasses = [utils.float(x) for x in chomp[3:]]

        line = next(inputfile)
        assert 'Force const' in line
        chomp = line.split()
        vibfconsts = [utils.float(x) for x in chomp[3:]]

        line = next(inputfile)
        assert 'Turning point' in line

        line = next(inputfile)
        assert 'RMS dev' in line

        line = next(inputfile)
        assert 'Char temp' in line

        line = next(inputfile)
        assert '---' in line

        line = next(inputfile)
        vibdisps = [ [] for i in range(n)]
        while len(line.strip()) > 0:
            chomp = line.split()
            for i in range(n):
                start = len(chomp) - (n - i) * 3
                stop = start + 3
                mode_disps = [float(c) for c in chomp[start:stop]]
                vibdisps[i].append(mode_disps)

            line = next(inputfile)

        return vibfreqs, vibsyms, vibdisps, vibrmasses, vibfconsts

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
