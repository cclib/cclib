# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for MOLDEN format files."""

import os.path
import math
import decimal
import numpy 

from cclib.parser import utils
from cclib.io import filewriter


def round_molden(num, p=6):
    """Molden style number rounding in [Atoms] section."""
    # Digit at pth position after dot.
    p_digit = math.floor(abs(num) * 10 ** p) % 10
    # If the 6th digit after dot is greater than 5, but is not 7,
    # round the number upto 6th place.
    # Else truncate at 6th digit after dot.
    if p_digit > 5 and p_digit != 7:
        return round(num, p)
    if num >= 0:
        return math.floor(num * 10 ** p) / 10 ** p
    else:
        return math.ceil(num * 10 ** p) / 10 ** p


class MOLDEN(filewriter.Writer):
    """A writer for MOLDEN files."""

    required_attrs = ('atomcoords', 'atomnos', 'natom')

    def _title(self, path):
        """Return filename without extension to be used as title."""
        title = os.path.basename(os.path.splitext(path)[0])
        return title

    def _coords_from_ccdata(self, index):
        """Create [Atoms] section using geometry at the given index."""
        elements = [self.pt.element[Z] for Z in self.ccdata.atomnos]
        if self.ghost is not None:
            elements = [self.ghost if e is None else e for e in elements]
        elif None in elements:
            raise ValueError(
                f"It seems that there is at least one ghost atom in these elements. Please use the ghost flag to specify a label for the ghost atoms."
            )
        atomcoords = self.ccdata.atomcoords[index]
        atomnos = self.ccdata.atomnos
        nos = range(self.ccdata.natom)

        # element_name number atomic_number x y z
        atom_template = "{:2s} {:5d} {:2d} {:12.6f} {:12.6f} {:12.6f}"
        lines = []
        for element, no, atomno, coord in zip(elements, nos, atomnos,
                                              atomcoords):
            x, y, z = map(round_molden, coord)
            lines.append(atom_template.format(element, no + 1, atomno, x, y, z))

        return lines

    def _gto_from_ccdata(self):
        """Create [GTO] section using gbasis.

        atom_sequence_number1 0
        shell_label number_of_primitives 1.00
        exponent_primitive_1 contraction_coeff_1 (contraction_coeff_1)
        ...
        empty line
        atom_sequence__number2 0
        """

        gbasis = self.ccdata.gbasis
        label_template = "{:s} {:5d} 1.00"
        basis_template = "{:15.9e} {:15.9e}"
        lines = []

        for no, basis in enumerate(gbasis):
            lines.append(f"{no + 1:3d} 0")
            for prims in basis:
                lines.append(label_template.format(prims[0].lower(),
                                                   len(prims[1])))
                for prim in prims[1]:
                    lines.append(basis_template.format(prim[0], prim[1]))
            lines.append('')
        lines.append('')
        return lines

    def _scfconv_from_ccdata(self):
        """Create [SCFCONV] section using gbasis.

        scf-first    1 THROUGH   12
           -672.634394
           ...
           -673.590571
           -673.590571
        """

        lines = [f"scf-first    1 THROUGH   {len(self.ccdata.scfenergies)}"]

        for scfenergy in self.ccdata.scfenergies:
            lines.append(f"{scfenergy:15.6f}")

        return lines

    def _rearrange_mocoeffs(self, mocoeffs):
        """Rearrange cartesian F functions in mocoeffs.

        Molden's order:
        xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
        cclib's order:
        XXX, YYY, ZZZ, XXY, XXZ, YYX, YYZ, ZZX, ZZY, XYZ
        cclib's order can be converted by:
        moving YYX two indexes ahead, and
        moving YYZ two indexes back.
        """
        aonames = self.ccdata.aonames
        mocoeffs = mocoeffs.tolist()

        pos_yyx = [key for key, val in enumerate(aonames) if '_YYX' in val]
        pos_yyz = [key for key, val in enumerate(aonames) if '_YYZ' in val]

        if pos_yyx:
            for pos in pos_yyx:
                mocoeffs.insert(pos-2, mocoeffs.pop(pos))
        if pos_yyz:
            for pos in pos_yyz:
                mocoeffs.insert(pos+2, mocoeffs.pop(pos))

        return mocoeffs


    def _mo_from_ccdata(self):
        """Create [MO] section.

        Sym= symmetry_label_1
        Ene= mo_energy_1
        Spin= (Alpha|Beta)
        Occup= mo_occupation_number_1
        ao_number_1 mo_coefficient_1
        ...
        ao_number_n mo_coefficient_n
        ...
        """

        moenergies = self.ccdata.moenergies
        mocoeffs = self.ccdata.mocoeffs
        homos = self.ccdata.homos
        mult = self.ccdata.mult

        has_syms = False
        lines = []

        # Sym attribute is optional in [MO] section.
        if hasattr(self.ccdata, 'mosyms'):
            has_syms = True
            syms = self.ccdata.mosyms
        else:
            syms = numpy.full_like(moenergies, 'A', dtype=str)
        unres = len(moenergies) > 1
        openshell = len(homos) > 1

        spin = 'Alpha'
        for i in range(len(moenergies)):
            for j in range(len(moenergies[i])):
                lines.append(f" Sym= {syms[i][j]}")
                moenergy = utils.convertor(moenergies[i][j], "eV", "hartree")
                lines.append(f" Ene= {moenergy:10.4f}")
                lines.append(f" Spin= {spin}")
                if unres and openshell:
                    if j <= homos[i]:
                        lines.append(f" Occup= {1.0:10.6f}")
                    else:
                        lines.append(f" Occup= {0.0:10.6f}")
                elif not unres and openshell:
                    occ = numpy.sum(j <= homos)
                    if j <= homos[i]:
                        lines.append(f" Occup= {occ:10.6f}")
                    else:
                        lines.append(f" Occup= {0.0:10.6f}")
                else:
                    if j <= homos[i]:
                        lines.append(f" Occup= {2.0:10.6f}")
                    else:
                        lines.append(f" Occup= {0.0:10.6f}")
                # Rearrange mocoeffs according to Molden's lexicographical order.
                mocoeffs[i][j] = self._rearrange_mocoeffs(mocoeffs[i][j])
                for k, mocoeff in enumerate(mocoeffs[i][j]):
                    lines.append(f"{k + 1:4d}  {mocoeff:10.6f}")

            spin = 'Beta'

        return lines

    def generate_repr(self):
        """Generate the MOLDEN representation of the logfile data."""

        molden_lines = ['[Molden Format]']

        # Title of file.
        if self.jobfilename is not None:
            molden_lines.append('[Title]')
            molden_lines.append(self._title(self.jobfilename))

        # Coordinates for the Electron Density/Molecular orbitals.
        # [Atoms] (Angs|AU)
        unit = "Angs"
        molden_lines.append(f"[Atoms] {unit}")
        # Last set of coordinates for geometry optimization runs.
        index = -1
        molden_lines.extend(self._coords_from_ccdata(index))

        # Either both [GTO] and [MO] should be present or none of them.
        if hasattr(self.ccdata, 'gbasis') and hasattr(self.ccdata, 'mocoeffs')\
                and hasattr(self.ccdata, 'moenergies'):

            molden_lines.append('[GTO]')
            molden_lines.extend(self._gto_from_ccdata())

            molden_lines.append('[MO]')
            molden_lines.extend(self._mo_from_ccdata())

        # Omitting until issue #390 is resolved.
        # https://github.com/cclib/cclib/issues/390
        # if hasattr(self.ccdata, 'scfenergies'):
        #     if len(self.ccdata.scfenergies) > 1:
        #         molden_lines.append('[SCFCONV]')
        #         molden_lines.extend(self._scfconv_from_ccdata())

        # molden_lines.append('')

        return '\n'.join(molden_lines)


class MoldenReformatter:
    """Reformat Molden output files."""

    def __init__(self, filestring):
        self.filestring = filestring

    def scinotation(self, num):
        """Convert Molden style number formatting to scientific notation.
        0.9910616900D+02 --> 9.910617e+01
        """
        num = num.replace("D", "e")
        return str(f"{decimal.Decimal(num):.9e}")

    def reformat(self):
        """Reformat Molden output file to:
        - use scientific notation,
        - split sp molecular orbitals to s and p, and
        - replace multiple spaces with single."""
        filelines = iter(self.filestring.split("\n"))
        lines = []

        for line in filelines:
            line = line.replace('\n', '')
            # Replace multiple spaces with single spaces.
            line = ' '.join(line.split())

            # Check for [Title] section.
            if '[title]' in line.lower():
                # skip the title
                line = next(filelines)
                line = next(filelines)

            # Exclude SCFCONV section until issue #390 is resolved.
            # https://github.com/cclib/cclib/issues/390
            if '[scfconv]' in line.lower():
                break

            # Although Molden format specifies Sym in [MO] section,
            # the Molden program does not print it.
            if 'sym' in line.lower():
                continue

            # Convert D notation to scientific notation.
            if 'D' in line:
                vals = line.split()
                vals = [self.scinotation(i) for i in vals]
                lines.append(' '.join(vals))

            # Convert sp to s and p orbitals.
            elif 'sp' in line:
                n_prim = int(line.split()[1])
                new_s = [f"s {str(n_prim)} 1.00"]
                new_p = [f"p {str(n_prim)} 1.00"]
                while n_prim > 0:
                    n_prim -= 1
                    line = next(filelines).split()
                    new_s.append(
                        f"{self.scinotation(line[0])} {self.scinotation(line[1])}"
                    )
                    new_p.append(
                        f"{self.scinotation(line[0])} {self.scinotation(line[2])}"
                    )
                lines.extend(new_s)
                lines.extend(new_p)
            else:
                lines.append(line)

        return '\n'.join(lines)
