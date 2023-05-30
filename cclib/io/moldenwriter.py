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
from itertools import zip_longest

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
    """A writer for MOLDEN files.

    Documentation of the format is located at
    https://www3.cmbi.umcn.nl/molden/molden_format.html.
    """

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

        lines = []
        for element, no, atomno, coord in zip(elements, nos, atomnos,
                                              atomcoords):
            x, y, z = map(round_molden, coord)
            lines.append(f"{element:2s} {no + 1:5d} {atomno:2d} {x:12.6f} {y:12.6f} {z:12.6f}")

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
        lines = []

        for no, basis in enumerate(gbasis):
            lines.append(f"{no + 1:3d} 0")
            for prims in basis:
                lines.append(f"{prims[0].lower():s} {len(prims[1]):5d} 1.00")
                for prim in prims[1]:
                    lines.append(f"{prim[0]:15.9e} {prim[1]:15.9e}")
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

    def _syms_energies_occs_coeffs_from_ccdata_for_moldenwriter(self, data=None):
        syms = None
        energies = None
        occs = None
        coeffs = None

        if data is None:
            data = self.ccdata

        if not self.naturalorbitals and hasattr(data, 'moenergies') and hasattr(data, 'mocoeffs'):
            energies = data.moenergies
            coeffs = data.mocoeffs
            occs = numpy.zeros((len(data.homos),len(energies[0])))
            occval = 2 // len(data.homos)
            for i in range(len(data.homos)):
                occs[i][0:data.homos[i]+1] = occval
        elif self.naturalorbitals and hasattr(data, 'nooccnos') and hasattr(data, "nocoeffs"):
            energies = numpy.array([data.nooccnos])
            coeffs = numpy.array([data.nocoeffs])
            occs = numpy.array([data.nooccnos])
        elif self.naturalorbitals and hasattr(data, 'nsooccnos') and hasattr(data, "nsocoeffs"):
            energies = data.nsooccnos
            coeffs = data.nsocoeffs
            occs = data.nsooccnos

        if hasattr(data, 'mosyms') and not self.naturalorbitals:
            syms = data.mosyms
        else:
            syms = numpy.full_like(energies, 'A', dtype=str)


        return syms, energies, occs, coeffs

    def _mo_from_ccdata(self, mosyms, moenergies, mooccs, mocoeffs):
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

        lines = []

        spin = 'Alpha'
        for i in range(len(mooccs)):
            for j in range(len(mooccs[i])):
                restricted_spin_idx = i % len(mocoeffs)
                lines.append(' Sym= {}'.format(mosyms[restricted_spin_idx][j]))
                moenergy = utils.convertor(moenergies[restricted_spin_idx][j], 'eV', 'hartree')
                lines.append(' Ene= {:10.4f}'.format(moenergy))
                lines.append(' Spin= {}'.format(spin))
                lines.append(' Occup= {:10.6f}'.format(mooccs[i][j]))
                # Rearrange mocoeffs according to Molden's lexicographical order.
                mocoeffs[restricted_spin_idx][j] = self._rearrange_mocoeffs(mocoeffs[restricted_spin_idx][j])
                for k, mocoeff in enumerate(mocoeffs[restricted_spin_idx][j]):
                    lines.append('{:4d}  {:10.6f}'.format(k + 1, mocoeff))
            spin = 'Beta'
        return lines

    def _freq_from_ccdata(self):

        lines = []

        if hasattr(self.ccdata, "vibfreqs"):
            vibfreqs = self.ccdata.vibfreqs
            vibfreqs_lines = ["[FREQ]"]
            vibfreqs_lines.extend([f"{vibfreq:16.8f}" for vibfreq in vibfreqs])
            lines.append("\n".join(vibfreqs_lines))

        if all(hasattr(self.ccdata, attrname) for attrname in ("atomcoords", "atomnos")) and \
           any(hasattr(self.ccdata, attrname) for attrname in ("vibfreqs", "vibdisps", "vibirs")):
            # Selecting the first set of coordinates works for when the
            # frequency calculation is done via finite difference, but not
            # with multi-part inputs, which there is currently no way of
            # detecting.
            atomcoords = utils.convertor(self.ccdata.atomcoords[0], "Angstrom", "bohr")
            atomsyms = (self.pt.element[atomno] for atomno in self.ccdata.atomnos)
            atomcoords_lines = ["[FR-COORD]"]
            for atomsym, atomcoord in zip(atomsyms, atomcoords):
                atomcoords_lines.append(
                    f"{atomsym:3s} {atomcoord[0]:15.8f} {atomcoord[1]:15.8f} {atomcoord[2]:15.8f}"
                )
            lines.append("\n".join(atomcoords_lines))

        if hasattr(self.ccdata, "vibdisps"):
            vibdisps = self.ccdata.vibdisps
            vibdisps_lines = ["[FR-NORM-COORD]"]
            for vibidx in range(vibdisps.shape[0]):
                vibdisps_lines.append(f"vibration {vibidx + 1}")
                for iatom in range(vibdisps.shape[1]):
                    vibdisp = vibdisps[vibidx, iatom]
                    vibdisps_lines.append(
                        f"{vibdisp[0]:15.8f} {vibdisp[1]:15.8f} {vibdisp[2]:15.8f}"
                    )
            lines.append("\n".join(vibdisps_lines))

        if hasattr(self.ccdata, "vibirs"):
            vibirs = self.ccdata.vibirs
            has_vibramans = hasattr(self.ccdata, "vibramans")
            vibramans = [] if not has_vibramans else self.ccdata.vibramans
            vibirs_lines = ["[INT]"]
            for vibir, vibraman in zip_longest(vibirs, vibramans):
                if not has_vibramans:
                    vibirs_lines.append(f"{vibir:12.6f}")
                else:
                    vibirs_lines.append(f"{vibir:12.6f} {vibraman:12.6f}")
            lines.append("\n".join(vibirs_lines))

        if lines:
            lines.append("")

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

        mosyms, moenergies, mooccs, mocoeffs = self._syms_energies_occs_coeffs_from_ccdata_for_moldenwriter()

        if hasattr(self.ccdata, 'gbasis'):
            molden_lines.append('[GTO]')
            molden_lines.extend(self._gto_from_ccdata())
        if all(attr is not None for attr in (mosyms, moenergies, mooccs)):
            molden_lines.append('[MO]')
            molden_lines.extend(self._mo_from_ccdata(mosyms, moenergies, mooccs, mocoeffs))

        # Omitting until issue #390 is resolved.
        # https://github.com/cclib/cclib/issues/390
        # if hasattr(self.ccdata, 'scfenergies'):
        #     if len(self.ccdata.scfenergies) > 1:
        #         molden_lines.append('[SCFCONV]')
        #         molden_lines.extend(self._scfconv_from_ccdata())

        # molden_lines.append('')

        molden_lines.extend(self._freq_from_ccdata())

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
        return f"{decimal.Decimal(num):.9e}"

    def reformat(self):
        """Reformat Molden output file to:
        - use scientific notation,
        - split sp molecular orbitals to s and p, and
        - replace multiple spaces with single."""
        filelines = iter(self.filestring.split("\n"))
        lines = []
        is_header = False

        for line in filelines:
            line = line.replace('\n', '')
            # Replace multiple spaces with single spaces.
            line = ' '.join(line.split())

            is_header = line and line[0] == "[" and line[-1] == "]"

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
            if not is_header and 'D' in line:
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
