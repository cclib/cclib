# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for MOLDEN format files."""

import os.path

from . import filewriter
from cclib.parser import utils


class MOLDEN(filewriter.Writer):
    """A writer for MOLDEN files."""

    def _title(self, path):
        """Return filename without extension to be used as title."""
        title = os.path.basename(os.path.splitext(path)[0])
        return title

    def _coords_from_ccdata(self, index):
        """Create [Atoms] section using geometry at the given index."""
        elements = [self.pt.element[Z] for Z in self.ccdata.atomnos]
        atomcoords = self.ccdata.atomcoords[index]
        atomnos = self.ccdata.atomnos
        nos = range(self.ccdata.natom)

        # element_name number atomic_number x y z
        atom_template = '{:2s} {:5d} {:2d} {:12.6f} {:12.6f} {:12.6f}'
        lines = []
        for element, no, atomno, (x, y, z) in zip(elements, nos, atomnos,
                                                  atomcoords):
            lines.append(atom_template.format(element, no + 1, atomno,
                                              x, y, z))

        return lines

    def _gto_from_ccdata(self):
        """Create [GTO] section using gbasis."""

        # atom_sequence_number1 0
        # shell_label number_of_primitives 1.00
        # exponent_primitive_1 contraction_coeff_1 (contraction_coeff_1)
        # ...
        # empty line
        # atom_sequence__number2 0
        gbasis = self.ccdata.gbasis
        label_template = '{:s} {:5d} 1.00'
        basis_template = '{:15.6e} {:15.6e}'
        lines = []

        for no, basis in enumerate(gbasis):
            lines.append('{:3d} 0'.format(no + 1))
            for prims in basis:
                lines.append(label_template.format(prims[0].lower(),
                                                   len(prims[1])))
                for prim in prims[1]:
                    lines.append(basis_template.format(prim[0], prim[1]))
            lines.append('')

        return lines

    def _scfconv_from_ccdata(self):
        """Create [SCFCONV] section using gbasis."""

        # scf-first    1 THROUGH   12
        #    -672.634394
        #    ...
        #    -673.590571
        #    -673.590571
        lines = ["scf-first    1 THROUGH   %d" % len(self.ccdata.scfenergies)]

        for scfenergy in self.ccdata.scfenergies:
            lines.append('{:15.6f}'.format(scfenergy))

        return lines

    def _mo_from_ccdata(self):
        """Create [MO] section."""

        # Sym= symmetry_label_1
        # Ene= mo_energy_1
        # Spin= (Alpha|Beta)
        # Occup= mo_occupation_number_1
        # ao_number_1 mo_coefficient_1
        # ...
        # ao_number_n mo_coefficient_n
        # ...
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

        for i in range(mult):
            spin = 'Alpha'
            for j in range(len(moenergies[i])):
                if has_syms:
                    lines.append(' Sym= %s' % syms[i][j])
                moenergy = utils.convertor(moenergies[i][j], 'eV', 'hartree')
                lines.append(' Ene= {:10.4f}'.format(moenergy))
                lines.append(' Spin= %s' % spin)
                if j <= homos[i]:
                    lines.append(' Occup= {:10.6f}'.format(2.0 / mult))
                else:
                    lines.append(' Occup= {:10.6f}'.format(0.0))
                for k, mocoeff in enumerate(mocoeffs[i][j]):
                    lines.append('{:4d}  {:10.6f}'.format(k + 1, mocoeff))

            spin = 'Beta'

        return lines

    def generate_repr(self):
        """Generate the MOLDEN representation of the logfile data."""

        molden_lines = ['[Molden Format]']

        # Title of file.
        molden_lines.append('[Title]')
        molden_lines.append(self._title(self.jobfilename))

        # Coordinates for the Electron Density/Molecular orbitals.
        # [Atoms] (Angs|AU)
        unit = "Angs"
        molden_lines.append('[Atoms] %s' % unit)
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

        if hasattr(self.ccdata, 'scfenergies'):
            if len(self.ccdata.scfenergies) > 1:
                molden_lines.append('[SCFCONV]')
                molden_lines.extend(self._scfconv_from_ccdata())

        molden_lines.append('')

        return '\n'.join(molden_lines)


if __name__ == "__main__":
    pass
