# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for wfx format files."""

import os.path

from . import filewriter

class WFXWriter(filewriter.Writer):
    """A writer for wfx files."""

    required_attrs = ['atomcoords', 'atomnos', 'gbasis', 'charge', 'homos']
    keyword = 'GTO'

    @staticmethod
    def _section(section_name, section_data):
        """Add opening/closing section_name tags to data."""
        opening_tag = ['<' + section_name + '>']
        closing_tag = ['</' + section_name + '>']
        
        if isinstance(section_data, list):
            return opening_tag + section_data + closing_tag
        elif isinstance(section_data, str):
            return opening_tag + section_data.split() + closing_tag
        elif isinstance(section_data, int) or isinstance(section_data, float):
            return opening_tag + [str(section_data)] + closing_tag
        else:
            return opening_tag + closing_tag

    def _title(self):
        """Return filename without extension to be used as title."""
        if self.jobfilename is not None:
            title = os.path.basename(os.path.splitext(self.jobfilename)[0])
            return title

    def _keywords(self):
        """Return one of GTO, GIAO, CSGT keyword."""
        # AIMALL currently supports on GTO keyword.
        return self.keyword

    def _no_of_nuclei(self):
        """Number of nuclei in the molecule."""
        no_nuclei = len(self.ccdata.atomnos)
        return no_nuclei

    def _no_of_prims(self):
        """Number of Primitives."""
        nprims = 0
        for atom in self.ccdata.gbasis:
            for prims in atom:
                nprims += len(prims)
        return nprims

    def _no_of_mos(self):
        """Number of occupied MOs."""
        return int(max(self.ccdata.homos))

    def _no_of_perturbations(self):
        """Number of Perturbation."""
        # This is usually zero.  For GIAO it should be 3 (corresponding to Lx, Ly and Lz) and for
        # CSGT it should be 6 (corresponding to Lx, Ly, Lz, Px, Py and Pz).
        nperturb = 0
        if self.keyword == 'GIAO':
            nperturb = 3
        elif self.keyword == 'CSGT':
            nperturb = 6
        return nperturb

    def _nuclear_names(self):
        """Names of nuclei present in the molecule."""
        # O1 H2 H3 etc.
        return [self.pt.element[Z]+str(i) for i, Z in enumerate(self.ccdata.atomnos)]

    def _atomic_nos(self):
        """Atomic numbers of elements."""
        return [str(Z) for Z in self.ccdata.atomnos]

    def _nuclear_charges(self):
        """Nuclear charges."""
        return ['{:.14e}'.format(Z) for Z in self.ccdata.atomnos]

    def _nuclear_coords(self):
        """Nuclear charges."""
        coord_template = '%.14E %.14E %.14E'        
        nuc_coords = [coord_template % tuple(coord) for coord in self.ccdata.atomcoords[-1]]
        return nuc_coords

    def _net_charge(self):
        """Net charge on molecule."""
        return '%.14E' % self.ccdata.charge

    def _no_electrons(self):
        """Number of electrons in molecule."""
        homos = self.ccdata.homos
        charge = self.ccdata.charge
        print (sum(homos), charge)
        if len(homos) > 1:
            return int(sum(homos))
        elif charge > 0:
            return int(homos[0]*2 - (charge % 2))

    def generate_repr(self):
        """Generate the wfx representation of the logfile data."""

        sections = [
            (self._title, "Title"),
            (self._keywords, "Keywords"),
            (self._no_of_prims, "Number of Primitives"),
            (self._no_of_mos, "Number of Occupied Molecular Orbitals"),
            (self._no_of_perturbations, "Number of Perturbations"),
            (self._nuclear_names, "Nuclear Names"),
            (self._atomic_nos, "Atomic Numbers"),
            (self._nuclear_charges, "Nuclear Charges"),
            (self._nuclear_coords, "Nuclear Cartesian Coordinates"),
            (self._net_charge, "Net Charge"),
            (self._no_electrons, "Number of Electrons"),
            # (self._no_alpha_electrons, "Number of Alpha Electrons"),
            # (self._no_beta_electrons, "Number of Beta Electrons"),
            # (self._spin_mult, "Electronic Spin Multiplicity"),
            # (self._model, "Model"),
            # (self._prim_centers, "Primitive Centers"),
            # (self._prim_types, "Primitive Types"),
            # (self._prim_exps, "Primitive Exponents"),
            # (self._mo_occup_nos, "Molecular Orbital Occupation Numbers"),
            # (self._mo_energies, "Molecular Orbital Energies"),
            # (self._mo_spin_types, "Molecular Orbital Spin Types"),
            # (self._mo_prim_coeffs, "Molecular Orbital Primitive Coefficients"),
            # (self._energy, "Energy"),
            # (self._nuc_energy_gradients, "Nuclear Cartesian Energy Gradients"),
            # (self._nuc_virial, "Nuclear Virial of Energy-Gradient-Based Forces on Nuclei, W"),
            # (self._virial_ratio, "Full Virial Ratio"),
        ]

        wfx_lines = []

        for section_module, section_name in sections:
            try:
                section_data = section_module()
            except:
                section_data = None
            wfx_lines.extend(self._section(section_name, section_data))

        wfx_lines.append('')
        return '\n'.join(wfx_lines)


if __name__ == "__main__":
    pass
