# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for wfx format files."""

import os.path
import numpy
from cclib.parser import utils

from . import filewriter

class WFXWriter(filewriter.Writer):
    """A writer for wfx files."""

    required_attrs = ('atomcoords', 'atomnos', 'gbasis', 'charge', 'homos',
                      'mult')
    keyword = 'GTO'
    # Number of orbitals of type key.
    num_orb = {'s':1, 'p':3, 'd':6, 'f':10, 'g':15, 'h':21,
               'S':1, 'P':3, 'D':6, 'F':10, 'G':15, 'H':21}
    # Index of first orbital of type key in a list of orbitals.
    idx_orb = {'s':1, 'p':2, 'd':5, 'f':11, 'g':21, 'h':36,
               'S':1, 'P':2, 'D':5, 'F':11, 'G':21, 'H':36}

    @staticmethod
    def _section(section_name, section_data):
        """Add opening/closing section_name tags to data."""
        opening_tag = ['<' + section_name + '>']
        closing_tag = ['</' + section_name + '>']

        if isinstance(section_data, list):
            return opening_tag + section_data + closing_tag
        elif isinstance(section_data, str):
            return opening_tag + (' ' + section_data).split('\n') + closing_tag
        elif isinstance(section_data, int) or isinstance(section_data, float):
            return opening_tag + [' ' + str(section_data)] + closing_tag


    @staticmethod
    def _list_format(data, per_line, style='%22.11E'):
        """Format lists for pretty print."""
        template = style * per_line
        leftover = len(data) % per_line
        # Template for last line.
        last_template = style * leftover

        pretty_list = [template % tuple(data[i:i+per_line])
                       for i in range(0, len(data) - leftover, per_line)]
        if leftover:
            return pretty_list + [last_template % tuple(data[-1*leftover:])]
        return  pretty_list

    def _title(self):
        """Return filename without extension to be used as title."""
        title = "Written by cclib-1.5."
        if self.jobfilename is not None:
            return os.path.basename(os.path.splitext(self.jobfilename)[0]) +\
                    '. ' + title
        return title

    def _keywords(self):
        """Return one of GTO, GIAO, CSGT keyword."""
        # AIMALL currently supports on GTO keyword.
        return self.keyword

    def _no_of_nuclei(self):
        """Number of nuclei in the molecule."""
        return len(self.ccdata.atomnos)

    def _no_of_prims(self):
        """Number of Primitives."""
        nprims = 0
        for atom in self.ccdata.gbasis:
            for prims in atom:
                nprims += self.num_orb[prims[0]] * len(prims[1])
        return nprims

    def _no_of_mos(self):
        """Number of occupied MOs."""
        return int(max(self.ccdata.homos)) + 1

    def _no_of_perturbations(self):
        """Number of Perturbation.

        This is usually zero.  For GIAO it should be 3
        (corresponding to Lx, Ly and Lz), and
        for CSGT it should be 6
        (corresponding to Lx, Ly, Lz, Px, Py and Pz).
        """
        if 'GIAO' in self._keywords():
            return 3
        elif 'CSGT' in self._keywords():
            return 6
        return 0

    def _nuclear_names(self):
        """Names of nuclei present in the molecule.

        O1
        H2
        H3
        """
        return [self.pt.element[Z]+str(i) for i, Z in
                enumerate(self.ccdata.atomnos, start=1)]

    def _atomic_nos(self):
        """Atomic numbers of elements."""
        return [str(Z) for Z in self.ccdata.atomnos]

    def _nuclear_charges(self):
        """Nuclear charges."""
        return ['{:22.11E}'.format(Z) for Z in self.ccdata.atomnos]

    def _nuclear_coords(self):
        """Nuclear coordinates in Bohr."""
        coord_template = '%22.11E %.12E %.12E'
        to_bohr = lambda x: utils.convertor(x, 'Angstrom', 'bohr')
        nuc_coords = [coord_template % tuple(to_bohr(coord))
                      for coord in self.ccdata.atomcoords[-1]]
        return nuc_coords

    def _net_charge(self):
        """Net charge on molecule."""
        return '%22.11E' % self.ccdata.charge

    def _no_electrons(self):
        """Number of electrons in molecule."""
        homos = self.ccdata.homos
        charge = self.ccdata.charge if hasattr(self.ccdata, 'charge') else 0.0
        if len(homos) > 1:
            return int(sum(homos)) + 2
        else:
            return int((homos[0] + 1)*2 - (charge % 2))

    def _no_alpha_electrons(self):
        """Number of alpha electrons in molecule."""
        return int(numpy.ceil(self._no_electrons() / 2.0))

    def _no_beta_electrons(self):
        """Number of beta electrons in molecule."""
        return self._no_electrons() - self._no_alpha_electrons()

    def _spin_mult(self):
        """Electronic Spin Multiplicity"""
        return self.ccdata.mult

    def _prim_centers(self):
        """List of nuclear numbers upon which the primitive basis functions
        are centered."""
        prim_centers = []
        for nuc_num, atom in enumerate(self.ccdata.gbasis, start=1):
            for prims in atom:
                prim_centers += [nuc_num] * self.num_orb[prims[0]]\
                                * len(prims[1])

        return self._list_format(prim_centers, 10, '%d ')

    def _rearrange_modata(self, data):
        """Rearrange MO related data."""
        prim_types = self._prim_types()
        if isinstance(data, numpy.ndarray):
            data = data.tolist()

        pos_yyx = [key for key, val in enumerate(prim_types)
                   if val == 17]
        pos_yyz = [key for key, val in enumerate(prim_types)
                   if val == 16]

        if pos_yyx:
            for pos in pos_yyx:
                data.insert(pos-3, data.pop(pos))
        if pos_yyz:
            for pos in pos_yyz:
                data.insert(pos+3, data.pop(pos + 1))

        return data


    def _prim_types(self):
        """Primitive type numbering starts at 1.
        The definition of the Cartesian Gaussian primitive types is:
        1 S, 2 PX, 3 PY, 4 PZ, 5 DXX, 6 DYY, 7 DZZ, 8 DXY, 9 DXZ, 10 DYZ,
        11 FXXX, 12 FYYY, 13 FZZZ, 14 FXXY, 15 FXXZ, 16 FYYZ, 17 FXYY,
        18 FXZZ, 19 FYZZ, 20 FXYZ,
        21 GXXXX, 22 GYYYY, 23 GZZZZ, 24 GXXXY, 25 GXXXZ, 26 GXYYY,
        27 GYYYZ, 28 GXZZZ,
        29 GYZZZ, 30 GXXYY, 31 GXXZZ, 32 GYYZZ, 33 GXXYZ, 34 GXYYZ,
        35 GXYZZ,
        36 HZZZZZ, 37 HYZZZZ, 38 HYYZZZ, 39 HYYYZZ,
        40 HYYYYZ, 41 HYYYYY, 42 HXZZZZ, 43 HXYZZZ, 44 HXYYZZ,
        45 HXYYYZ, 46 HXYYYY, 47 HXXZZZ, 48 HXXYZZ, 49 HXXYYZ,
        50 HXXYYY, 51 HXXXZZ, 52 HXXXYZ, 53 HXXXYY, 54 HXXXXZ, 55 HXXXXY,
        56 HXXXXX
        """

        prim_types = []
        for atom in self.ccdata.gbasis:
            for prims in atom:
                prim_orb = []
                for i in range(self.num_orb[prims[0]]):
                    prim_orb += [(self.idx_orb[prims[0]]  + i)]\
                                * len(prims[1])
                prim_types += prim_orb
        return prim_types

    def _prim_types_section(self):
        prim_types = self._prim_types()
        # GAMESS specific reordering.
        if self.ccdata.metadata['package'] == 'GAMESS':
            prim_types = self._rearrange_modata(prim_types)
        return self._list_format(prim_types, 10, '%d ')

    def _prim_exps(self):
        """Space-separated list of primitive exponents."""
        prim_exps = []
        for atom in self.ccdata.gbasis:
            for prims in atom:
                prim_exps += [prim[0] for prim in prims[1]]\
                            * self.num_orb[prims[0]]
        return self._list_format(prim_exps, 5)

    def _mo_occup_nos(self):
        """MO occupation numbers."""
        occup = []
        electrons = self._no_electrons()
        alpha = self._no_alpha_electrons()
        beta = self._no_beta_electrons()
        if len(self.ccdata.homos) == 1:
            occup += ['{:22.11E}'.format(2)] * int(electrons / 2) +\
                        ['{:22.11E}'.format(1)] * (electrons % 2)
        else:
            occup += ['{:22.11E}'.format(1)] *  + alpha +\
                        ['{:22.11E}'.format(1)] * beta
        return occup

    def _mo_energies(self):
        """List of MO Energies."""
        mo_energies = []
        alpha_elctrons = self._no_alpha_electrons()
        beta_electrons = self._no_beta_electrons()
        for mo_energy in self.ccdata.moenergies[0][:alpha_elctrons]:
            mo_energies.append('{:22.11E}'.format(
                utils.convertor(mo_energy, 'eV', 'hartree')))
        if self.ccdata.mult > 1:
            for mo_energy in self.ccdata.moenergies[1][:beta_electrons]:
                mo_energies.append('{:22.11E}'.format(
                    utils.convertor(mo_energy, 'eV', 'hartree')))
        return mo_energies

    def _mo_spin_types(self):
        """MO Spin Types."""
        spin_types = []
        electrons = self._no_electrons()
        alpha = self._no_alpha_electrons()
        beta = self._no_beta_electrons()
        if len(self.ccdata.homos) == 1:
            spin_types += ['Alpha and Beta'] * int(electrons / 2) +\
                            ['Alpha'] * (electrons % 2)
        else:
            spin_types += ['Alpha'] * alpha +\
                            ['Beta'] * beta
        return spin_types

    pi3 = (2.0 / numpy.pi) ** 3

    # Precomputed values for l+m+n/
    _L = dict(
        [(prim_type, 0) for prim_type in range(1, 2)] +   # s
        [(prim_type, 1) for prim_type in range(2, 5)] +   # p
        [(prim_type, 2) for prim_type in range(5, 11)] +  # d
        [(prim_type, 3) for prim_type in range(11, 21)] + # f
        [(prim_type, 4) for prim_type in range(21, 36)]   # g
    )

    # Precomputed values for ((2l-1)!! * (2m-1)!! * (2n-1)!!).
    # Ref. Molden2AIM
    _M = dict(
        [(L, 1) for L in range(1, 5)] +
        [(L, 9) for L in range(5, 8)] +
        [(L, 1) for L in range(8, 11)] +
        [(L, 225) for L in range(11, 14)] +
        [(L, 9) for L in range(14, 20)] +
        [(L, 1) for L in range(20, 21)] +
        [(L, 11025) for L in range(21, 24)] +
        [(L, 225) for L in range(24, 30)] +
        [(L, 81) for L in range(30, 33)] +
        [(L, 9) for L in range(33, 36)]
    )

    def _normalize(self, prim_type, alpha=1.0):
        """Normalization factor for Cartesian Gaussian Functions.

        N**4 = (2/pi)**3 * 2**(l+m+n) * alpha**(3 + 2(l+m+n)) /
                            ((2l-1)!! * (2m-1)!! * (2n-1)!!)**2
            = (2/pi)**3 * 2**(L) * alpha**(3 + 2L) /
                            M**2,
        L = l+m+n,
        M = ((2l-1)!! * (2m-1)!! * (2n-1)!!)
        """
        L = self._L[prim_type]
        M = self._M[prim_type]
        norm_four = self.pi3 * 2**(4*L) * alpha**(3+2*L) / M
        norm = numpy.power(norm_four, 1/4.0)
        return norm

    def _rearrange_mocoeffs(self, mocoeffs):
        """Rearrange cartesian F functions in mocoeffs.
        Expected order:
        xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
        cclib's order:
        XXX, YYY, ZZZ, XXY, XXZ, YYX, YYZ, ZZX, ZZY, XYZ
        cclib's order can be converted by:
        moving YYX two indexes ahead, and
        moving YYZ two indexes back.
        """

        aonames = self.ccdata.aonames
        mocoeffs = mocoeffs.tolist()

        pos_yyx = [key for key, val in enumerate(aonames)
                   if '_YYX' in val]
        pos_yyz = [key for key, val in enumerate(aonames)
                   if '_YYZ' in val]

        if pos_yyx:
            for pos in pos_yyx:
                mocoeffs.insert(pos-2, mocoeffs.pop(pos))
        if pos_yyz:
            for pos in pos_yyz:
                mocoeffs.insert(pos+2, mocoeffs.pop(pos))

        return mocoeffs

    def _normalized_mocoeffs(self):
        """Raw-Primitive Expansion coefficients for each normalized MO."""
        mo_count = []
        alpha = []
        prim_coeff = []
        prim_mocoeff = []

        for atom in self.ccdata.gbasis:
            for prims in atom:
                prim_orb = []
                mo_count += [len(prims[1])] * self.num_orb[prims[0]]
                for i in range(self.num_orb[prims[0]]):
                    norb = self.idx_orb[prims[0]]
                    prim_orb += [norb + i]
                    alpha += [prim[0] for prim in prims[1]]
                    prim_coeff += [prim[1] for prim in prims[1]]

        # GAMESS specific reordering.
        if self.ccdata.metadata['package'] == 'GAMESS':
            prim_type = self._rearrange_modata(self._prim_types())
            alpha = self._rearrange_modata(alpha)
            prim_coeff = self._rearrange_modata(prim_coeff)

        # Normalization Matrix.
        norm_mat = [self._normalize(prim_type[i], alpha[i]) * prim_coeff[i]
                    for i in range(len(prim_coeff))]

        # Number of molecular orbitals.
        nmos = self._no_electrons()\
                if self.ccdata.mult > 1\
                else self._no_of_mos()

        for i in range(self.ccdata.mult):
            for j in range(nmos):
                mocoeffs = self.ccdata.mocoeffs[i][j]
                if self.ccdata.metadata['package'] == 'GAMESS':
                    mocoeffs = self._rearrange_mocoeffs(self.ccdata.mocoeffs[i][j])
                for k, mocoeff in enumerate(mocoeffs):
                    prim_mocoeff += [mocoeff] * mo_count[k]

        norm_mocoeffs = []
        for mo_num in range(nmos):
            norm_mocoeffs.append([norm_mat[i] *
                                  prim_mocoeff[i + mo_num * len(prim_coeff)]
                                  for i in range(len(prim_coeff))])

        return norm_mocoeffs

    def _mo_prim_coeffs(self):
        # Normalized MO Coeffs.
        norm_mocoeffs = self._normalized_mocoeffs()
        mocoeffs_section = []

        for mo_num, mocoeffs in enumerate(norm_mocoeffs):
            mocoeffs_section.extend(self._section('MO Number', mo_num + 1))
            mocoeffs_section.extend(self._list_format
                                    (mocoeffs, 5))
        return mocoeffs_section

    def _energy(self):
        """The total energy of the molecule.
        HF and KSDFT: SCF energy        (scfenergies),
        MP2         : MP2 total energy  (mpenergies),
        CCSD        : CCSD total energy (ccenergies).
        """
        energy = 0
        if hasattr(self.ccdata, 'ccenergies'):
            energy = self.ccdata.ccenergies[-1]
        elif hasattr(self.ccdata, 'mpenergies'):
            energy = self.ccdata.mpenergies[-1][-1]
        elif hasattr(self.ccdata, 'scfenergies'):
            energy = self.ccdata.scfenergies[-1]
        else:
            raise filewriter.MissingAttributeError(
                'scfenergies/mpenergies/ccenergies')
        return '{:22.11E}'.format(utils.convertor(energy, 'eV', 'hartree'))

    def _virial_ratio(self):
        """Ratio of kinetic energy to potential energy."""
        # Hardcoding expected value for Required Field.
        return '{:22.11E}'.format(2.0)

    def generate_repr(self):
        """Generate the wfx representation of the logfile data."""

        sections = [
            (self._title, "Title", True),
            (self._keywords, "Keywords", True),
            (self._no_of_nuclei, "Number of Nuclei", True),
            (self._no_of_prims, "Number of Primitives", True),
            (self._no_of_mos, "Number of Occupied Molecular Orbitals", True),
            (self._no_of_perturbations, "Number of Perturbations", True),
            (self._nuclear_names, "Nuclear Names", True),
            (self._atomic_nos, "Atomic Numbers", True),
            (self._nuclear_charges, "Nuclear Charges", True),
            (self._nuclear_coords, "Nuclear Cartesian Coordinates", True),
            (self._net_charge, "Net Charge", True),
            (self._no_electrons, "Number of Electrons", True),
            (self._no_alpha_electrons, "Number of Alpha Electrons", True),
            (self._no_beta_electrons, "Number of Beta Electrons", True),
            (self._spin_mult, "Electronic Spin Multiplicity", False),
            # (self._model, "Model", False),
            (self._prim_centers, "Primitive Centers", True),
            (self._prim_types_section, "Primitive Types", True),
            (self._prim_exps, "Primitive Exponents", True),
            (self._mo_occup_nos,
             "Molecular Orbital Occupation Numbers", True),
            (self._mo_energies, "Molecular Orbital Energies", True),
            (self._mo_spin_types, "Molecular Orbital Spin Types", True),
            (self._mo_prim_coeffs,
             "Molecular Orbital Primitive Coefficients", True),
            (self._energy, "Energy = T + Vne + Vee + Vnn", True),
            # (self._nuc_energy_gradients,
            #  "Nuclear Cartesian Energy Gradients", False),
            # (self._nuc_virial,
            #  "Nuclear Virial of Energy-Gradient-Based Forces on Nuclei, W",
            #  False),
            (self._virial_ratio, "Virial Ratio (-V/T)", True),
        ]

        wfx_lines = []

        for section_module, section_name, section_required in sections:
            try:
                section_data = section_module()
                wfx_lines.extend(self._section(section_name, section_data))
            except:
                if section_required:
                    raise filewriter.MissingAttributeError(
                        'Unable to write required wfx section: '
                        + section_name)

        wfx_lines.append('')
        return '\n'.join(wfx_lines)


if __name__ == "__main__":
    pass
