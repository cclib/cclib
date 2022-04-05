# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A writer for wfx format files."""

import os.path
import numpy

from cclib.io import filewriter
from cclib.parser import utils


# Number of orbitals of type key.
# There are 3 p type, 6 d type orbitals etc.
ORBITAL_COUNT = {'S':1, 'P':3, 'D':6, 'F':10, 'G':15, 'H':21}

# Index of first orbital of type key in a list of orbitals.
# The first s orbital has index 1, first p orbital has index 2, and first d
# has index 5.
ORBITAL_INDICES = {'S': 1}
ORBITAL_NAMES = 'SPDFGH'
for idx, name in enumerate(ORBITAL_NAMES[1:], start=1):
    prev_orbital_name = ORBITAL_NAMES[idx - 1]
    prev_orbital_count = ORBITAL_COUNT[prev_orbital_name]
    prev_orbital_index = ORBITAL_INDICES[prev_orbital_name]
    ORBITAL_INDICES[name] = prev_orbital_count + prev_orbital_index

PI_CUBE_INV = (2.0 / numpy.pi) ** 3

# Float formatting template.
WFX_FIELD_FMT = '%22.11E'

# Precomputed values for l+m+n to be used in MO normalization.
_L = dict(
    [(prim_type, 0) for prim_type in range(1, 2)] +   # s
    [(prim_type, 1) for prim_type in range(2, 5)] +   # p
    [(prim_type, 2) for prim_type in range(5, 11)] +  # d
    [(prim_type, 3) for prim_type in range(11, 21)] + # f
    [(prim_type, 4) for prim_type in range(21, 36)]   # g
)

# Precomputed values for ((2l-1)!! * (2m-1)!! * (2n-1)!!).
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


def _section(section_name, section_data):
    """Add opening/closing section_name tags to data."""
    opening_tag = [f"<{section_name}>"]
    closing_tag = [f"</{section_name}>"]

    section = None
    if isinstance(section_data, list):
        section = opening_tag + section_data + closing_tag
    elif isinstance(section_data, str):
        section = opening_tag + f" {section_data}".split("\n") + closing_tag
    elif isinstance(section_data, int) or isinstance(section_data, float):
        section = opening_tag + [f" {str(section_data)}"] + closing_tag
    return section


def _list_format(data, per_line, style=WFX_FIELD_FMT):
    """Format lists for pretty print."""
    template = style * per_line
    leftover = len(data) % per_line
    # Template for last line.
    last_template = style * leftover

    pretty_list = [
        template % tuple(data[i : i + per_line])
        for i in range(0, len(data) - leftover, per_line)
    ]
    if leftover:
        return pretty_list + [last_template % tuple(data[-1 * leftover :])]
    return pretty_list


class WFXWriter(filewriter.Writer):
    """A writer for wfx files."""

    required_attrs = ('natom', 'atomcoords', 'atomnos', 'gbasis', 'charge',
                       'homos', 'mult', 'mocoeffs')

    def _title(self):
        """Section: Title
        Return filename without extension to be used as title."""
        title = "Written by cclib."
        if self.jobfilename is not None:
            return f"{os.path.basename(os.path.splitext(self.jobfilename)[0])}. {title}"
        return title

    def _keywords(self):
        """Section: Keywords.
        Return one of GTO, GIAO, CSGT keyword."""
        # Currently only GTO is supported.
        return 'GTO'

    def _no_of_nuclei(self):
        """Section: Number of Nuclei."""
        return self.ccdata.natom

    def _no_of_prims(self):
        """Section: Number of Primitives."""
        nprims = 0
        for atom in self.ccdata.gbasis:
            for prims in atom:
                nprims += ORBITAL_COUNT[prims[0]] * len(prims[1])
        return nprims

    def _no_of_mos(self):
        """Section: Number of Occupied Molecular Orbitals."""
        return int(max(self.ccdata.homos)) + 1

    def _no_of_perturbations(self):
        """Section: Number of Perturbation.

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
        """Section: Nuclear Names.
        Names of nuclei present in the molecule.

        O1
        H2
        H3
        """
        return [self.pt.element[Z]+str(i) for i, Z in
                enumerate(self.ccdata.atomnos, start=1)]

    def _atomic_nos(self):
        """Section: Atomic Numbers."""
        return [str(Z) for Z in self.ccdata.atomnos]

    def _nuclear_charges(self):
        """Section: Nuclear Charges."""
        nuclear_charge = [WFX_FIELD_FMT % Z for Z in self.ccdata.atomnos]
        if hasattr(self.ccdata, "coreelectrons"):
            nuclear_charge = [
                WFX_FIELD_FMT % Z
                for Z in self.ccdata.atomnos - self.ccdata.coreelectrons
            ]
        return nuclear_charge

    def _nuclear_coords(self):
        """Section: Nuclear Cartesian Coordinates.
        Nuclear coordinates in Bohr."""
        coord_template = WFX_FIELD_FMT * 3
        to_bohr = lambda x: utils.convertor(x, "Angstrom", "bohr")
        nuc_coords = [
            coord_template % tuple(to_bohr(coord))
            for coord in self.ccdata.atomcoords[-1]
        ]
        return nuc_coords

    def _net_charge(self):
        """Section: Net Charge.
        Net charge on molecule."""
        return WFX_FIELD_FMT % self.ccdata.charge

    def _no_electrons(self):
        """Section: Number of Electrons."""
        return int(self.ccdata.nelectrons)

    def _no_alpha_electrons(self):
        """Section: Number of Alpha Electrons."""
        no_electrons = numpy.sum(self.ccdata.atomnos - self.ccdata.coreelectrons) - self.ccdata.charge
        no_alpha = (no_electrons + (self.ccdata.mult - 1))//2
        return int(no_alpha)

    def _no_beta_electrons(self):
        """Section: Number of Beta Electrons."""
        return int(self.ccdata.nelectrons - self._no_alpha_electrons())

    def _spin_mult(self):
        """Section: Electronic Spin Multiplicity"""
        return self.ccdata.mult

    def _prim_centers(self):
        """Section: Primitive Centers.
        List of nuclear numbers upon which the primitive basis functions
        are centered."""
        prim_centers = []
        for nuc_num, atom in enumerate(self.ccdata.gbasis, start=1):
            for prims in atom:
                prim_centers += [nuc_num] * ORBITAL_COUNT[prims[0]]\
                                * len(prims[1])

        return _list_format(prim_centers, 10, "%d ")

    def _rearrange_modata(self, data):
        """Rearranges MO related data according the expected order of
        Cartesian gaussian primitive types in wfx format.
        cclib parses mocoeffs in the order they occur in output files.
        """
        prim_types = self._get_prim_types()
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


    def _get_prim_types(self):
        """List of primitive types.
        Definition of the Cartesian Gaussian primitive types is as follows:
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
        Spherical basis are not currently supported by the writer.
        """
        prim_types = []
        for atom in self.ccdata.gbasis:
            for prims in atom:
                prim_orb = []
                for i in range(ORBITAL_COUNT[prims[0]]):
                    prim_orb += [(ORBITAL_INDICES[prims[0]]  + i)]\
                                * len(prims[1])
                prim_types += prim_orb
        return prim_types

    def _prim_types(self):
        """Section: Primitive Types."""
        prim_types = self._get_prim_types()
        # GAMESS specific reordering.
        if self.ccdata.metadata['package'] == 'GAMESS':
            prim_types = self._rearrange_modata(prim_types)
        return _list_format(prim_types, 10, "%d ")

    def _prim_exps(self):
        """Section: Primitive Exponents.
        Space-separated list of primitive exponents."""
        prim_exps = []
        for atom in self.ccdata.gbasis:
            for prims in atom:
                prim_exps += [prim[0] for prim in prims[1]]\
                            * ORBITAL_COUNT[prims[0]]
        return _list_format(prim_exps, 5)

    def _mo_occup_nos(self):
        """Section: Molecular Orbital Occupation Numbers."""
        occup = []
        electrons = self._no_electrons()
        alpha = self._no_alpha_electrons()
        beta = self._no_beta_electrons()
        if len(self.ccdata.homos) == 1:
            occup += [WFX_FIELD_FMT % (2)] * int(electrons / 2) + [
                WFX_FIELD_FMT % (1)
            ] * (electrons % 2)
        else:
            occup += [WFX_FIELD_FMT % (1)] * +alpha + [WFX_FIELD_FMT % (1)] * beta
        return occup

    def _mo_energies(self):
        """Section: Molecular Orbital Energies."""
        mo_energies = []
        alpha_elctrons = self._no_alpha_electrons()
        beta_electrons = self._no_beta_electrons()
        for mo_energy in self.ccdata.moenergies[0][:alpha_elctrons]:
            mo_energies.append(
                WFX_FIELD_FMT % (utils.convertor(mo_energy, "eV", "hartree"))
            )
        if self.ccdata.mult > 1:
            for mo_energy in self.ccdata.moenergies[1][:beta_electrons]:
                mo_energies.append(
                    WFX_FIELD_FMT % (utils.convertor(mo_energy, "eV", "hartree"))
                )
        return mo_energies

    def _mo_spin_types(self):
        """Section: Molecular Orbital Spin Types."""
        spin_types = []
        electrons = self._no_electrons()
        alpha = self._no_alpha_electrons()
        beta = self._no_beta_electrons()
        if len(self.ccdata.homos) == 1:
            spin_types += ["Alpha and Beta"] * int(electrons / 2) + ["Alpha"] * (
                electrons % 2
            )
        else:
            spin_types += ['Alpha'] * alpha +\
                            ['Beta'] * beta
        return spin_types

    def _normalize(self, prim_type, alpha=1.0):
        """Normalization factor for Cartesian Gaussian Functions.

        N**4 = (2/pi)**3 * 2**(l+m+n) * alpha**(3 + 2(l+m+n)) /
                            ((2l-1)!! * (2m-1)!! * (2n-1)!!)**2
            = (2/pi)**3 * 2**(L) * alpha**(3 + 2L) /
                            M**2,
        L = l+m+n,
        M = ((2l-1)!! * (2m-1)!! * (2n-1)!!)
        """
        L = _L[prim_type]
        M = _M[prim_type]
        norm_four = PI_CUBE_INV * 2**(4*L) * alpha**(3+2*L) / M
        norm = numpy.power(norm_four, 1/4.0)
        return norm

    def _rearrange_mocoeffs(self, mocoeffs):
        """Rearrange cartesian F functions in mocoeffs.
        Expected order:
        xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
        cclib's order for GAMESS:
        XXX, YYY, ZZZ, XXY, XXZ, YYX, YYZ, ZZX, ZZY, XYZ
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

    def _norm_mat(self):
        """Calculate normalization matrix for normalizing MOcoeffs."""
        alpha = []
        prim_coeff = []
        mo_count = []
        prim_type = self._get_prim_types()

        for atom in self.ccdata.gbasis:
            for prims in atom:
                prim_orb = []
                mo_count += [len(prims[1])] * ORBITAL_COUNT[prims[0]]
                for i in range(ORBITAL_COUNT[prims[0]]):
                    norb = ORBITAL_INDICES[prims[0]]
                    prim_orb += [norb + i]
                    alpha += [prim[0] for prim in prims[1]]
                    prim_coeff += [prim[1] for prim in prims[1]]

        # GAMESS specific reordering.
        if self.ccdata.metadata['package'] == 'GAMESS':
            prim_type = self._rearrange_modata(self._get_prim_types())
            alpha = self._rearrange_modata(alpha)
            prim_coeff = self._rearrange_modata(prim_coeff)

        norm_mat = [self._normalize(prim_type[i], alpha[i]) * prim_coeff[i]
                    for i in range(len(prim_coeff))]

        return (norm_mat, mo_count, prim_coeff)

    def _nmos(self):
        """Return number of molecular orbitals to be printed."""

        return self.ccdata.nelectrons if self.ccdata.mult > 1\
                else self._no_of_mos()

    def _prim_mocoeff(self, mo_count):
        """Return primitve mocoeffs array."""
        prim_mocoeff = []

        for i in range(len(self.ccdata.mocoeffs)):
            for j in range(self._nmos()):
                mocoeffs = self.ccdata.mocoeffs[i][j]
                if self.ccdata.metadata['package'] == 'GAMESS':
                    mocoeffs = self._rearrange_mocoeffs(self.ccdata.mocoeffs[i][j])
                for k, mocoeff in enumerate(mocoeffs):
                    prim_mocoeff += [mocoeff] * mo_count[k]

        return prim_mocoeff

    def _normalized_mocoeffs(self):
        """Raw-Primitive Expansion coefficients for each normalized MO."""
        # Normalization Matrix.
        norm_mat, mo_count, prim_coeff = self._norm_mat()

        prim_mocoeff = self._prim_mocoeff(mo_count)

        norm_mocoeffs = []
        for mo_num in range(self._nmos()):
            norm_mocoeffs.append([norm_mat[i] *
                                  prim_mocoeff[i + mo_num * len(prim_coeff)]
                                  for i in range(len(prim_coeff))])

        return norm_mocoeffs

    def _mo_prim_coeffs(self):
        """Section: Molecular Orbital Primitive Coefficients."""
        # Normalized MO Coeffs.
        norm_mocoeffs = self._normalized_mocoeffs()
        mocoeffs_section = []

        for mo_num, mocoeffs in enumerate(norm_mocoeffs):
            mocoeffs_section.extend(_section('MO Number', mo_num + 1))
            mocoeffs_section.extend(_list_format
                                    (mocoeffs, 5))
        return mocoeffs_section

    def _energy(self):
        """Section: Energy = T + Vne + Vee + Vnn.
        The total energy of the molecule.
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
            raise filewriter.MissingAttributeError("scfenergies/mpenergies/ccenergies")
        return WFX_FIELD_FMT % (utils.convertor(energy, "eV", "hartree"))

    def _virial_ratio(self):
        """Ratio of kinetic energy to potential energy."""
        # Hardcoding expected value for Required Field.
        return WFX_FIELD_FMT % (2.0)

    def generate_repr(self):
        """Generate the wfx representation of the logfile data."""

        # sections:(Function returning data for section,
        #           Section heading,
        #           Required)

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
            (self._prim_types, "Primitive Types", True),
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
                wfx_lines.extend(_section(section_name, section_data))
            except:
                if section_required:
                    raise filewriter.MissingAttributeError(
                        f"Unable to write required wfx section: {section_name}"
                    )

        wfx_lines.append('')
        return '\n'.join(wfx_lines)
