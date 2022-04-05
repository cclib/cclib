# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Formatted Checkpoint files"""


from __future__ import print_function

import re

import numpy

from cclib.parser import data
from cclib.parser import logfileparser
from cclib.parser import utils

SHELL_ORBITALS = {
    0: ['S'],
    1: ['PX', 'PY', 'PZ'],
    -1: ['S', 'PX', 'PY', 'PZ'],
    2: ['D1', 'D2', 'D3', 'D4', 'D5', 'D6'],
    -2: ['D1', 'D2', 'D3', 'D4', 'D5'],
    3:  ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10'],
    -3: ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7'],
    4: ['G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12','G13'],
    -4: ['G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9']

}

SHELL_START = {
    0: 1,
    1: 2,
    -1: 2,
    2: 3,
    -2: 3,
    3: 4,
    -3: 4
}


def _shell_to_orbitals(type, offset):
    """Convert a Fchk shell type and offset to a list of string representations.

    For example, shell type = -2 corresponds to d orbitals (spherical basis) with
    an offset = 1 would correspond to the 4d orbitals, so this function returns
    `['4D1', '4D2', '4D3', '4D4', '4D5']`.
    """

    return [f"{SHELL_START[type] + offset}{x}" for x in SHELL_ORBITALS[type]]


class FChk(logfileparser.Logfile):
    """A Formatted checkpoint file, which contains molecular and wavefunction information.

    These files are produced by Gaussian and Q-Chem.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(logname="FChk", *args, **kwargs)
        self.start = True

    def __str__(self):
        """Return a string representation of the object."""
        return f"Formatted checkpoint file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'FCHK("{self.filename}")'

    def normalisesym(self, symlabel):
        """Just return label"""
        return symlabel

    def extract(self, inputfile, line):

        # just opened file, skip first line to get basis
        if self.start:
            method = next(inputfile)
            self.metadata['basis_set'] = method.split()[-1]
            self.start = False

        if line[0:6] == 'Charge':
            self.set_attribute('charge', int(line.split()[-1]))

        if line[0:12] == 'Multiplicity':
            self.set_attribute('mult', int(line.split()[-1]))

        if line[0:14] == 'Atomic numbers':
            self.natom = int(line.split()[-1])
            atomnos = self._parse_block(inputfile, self.natom, int, 'Basic Information')
            self.set_attribute('atomnos', atomnos)

        if line[0:19] == 'Number of electrons':
            alpha = next(inputfile)
            alpha_homo = int(alpha.split()[-1]) - 1

            beta = next(inputfile)
            beta_homo = int(beta.split()[-1]) - 1

            self.set_attribute('homos', [alpha_homo, beta_homo])

        if line[0:29] == 'Current cartesian coordinates':
            count = int(line.split()[-1])
            assert count % 3 == 0

            coords = numpy.array(self._parse_block(inputfile, count, float, 'Coordinates'))
            coords.shape = (1, int(count / 3), 3)
            self.set_attribute('atomcoords', utils.convertor(coords, 'bohr', 'Angstrom'))

        if line[0:25] == 'Number of basis functions':
            self.set_attribute('nbasis', int(line.split()[-1]))

        if line[0:14] == 'Overlap Matrix':
            count = int(line.split()[-1])

            # triangle matrix, with number of elements in a row:
            # 1 + 2 + 3 + .... + self.nbasis
            assert count == (self.nbasis + 1) * self.nbasis / 2
            raw_overlaps = self._parse_block(inputfile, count, float, 'Overlap Matrix')

            # now turn into matrix
            overlaps = numpy.zeros((self.nbasis, self.nbasis))
            raw_index = 0
            for row in range(self.nbasis):
                for col in range(row + 1):
                    overlaps[row, col] = raw_overlaps[raw_index]
                    overlaps[col, row] = raw_overlaps[raw_index]
                    raw_index += 1

            self.set_attribute('aooverlaps', overlaps)

        if line[0:31] == 'Number of independent functions':
            self.set_attribute('nmo', int(line.split()[-1]))

        if line[0:21] == 'Alpha MO coefficients':
            count = int(line.split()[-1])
            assert count == self.nbasis * self.nmo

            coeffs = numpy.array(self._parse_block(inputfile, count, float, 'Alpha Coefficients'))
            coeffs.shape = (self.nmo, self.nbasis)
            self.set_attribute('mocoeffs', [coeffs])

        if line[0:22] == 'Alpha Orbital Energies':
            count = int(line.split()[-1])
            assert count == self.nmo

            energies = numpy.array(self._parse_block(inputfile, count, float, 'Alpha MO Energies'))
            self.set_attribute('moenergies', [utils.convertor(energies, 'hartree', 'eV')])

        if line[0:20] == 'Beta MO coefficients':
            count = int(line.split()[-1])
            assert count == self.nbasis * self.nmo

            coeffs = numpy.array(self._parse_block(inputfile, count, float, 'Beta Coefficients'))
            coeffs.shape = (self.nmo, self.nbasis)
            self.append_attribute('mocoeffs', coeffs)

        if line[0:21] == 'Beta Orbital Energies':
            count = int(line.split()[-1])
            assert count == self.nmo

            energies = numpy.array(self._parse_block(inputfile, count, float, 'Alpha MO Energies'))
            self.append_attribute('moenergies', utils.convertor(energies, 'hartree', 'eV'))

        if line[0:11] == 'Shell types':
            self.parse_aonames(line, inputfile)

        if line[0:19] == 'Real atomic weights':
            count = int(line.split()[-1])
            assert count == self.natom

            atommasses = numpy.array(self._parse_block(inputfile, count, float, 'Atomic Masses'))

            self.set_attribute('atommasses', atommasses)

        if line[0:10] == 'SCF Energy':
            self.scfenergy = float(line.split()[-1])

            self.set_attribute('scfenergies', [utils.convertor(self.scfenergy,'hartree','eV')])

        if line[0:18] == 'Cartesian Gradient':
            count = int(line.split()[-1])
            assert count == self.natom*3

            gradient = numpy.array(self._parse_block(inputfile, count, float, 'Gradient'))

            self.set_attribute('grads', gradient)

        if line[0:25] == 'Cartesian Force Constants':
            count = int(line.split()[-1])
            assert count == (3*self.natom*(3*self.natom+1))/2

            hessian = numpy.array(self._parse_block(inputfile, count, float, 'Hessian'))

            self.set_attribute('hessian', hessian)

        if line[0:13] == 'ETran scalars':
            count = int(line.split()[-1])

            etscalars = self._parse_block(inputfile, count, int, 'ET Scalars')

            # Set attribute: self.netroot (number of excited estates)
            self.netroot = etscalars[4]

        if line[0:10] == 'ETran spin':
            count = int(line.split()[-1])

            etspin = self._parse_block(inputfile, count, int, 'ET Spin')

            spin_labels = { 0:'Singlet',
                            2:'Triplet',
                           -1:'Unknown'}
            etsyms = []
            for i in etspin:
                if i in spin_labels:
                    etsyms.append(spin_labels[i])
                else:
                    etsyms.append(spin_labels[-1])

            # The extracted property does not contain the actual irrep label
            # (contrarily to that extracted from the Gaussian log)
            # After this, 'Etran sym' appears (and would need to be parsed), 
            # but at least in Gaussian this contains only zeroes regardless of the irrep.

            self.set_attribute('etsyms', etsyms)

        if line[0:18] == 'ETran state values':
            # This section is organized as follows:
            # ·First the properties of each excited state (up to net):
            # E, {muNx,muNy,muNz,muvelNx,muvelNy,muvelNz,mmagNx,mmagNy,mmagNz,unkX,unkY,unkZ,unkX,unkY,unkZ}_N=1,net
            # ·Then come 48 items (only if Freq is requested)
            # They were all 0.000 in G09, but get an actual value in G16
            # ·Then, the derivates of each property with respect to Cartesian coordiates only for target state (netroot)
            # For each Cartesian coordiate, all derivatives wrt to it are included:
            #  dE/dx1 dmux/dx1 dmuy/dx1 ... unkZ/dx1
            #  dE/dy1 dmux/dy1 dmuy/dy1 ... unkZ/dy1
            #  ...
            #  dE/dzN dmux/dzN dmuy/dzN ... unkZ/dzN
            # The number of items is therefore:
            ### 16*net (no Freq jobs)
            ### 16*net + 48 + 3*self.natom*16 (Freq jobs)
            count = int(line.split()[-1])
            if hasattr(self,'etsyms'):
                net = len(self.etsyms)
            else:
                net = 0 # This forces an AssertionError below
            assert count in [16*net, 16*net+48+3*self.natom*16]

            etvalues = self._parse_block(inputfile, count, float, 'ET Values')

            # ETr energies (1/cm)
            etenergies_au = [ e_es-self.scfenergy for e_es in etvalues[0:net*16:16] ]
            etenergies = [ utils.convertor(etr,'hartree','wavenumber') for etr in etenergies_au ]
            self.set_attribute('etenergies', etenergies)

            # ETr dipoles (length-gauge)
            etdips = []
            for k in range(1,16*net,16):
                etdips.append(etvalues[k:k+3])
            self.set_attribute('etdips',etdips)

            # Osc. Strength from Etr dipoles
            # oscs = 2/3 * Etr(au) * dip²(au)
            etoscs = [ 2.0/3.0*e*numpy.linalg.norm(numpy.array(dip))**2 for e,dip in zip(etenergies_au,etdips) ]
            self.set_attribute('etoscs', etoscs)

            # ETr dipoles (velocity-gauge)
            etveldips = []
            for k in range(4,16*net,16):
                etveldips.append(etvalues[k:k+3])
            self.set_attribute('etveldips',etveldips)

            # ETr magnetic dipoles
            etmagdips = []
            for k in range(7,16*net,16):
                etmagdips.append(etvalues[k:k+3])
            self.set_attribute('etmagdips',etmagdips)

    def parse_aonames(self, line, inputfile):
        # e.g.: Shell types                                I   N=          28
        count = int(line.split()[-1])
        shell_types = self._parse_block(inputfile, count, int, 'Atomic Orbital Names')

        # e.g.: Number of primitives per shell             I   N=          28
        next(inputfile)
        self._parse_block(inputfile, count, int, 'Atomic Orbital Names')

        # e.g. Shell to atom map                          I   N=          28
        next(inputfile)
        shell_map = self._parse_block(inputfile, count, int, 'Atomic Orbital Names')

        elements = (self.table.element[x] for x in self.atomnos)
        atom_labels = [f"{y}{x}" for x, y in enumerate(elements, 1)]

        # get orbitals for first atom and start aonames and atombasis lists
        atom = shell_map[0] - 1
        shell_offset = 0
        orbitals = _shell_to_orbitals(shell_types[0], shell_offset)
        aonames = [f"{atom_labels[atom]}_{x}" for x in orbitals]
        atombasis = [list(range(len(orbitals)))]

        # get rest
        for i in range(1, len(shell_types)):
            _type = shell_types[i]
            atom = shell_map[i] - 1
            shell_offset += 1
            basis_offset = atombasis[-1][-1] + 1 # atombasis is increasing numbers, so just grab last

            # if we've move to next atom, need to update offset of shells (e.g. start at 1S)
            # and start new list for atom basis
            if atom != shell_map[i - 1] - 1:
                shell_offset = 0
                atombasis.append([])

            # determine if we've changed shell type (e.g. from S to P)
            if _type != shell_types[i - 1]:
                shell_offset = 0

            orbitals = _shell_to_orbitals(_type, shell_offset)
            aonames.extend([f"{atom_labels[atom]}_{x}" for x in orbitals])
            atombasis[-1].extend(list(range(basis_offset, basis_offset + len(orbitals))))

        assert (
            len(aonames) == self.nbasis
        ), f"Length of aonames != nbasis: {len(aonames)} != {self.nbasis}"
        self.set_attribute("aonames", aonames)

        assert (
            len(atombasis) == self.natom
        ), f"Length of atombasis != natom: {len(atombasis)} != {self.natom}"
        self.set_attribute("atombasis", atombasis)

    def after_parsing(self):
        """Correct data or do parser-specific validation after parsing is finished."""

        # If restricted calculation, need to remove beta homo
        if len(self.moenergies) == len(self.homos) - 1:
            self.homos.pop()

    def _parse_block(self, inputfile, count, type, msg):
        atomnos = []
        while len(atomnos) < count :
            self.updateprogress(inputfile, msg, self.fupdate)
            line = next(inputfile)
            atomnos.extend([ type(x) for x in line.split()])
        return atomnos
