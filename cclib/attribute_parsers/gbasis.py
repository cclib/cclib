# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np
from typing import Optional


class gbasis(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["psi4", "qchem"]

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        dependency_list = ["natom", "atomnos"]
        line = file_handler.last_line
        if line.strip() == "==> AO Basis Functions <==":
            if base_parser.check_dependencies(dependency_list, ccdata, "gbasis"):

                def get_symmetry_atom_basis(parsed_gbasis):
                    """Get symmetry atom by replicating the last atom in gbasis of the same element."""
                    missing_index = len(parsed_gbasis)
                    missing_atomno = ccdata.atomnos[missing_index]
                    ngbasis = len(parsed_gbasis)
                    last_same = ngbasis - ccdata.atomnos[:ngbasis][::-1].index(missing_atomno) - 1
                    return parsed_gbasis[last_same]

                dfact = lambda n: (n <= 0) or n * dfact(n - 2)

                # Early beta versions of Psi4 normalize basis function
                # coefficients when printing.
                version_4_beta = False  # don't know what to do with this yet. will have to add to ccdata maybe (psi4 specific metadata?).
                if version_4_beta:

                    def get_normalization_factor(exp, lx, ly, lz):
                        norm_s = (2 * exp / numpy.pi) ** 0.75
                        if lx + ly + lz > 0:
                            nom = (4 * exp) ** ((lx + ly + lz) / 2.0)
                            den = numpy.sqrt(
                                dfact(2 * lx - 1) * dfact(2 * ly - 1) * dfact(2 * lz - 1)
                            )
                            return norm_s * nom / den
                        else:
                            return norm_s
                else:
                    get_normalization_factor = lambda exp, lx, ly, lz: 1

                file_handler.skip_lines(["b", "basisname"], virtual=True)
                line = file_handler.virtual_next()
                parsed_gbasis = []
                file_handler.skip_lines(["stars"], virtual=True)
                line = file_handler.virtual_next()
                while line.strip():
                    element, index = line.split()
                    if len(element) > 1:
                        element = element[0] + element[1:].lower()
                    index = int(index)
                    # This is the code that adds missing atoms when symmetry atoms are excluded
                    # from the basis set printout. Again, this will work only if all atoms of
                    # the same element use the same basis set.
                    while index > len(parsed_gbasis) + 1:
                        parsed_gbasis.append(get_symmetry_atom_basis(parsed_gbasis))

                    parsed_gbasis.append([])
                    line = file_handler.virtual_next()
                    while line.find("*") == -1:
                        # The shell type and primitive count is in the first line.
                        shell_type, nprimitives, _ = line.split()
                        nprimitives = int(nprimitives)

                        # Get the angular momentum for this shell type.
                        momentum = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4, "H": 5, "I": 6}[
                            shell_type.upper()
                        ]

                        # Read in the primitives.
                        primitives_lines = [file_handler.virtual_next() for i in range(nprimitives)]
                        primitives = [list(map(float, pl.split())) for pl in primitives_lines]

                        # Un-normalize the coefficients. Psi prints the normalized coefficient
                        # of the highest polynomial, namely XX for D orbitals, XXX for F, and so on.
                        for iprim, prim in enumerate(primitives):
                            exp, coef = prim
                            coef = coef / get_normalization_factor(exp, momentum, 0, 0)
                            primitives[iprim] = [exp, coef]

                        primitives = [tuple(p) for p in primitives]
                        shell = [shell_type, primitives]
                        parsed_gbasis[-1].append(shell)

                        line = file_handler.virtual_next()

                    line = file_handler.virtual_next()

                # We will also need to add symmetry atoms that are missing from the input
                # at the end of this block, if the symmetry atoms are last.
                while len(parsed_gbasis) < ccdata.natom:
                    parsed_gbasis.append(get_symmetry_atom_basis(parsed_gbasis))

                constructed_data = {gbasis.__name__: parsed_gbasis}
                return constructed_data
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> dict | None:
        constructed_data = None
        line = file_handler.last_line
        parsed_gbasis = []
        # Parse the general basis for `gbasis`, in the style used by Gaussian.
        if "Basis set in general basis input format:" in line:
            file_handler.skip_lines(["d", "$basis"], virtual=True)
            line = file_handler.virtual_next()
            # The end of the general basis block.
            while "$end" not in line:
                atom = []
                # 1. Contains element symbol and atomic index of
                # basis functions; if 0, applies to all atoms of
                # same element.
                assert len(line.split()) == 2
                line = file_handler.virtual_next()
                # The end of each atomic block.
                while "****" not in line:
                    # 2. Contains the type of basis function {S, SP,
                    # P, D, F, G, H, ...}, the number of primitives,
                    # and the weight of the final contracted function.
                    bfsplitline = line.split()
                    assert len(bfsplitline) == 3
                    bftype = bfsplitline[0]
                    nprim = int(bfsplitline[1])
                    line = file_handler.virtual_next()
                    # 3. The primitive basis functions that compose
                    # the contracted basis function; there are `nprim`
                    # of them. The first value is the exponent, and
                    # the second value is the contraction
                    # coefficient. If `bftype == 'SP'`, the primitives
                    # are for both S- and P-type basis functions but
                    # with separate contraction coefficients,
                    # resulting in three columns.
                    if bftype == "SP":
                        primitives_S = []
                        primitives_P = []
                    else:
                        primitives = []
                    # For each primitive in the contracted basis
                    # function...
                    for iprim in range(nprim):
                        primsplitline = line.split()
                        exponent = float(primsplitline[0])
                        if bftype == "SP":
                            assert len(primsplitline) == 3
                            coefficient_S = float(primsplitline[1])
                            coefficient_P = float(primsplitline[2])
                            primitives_S.append((exponent, coefficient_S))
                            primitives_P.append((exponent, coefficient_P))
                        else:
                            assert len(primsplitline) == 2
                            coefficient = float(primsplitline[1])
                            primitives.append((exponent, coefficient))
                        line = file_handler.virtual_next()
                    if bftype == "SP":
                        bf_S = ("S", primitives_S)
                        bf_P = ("P", primitives_P)
                        atom.append(bf_S)
                        atom.append(bf_P)
                    else:
                        bf = (bftype, primitives)
                        atom.append(bf)
                    # Move to the next contracted basis function
                    # as long as we don't hit the '****' atom
                    # delimiter.
                parsed_gbasis.append(atom)
                line = file_handler.virtual_next()
        if len(parsed_gbasis) > 0:
            constructed_data = {gbasis.__name__: parsed_gbasis}
        return constructed_data

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in gbasis.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(gbasis, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
