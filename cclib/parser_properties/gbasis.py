from cclib.parser_properties import utils
from cclib.parser_properties.base_parser import base_parser

import numpy as np


class gbasis(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["psi4"]

    @staticmethod
    def psi4(file_handler, ccdata) -> list | None:
        dependency_list = ["natom", "atomnos"]
        line = file_handler.last_line
        if line.strip() == "==> AO Basis Functions <==":
            if base_parser.check_dependencies(dependency_list, ccdata, "gbasis"):

                def get_symmetry_atom_basis(gbasis):
                    """Get symmetry atom by replicating the last atom in gbasis of the same element."""

                    missing_index = len(gbasis)
                    missing_atomno = ccdata.atomnos[missing_index]

                    ngbasis = len(gbasis)
                    last_same = ngbasis - ccdata.atomnos[:ngbasis][::-1].index(missing_atomno) - 1
                    return gbasis[last_same]

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
                gbasis = []
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
                    while index > len(gbasis) + 1:
                        gbasis.append(get_symmetry_atom_basis(gbasis))

                    gbasis.append([])
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
                        gbasis[-1].append(shell)

                        line = file_handler.virtual_next()

                    line = file_handler.virtual_next()

                # We will also need to add symmetry atoms that are missing from the input
                # at the end of this block, if the symmetry atoms are last.
                while len(gbasis) < ccdata.natom:
                    gbasis.append(get_symmetry_atom_basis(gbasis))

                return gbasis
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> list | None:
        constructed_data = None
        if program in gbasis.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(gbasis, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
