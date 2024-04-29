# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
import re

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class atombasis(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "psi4", "qchem"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> dict | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        dependency_list = ["nmo", "nbasis"]
        line = file_handler.last_line
        constructed_data = None
        if (
            line[5:35] == "Molecular Orbital Coefficients"
            or line[5:41] == "Alpha Molecular Orbital Coefficients"
            or line[5:40] == "Beta Molecular Orbital Coefficients"
        ):
            constructed_atombasis = []
            if not base_parser.check_dependencies(dependency_list, ccdata, "atombasis"):
                return None
            beta = False
            if line[5:40] == "Beta Molecular Orbital Coefficients":
                beta = True
            colNames = file_handler.virtual_next()
            symmetries = file_handler.virtual_next()
            eigenvalues = file_handler.virtual_next()
            base = 0
            curr_atombasis = []
            for base in range(0, ccdata.nmo, 5):
                for i in range(ccdata.nbasis):
                    line = file_handler.virtual_next()
                    if i == 0:
                        # Find location of the start of the basis function name
                        start_of_basis_fn_name = line.find(line.split()[3]) - 1
                    if base == 0 and not beta:  # Just do this the first time 'round
                        parts = line[:start_of_basis_fn_name].split()
                        if len(parts) > 1:  # New atom
                            if i > 0:
                                constructed_atombasis.append(curr_atombasis)
                            curr_atombasis = []
                        curr_atombasis.append(i)
                    curr_atombasis.append(i)
            constructed_data = {atombasis.__name__: constructed_atombasis}
            return constructed_data
        return None

    @staticmethod
    def psi4(file_handler, ccdata) -> dict | None:
        dependency_list = ["nmo", "nbasis"]
        if getattr(ccdata, "atombasis") == None:
            line = file_handler.last_line
            if line.strip() == "-Contraction Scheme:":
                file_handler.skip_lines(["headers", "d"], virtual=True)
                line = file_handler.virtual_next()
                constructed_atombasis = []
                atombasis_pos = 0
                while line.strip():
                    ao_count = 0
                    shells = line.split("//")[1].split()
                    for s in shells:
                        count, basistype = s
                        multiplier = 3 * (basistype == "p") or 1
                        ao_count += multiplier * int(count)
                    if len(constructed_data) > 0:
                        atombasis_pos = constructed_data[-1][-1] + 1
                    constructed_atombasis.append(
                        list(range(atombasis_pos, atombasis_pos + ao_count))
                    )
                    line = file_handler.virtual_next()
                constructed_data = {atombasis.__name__: constructed_atombasis}
                return constructed_data
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> dict | None:
        # This block comes from `print_orbitals = true/{int}`. Less
        # precision than `scf_final_print >= 2` for `mocoeffs`, but
        # important for `aonames` and `atombasis`.
        alpha_mo_coefficient_headers = []
        qchem_re_atomindex = re.compile(r"(\d+)_")

        dependency_list = ["aonames", "atomnos", "nbasis", "natom"]
        line = file_handler.last_line
        if any(header in line for header in alpha_mo_coefficient_headers):
            if base_parser.check_dependencies(dependency_list, ccdata, "atombasis"):
                constructed_atombasis = []
                for n in range(ccdata.natom):
                    constructed_atombasis.append([])
                # Go back through `aonames` to create `atombasis`.
                assert len(ccdata.aonames) == ccdata.nbasis
                for aoindex, aoname in enumerate(ccdata.aonames):
                    atomindex = int(qchem_re_atomindex.search(aoname).groups()[0]) - 1
                    constructed_atombasis.atombasis[atomindex].append(aoindex)
                assert len(self.atombasis) == len(self.atomnos)
                if not hasattr(ccdata, "atombasis"):
                    constructed_data = {atombasis.__name__: constructed_atombasis}
                    return constructed_data
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> dict | None:
        constructed_data = None
        if program in atombasis.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(atombasis, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
