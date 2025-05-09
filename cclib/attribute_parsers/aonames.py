# Copyright (c) 2024-2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class aonames(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        dependency_list = ["nmo", "nbasis"]
        line = file_handler.last_line
        constructed_data = None
        if (
            line[5:35] == "Molecular Orbital Coefficients"
            or line[5:41] == "Alpha Molecular Orbital Coefficients"
            or line[5:40] == "Beta Molecular Orbital Coefficients"
        ):
            constructed_aonames = []
            if not base_parser.check_dependencies(dependency_list, ccdata, "atombasis"):
                return None
            beta = False
            if line[5:40] == "Beta Molecular Orbital Coefficients":
                beta = True
            base = 0
            curr_atombasis = []
            for base in range(0, ccdata.nmo, 5):
                colNames = file_handler.virtual_next()  # noqa: F841
                symmetries = file_handler.virtual_next()  # noqa: F841
                eigenvalues = file_handler.virtual_next()  # noqa: F841
                for i in range(ccdata.nbasis):
                    line = file_handler.virtual_next()
                    if i == 0:
                        # Find location of the start of the basis function name
                        start_of_basis_fn_name = line.find(line.split()[3]) - 1
                    if base == 0 and not beta:  # Just do this the first time 'round
                        parts = line[:start_of_basis_fn_name].split()
                        if len(parts) > 1:  # New atom
                            atomname = f"{parts[2]}{parts[1]}"
                        orbital = line[start_of_basis_fn_name:20].strip()
                        constructed_aonames.append(f"{atomname}_{orbital}")
            constructed_data = {aonames.__name__: constructed_aonames}
            return constructed_data

        # Natural orbital coefficients (nocoeffs) and occupation numbers (nooccnos),
        # which are respectively define the eigenvectors and eigenvalues of the
        # diagonalized one-electron density matrix. These orbitals are formed after
        # configuration interaction (CI) calculations, but not only. Similarly to mocoeffs,
        # we can parse and check aonames and atombasis here.
        #
        #     Natural Orbital Coefficients:
        #                           1         2         3         4         5
        #     Eigenvalues --     2.01580   2.00363   2.00000   2.00000   1.00000
        #   1 1   O  1S          0.00000  -0.15731  -0.28062   0.97330   0.00000
        #   2        2S          0.00000   0.75440   0.57746   0.07245   0.00000
        # ...
        #
        def natural_orbital_single_spin_parsing(fh):
            coeffs = np.zeros((ccdata.nmo, ccdata.nbasis), "d")  # noqa: F841
            this_atombasis = []
            for base in range(0, ccdata.nmo, 5):
                colmNames = fh.virtual_next()  # noqa: F841
                eigenvalues = fh.virtual_next()  # noqa: F841
                for i in range(ccdata.nbasis):
                    line = fh.virtual_next()
                    # Just do this the first time 'round.
                    if base == 0:
                        parts = line[:12].split()
                        # New atom.
                        if len(parts) > 1:
                            if i > 0:
                                this_atombasis.append(basisonatom)  # noqa: F821
                            basisonatom = []
                        orbital = line[11:20].strip()  # noqa: F841
                        basisonatom.append(i)
                    part = line[21:].replace("D", "E").rstrip()
                    temp = []
                    for j in range(0, len(part), 10):
                        temp.append(float(part[j : j + 10]))
                # Do the last update of atombasis.
                if base == 0:
                    this_atombasis.append(basisonatom)
            return this_atombasis

        if line[5:33] == "Natural Orbital Coefficients":
            parsed_atombasis = natural_orbital_single_spin_parsing(file_handler)
            return {atombasis.__name__: parsed_atombasis}

        if line[5:39] == "Alpha Natural Orbital Coefficients":
            parsed_atombasis = natural_orbital_single_spin_parsing(file_handler)
            return {atombasis.__name__: parsed_atombasis}
        if line[5:38] == "Beta Natural Orbital Coefficients":
            parsed_atombasis = natural_orbital_single_spin_parsing(file_handler)
            return {atombasis.__name__: parsed_atombasis}

        return None

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        # Molecular orbital overlap matrix.
        # Has to deal with lines such as:
        #   *** Overlap ***
        #   ****** Overlap ******
        # Note that Gaussian sometimes drops basis functions,
        #  causing the overlap matrix as parsed below to not be
        line = file_handler.last_line
        if line[1:4] == "***" and (line[5:12] == "Overlap" or line[8:15] == "Overlap"):
            # Ensure that this is the main calc and not a fragment
            # if ccdata.counterpoise != 0:
            #   return
            constructed_aonames = np.zeros((ccdata.nbasis, ccdata.nbasis), "d")
            # Overlap integrals for basis fn#1 are in aonames[0]
            base = 0
            colmNames = file_handler.virtual_next()  # noqa: F841
            while base < ccdata.nbasis:
                for i in range(ccdata.nbasis - base):  # Fewer lines this time
                    line = file_handler.virtual_next()
                    parts = line.split()
                    for j in range(len(parts) - 1):  # Some lines are longer than others
                        k = float(parts[j + 1].replace("D", "E"))
                        constructed_aonames[base + j, i + base] = k
                        constructed_aonames[i + base, base + j] = k
                base += 5
                colmNames = file_handler.virtual_next()  # noqa: F841
            return {aonames.__name__: constructed_aonames}

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in aonames.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(aonames, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
