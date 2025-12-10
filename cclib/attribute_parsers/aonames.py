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

    known_codes = ["gaussian", "ORCA"]

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
            this_aonames = []
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
                                this_aonames.append(basisonatom)
                            basisonatom = []
                        orbital = line[11:20].strip()  # noqa: F841
                        basisonatom.append(i)
                    part = line[21:].replace("D", "E").rstrip()
                    temp = []
                    for j in range(0, len(part), 10):
                        temp.append(float(part[j : j + 10]))
                # Do the last update of aonames.
                if base == 0:
                    this_aonames.append(basisonatom)
            return this_aonames

        if "Natural Orbital Coefficients" in line:
            parsed_aonames = natural_orbital_single_spin_parsing(file_handler)
            return {aonames.__name__: parsed_aonames}

        return None

    @staticmethod
    def ORCA(file_handler, ccdata) -> Optional[dict]:
        # Molecular orbital coefficients are parsed here, but also related things
        # like atombasis and aonames if possible.
        #
        # Normally the output is easy to parse like this:
        # ------------------
        # MOLECULAR ORBITALS
        # ------------------
        #                       0         1         2         3         4         5
        #                  -19.28527 -19.26828 -19.26356 -19.25801 -19.25765 -19.21471
        #                    2.00000   2.00000   2.00000   2.00000   2.00000   2.00000
        #                   --------  --------  --------  --------  --------  --------
        #   0C   1s         0.000002 -0.000001  0.000000  0.000000 -0.000000  0.000001
        #   0C   2s        -0.000007  0.000006 -0.000002 -0.000000  0.000001 -0.000003
        #   0C   3s        -0.000086 -0.000061  0.000058 -0.000033 -0.000027 -0.000058
        # ...
        #
        # But when the numbers get big, things get yucky since ORCA does not use
        # fixed width formatting for the floats, and does not insert extra spaces
        # when the numbers get wider. So things get stuck together overflowing columns,
        # like this:
        #   12C   6s       -11.608845-53.775398161.302640-76.633779 29.914985 22.083999
        #
        # One assumption that seems to hold is that there are always six significant
        # digits in the coefficients, so we can try to use that to delineate numbers
        # when the parsing gets rough. This is what we do below with a regex, and a case
        # like this is tested in regression ORCA/ORCA4.0/invalid-literal-for-float.out
        # which was reported in https://github.com/cclib/cclib/issues/629
        line = file_handler.last_line
        if line[0:18] == "MOLECULAR ORBITALS":
            line = file_handler.virtual_next() # dashes

            parsed_aonames = []
            for spin in range(len(self.moenergies)):
                if spin == 1:
                    line = file_handler.virtual_next() # blank
                for i in range(0, self.nbasis, 6):
                    file_handler.skip_lines(["numbers", "energies, "occs", "d"], virtual=True)
                    for j in range(self.nbasis):
                        line = file_handler.virtual_next()
                        # Only need this in the first iteration.
                        if spin == 0 and i == 0:
                            atomname = line[3:5].split()[0]
                            num = int(line[0:3])
                            orbital = line.split()[1].upper()
                            parsed_aonames.append(f"{atomname}{int(num + 1)}_{orbital}")
            return {aonames.__name__: parsed_aonames}
        return None


    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in aonames.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(aonames, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
