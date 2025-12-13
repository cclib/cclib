# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers.base_parser import base_parser


class mocoeffs(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["ORCA", "psi4"]

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if "Molecular Orbitals" in line:
            file_handler.skip_lines(["b"], virtual=True)
            indices = file_handler.virtual_next()
            mocoeffs = []
            while indices.strip():
                if indices[:3] == "***":
                    break

                indices = [int(i) for i in indices.split()]

                if len(mocoeffs) < indices[-1]:
                    for i in range(len(indices)):
                        mocoeffs.append([])
                else:
                    assert len(mocoeffs) == indices[-1]

                file_handler.skip_lines(["b"], virtual=True)

                n = len(indices)
                line = file_handler.virtual_next()
                while line.strip():
                    chomp = line.split()
                    m = len(chomp)
                    iao = int(chomp[0])  # noqa: F841
                    coeffs = [float(c) for c in chomp[m - n :]]
                    for i, c in enumerate(coeffs):
                        mocoeffs[indices[i] - 1].append(c)
                    line = file_handler.virtual_next()

                line = file_handler.virtual_next()
                line = file_handler.virtual_next()
                line = file_handler.virtual_next()
                file_handler.skip_lines(["b", "b"], virtual=True)
                indices = file_handler.virtual_next()

            if hasattr(ccdata, "mocoeffs") and getattr(ccdata, "mocoeffs") is not None:
                extended_mocoeffs = [ccdata.mocoeffs, mocoeffs]
                return {mocoeffs.__name__: extended_mocoeffs}
            else:
                return {mocoeffs.__name__: mocoeffs}
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
        dependency_list = ["natom", "nbasis", "moenergies"]
        if not base_parser.check_dependencies(dependency_list, ccdata, "atombasis"):
            return None
        if line[0:18] == "MOLECULAR ORBITALS":
            file_handler.skip_lines(["dashes"], virtual=True)
            constructed_mocoeffs = [numpy.zeros((ccdata.nbasis, ccdata.nbasis), "d")]
            for spin in range(len(ccdata.moenergies)):
                if spin == 1:
                    file_handler.skip_lines(["blank"], virtual=True)
                    constructed_mocoeffs.append(numpy.zeros((ccdata.nbasis, ccdata.nbasis), "d"))
                for i in range(0, ccdata.nbasis, 6):
                    #self.updateprogress(inputfile, "Coefficients")
                    file_handler.skip_lines(["numbers", "energies", "occs", "d"], virtual=True)
                    for j in range(ccdata.nbasis):
                        line = file_handler.virtual_next()
                        # This regex will tease out all number with exactly
                        # six digits after the decimal point.
                        coeffs = re.findall(r"-?\d+\.\d{6}", line)
                        # Something is very wrong if this does not hold.
                        assert len(coeffs) <= 6
                        constructed_mocoeffs[spin][i : i + len(coeffs), j] = [float(c) for c in coeffs]
            constructed_data = {mocoeffs.__name__: constructed_mocoeffs}
            return constructed_data
        return None


    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in mocoeffs.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(mocoeffs, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
