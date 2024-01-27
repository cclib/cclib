from cclib.parser_properties import utils
from cclib.parser_properties.base_parser import base_parser

import numpy as np


class mocoeffs(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["psi4"]

    @staticmethod
    def psi4(file_handler, ccdata) -> list | None:
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
                    iao = int(chomp[0])
                    coeffs = [float(c) for c in chomp[m - n :]]
                    for i, c in enumerate(coeffs):
                        mocoeffs[indices[i] - 1].append(c)
                    line = file_handler.virtual_next()

                line = file_handler.virtual_next()
                line = file_handler.virtual_next()
                line = file_handler.virtual_next()
                file_handler.skip_lines(["b", "b"], virtual=True)
                indices = file_handler.virtual_next()

            if getattr(ccdata, "mocoeffs") != None:
                extended_mocoeffs = [ccdata.mocoeffs, mocoeffs]
                return extended_mocoeffs
            else:
                return mocoeffs
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> list | None:
        constructed_data = None
        if program in mocoeffs.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(mocoeffs, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
