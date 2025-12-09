# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class aooverlaps(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "ORCA"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        # Molecular orbital overlap matrix.
        # Has to deal with lines such as:
        #   *** Overlap ***
        #   ****** Overlap ******
        # Note that Gaussian sometimes drops basis functions,
        #  causing the overlap matrix as parsed below to not be
        dependency_list = ["nbasis"]
        if base_parser.check_dependencies(dependency_list, ccdata, "aooverlaps"):
            line = file_handler.last_line
            if line[1:4] == "***" and (line[5:12] == "Overlap" or line[8:15] == "Overlap"):
                # Ensure that this is the main calc and not a fragment
                # if ccdata.counterpoise != 0:
                #   return
                constructed_aooverlaps = np.zeros((ccdata.nbasis, ccdata.nbasis), "d")
                # Overlap integrals for basis fn#1 are in aooverlaps[0]
                base = 0
                colmNames = file_handler.virtual_next()  # noqa: F841
                while base < ccdata.nbasis:
                    for i in range(ccdata.nbasis - base):  # Fewer lines this time
                        line = file_handler.virtual_next()
                        parts = line.split()
                        for j in range(len(parts) - 1):  # Some lines are longer than others
                            k = float(parts[j + 1].replace("D", "E"))
                            constructed_aooverlaps[base + j, i + base] = k
                            constructed_aooverlaps[i + base, base + j] = k
                    base += 5
                    colmNames = file_handler.virtual_next()  # noqa: F841
                return {aooverlaps.__name__: constructed_aooverlaps}
        return None

    @staticmethod
    def ORCA(file_handler, ccdata) -> Optional[dict]:
        dependency_list = ["nbasis"]
        if base_parser.check_dependencies(dependency_list, ccdata, "aooverlaps"):
            line = file_handler.last_line
            if line[0:14] == "OVERLAP MATRIX":
                dashes = file_handler.virtual_next()  # noqa: F841
                constructed_aooverlaps = np.zeros((ccdata.nbasis, ccdata.nbasis), "d")
                for i in range(0, ccdata.nbasis, 6):
                    #self.updateprogress(inputfile, "Overlap")
                    header = file_handler.virtual_next()
                    size = len(header.split())
                    for j in range(ccdata.nbasis):
                        line = file_handler.virtual_next()
                        broken = line.split()
                        constructed_aooverlaps[j, i : i + size] = list(map(float, broken[1 : size + 1]))
                return {aooverlaps.__name__: constructed_aooverlaps}
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in aooverlaps.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(aooverlaps, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
