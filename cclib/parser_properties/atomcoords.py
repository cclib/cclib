# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from cclib.parser_properties import utils
from cclib.parser_properties.base_parser import base_parser

import numpy as np


class atomcoords(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata):
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_data = None
        if line.strip() == "Standard orientation:":
            line = file_handler.skip_lines(["d", "Center", "Number", "d"], virtual=True)
            line = file_handler.next()
            constructed_data = []
            atomcoords = []
            while list(set(line.strip())) != ["-"]:
                broken = line.split()
                atomcoords.append(list(map(float, broken[-3:])))
                line = file_handler.virtual_next()
            constructed_data.append(atomcoords)
            return constructed_data
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> np.ndarray:
        constructed_data = None
        if program in atomcoords.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(atomcoords, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
