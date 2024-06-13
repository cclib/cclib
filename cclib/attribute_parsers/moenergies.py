# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class moenergies(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if line[1:6] == "Alpha" and line.find("eigenvalues") >= 0:
            # For counterpoise fragments, skip these lines.
            if getattr(ccdata, "moenergies") is None:
                constructed_moenergies = [[]]
                while line.find("Alpha") == 1:
                    part = line[28:]
                    i = 0
                    while i * 10 + 4 < len(part):
                        s = part[i * 10 : (i + 1) * 10]
                        try:
                            x = utils.float(s)
                        except ValueError:
                            x = np.nan
                        constructed_moenergies[0].append(x)
                        i += 1
                    line = file_handler.virtual_next()
                if line.find("Beta") == 2:
                    constructed_moenergies.append([])

                while line.find("Beta") == 2:
                    part = line[28:]
                    i = 0
                    while i * 10 + 4 < len(part):
                        x = part[i * 10 : (i + 1) * 10]
                        constructed_moenergies[1].append(utils.float(x))
                        i += 1
                    line = file_handler.virtual_next()

                constructed_moenergies = [np.array(x, "d") for x in constructed_moenergies]
                return {moenergies.__name__: constructed_moenergies}
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in moenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(moenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
