# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class coreelectrons(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> dict | None:
        dependency_list = ["natom"]
        line = file_handler.last_line
        if base_parser.check_dependencies(dependency_list, ccdata, "coreelectrons"):
            constructed_coreelectrons = np.zeros(ccdata.natom, "i")
            if line.find("Pseudopotential Parameters") > -1:
                file_handler.skip_lines(["e", "label1", "label2", "e"])

                line = file_handler.virtual_next()
                if line.find("Centers:") < 0:
                    return
                    # This was continue before parser refactoring.
                 # continue

                # Needs to handle code like the following:
                #
                #  Center     Atomic      Valence      Angular      Power                                                       Coordinates
                #  Number     Number     Electrons     Momentum     of R      Exponent        Coefficient                X           Y           Z
                # ===================================================================================================================================
                # Centers:   1
                # Centers:  16
                # Centers:  21 24
                # Centers:  99100101102
                #    1         44           16                                                                      -4.012684 -0.696698  0.006750
                #                                      F and up
                #                                                     0      554.3796303       -0.05152700
                centers = []
                while line.find("Centers:") >= 0:
                    temp = line[10:]
                    for i in range(0, len(temp) - 3, 3):
                       centers.append(int(temp[i : i + 3]))
                    line = file_handler.virtual_next()
                centers.sort()  # Not always in increasing order

                constructed_coreelectrons = np.zeros(ccdata.natom, "i")

                for center in centers:
                    front = line[:10].strip()
                    while not (front and int(front) == center):
                        line = file_handler.virtual_next()
                        front = line[:10].strip()
                    info = line.split()
                    constructed_coreelectrons[center - 1] = int(info[1]) - int(info[2])
                    line = file_handler.virtual_next()
            return {coreelectrons.__name__: constructed_coreelectrons}

        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> dict | None:
        constructed_data = None
        if program in coreelectrons.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(coreelectrons, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
