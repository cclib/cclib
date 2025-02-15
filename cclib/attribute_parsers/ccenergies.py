# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class ccenergies(base_parser):
    """
    Docstring? Units?
    """

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        ccsd_trigger = "* CCSD total energy"
        ccsd_t_trigger = "* CCSD(T) total energy"
        line = file_handler.last_line
        if getattr(ccdata, "ccenergies") is None:
            this_ccenergies = []
        else:
            this_ccenergies = ccdata.ccenergies
        if line.lstrip().startswith(ccsd_trigger):
            this_ccenergies.append(float(line.split()[-1]))
            return {ccenergies.__name__: np.array(this_ccenergies)}
        if line.strip().startswith(ccsd_t_trigger):
            # Not sure how to deal with metadata  yet
            # assert ccdata.metadata["methods"][-1] == "CCSD"
            # self.metadata["methods"].append("CCSD(T)")
            this_ccenergies[-1] = float(line.split()[-1])
            return {ccenergies.__name__: np.array(this_ccenergies)}
        return None

        # The geometry convergence targets and values are printed in a table, with the legends

    known_codes = ["psi4"]

    @staticmethod
    def parse(file_handler, program, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in ccenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(ccenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
