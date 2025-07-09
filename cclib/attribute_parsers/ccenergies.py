# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class ccenergies(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["psi4", "gaussian"]

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

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if getattr(ccdata, "ccenergies") is None:
            this_ccenergies = []
        else:
            this_ccenergies = []  # note we only save the last ccenergy
        # Total energies after Coupled Cluster corrections.
        # Second order MBPT energies (MP2) are also calculated for these runs,
        # but the output is the same as when parsing for mpenergies.
        # Read the consecutive correlated energies
        # but append only the last one to ccenergies.
        # Only the highest level energy is appended - ex. CCSD(T), not CCSD.
        parsed_ccenergy = None
        if line[1:10] == "DE(Corr)=" and line[27:35] == "E(CORR)=":
            parsed_ccenergy = utils.float(line.split()[3])
        if line[1:10] == "T5(CCSD)=":
            line = file_handler.virtual_next()
            if line[1:9] == "CCSD(T)=":
                parsed_ccenergy = utils.float(line.split()[1])

        if parsed_ccenergy:
            this_ccenergies.append(parsed_ccenergy)
            return {ccenergies.__name__: np.array(this_ccenergies)}
        return None

    @staticmethod
    def parse(file_handler, program, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in ccenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(ccenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
