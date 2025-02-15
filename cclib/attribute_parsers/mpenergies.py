# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class mpenergies(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["psi4"]

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        # This is for the older conventional MP2 code in 4.0b5.
        mp_trigger = "MP2 Total Energy (a.u.)"
        if line.strip()[: len(mp_trigger)] == mp_trigger:
            mpenergy = utils.convertor(float(line.split()[-1]), "hartree", "eV")
            if getattr(ccdata, "mpenergies") is None:
                this_mpenergies = []
            this_mpenergies.append([mpenergy])
            return {mpenergies.__name__: this_mpenergies}
        # This is for the newer DF-MP2 code in 4.0.
        if "DF-MP2 Energies" in line:
            while "Total Energy" not in line:
                line = file_handler.virtual_next()
            mpenergy = utils.convertor(float(line.split()[3]), "hartree", "eV")
            if getattr(ccdata, "mpenergies") is None:
                this_mpenergies = []
            this_mpenergies.append([mpenergy])
            return {mpenergies.__name__: np.array(this_mpenergies)}
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in mpenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(mpenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
