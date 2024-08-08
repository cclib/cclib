# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers.base_parser import base_parser
from cclib.parser import utils

import numpy as np


class dispersionenergies(base_parser):
    """
    Docstring? Units?
    """

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        ccsd_trigger = "* CCSD total energy"
        ccsd_t_trigger = "* CCSD(T) total energy"
        line = file_handler.last_line
        if getattr(ccdata, "dispersionenergies") is None:
            this_dispersionenergies = []
        else:
            this_dispersionenergies = ccdata.dispoersionenergies
        if "Empirical Dispersion Energy" in line:
            dispersion = utils.convertor(float(line.split()[-1]), "hartree", "eV")
            this_dispersionenergies.append(dispersion)
            return {dispersionenergies.__name__: np.array(this_dispersionenergies)}
        return None

        # The geometry convergence targets and values are printed in a table, with the legends

    known_codes = ["psi4"]

    @staticmethod
    def parse(file_handler, program, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in dispersionenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(dispersionenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
