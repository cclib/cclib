# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class mult(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "psi4", "qchem"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        if line.find("Multiplicity") > 0:
            constructed_data = int(line.split()[5])
            return {mult.__name__: constructed_data}
        return None

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        if line[2:16].lower() == "multiplicity =":
            constructed_data = int(line.split()[-1])
            return {mult.__name__: constructed_data}
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> dict | None:
        # TODO: ecp charge may be different!
        line = file_handler.last_line
        constructed_data = None
        # Number of electrons.
        # Useful for determining the number of occupied/virtual orbitals.
        if "Nuclear Repulsion Energy" in line:
            line = file_handler.virtual_next()
            nelec_re_string = r"There are(\s+[0-9]+) alpha and(\s+[0-9]+) beta electrons"
            match = re.findall(nelec_re_string, line.strip())
            nalpha_elec = int(match[0][0].strip())
            nbeta_elec = int(match[0][1].strip())
            # Calculate the spin multiplicity (2S + 1), where S is the
            # total spin of the system.
            S = (nalpha_elec - nbeta_elec) / 2
            constructed_mult = int(2 * S + 1)
            self.set_attribute("mult", mult)
            constructed_data = {mult.__name__: constructed_mult}
        return constructed_data

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in mult.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(mult, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
