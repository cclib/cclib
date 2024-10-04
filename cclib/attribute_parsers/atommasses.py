# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class atommasses(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "psi4", "qchem"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_data = None
        if line.strip() == "Isotopes and Nuclear Properties:":
            constructed_data = []
            line = file_handler.virtual_next()
            while line[1:16] != "Leave Link  101":
                if line[1:8] == "AtmWgt=":
                    constructed_data.extend(list(map(float, line.split()[1:])))
                line = file_handler.virtual_next()
            return {atommasses.__name__: constructed_data}
        return None

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if "Geometry (in" in line:
            this_atommasses = []
            file_handler.skip_lines(["blank"], virtual=True)
            line = file_handler.virtual_next()
            if line.split()[0] == "Center":
                file_handler.skip_lines(["dashes"], virtual=True)
                line = file_handler.virtual_next()
            while line.strip():
                chomp = line.split()
                # Newer versions of Psi4 print atomic masses.
                if len(chomp) == 5:
                    this_atommasses.append(float(chomp[4]))
                line = file_handler.virtual_next()
            return {atommasses.__name__: this_atommasses}
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        if "STANDARD THERMODYNAMIC QUANTITIES AT" in line:
            file_handler.skip_lines(["blank"], virtual=True)
            line = file_handler.virtual_next()
            if ccdata.natom == 1:
                assert "Translational Enthalpy" in line
            else:
                assert "Imaginary Frequencies" in line
                line = file_handler.virtual_next()
                while "Translational Enthalpy" not in line:
                    if "Has Mass" in line:
                        atommass = float(line.split()[6])
                        atommasses.append(atommass)
                    line = file_handler.virtual_next()
                if not hasattr(self, "atommasses"):  # noqa: F821
                    return {atommasses.__name__: np.array(atommasses)}
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in atommasses.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(atommasses, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
