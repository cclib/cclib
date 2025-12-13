# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class parser_state(base_parser):
    """
    A temporary variable
    """

    known_codes = ["psi4", "ORCA"]

    @staticmethod
    def psi4(file_handler, ccdata):
        line = file_handler.last_line
        if "Properties will be evaluated at" in line.strip():
            if getattr(ccdata, "parser_state"):
                this_metadata = ccdata.parser_state
            else:
                this_metadata = {}
            tokens = line.split()
            assert tokens[-1] in ["Bohr", "[a0]"]
            this_metadata["origin"] = utils.convertor(
                np.array([float(x.strip(",")) for x in line.split()[-4:-1]]), "bohr", "Angstrom"
            )
            return {parser_state.__name__: this_metadata}
        return None

    @staticmethod
    def ORCA(file_handler, ccdata):
        line = file_handler.last_line
        # Needed for GBASIS
        # Basis set information
        # ORCA prints this out in a somewhat indirect fashion.
        # Therefore, parsing occurs in several steps:
        # 1. read which atom belongs to which basis set group
        this_metadata = {}
        if getattr(ccdata, "parser_state"):
            this_metadata = ccdata.parser_state

        if line[0:21] == "BASIS SET INFORMATION":
            gbasis_tmp_atnames = []  # temporary attribute, needed later
            if "gbasis_tmp_atnames" in ccdata.parser_state:
                gbasis_tmp_atnames = ccdata.parser_state["gbasis_tmp_atnames"]
            line = file_handler.virtual_next()
            line = file_handler.virtual_next()
            while not line[0:5] == "-----":
                if line[0:4] == "Atom":
                    gbasis_tmp_atnames.append(line[8:12].strip())
                line = file_handler.virtual_next()
            this_metadata["gbasis_tmp_atnames"] = gbasis_tmp_atnames
            return {parser_state.__name__: this_metadata}
        # do we use symmetry
        if line[1:18] == "Symmetry handling":
            this_metadata["uses_symmetry"] = True
            return {parser_state.__name__: this_metadata}
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in parser_state.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(parser_state, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
