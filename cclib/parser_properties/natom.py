# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from cclib.parser_properties import utils
from cclib.parser_properties.base_parser import base_parser

import numpy as np


class natom(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["psi4"]

    @staticmethod
    def psi4(file_handler, ccdata) -> int | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        dependency_list = ["atomnos"]
        if base_parser.check_dependencies(dependency_list, ccdata, "natom"):
            return len(ccdata.atomnos)

        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> int | None:
        constructed_data = None
        if program in natom.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(natom, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
