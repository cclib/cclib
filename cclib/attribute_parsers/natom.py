# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class natom(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "psi4", "qchem"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        # The remaining part will allow us to get the atom count.
        # When coordinates are given, there is a blank line at the end, but if
        # there is a Z-matrix here, there will also be variables and we need to
        # stop at those to get the right atom count.
        # Also, in older versions there is bo blank line (G98 regressions),
        # so we need to watch out for leaving the link.
        dependency_list = ["atomnos"]
        line = file_handler.last_line
        if base_parser.check_dependencies(dependency_list, ccdata, "natom"):
            return {natom.__name__: len(ccdata.atomnos)}

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        dependency_list = ["atomnos"]
        if base_parser.check_dependencies(dependency_list, ccdata, "natom"):
            return {natom.__name__: len(ccdata.atomnos)}
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> dict | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        dependency_list = ["atomnos"]
        if base_parser.check_dependencies(dependency_list, ccdata, "natom"):
            return {natom.__name__: len(ccdata.atomnos)}
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in natom.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(natom, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
