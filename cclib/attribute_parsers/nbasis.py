# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class nbasis(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "psi4", "qchem"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if line[1:7] == "NBasis" or line[4:10] == "NBasis":
            # For counterpoise fragment, skip these lines.
            #
            # AED: not sure how to handle these cases at the moment
            # if self.counterpoise != 0:
            #   return
            # For ONIOM calcs, ignore this section in order to bypass assertion failure.
            # if self.oniom:
            #   return
            # If nbasis was already parsed, check if it changed. If it did, issue a warning.
            # In the future, we will probably want to have nbasis, as well as nmo below,
            # as a list so that we don't need to pick one value when it changes.
            constructed_data = int(line.split("=")[1].split()[0])
            if hasattr(ccdata, "nbasis"):
                try:
                    assert constructed_data == ccdata.nbasis
                except AssertionError:
                    Warning(
                        f"Number of basis functions (nbasis) is different from that stores in ccdata. Changing from {ccdata.nbasis} to {constructed_data}"
                    )
            return {nbasis.__name__: constructed_data}
        return None

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        if "Primary Basis" in line:
            file_handler.skip_lines(["Basis Set", "Blend", "Number of", "Number"], virtual=True)
            line = file_handler.virtual_next()
            constructed_data = int(line.split()[-1])
            return {nbasis.__name__: constructed_data}
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_nbasis = None
        constructed_data = None
        if "basis functions" in line:
            if not hasattr(self, "nbasis"):
                constructed_nbasis = int(line.split()[-3])
        if constructed_nbasis is not None:
            constructed_data = {nbasis.__name__: constructed_nbasis}
        return constructed_data

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in nbasis.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(nbasis, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
