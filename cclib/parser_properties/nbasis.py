from cclib.parser_properties import utils
from cclib.parser_properties.base_parser import base_parser

import numpy as np


class nbasis(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "psi4"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> int | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
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
                    assert nbasis == ccdata.nbasis
                except AssertionError:
                    Warning(
                        f"Number of basis functions (nbasis) is different from that stores in ccdata. Changing from {ccdata.nbasis} to {constructed_data}"
                    )
            return constructed_data
        return None

    def psi4(file_handler, ccdata) -> int | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        if "Primary Basis" in line:
            file_handler.skip_lines(["Basis Set", "Blend", "Number of", "Number"], virtual=True)
            line = file_handler.virtual_next()
            nbasis = int(line.split()[-1])
            constructed_data = nbasis
            return constructed_data
        return None

    @staticmethod
    def psi4(file_handler, ccdata) -> int | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        if "Primary Basis" in line:
            file_handler.skip_lines(["Basis Set", "Blend", "Number of", "Number"], virtual=True)
            line = file_handler.virtual_next()
            nbasis = int(line.split()[-1])
            constructed_data = nbasis
            return constructed_data
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> int | None:
        constructed_data = None
        if program in nbasis.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(nbasis, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
