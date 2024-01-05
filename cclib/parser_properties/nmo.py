from cclib.parser_properties import utils
from cclib.parser_properties.base_parser import base_parser

import numpy as np


class nmo(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> int | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_data = None
        if line[1:7] == "NBsUse":
            # For counterpoise fragment, skip these lines.
            #
            # AED: not sure how to handle these cases at the moment
            # if ccdata.counterpoise != 1:
            #   return
            # For ONIOM calcs, ignore this section in order to bypass assertion failure.
            # if ccdata.oniom:
            #   return
            constructed_data = int(line.split("=")[1].split()[0])
            return constructed_data
        return None
  
    @staticmethod
    def parse(file_handler, program: str, ccdata) -> int | None:
        constructed_data = None
        if program in nmo.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(nmo, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
