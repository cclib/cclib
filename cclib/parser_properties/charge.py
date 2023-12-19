from cclib.parser_properties import utils
from cclib.parser_properties.base_parser import base_parser

import numpy as np


class charge(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> list | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_data = None
        if line[1:7] == "Charge":
            constructed_data = line.split()[2]
            return constructed_data
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> list | None:
        constructed_data = None
        if program in charge.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(charge, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
