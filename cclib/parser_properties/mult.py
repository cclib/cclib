from cclib.parser_properties.base_parser import base_parser
from cclib.parser_properties import utils
import numpy as np

class mult(base_parser):
    """
    Docstring? Units?
    """
    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> list | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_data = None
        if line.find("Multiplicity") > 0:
            constructed_data = line.split()[5]
            return constructed_data
        return None

    @staticmethod
    def parse(file_handler, program:str, ccdata) -> list | None:
        constructed_data = None
        if program in mult.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(mult, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data


