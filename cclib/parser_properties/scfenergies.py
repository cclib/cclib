from cclib.parser_properties.base_parser import base_parser
from cclib.parser_properties import utils

class scfenergies(base_parser):
    """
    Docstring? Units?
    """
    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata):
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_data = None
        if line[1:9] == 'SCF Done':
            constructed_data = utils.float(line.split()[4])
            return constructed_data
        return None

    @staticmethod
    def parse(file_handler, program, ccdata):
        constructed_data = None
        if program in scfenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(scfenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
