from cclib.parser_properties.base_parser import base_parser
from cclib.parser_properties import utils
import numpy as np

class atomnos(base_parser):
    """
    Docstring? Units?
    """
    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> list | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_data = None
        print(f'once in function the line is {line}')
        if line.strip() == "Standard orientation:":
            print('found the line for atom nos')
            line = file_handler.skip_lines(['d','Center','Number','d'],virtual=True)
            line = file_handler.skip_lines(['d'],virtual=True)
            constructed_data=[]
            while list(set(line.strip())) != ["-"]:
                broken = line.split()
                constructed_data.append(int(broken[1]))
                line = file_handler.virtual_next()
            return constructed_data
        return None

    @staticmethod
    def parse(file_handler, program:str, ccdata) -> list | None:
        constructed_data = None
        if program in atomnos.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(atomnos, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data


