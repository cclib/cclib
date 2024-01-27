from cclib.parser_properties import utils
from cclib.parser_properties.base_parser import base_parser

import numpy as np


class atomcoords(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "psi4"]

    @staticmethod
    def gaussian(file_handler, ccdata):
        line = file_handler.last_line
        constructed_data = None
        if line.strip() == "Standard orientation:":
            line = file_handler.skip_lines(["d", "Center", "Number", "d"], virtual=True)
            line = file_handler.virtual_next()
            constructed_data = []
            atomcoords = []
            while list(set(line.strip())) != ["-"]:
                broken = line.split()
                atomcoords.append(np.asarray(map(float, broken[-3:])))
                line = file_handler.virtual_next()
            constructed_data.append(atomcoords)
            return constructed_data
        return None

    @staticmethod
    def psi4(file_handler, ccdata):
        line = file_handler.last_line
        if "Geometry (in " in line:
            ## I am not handling units here, i think this should be done on the ccdata object and not the parser as mentioned in #1124
            ## but will need a way to indicate which units are stored i guess. hm.
            tokens = line.split()
            units = tokens[2][:-2]
            assert units in ("Angstrom", "Bohr")
            file_handler.skip_lines(["blank"], virtual=True)
            line = file_handler.virtual_next()
            if line.split()[0] == "Center":
                file_handler.skip_lines(["dashes"], virtual=True)
                line = file_handler.virtual_next()
            constructed_data = []
            while line.strip():
                chomp = line.split()
                _el, x, y, z = chomp[:4]
                constructed_data.append(np.asarray([float(x), float(y), float(z)]))
                line = file_handler.virtual_next()
            return constructed_data
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> np.ndarray:
        constructed_data = None
        if program in atomcoords.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(atomcoords, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
