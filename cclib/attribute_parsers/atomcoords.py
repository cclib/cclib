# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class atomcoords(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "psi4", "qchem"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if line.strip() == "Standard orientation:":
            line = file_handler.skip_lines(["d", "Center", "Number", "d"], virtual=True)
            line = file_handler.virtual_next()
            constructed_atomcoords = []
            curr_atomcoords = []
            while list(set(line.strip())) != ["-"]:
                broken = line.split()
                curr_atomcoords.append(list(map(float, broken[-3:])))
                line = file_handler.virtual_next()
            constructed_atomcoords.append(curr_atomcoords)
            return {atomcoords.__name__: np.array(constructed_atomcoords)}
        return None

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if "Geometry (in " in line:
            if getattr(ccdata, "atomcoords") is None:
                constructed_data = []
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
            constructed_atomcoords = []
            while line.strip():
                chomp = line.split()
                _el, x, y, z = chomp[:4]
                constructed_atomcoords.append(np.asarray([float(x), float(y), float(z)]))
                line = file_handler.virtual_next()
            constructed_data.append(constructed_atomcoords)
            return {atomcoords.__name__: np.array(constructed_data)}
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> Optional[dict]:
        # Extract the atomic numbers and coordinates of the atoms.
        # TODO: afterparsing for geometries
        line = file_handler.last_line
        if "Standard Nuclear Orientation" in line:
            if "Angstroms" in line:

                def convertor(x: float) -> float:
                    return x
            elif "Bohr" in line:

                def convertor(x: float) -> float:
                    return utils.convertor(x, "bohr", "Angstrom")
            else:
                raise ValueError(f"Unknown units in coordinate header: {line}")
            file_handler.skip_lines(["cols", "dashes"], virtual=True)
            atomelements = []
            atomcoords = []
            line = next(inputfile)  # noqa: F821
            while list(set(line.strip())) != ["-"]:
                entry = line.split()
                atomelements.append(entry[1])
                atomcoords.append([convertor(float(value)) for value in entry[2:]])
                line = file_handler.virtual_next()
            constructed_data = {atomcoords.__name__: constructed_atomcoords}  # noqa: F821
            return constructed_data

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in atomcoords.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(atomcoords, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
