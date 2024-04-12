# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser
from typing import Optional
import numpy as np


class atomnos(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "psi4", "qchem"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        if line.strip() == "Standard orientation:":
            file_handler.skip_lines(["d", "Center", "Number", "d"], virtual=True)
            line = file_handler.virtual_next()
            constructed_data = []
            # print('line in file is')
            # print(line)
            while list(set(line.strip())) != ["-"]:
                broken = line.split()
                constructed_data.append(int(broken[1]))
                line = file_handler.virtual_next()
            return {atomnos.__name__: constructed_data}
        return None

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        table = utils.PeriodicTable()
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        if line.strip() == "-Contraction Scheme:":
            file_handler.skip_lines(["headers", "d"], virtual=True)
            line = file_handler.virtual_next()
            constructed_data = []
            while line.strip():
                element = line.split()[1]
                if len(element) > 1:
                    element = element[0] + element[1:].lower()
                constructed_data.append(table.number[element])
                line = file_handler.virtual_next()
            return {atomnos.__name__: constructed_data}
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> dict | None:
        table = utils.PeriodicTable()
        constructed_data = None
        # Extract the atomic numbers
        line = file_handler.last_line
        if "Standard Nuclear Orientation" in line:
            file_handler.skip_lines(["cols", "dashes"], virtual=True)
            atomelements = []
            line = next(inputfile)
            while list(set(line.strip())) != ["-"]:
                entry = line.split()
                atomelements.append(entry[1])
                line = file_handler.virtual_next()
            # We calculate and handle atomnos no matter what, since in
            # the case of fragment calculations the atoms may change,
            # along with the charge and spin multiplicity.
            constructed_atom_nos = []
            for atomelement in atomelements:
                self.atomelements.append(atomelement)
                if atomelement == "GH":
                    constructed_atomnos.append(0)
                else:
                    constructed_atomnos.append(table.number[atomelement])
        if constructued_atomnos:
            constructed_data = {atomnos.__name__: constructed_atomnos}
        return constructed_data

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in atomnos.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(atomnos, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
