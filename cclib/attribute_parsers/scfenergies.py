# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser
from cclib import ureg


class scfenergies(base_parser):
    """
    molecular electronic energies after SCF (Hartree-Fock, DFT)
    Units: Hartree
    """

    known_codes = ["gaussian", "psi4", "qchem"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        if line[1:9] == "SCF Done":
            constructed_data = utils.float(line.split()[4])
            print(constructed_data)
            print(ureg.hartree)
            constructed_data = [constructed_data] * ureg.hartree
            print(constructed_data)
            return {scfenergies.__name__: constructed_data}
        return None

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if "Final Energy" in line:
            constructed_data = float(line.split()[3])
            return {scfenergies.__name__: [constructed_data] * ureg.hartree}
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        constructed_scfenergies = None
        constructed_data = None
        if "Total energy in the final basis set" in line:
            constructed_scfenergies = float(line.split()[-1])
            constructed_data *= ureg.hartree
        if constructed_scfenergies is not None:
            constructed_data = {scfenergies.__name__: constructed_data}
        return constructed_data

    @staticmethod
    def parse(file_handler, program, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in scfenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(scfenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
