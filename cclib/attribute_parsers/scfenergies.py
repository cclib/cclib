# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser
from cclib import ureg 
import numpy as np


class scfenergies(base_parser):
    """
    molecular electronic energies after SCF (Hartree-Fock, DFT)
    Units: Hartree
    """

    known_codes = ["gaussian", "psi4", "qchem"]
    cclib_unit = ureg.eV

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        gaussian_unit = ureg.hartree
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        existing = getattr(ccdata, "scfenergies")
        if existing is None:
            this_scfenergies = ureg.Quantity(np.array([]), scfenergies.cclib_unit)
        else:
            this_scfenergies=existing
        if line[1:9] == "SCF Done":
            constructed_data = utils.float(line.split()[4]) * gaussian_unit
            this_scfenergies = np.append(this_scfenergies,constructed_data)
            return {scfenergies.__name__: this_scfenergies}
        return None

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        psi4_unit = ureg.hartree
        line = file_handler.last_line
        existing = getattr(ccdata, "scfenergies")
        if existing is None:
            this_scfenergies = ureg.Quantity(np.array([]), scfenergies.cclib_unit)
        else:
            this_scfenergies=existing
        if "Final Energy" in line:
            constructed_data = float(line.split()[3]) * psi4_unit
            this_scfenergies = np.append(this_scfenergies,constructed_data)
            return {scfenergies.__name__: this_scfenergies}
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> Optional[dict]:
        qchem_unit = ureg.hartree
        line = file_handler.last_line
        existing = getattr(ccdata, "scfenergies")
        if existing is None:
            this_scfenergies = ureg.Quantity(np.array([]), scfenergies.cclib_unit)
        else:
            this_scfenergies=existing
        constructed_scfenergies = None
        if "Total energy in the final basis set" in line:
            constructed_scfenergies = float(line.split()[-1]) * qchem_unit
            if constructed_scfenergies is not None:
                this_scfenergies = np.append(this_scfenergies,constructed_scfenergies)
            return {scfenergies.__name__: this_scfenergies}


    @staticmethod
    def parse(file_handler, program, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in scfenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(scfenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
