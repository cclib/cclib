# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

from cclib.attribute_parsers.base_parser import base_parser
from cclib.units import ureg

import numpy as np


class dispersionenergies(base_parser):
    """
    This is a rank 1 array that contains the isolated dispersion energy for each geometry.
    This will be populated for empirical models, such as those from Grimme that only depend on relative atomic positions.

    Units: eV
    """

    known_codes = ["psi4", "gaussian"]
    cclib_unit = ureg.eV

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        psi4_unit = ureg.hartree
        line = file_handler.last_line
        existing = getattr(ccdata, "dispersionenergies", None)
        if existing is None:
            this_dispersionenergies = ureg.Quantity(np.array([]), dispersionenergies.cclib_unit)
        else:
            this_dispersionenergies = existing
        if "Empirical Dispersion Energy" in line:
            dispersion = float(line.split()[-1]) * psi4_unit
            this_dispersionenergies = np.append(this_dispersionenergies, dispersion)
            return {dispersionenergies.__name__: this_dispersionenergies}
        return None
        # The geometry convergence targets and values are printed in a table, with the legends

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        gaussian_unit = ureg.hartree
        line = file_handler.last_line
        existing = getattr(ccdata, "dispersionenergies", None)
        if existing is None:
            this_dispersionenergies = ureg.Quantity(np.array([]), dispersionenergies.cclib_unit)
        else:
            this_dispersionenergies = existing

        if "Dispersion energy=" in line:
            dispersion = float(line.split()[-2]) * gaussian_unit
            this_dispersionenergies = np.append(this_dispersionenergies, dispersion)
            return {dispersionenergies.__name__: this_dispersionenergies}
        return None

    @staticmethod
    def parse(file_handler, program, ccdata) -> dict | None:
        constructed_data = None
        if program in dispersionenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(dispersionenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
