# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser
from cclib.units import ureg

import numpy as np


class ccenergies(base_parser):
    """
    long: A one-dimensional array holds the total molecule energies including Coupled Cluster corrections. The array’s length is 1 for single point calculations and larger for optimisations. Only the highest theory level is parsed into this attribute (for example, CCSD energies as opposed to CCD energies, or CCSD(T) as opposed to CCSD energies).
    short: molecular energies with Coupled-Cluster corrections
    """

    known_codes = ["psi4", "gaussian"]
    cclib_unit = ureg.eV

    @staticmethod
    def psi4(file_handler, ccdata) -> dict | None:
        psi4_units = ureg.hartree
        ccsd_trigger = "* CCSD total energy"
        ccsd_t_trigger = "* CCSD(T) total energy"
        line = file_handler.last_line
        existing = getattr(ccdata, "ccenergies")
        if existing is None:
            this_ccenergies = ureg.Quantity(np.array([]), ccenergies.cclib_unit)
        else:
            this_ccenergies = existing
        if line.lstrip().startswith(ccsd_trigger):
            value = (float(line.split()[-1]) * psi4_units).to(ccenergies.cclib_unit)
            this_ccenergies = np.append(this_ccenergies, value)
            return {ccenergies.__name__: this_ccenergies}
        if line.strip().startswith(ccsd_t_trigger):
            # Not sure how to deal with metadata  yet
            # assert ccdata.metadata["methods"][-1] == "CCSD"
            # self.metadata["methods"].append("CCSD(T)")
            this_ccenergies[-1] = (float(line.split()[-1]) * psi4_units).to(ccenergies.cclib_unit)
            return {ccenergies.__name__: this_ccenergies}
        return None

    @staticmethod
    def gaussian(file_handler, ccdata) -> dict | None:
        gaussian_units = ureg.hartree
        line = file_handler.last_line
        this_ccenergies = ureg.Quantity(np.array([]), ccenergies.cclib_unit)
        # Total energies after Coupled Cluster corrections.
        # Second order MBPT energies (MP2) are also calculated for these runs,
        # but the output is the same as when parsing for mpenergies.
        # Read the consecutive correlated energies
        # but append only the last one to ccenergies.
        # Only the highest level energy is appended - ex. CCSD(T), not CCSD.
        parsed_ccenergy = None
        if line[1:10] == "DE(Corr)=" and line[27:35] == "E(CORR)=":
            parsed_ccenergy = (utils.float(line.split()[3]) * gaussian_units).to(
                ccenergies.cclib_unit
            )
        if line[1:10] == "T5(CCSD)=":
            line = file_handler.virtual_next()
            if line[1:9] == "CCSD(T)=":
                parsed_ccenergy = (utils.float(line.split()[1]) * gaussian_units).to(
                    ccenergies.cclib_unit
                )

        if parsed_ccenergy:
            this_ccenergies = np.append(this_ccenergies, parsed_ccenergy)
            return {ccenergies.__name__: this_ccenergies}
        return None

    @staticmethod
    def parse(file_handler, program, ccdata) -> dict | None:
        constructed_data = None
        if program in ccenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(ccenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
