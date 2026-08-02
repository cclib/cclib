# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser
from cclib.units import ureg

import numpy as np


class mpenergies(base_parser):
    """

    V2 note: not set up correctly for optimizations yet

    long: The attribute mpenergies holds the total molecule energies including Møller-Plesset correlation energy corrections in a two-dimensional array.
    The arrayi\’s shape is (n,L), where n is 1 for single point calculations and larger for optimisations, and L is the order at which the correction is truncated.

    short:

    Units: eV
    """

    known_codes = ["psi4", "gaussian"]
    cclib_units = ureg.eV

    def _append_row(ccdata, value_ev):
        existing = getattr(ccdata, "mpenergies", None)
        if existing is None:
            val = np.array([[value_ev]])
        else:
            val = np.vstack(existing.m_as(mpenergies.cclib_units), [value_ev])
        return ureg.Quantity(val, mpenergies.cclib_units)

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        psi4_units = ureg.hartree
        line = file_handler.last_line
        # This is for the older conventional MP2 code in 4.0b5.
        mp_trigger = "MP2 Total Energy (a.u.)"
        if line.strip()[: len(mp_trigger)] == mp_trigger:
            mpenergy = (float(line.split()[-1]) * psi4_units).to(mpenergies.cclib_units)
            this_mpenergies = mpenergies._append_row(ccdata, mpenergy.magnitude)
            return {mpenergies.__name__: this_mpenergies}
        # This is for the newer DF-MP2 code in 4.0.
        if "DF-MP2 Energies" in line:
            while "Total Energy" not in line:
                line = file_handler.virtual_next()
            mpenergy = (float(line.split()[3]) * psi4_units).to(mpenergies.cclib_units)
            this_mpenergies = mpenergies._append_row(ccdata, mpenergy.magnitude)
            return {mpenergies.__name__: this_mpenergies}
        return None

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        def get_row():
            existing = getattr(ccdata, "mpenergies", None)
            if existing is None:
                return [np.array([])]
            return [row.copy() for row in existing.m_as(mpenergies.cclib_units)]

        gaussian_units = ureg.hartree
        line = file_handler.last_line
        # Total energies after Moller-Plesset corrections.
        # Second order correction is always first, so its first occurance
        #   triggers creation of mpenergies (list of lists of energies).
        # Further MP2 corrections are appended as found.
        #
        # Example MP2 output line:
        #  E2 =    -0.9505918144D+00 EUMP2 =    -0.28670924198852D+03
        # Warning! this output line is subtly different for MP3/4/5 runs
        # Newer versions of gausian introduced a space between 'EUMP2' and '='...
        if "EUMP2 =" in line[27:36] or "EUMP2=" in line[27:35]:
            rows = get_row()
            x = (utils.float(line.split("=")[2]) * gaussian_units).to(mpenergies.cclib_units)
            rows[-1] = np.append(rows[-1], x.magnitude)
            this_mpenergies = ureg.Quantity(rows, mpenergies.cclib_units)
            return {mpenergies.__name__: this_mpenergies}

        # Example MP3 output line:
        #  E3=       -0.10518801D-01     EUMP3=      -0.75012800924D+02
        if line[34:40] == "EUMP3=":
            rows = get_row()
            x = (utils.float(line.split("=")[2]) * gaussian_units).to(mpenergies.cclib_units)
            rows[-1] = np.append(rows[-1], x.magnitude)
            this_mpenergies = ureg.Quantity(rows, mpenergies.cclib_units)
            return {mpenergies.__name__: this_mpenergies}

        # Example MP4 output lines:
        #  E4(DQ)=   -0.31002157D-02        UMP4(DQ)=   -0.75015901139D+02
        #  E4(SDQ)=  -0.32127241D-02        UMP4(SDQ)=  -0.75016013648D+02
        #  E4(SDTQ)= -0.32671209D-02        UMP4(SDTQ)= -0.75016068045D+02
        # Energy for most substitutions is used only (SDTQ by default)
        if line[34:43] == "UMP4(DQ)=":
            mp4energy = utils.float(line.split("=")[2]) * gaussian_units
            line = file_handler.virtual_next()
            if line[34:44] == "UMP4(SDQ)=":
                mp4energy = utils.float(line.split("=")[2]) * gaussian_units
                line = file_handler.virtual_next()
                if line[34:45] == "UMP4(SDTQ)=":
                    mp4energy = utils.float(line.split("=")[2]) * gaussian_units
            mp4energy = mp4energy.to(mpenergies.cclib_units)
            rows = get_row()
            rows[-1] = np.append(rows[-1], mp4energy.magnitude)
            this_mpenergies = ureg.Quantity(rows, mpenergies.cclib_units)
            return {mpenergies.__name__: this_mpenergies}

        # Example MP5 output line:
        #  DEMP5 =  -0.11048812312D-02 MP5 =  -0.75017172926D+02
        if line[29:34] == "MP5 =":
            mp5energy = (utils.float(line.split("=")[2]) * gaussian_units).to(
                mpenergies.cclib_units
            )
            rows = get_row()
            rows[-1] = np.append(rows[-1], mp5energy.magnitude)
            this_mpenergies = ureg.Quantity(rows, mpenergies.cclib_units)
            return {mpenergies.__name__: this_mpenergies}
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> dict | None:
        constructed_data = None
        if program in mpenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(mpenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
