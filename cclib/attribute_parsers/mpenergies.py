# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class mpenergies(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["psi4", "gaussian"]

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        # This is for the older conventional MP2 code in 4.0b5.
        mp_trigger = "MP2 Total Energy (a.u.)"
        if line.strip()[: len(mp_trigger)] == mp_trigger:
            mpenergy = utils.convertor(float(line.split()[-1]), "hartree", "eV")
            if getattr(ccdata, "mpenergies") is None:
                this_mpenergies = []
            this_mpenergies.append([mpenergy])
            return {mpenergies.__name__: this_mpenergies}
        # This is for the newer DF-MP2 code in 4.0.
        if "DF-MP2 Energies" in line:
            while "Total Energy" not in line:
                line = file_handler.virtual_next()
            mpenergy = utils.convertor(float(line.split()[3]), "hartree", "eV")
            if getattr(ccdata, "mpenergies") is None:
                this_mpenergies = []
            this_mpenergies.append([mpenergy])
            return {mpenergies.__name__: np.array(this_mpenergies)}
        return None

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if getattr(ccdata, "mpenergies") is None:
            this_mpenergies = []
        else:
            this_mpenergies = ccdata.mpenergies.tolist()
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
            this_mpenergies.append([utils.float(line.split("=")[2])])
            return {mpenergies.__name__: np.array(this_mpenergies)}

        # Example MP3 output line:
        #  E3=       -0.10518801D-01     EUMP3=      -0.75012800924D+02
        if line[34:40] == "EUMP3=":
            this_mpenergies[-1].append(utils.float(line.split("=")[2]))
            return {mpenergies.__name__: np.array(this_mpenergies)}

        # Example MP4 output lines:
        #  E4(DQ)=   -0.31002157D-02        UMP4(DQ)=   -0.75015901139D+02
        #  E4(SDQ)=  -0.32127241D-02        UMP4(SDQ)=  -0.75016013648D+02
        #  E4(SDTQ)= -0.32671209D-02        UMP4(SDTQ)= -0.75016068045D+02
        # Energy for most substitutions is used only (SDTQ by default)
        if line[34:43] == "UMP4(DQ)=":
            mp4energy = utils.float(line.split("=")[2])
            line = file_handler.virtual_next()
            if line[34:44] == "UMP4(SDQ)=":
                mp4energy = utils.float(line.split("=")[2])
                line = file_handler.virtual_next()
                if line[34:45] == "UMP4(SDTQ)=":
                    mp4energy = utils.float(line.split("=")[2])
            this_mpenergies[-1].append(mp4energy)
            return {mpenergies.__name__: np.array(this_mpenergies)}

        # Example MP5 output line:
        #  DEMP5 =  -0.11048812312D-02 MP5 =  -0.75017172926D+02
        if line[29:34] == "MP5 =":
            this_mpenergies[-1].append(utils.float(line.split("=")[2]))
            return {mpenergies.__name__: np.array(this_mpenergies)}
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in mpenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(mpenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
            if constructed_data:
                print(constructed_data)
        return constructed_data
