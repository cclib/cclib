# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class moenergies(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "ORCA"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if line[1:6] == "Alpha" and line.find("eigenvalues") >= 0:
            # For counterpoise fragments, skip these lines.
            if getattr(ccdata, "moenergies") is None:
                constructed_moenergies = [[]]
                while line.find("Alpha") == 1:
                    part = line[28:]
                    i = 0
                    while i * 10 + 4 < len(part):
                        s = part[i * 10 : (i + 1) * 10]
                        try:
                            x = utils.float(s)
                        except ValueError:
                            x = np.nan
                        constructed_moenergies[0].append(x)
                        i += 1
                    line = file_handler.virtual_next()
                if line.find("Beta") == 2:
                    constructed_moenergies.append([])

                while line.find("Beta") == 2:
                    part = line[28:]
                    i = 0
                    while i * 10 + 4 < len(part):
                        x = part[i * 10 : (i + 1) * 10]
                        constructed_moenergies[1].append(utils.float(x))
                        i += 1
                    line = file_handler.virtual_next()

                constructed_moenergies = [np.array(x, "d") for x in constructed_moenergies]
                return {moenergies.__name__: constructed_moenergies}
        return None

    @staticmethod
    def ORCA(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if line[0:16] == "ORBITAL ENERGIES":
            file_handler.skip_lines(["d", "text", "text"], virtual=True)
            #constructed_mooccnos = [[]]
            constructed_moenergies = [[]]
            #constructed_mosyms = [[]]
            line = file_handler.virtual_next()
            while len(line) > 20:  # restricted calcs are terminated by ------
                info = line.split()
                #mooccno = int(float(info[1]))
                moenergy = float(info[2])
                #mosym = "A"
                #if self.uses_symmetry:
                #    mosym = self.normalisesym(info[4].split("-")[1])
                #self.mooccnos[0].append(mooccno)
                constructed_moenergies[0].append(moenergy)
                #self.mosyms[0].append(mosym)
                line = file_handler.virtual_next()
            line = file_handler.virtual_next()
            # handle beta orbitals for UHF
            if line[17:35] == "SPIN DOWN ORBITALS":
                file_handler.skip_lines(["text"], virtual=True)
                #constructed_mooccnos.append([])
                constructed_moenergies.append([])
                #constructed_mosyms.append([])
                line = file_handler.virtual_next()
                while len(line) > 20:  # actually terminated by ------
                    info = line.split()
                    #mooccno = int(float(info[1]))
                    moenergy = float(info[2])
                    #mosym = "A"
                    #if self.uses_symmetry:
                    #    mosym = self.normalisesym(info[4].split("-")[1])
                    #constructed_mooccnos[1].append(mooccno)
                    constructed_moenergies[1].append(moenergy)
                    #constructed_mosyms[1].append(mosym)
                    line = file_handler.virtual_next()
            constructed_moenergies = [np.array(x, "d") for x in constructed_moenergies]
            return {moenergies.__name__: constructed_moenergies}
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in moenergies.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(moenergies, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
