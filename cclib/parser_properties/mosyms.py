# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from cclib.parser_properties import utils
from cclib.parser_properties.base_parser import base_parser

import numpy as np


def gaussian_normalizesym(label):
    """Use standard symmetry labels instead of Gaussian labels.

    To normalise:
    (1) If label is one of [SG, PI, PHI, DLTA], replace by [sigma, pi, phi, delta]
    (2) replace any G or U by their lowercase equivalent
    """
    # note: DLT must come after DLTA
    greek = [("SG", "sigma"), ("PI", "pi"), ("PHI", "phi"), ("DLTA", "delta"), ("DLT", "delta")]
    for k, v in greek:
        if label.startswith(k):
            tmp = label[len(k) :]
            label = v
            if tmp:
                label = f"{v}.{tmp}"

    ans = label.replace("U", "u").replace("G", "g")
    return ans


class mosyms(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> dict | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_data = None
        if line[1:20] == "Orbital symmetries:":
            line = file_handler.virtual_next()
            constructed_data = [[]]
            unres = False
            if line.find("Alpha Orbitals") == 1:
                unres = True
                line = file_handler.virtual_next()
            i = 0
            while len(line) > 18 and line[17] == "(":
                parts = line[17:].split()
                for x in parts:
                    constructed_data[0].append(gaussian_normalizesym(x.strip("()")))
                    i += 1
                line = file_handler.virtual_next()
            if unres:
                line = file_handler.virtual_next()
                # Repeat with beta orbital information
                i = 0
                constructed_data.append([])
                while len(line) > 18 and line[17] == "(":
                    parts = line[17:].split()
                    for x in parts:
                        constructed_data[1].append(gaussian_normalisesym(x.strip("()")))
                        i += 1
                    line = file_handler.virtual_next()
            # Some calculations won't explicitly print the number of basis sets used,
            # and will occasionally drop some without warning. We can infer the number,
            # however, from the MO symmetries printed here. Specifically, this fixes
            # regression Gaussian/Gaussian09/dvb_sp_terse.log (#23 on github).
            return {mosyms.__name__: constructed_data}
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> dict | None:
        constructed_data = None
        if program in mosyms.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(mosyms, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
