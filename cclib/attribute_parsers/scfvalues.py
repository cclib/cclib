# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import re
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


# Extract only well-formed numbers in scientific notation.
class scfvalues(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        re_scinot = re.compile(r"(\w*)=\s*(-?\d\.\d{2}D[+-]\d{2})")
        line = file_handler.last_line
        # Extract SCF convergence information (QM calcs).
        if line[1:10] == "Cycle   1":
            if getattr(ccdata, "scfvalues") is not None:
                constructed_scfvalues = ccdata.scfvalues
            else:
                constructed_scfvalues = []

            this_scfvalue = []
            line = file_handler.virtual_next()
            while line.find("SCF Done") == -1:
                #  RMSDP=3.74D-06 MaxDP=7.27D-05 DE=-1.73D-07 OVMax= 3.67D-05
                # or
                #  RMSDP=1.13D-05 MaxDP=1.08D-04              OVMax= 1.66D-04
                if line.find(" RMSDP") == 0:
                    # Fields of interest:
                    # RMSDP
                    # MaxDP
                    # (DE) -> Only add the energy if it's a target criteria

                    matches = re_scinot.findall(line)
                    matches = {match[0]: utils.float(match[1]) for match in matches}
                    scfvalues_step = [matches.get("RMSDP", np.nan), matches.get("MaxDP", np.nan)]
                    if (ccdata.scftargets is not None) and len(ccdata.scftargets[0]) == 3:
                        scfvalues_step.append(matches.get("DE", np.nan))
                    this_scfvalue.append(scfvalues_step)
                try:
                    line = file_handler.virtual_next()
                # May be interupted by EOF.
                except StopIteration:
                    break
            constructed_scfvalues.append(np.array(this_scfvalue))

            return {scfvalues.__name__: constructed_scfvalues}
        return None

    @staticmethod
    def parse(file_handler, program, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in scfvalues.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(scfvalues, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
