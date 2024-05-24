# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


class geotargets(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["psi4"]


    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if ((line.strip() == "==> Convergence Check <==")
            and (getattr(ccdata, "finite_difference") is None)
        ):
            this_geotargets = []
            file_handler.skip_lines(["b", "units"])
            if ccdata.package_version.minor >= 7:
                file_handler.self.skip_line(["b"])
            file_handler.skip_lines(["comment"])
            if ccdata.package_version.minor >= 7:
                file_handler.skip_lines(["b"])
            file_handler.skip_lines(["dash+tilde", "header", "dash+tilde"])

            # These are the position in the line at which numbers should start.
            starts = [27, 41, 55, 69, 83]

            criteria = file_next.virtual_next()
            for istart in starts:
                if criteria[istart : istart + 9].strip():
                    this_geotargets.append(float(criteria[istart : istart + 9]))
                else:
                    this_geotargets.append(numpy.inf)

            file_handler.skip_line(, "dashes")

            values = file_handler.virtual_next()
            step = int(values.split()[0])
            geovalues = []
            for istart in starts:
                if values[istart : istart + 9].strip():
                    geovalues.append(float(values[istart : istart + 9]))

            if step == 1:
                self.optstatus[-1] += data.ccData.OPT_NEW

            # This assertion may be too restrictive, but we haven't seen the geotargets change.
            # If such an example comes up, update the value since we're interested in the last ones.
            if not hasattr(self, "geotargets"):
                self.geotargets = geotargets
            else:
                assert self.geotargets == geotargets

            if not hasattr(self, "geovalues"):
                self.geovalues = []
            self.geovalues.append(geovalues)

            line = file_handler.virtual_next()
            while line.strip():
                if "Energy threshold" in line:
                    etarget = float(line.split()[-1])
                if "Density threshold" in line:
                    dtarget = float(line.split()[-1])
                line = file_handler.virtual_next()

            if getattr(ccdata, "scftargets") is None:
                this_scftargets = []
            else:
                this_scftargers = ccdata.scftargets
            this_scftargets.append([etarget, dtarget])

            return {geotargets.__name__:this_scftargets}
        return None

        # This section prints contraction information before the atomic basis set functions and
    @staticmethod
    def parse(file_handler, program, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in scftargets.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(scftargets, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
