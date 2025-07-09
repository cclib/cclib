# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy


class geotargets(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["psi4", "gaussian"]

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if (line.strip() == "==> Convergence Check <==") and (
            getattr(ccdata, "finite_difference") is None
        ):
            this_geotargets = []
            file_handler.skip_lines(["b", "units"])
            if ccdata.package_version.minor >= 7:
                file_handler.skip_lines(["b"])
            file_handler.skip_lines(["comment"])
            if ccdata.package_version.minor >= 7:
                file_handler.skip_lines(["b"])
            file_handler.skip_lines(["dash+tilde", "header", "dash+tilde"])

            # These are the position in the line at which numbers should start.
            starts = [27, 41, 55, 69, 83]

            criteria = file_handler.virtual_next()
            for istart in starts:
                if criteria[istart : istart + 9].strip():
                    this_geotargets.append(float(criteria[istart : istart + 9]))
                else:
                    this_geotargets.append(numpy.inf)

            file_handler.skip_lines(["dashes"])

            values = file_handler.virtual_next()
            step = int(values.split()[0])
            geovalues = []
            for istart in starts:
                if values[istart : istart + 9].strip():
                    geovalues.append(float(values[istart : istart + 9]))

            if step == 1:
                self.optstatus[-1] += data.ccData.OPT_NEW  # noqa: F821

            # This assertion may be too restrictive, but we haven't seen the geotargets change.
            # If such an example comes up, update the value since we're interested in the last ones.
            if not hasattr(self, "geotargets"):  # noqa: F821
                self.geotargets = geotargets  # noqa: F821
            else:
                assert self.geotargets == geotargets  # noqa: F821

            self.append_attribute("geovalues", geovalues)  # noqa: F821

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
                this_scftargets = ccdata.scftargets
            this_scftargets.append([etarget, dtarget])

            return {geotargets.__name__: this_scftargets}
        return None

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        # Geometry convergence information.
        if line[49:59] == "Converged?":
            parsed_geotargets = []
            if not hasattr(ccdata, "geotargets"):
                parsed_geotargets = numpy.array([0.0, 0.0, 0.0, 0.0], "d")
            allconverged = True
            newlist = [0] * 4
            for i in range(4):
                line = file_handler.virtual_next()
                parts = line.split()
                if "NO" in parts[-1]:
                    allconverged = False  # noqa: F841
                try:
                    value = utils.float(parts[2])
                except ValueError:
                    pass  # todo logging
                    # self.logger.error(
                    #    "Problem parsing the value for geometry optimisation: %s is not a number.",
                    #    parts[2],
                    # )
                else:
                    newlist[i] = value
                parsed_geotargets[i] = utils.float(parts[3])
            return {geotargets.__name__: parsed_geotargets}
        return None

    @staticmethod
    def parse(file_handler, program, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in geotargets.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(geotargets, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
