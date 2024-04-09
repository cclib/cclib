# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser


class scftargets(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> dict | None:
        line = file_handler.last_line
        if not ccdata.BOMD and line[1:44] == "Requested convergence on RMS density matrix":
            constructed_scftargets = ccdata.scftargets
            if constructed_scftargets is None:
                constructed_scftargets = []
            elif type(ccdata.scftargets) == type(numpy.array([])):
                # This case can happen with ONIOM which are mixed SCF
                # and semi-empirical
                constructed_scftargets = []

            curr_scftargets = []

            # The RMS density matrix.
            curr_scftargets.append(utils.float(line.split("=")[1].split()[0]))
            line = file_handler.virtual_next()
            # The MAX density matrix.
            curr_scftargets.append(utils.float(line.strip().split("=")[1][:-1]))
            line = file_handler.virtual_next()

            # For G03, there's also the energy (not for G98).
            if line[1:10] == "Requested":
                curr_scftargets.append(utils.float(line.strip().split("=")[1][:-1]))

            constructed_scftargets.append(curr_scftargets)
            return {scftargets.__name__: constructed_scftargets}

        # Extract SCF convergence information (QM calcs).
        if line[1:10] == "Cycle   1":
            constructed_scftargets = ccdata.scftargets
            if constructed_scftargets is None:
                constructed_scftargets = []

        if line[1:4] == "It=":
            curr_scftargets = numpy.array([1e-7], "d")  # This is the target value for the rms

            if line.find(" Energy") == 0:
                constructed_scftargets = curr_scftargets

            return {scftargets.__name__: constructed_scftargets}
        return None

    @staticmethod
    def parse(file_handler, program, ccdata) -> dict | None:
        constructed_data = None
        if program in scftargets.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(scftargets, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
