# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers.base_parser import base_parser


class nmo(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "qchem", "ORCA"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_data = None
        if line[1:7] == "NBsUse":
            # For counterpoise fragment, skip these lines.
            #
            # AED: not sure how to handle these cases at the moment
            # if ccdata.counterpoise != 1:
            #   return
            # For ONIOM calcs, ignore this section in order to bypass assertion failure.
            # if ccdata.oniom:
            #   return
            constructed_data = int(line.split("=")[1].split()[0])
            return {nmo.__name__: constructed_data}
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> Optional[dict]:
        # Don't reset this all the time
        dependency_list = ["nmo"]
        if base_parser.check_dependencies(dependency_list, ccdata, "nmo"):
            return None

        dependency_list = ["moenergies"]
        if base_parser.check_dependencies(dependency_list, ccdata, "nmo"):
            return {nmo.__name__: len(ccdata.monergies[0])}
        return None

    @staticmethod
    def ORCA(file_handler, ccdata) -> Optional[dict]:
        # Don't reset this all the time
        dependency_list = ["nmo"]
        if base_parser.check_dependencies(dependency_list, ccdata, "nmo"):
            return None
        dependency_list = ["nbasis"]

        if base_parser.check_dependencies(dependency_list, ccdata, "nmo"):
            return {nmo.__name__: ccdata.nbasis}
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in nmo.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(nmo, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
