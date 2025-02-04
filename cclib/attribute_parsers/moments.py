# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy


class moments(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["psi4"]

    @staticmethod
    def psi4(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        if line.strip() == "Dipole Moment: (a.u.)":
            line = file_handler.virtual_next()
            tokens = line.split()
            dipole = utils.convertor(
                numpy.array([float(tokens[1]), float(tokens[3]), float(tokens[5])]),
                "ebohr",
                "Debye",
            )

            # AED: I don't know how to handle this in version 2 yet
            # this line below is an ugly fix to always assume the origin is at zero until i figure it out.

            origin = numpy.array([0.0, 0.0, 0.0])
            if getattr(ccdata, "moments") is None:
                # Old versions of Psi4 don't print the origin; assume
                # it's at zero.
                if getattr(ccdata, "origin") is None:
                    # AED: I don't know how to handle this in version 2 yet
                    origin = numpy.array([0.0, 0.0, 0.0])
                return {moments.__name__: [origin, dipole]}
            else:
                try:
                    assert numpy.all(ccdata.moments[1] == dipole)
                except AssertionError:
                    return {moments.__name__: [origin, dipole]}
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in moments.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(moments, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
