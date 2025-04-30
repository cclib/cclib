# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


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
                np.array([float(tokens[1]), float(tokens[3]), float(tokens[5])]), "ebohr", "Debye"
            )

            if getattr(ccdata, "moments") is None:
                # Old versions of Psi4 don't print the origin; assume
                # it's at zero.
                if "origin" in ccdata.parser_state.keys():
                    origin = ccdata.parser_state["origin"]
                else:
                    origin = np.array([0.0, 0.0, 0.0])

                return {moments.__name__: [origin, dipole]}
            else:
                try:
                    assert np.all(ccdata.moments[1] == dipole)
                except AssertionError:
                    return {moments.__name__: [origin, dipole]}

        if line.strip() == "Multipole Moments:":
            if "origin" in ccdata.parser_state.keys():
                origin = ccdata.parser_state["origin"]
            else:
                origin = np.array([0.0, 0.0, 0.0])
            file_handler.skip_lines(["b", "d", "header", "d", "b"])

            # The reference used here should have been printed somewhere
            # before the properties and parsed above.
            this_moments = [origin]

            line = file_handler.virtual_next()
            while "----------" not in line.strip():
                rank = int(line.split()[2].strip("."))

                multipole = []
                line = file_handler.virtual_next()
                while line.strip():
                    tokens = line.split()
                    if tokens[0] in ("Magnitude", "Traceless"):
                        line = file_handler.virtual_next()
                        continue
                    value = float(tokens[-1])
                    fromunits = f"ebohr{(rank > 1) * f'{int(rank)}'}"
                    tounits = f"Debye{(rank > 1) * '.ang'}{(rank > 2) * f'{int(rank - 1)}'}"
                    value = utils.convertor(value, fromunits, tounits)
                    multipole.append(value)

                    line = file_handler.virtual_next()
                this_moments.append(np.array(multipole))
                line = file_handler.virtual_next()
            if getattr(ccdata, "moments") is None:
                return {moments.__name__: this_moments}
            else:
                existing_moments_list = getattr(ccdata, "moments")
                for m_idx, m in enumerate(this_moments):
                    if len(ccdata.moments) <= m_idx:
                        existing_moments_list.append(m)
                    else:
                        assert np.allclose(existing_moments_list[m_idx], m, atol=1.0e4)
                return {moments.__name__: existing_moments_list}
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
