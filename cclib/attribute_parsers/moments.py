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

    known_codes = ["psi4", "gaussian"]

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
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        # Dipole moment
        # e.g. from G09
        #  Dipole moment (field-independent basis, Debye):
        #    X=              0.0000    Y=              0.0000    Z=              0.0930
        # e.g. from G03
        #     X=     0.0000    Y=     0.0000    Z=    -1.6735  Tot=     1.6735
        # need the "field independent" part - ONIOM and other calc use diff formats
        if line[1:39] == "Dipole moment (field-independent basis":
            reference = [0.0, 0.0, 0.0]
            parsed_moments = []

            line = file_handler.virtual_next()
            tokens = line.split()
            # split - dipole would need to be *huge* to fail a split
            # and G03 and G09 use different spacing
            if len(tokens) >= 6:
                dipole = (float(tokens[1]), float(tokens[3]), float(tokens[5]))
                if getattr(ccdata, "moments") is None:
                    parsed_moments = [reference, dipole]
                else:
                    parsed_moments = ccdata.moments
                    parsed_moments.append(dipole)
                return {moments.__name__: parsed_moments}

        if line[1:43] == "Quadrupole moment (field-independent basis":
            # e.g. (g09)
            # Quadrupole moment (field-independent basis, Debye-Ang):
            #   XX=             -6.1213   YY=             -4.2950   ZZ=             -5.4175
            #   XY=              0.0000   XZ=              0.0000   YZ=              0.0000
            # or from g03
            #   XX=    -6.1213   YY=    -4.2950   ZZ=    -5.4175
            quadrupole = {}
            for j in range(2):  # two rows
                line = file_handler.virtual_next()
                if line[22] == "=":  # g03 file
                    for i in (1, 18, 35):
                        quadrupole[line[i : i + 4]] = float(line[i + 5 : i + 16])
                else:
                    for i in (1, 27, 53):
                        quadrupole[line[i : i + 4]] = float(line[i + 5 : i + 25])

            lex = sorted(quadrupole.keys())
            quadrupole = [quadrupole[key] for key in lex]
            parsed_moments = []
            if getattr(ccdata, "moments") is None or len(ccdata.moments) < 2:
                # logger.warning("Found quadrupole moments but no previous dipole")
                reference = [0.0, 0.0, 0.0]
                parsed_moments = [reference, None, quadrupole]
                return {moments.__name__: parsed_moments}
            else:
                parsed_moments = ccdata.moments
                if len(parsed_moments) == 2:
                    parsed_moments.append(quadrupole)
                    return {moments.__name__: parsed_moments}
                else:
                    assert parsed_moments[2] == quadrupole

        if line[1:41] == "Octapole moment (field-independent basis":
            # e.g.
            # Octapole moment (field-independent basis, Debye-Ang**2):
            #  XXX=              0.0000  YYY=              0.0000  ZZZ=             -0.1457  XYY=              0.0000
            #  XXY=              0.0000  XXZ=              0.0136  XZZ=              0.0000  YZZ=              0.0000
            #  YYZ=             -0.5848  XYZ=              0.0000
            octapole = {}
            for j in range(2):  # two rows
                line = file_handler.virtual_next()
                if line[22] == "=":  # g03 file
                    for i in (1, 18, 35, 52):
                        octapole[line[i : i + 4]] = float(line[i + 5 : i + 16])
                else:
                    for i in (1, 27, 53, 79):
                        octapole[line[i : i + 4]] = float(line[i + 5 : i + 25])

            # last line only 2 moments
            line = file_handler.virtual_next()
            if line[22] == "=":  # g03 file
                for i in (1, 18):
                    octapole[line[i : i + 4]] = float(line[i + 5 : i + 16])
            else:
                for i in (1, 27):
                    octapole[line[i : i + 4]] = float(line[i + 5 : i + 25])

            lex = sorted(octapole.keys())
            octapole = [octapole[key] for key in lex]

            parsed_moments = []
            if getattr(ccdata, "moments") is None or len(ccdata.moments) < 3:
                # logger.warning("Found quadrupole moments but no previous dipole")
                reference = [0.0, 0.0, 0.0]
                parsed_moments = [reference, None, None, octapole]
                return {moments.__name__: parsed_moments}
            else:
                parsed_moments = ccdata.moments
                if len(parsed_moments) == 3:
                    parsed_moments.append(octapole)
                    return {moments.__name__: parsed_moments}
                else:
                    assert parsed_moments[3] == octapole

        if line[1:20] == "Hexadecapole moment":
            # e.g.
            # Hexadecapole moment (field-independent basis, Debye-Ang**3):
            # XXXX=             -3.2614 YYYY=             -6.8264 ZZZZ=             -4.9965 XXXY=              0.0000
            # XXXZ=              0.0000 YYYX=              0.0000 YYYZ=              0.0000 ZZZX=              0.0000
            # ZZZY=              0.0000 XXYY=             -1.8585 XXZZ=             -1.4123 YYZZ=             -1.7504
            # XXYZ=              0.0000 YYXZ=              0.0000 ZZXY=              0.0000
            hexadecapole = {}
            # read three lines worth of 4 moments per line
            for j in range(3):
                line = file_handler.virtual_next()
                if line[22] == "=":  # g03 file
                    for i in (1, 18, 35, 52):
                        hexadecapole[line[i : i + 4]] = float(line[i + 5 : i + 16])
                else:
                    for i in (1, 27, 53, 79):
                        hexadecapole[line[i : i + 4]] = float(line[i + 5 : i + 25])

            # last line only 3 moments
            line = file_handler.virtual_next()
            if line[22] == "=":  # g03 file
                for i in (1, 18, 35):
                    hexadecapole[line[i : i + 4]] = float(line[i + 5 : i + 16])
            else:
                for i in (1, 27, 53):
                    hexadecapole[line[i : i + 4]] = float(line[i + 5 : i + 25])

            lex = sorted(hexadecapole.keys())
            hexadecapole = [hexadecapole[key] for key in lex]

            parsed_moments = []
            if getattr(ccdata, "moments") is None or len(ccdata.moments) < 4:
                # logger.warning("Found quadrupole moments but no previous dipole")
                reference = [0.0, 0.0, 0.0]
                parsed_moments = [reference, None, None, None, hexadecapole]
                return {moments.__name__: parsed_moments}
            else:
                parsed_moments = ccdata.moments
                if len(parsed_moments) == 4:
                    parsed_moments.append(hexadecapole)
                    return {moments.__name__: parsed_moments}
                else:
                    assert parsed_moments[4] == hexadecapole
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
