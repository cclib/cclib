# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from cclib.parser_properties import utils
from cclib.parser_properties.base_parser import base_parser

import numpy as np


def orca_parse_charge_section(file_handler, chargestype):
    """Parse a charge section

    Parameters
    ----------
    file_handler :
      generates lines
    chargestype : str
      what type of charge we're dealing with, must be one of
      'mulliken', 'lowdin', 'chelpg' or 'hirshfeld'
    """
    line = file_handler.last_line
    has_spins = "AND SPIN POPULATIONS" in line

    self.skip_line(inputfile, "dashes")

    # depending on chargestype, decide when to stop parsing lines
    # start, stop - indices for slicing lines and grabbing values
    # should_stop: when to stop parsing
    if chargestype == "mulliken":
        should_stop = lambda x: x.startswith("Sum of atomic charges")
        start, stop = 8, 20
    elif chargestype == "lowdin":
        should_stop = lambda x: not bool(x.strip())
        start, stop = 8, 20
    elif chargestype == "chelpg":
        should_stop = lambda x: x.startswith("---")
        start, stop = 11, 26
    elif chargestype == "hirshfeld":
        should_stop = lambda x: not bool(x.strip())
        start, stop = 9, 18
        self.skip_lines(
            inputfile,
            ["d", "b", "Total integrated alpha density", "Total integrated beta density", "header"],
        )
    else:
        raise RuntimeError(f"unknown chargestype: {chargestype}")

    charges = []
    spins = []

    line = next(inputfile)
    while not should_stop(line):
        # Don't add point charges or embedding potentials.
        if "Q :" not in line:
            charges.append(float(line[start:stop]))
            if has_spins:
                spins.append(float(line[stop:]))
        line = next(inputfile)

    self.atomcharges[chargestype] = charges
    if has_spins:
        self.atomspins[chargestype] = spins


class atomcharges(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["orca", "NBO"]

    @staticmethod
    def orca(file_handler, ccdata) -> list | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_data = None
        if line.strip() == "Isotopes and Nuclear Properties:":
            constructed_data = []
            line = file_handler.virtual_next()
            while line[1:16] != "Leave Link  101":
                if line[1:8] == "AtmWgt=":
                    constructed_data.extend(list(map(float, line.split()[1:])))
                line = file_handler.virtual_next()
            return constructed_data
        return None

    @staticmethod
    def NBO(file_handler, ccdata) -> list | None:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        charges = None
        line = file_handler.last_line
        if "  Atom No    Charge" in line:
            if not hasattr(self, "atomcharges"):
                self.atomcharges = dict()

            line = file_handler.virtual_next()
            line = file_handler.virtual_next()

            parsed_charges = []

            while "==============" not in line:
                population_analysis = line.split()

                atom = population_analysis[0]
                no = int(population_analysis[1])
                natural_charge = float(population_analysis[2])
                core = float(population_analysis[3])
                valence = float(population_analysis[4])
                rydberg = float(population_analysis[5])
                total = float(population_analysis[6])
                parsed_charges.append(natural_charge)
                line = file_handler.virtual_next()

            self.atomcharges["nbo"] = parsed_charges
        if parsed_charges != []:
            charges = parsed_charges
        return charges

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> list | None:
        constructed_data = None
        if program in atomcharges.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(atomcharges, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
