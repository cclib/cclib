# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

from typing import Optional
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
    atomcharges = dict()
    atomspins = dict()
    line = file_handler.last_line
    has_spins = "AND SPIN POPULATIONS" in line

    file_handler.skip_lines(["dashes"], virtual=True)

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
        file_handler.skip_lines(
            ["d", "b", "Total integrated alpha density", "Total integrated beta density", "header"],
            virtual=True,
        )
    else:
        raise RuntimeError(f"unknown chargestype: {chargestype}")

    charges = []
    spins = []

    line = file_handler.virtual_next()
    while not should_stop(line):
        # Don't add point charges or embedding potentials.
        if "Q :" not in line:
            charges.append(float(line[start:stop]))
            if has_spins:
                spins.append(float(line[stop:]))
        line = file_handler.virtual_next()

    atomcharges[chargestype] = charges
    if has_spins:
        atomspins[chargestype] = spins
    return atomcharges, atomspins


def qchem_parse_charge_section(file_handler, chargetype):
    """Parse the population analysis charge block."""
    atomcharges = dict()
    atomspins = dict()
    charges = []
    spins = []
    line = file_handler.last_line
    file_handler.skip_lines(["blank"], virtual=True)
    line = file_handler.virtual_next()
    has_spins = False
    if "Spin" in line:
        has_spins = True
    file_handler.skip_lines(["dashes"], virtual=True)
    line = file_handler.virtual_next()

    while list(set(line.strip())) != ["-"]:
        elements = line.split()
        charge = utils.float(elements[2])
        charges.append(charge)
        if has_spins:
            spin = utils.float(elements[3])
            spins.append(spin)
        line = file_handler.virtual_next()

    atomcharges[chargetype] = charges
    if has_spins:
        atomspins[chargestype] = spins
    return atomcharges, atomspins


class atomcharges(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["ORCA", "NBO", "qchem"]

    @staticmethod
    def ORCA(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_charge_data = None
        constructed_spin_data = None

        # ORCA will print atomic charges along with the spin populations,
        #   so care must be taken about choosing the proper column.
        # Population analyses are performed usually only at the end
        #   of a geometry optimization or other run, so we want to
        #   leave just the final atom charges.
        # Here is an example for Mulliken charges:
        # --------------------------------------------
        # MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS
        # --------------------------------------------
        #    0 H :    0.126447    0.002622
        #    1 C :   -0.613018   -0.029484
        #    2 H :    0.189146    0.015452
        #    3 H :    0.320041    0.037434
        # ...
        # Sum of atomic charges         :   -0.0000000
        # Sum of atomic spin populations:    1.0000000
        if line[:23] == "MULLIKEN ATOMIC CHARGES":
            constructed_charge_data, constructed_spin_data = orca_parse_charge_section(
                file_handler, "mulliken"
            )
        # Things are the same for Lowdin populations, except that the sums
        #   are not printed (there is a blank line at the end).
        if line[:22] == "LOEWDIN ATOMIC CHARGES":
            constructed_charge_data, constructed_spin_data = orca_parse_charge_section(
                file_handler, "lowdin"
            )
        # ------------------
        # HIRSHFELD ANALYSIS
        # ------------------
        #
        # Total integrated alpha density =    142.999988722
        # Total integrated beta density  =    142.999988722
        #
        #   ATOM     CHARGE      SPIN
        #    0 H    0.157924    0.000000
        #    1 O   -0.209542    0.000000
        #    2 C    0.030659    0.000000
        # ...
        #   TOTAL  -0.999977    0.000000
        if line[:18] == "HIRSHFELD ANALYSIS":
            constructed_charge_data, constructed_spin_data = orca_parse_charge_section(
                file_handler, "hirshfeld"
            )
        # CHELPG Charges
        # --------------------------------
        #  0   C   :       0.363939
        #  1   H   :       0.025695
        # ...
        # --------------------------------
        # Total charge:    -0.000000
        # --------------------------------
        if line.startswith("CHELPG Charges"):
            constructed_charge_data, constructed_spin_data = orca_parse_charge_section(
                file_handler, "chelpg"
            )
        # TODO handle atomspins
        constructed_data = dict()
        if constructed_charge_data:
            if ccdata.atomcharges:
                constructed_data["atomcharges"] = {**ccdata.atomcharges, **constructed_charge_data}
            else:
                constructed_data["atomcharges"] = {**constructed_charge_data}
        if constructed_spin_data:
            if ccdata.atomspins:
                constructed_data["atomspins"] = {**ccdata.atomspins, **constructed_spin_data}
            else:
                constructed_data["atomspins"] = {**constructed_spin_data}
        if constructed_data:
            return constructed_data
        return None

    @staticmethod
    def NBO(file_handler, ccdata) -> Optional[dict]:
        atomcharges = dict()
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        charges = None
        line = file_handler.last_line
        if "  Atom No    Charge" in line:
            parsed_charges = []
            line = file_handler.virtual_next()
            line = file_handler.virtual_next()
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
            atomcharges["nbo"] = parsed_charges
        constructed_data = dict()
        if atomcharges != dict():
            if ccdata.atomcharges:
                constructed_data["atomcharges"] = {**ccdata.atomcharges, **atomcharges}
            else:
                constructed_data["atomcharges"] = {**atomcharges}
            return constructed_data
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> list | None:
        line = file_handler.last_line
        constructed_charge_data = None
        constructed_spin_data = None
        constructed_data = dict()
        if "Ground-State Mulliken Net Atomic Charges" in line:
            constructed_charge_data, constructed_spin_data = qchem_parse_charge_section(
                inputfile, "mulliken"
            )
        if "Hirshfeld Atomic Charges" in line:
            constructed_charge_data, constructed_spin_data = qchem_parse_charge_section(
                inputfile, "hirshfeld"
            )
        if "Charge Model 5" in line:
            constructed_charge_data, constructed_spin_data = qchem_parse_charge_section(
                inputfile, "cm5"
            )
        if "Ground-State ChElPG Net Atomic Charges" in line:
            constructed_charge_data, constructed_spin_data = qchem_parse_charge_section(
                inputfile, "chelpg"
            )
        if "Merz-Kollman ESP Net Atomic Charges" in line:
            constructed_charge_data, constructed_spin_data = qchem_parse_charge_section(
                inputfile, "esp"
            )
        if "Merz-Kollman RESP Net Atomic Charges" in line:
            constructed_charge_data, constructed_spin_data = qchem_parse_charge_section(
                inputfile, "resp"
            )
        constructed_data = dict()
        if constructed_charge_data:
            if ccdata.atomcharges:
                constructed_data["atomcharges"] = {**ccdata.atomcharges, **constructed_charge_data}
            else:
                constructed_data["atomcharges"] = {**constructed_charge_data}
        if constructed_spin_data:
            if ccdata.atomspins:
                constructed_data["atomspins"] = {**ccdata.atomspins, **constructed_spin_data}
            else:
                constructed_data["atomspins"] = {**constructed_spin_data}
        if constructed_data:
            return constructed_data
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in atomcharges.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(atomcharges, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
