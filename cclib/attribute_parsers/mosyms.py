# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
import re
from typing import Optional

from cclib.attribute_parsers import utils
from cclib.attribute_parsers.base_parser import base_parser

import numpy as np


def gaussian_normalizesym(label):
    """Use standard symmetry labels instead of Gaussian labels.

    To normalise:
    (1) If label is one of [SG, PI, PHI, DLTA], replace by [sigma, pi, phi, delta]
    (2) replace any G or U by their lowercase equivalent
    """
    # note: DLT must come after DLTA
    greek = [("SG", "sigma"), ("PI", "pi"), ("PHI", "phi"), ("DLTA", "delta"), ("DLT", "delta")]
    for k, v in greek:
        if label.startswith(k):
            tmp = label[len(k) :]
            label = v
            if tmp:
                label = f"{v}.{tmp}"

    ans = label.replace("U", "u").replace("G", "g")
    return ans

def orca_normalizesym(label):
        """ORCA does not require normalizing symmetry labels."""
        return label

 

def qchem_parse_orbital_energies_and_symmetries(file_handler):
    """Parse the 'Orbital Energies (a.u.)' block appearing after SCF converges,
    which optionally includes MO symmetries. Based upon the
    Occupied/Virtual labeling, the HOMO is also parsed.
    """
    re_dashes_and_spaces = re.compile(r"^[\s-]+$")
    energies = []
    symbols = []

    line = next(inputfile)  # noqa: F821
    line = file_handler.virtual_next()
    # Sometimes Q-Chem gets a little confused...
    while "MOs" not in line:
        line = file_handler.virtual_next()
    line = file_handler.virtual_next()

    # The end of the block is either a blank line or only dashes.
    while not re_dashes_and_spaces.search(line) and "Warning : Irrep of orbital" not in line:
        if "Occupied" in line or "Virtual" in line:
            # A nice trick to find where the HOMO is.
            if "Virtual" in line:
                homo = len(energies) - 1
            line = file_handler.virtual_next()
        tokens = line.split()
        # If the line contains letters, it must be the MO
        # symmetries. Otherwise, it's the energies.
        if re.search("[a-zA-Z]", line):
            symbols.extend(tokens[1::2])
        else:
            for e in tokens:
                try:
                    energy = utils.convertor(utils.float(e), "hartree", "eV")
                except ValueError:
                    energy = np.nan
                energies.append(energy)
        line = file_handler.virtual_next()

    # MO symmetries are either not present or there is one for each MO
    # (energy).
    assert len(symbols) in (0, len(energies))
    return energies, symbols, homo


class mosyms(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian", "qchem", "ORCA"]

    @staticmethod
    def gaussian(file_handler, ccdata) -> Optional[dict]:
        # ccdata is "const" here and we don't need to modify it yet. The driver will set the attr
        line = file_handler.last_line
        constructed_data = None
        if line[1:20] == "Orbital symmetries:":
            line = file_handler.virtual_next()
            constructed_data = [[]]
            unres = False
            if line.find("Alpha Orbitals") == 1:
                unres = True
                line = file_handler.virtual_next()
            i = 0
            while len(line) > 18 and line[17] == "(":
                parts = line[17:].split()
                for x in parts:
                    constructed_data[0].append(gaussian_normalizesym(x.strip("()")))
                    i += 1
                line = file_handler.virtual_next()
            if unres:
                line = file_handler.virtual_next()
                # Repeat with beta orbital information
                i = 0
                constructed_data.append([])
                while len(line) > 18 and line[17] == "(":
                    parts = line[17:].split()
                    for x in parts:
                        constructed_data[1].append(gaussian_normalizesym(x.strip("()")))
                        i += 1
                    line = file_handler.virtual_next()
            # Some calculations won't explicitly print the number of basis sets used,
            # and will occasionally drop some without warning. We can infer the number,
            # however, from the MO symmetries printed here. Specifically, this fixes
            # regression Gaussian/Gaussian09/dvb_sp_terse.log (#23 on github).
            return {mosyms.__name__: constructed_data}
        return None

    @staticmethod
    def qchem(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        parsed_data = None
        # Molecular orbital energies and symmetries.
        if line.strip() == "Orbital Energies (a.u.) and Symmetries":
            #  --------------------------------------------------------------
            #              Orbital Energies (a.u.) and Symmetries
            #  --------------------------------------------------------------
            #
            #  Alpha MOs, Restricted
            #  -- Occupied --
            # -10.018 -10.018 -10.008 -10.008 -10.007 -10.007 -10.006 -10.005
            #   1 Bu    1 Ag    2 Bu    2 Ag    3 Bu    3 Ag    4 Bu    4 Ag
            #  -9.992  -9.992  -0.818  -0.755  -0.721  -0.704  -0.670  -0.585
            #   5 Ag    5 Bu    6 Ag    6 Bu    7 Ag    7 Bu    8 Bu    8 Ag
            #  -0.561  -0.532  -0.512  -0.462  -0.439  -0.410  -0.400  -0.397
            #   9 Ag    9 Bu   10 Ag   11 Ag   10 Bu   11 Bu   12 Bu   12 Ag
            #  -0.376  -0.358  -0.349  -0.330  -0.305  -0.295  -0.281  -0.263
            #  13 Bu   14 Bu   13 Ag    1 Au   15 Bu   14 Ag   15 Ag    1 Bg
            #  -0.216  -0.198  -0.160
            #   2 Au    2 Bg    3 Bg
            #  -- Virtual --
            #   0.050   0.091   0.116   0.181   0.280   0.319   0.330   0.365
            #   3 Au    4 Au    4 Bg    5 Au    5 Bg   16 Ag   16 Bu   17 Bu
            #   0.370   0.413   0.416   0.422   0.446   0.469   0.496   0.539
            #  17 Ag   18 Bu   18 Ag   19 Bu   19 Ag   20 Bu   20 Ag   21 Ag
            #   0.571   0.587   0.610   0.627   0.646   0.693   0.743   0.806
            #  21 Bu   22 Ag   22 Bu   23 Bu   23 Ag   24 Ag   24 Bu   25 Ag
            #   0.816
            #  25 Bu
            #
            #  Beta MOs, Restricted
            #  -- Occupied --
            # -10.018 -10.018 -10.008 -10.008 -10.007 -10.007 -10.006 -10.005
            #   1 Bu    1 Ag    2 Bu    2 Ag    3 Bu    3 Ag    4 Bu    4 Ag
            #  -9.992  -9.992  -0.818  -0.755  -0.721  -0.704  -0.670  -0.585
            #   5 Ag    5 Bu    6 Ag    6 Bu    7 Ag    7 Bu    8 Bu    8 Ag
            #  -0.561  -0.532  -0.512  -0.462  -0.439  -0.410  -0.400  -0.397
            #   9 Ag    9 Bu   10 Ag   11 Ag   10 Bu   11 Bu   12 Bu   12 Ag
            #  -0.376  -0.358  -0.349  -0.330  -0.305  -0.295  -0.281  -0.263
            #  13 Bu   14 Bu   13 Ag    1 Au   15 Bu   14 Ag   15 Ag    1 Bg
            #  -0.216  -0.198  -0.160
            #   2 Au    2 Bg    3 Bg
            #  -- Virtual --
            #   0.050   0.091   0.116   0.181   0.280   0.319   0.330   0.365
            #   3 Au    4 Au    4 Bg    5 Au    5 Bg   16 Ag   16 Bu   17 Bu
            #   0.370   0.413   0.416   0.422   0.446   0.469   0.496   0.539
            #  17 Ag   18 Bu   18 Ag   19 Bu   19 Ag   20 Bu   20 Ag   21 Ag
            #   0.571   0.587   0.610   0.627   0.646   0.693   0.743   0.806
            #  21 Bu   22 Ag   22 Bu   23 Bu   23 Ag   24 Ag   24 Bu   25 Ag
            #   0.816
            #  25 Bu
            #  --------------------------------------------------------------

            file_handler.skip_lines(["dashes"], virtual=True)
            line = file_handler.virtual_next()
            (energies_alpha, symbols_alpha, homo_alpha) = (
                qchem_parse_orbital_energies_and_symmetries(file_handler)
            )
            # Only look at the second block if doing an unrestricted calculation.
            # This might be a problem for ROHF/ROKS.
            # TODO check metadata if unrestricted and then parse beta
            # if self.unrestricted:
            #    (energies_beta, symbols_beta, homo_beta) = (
            #        qchem_parse_orbital_energies_and_symmetries(file_handler)
            #    )

            # For now, only keep the last set of MO energies, even though it is
            # printed at every step of geometry optimizations and fragment jobs.
            parsed_data = {mosyms.__name__: [symbols_alpha]}

            # if self.unrestricted:
            #    self.moenergies.append(np.array(energies_beta))
            #    self.homos.append(homo_beta)
            #    self.mosyms.append(symbols_beta)
        return parsed_data

    @staticmethod
    def ORCA(file_handler, ccdata) -> Optional[dict]:
        line = file_handler.last_line
        uses_symmetry = False
        if "uses_symmetry" in ccdata.parser_state:
            uses_symmetry = ccdata.parser_state["uses_symmetry"]
        if line[0:16] == "ORBITAL ENERGIES":
            file_handler.skip_lines(["d", "text", "text"], virtual=True)
            #constructed_mooccnos = [[]]
            #constructed_moenergies = [[]]
            constructed_mosyms = [[]]
            line = file_handler.virtual_next()
            while len(line) > 20:  # restricted calcs are terminated by ------
                info = line.split()
                #mooccno = int(float(info[1]))
                #moenergy = float(info[2])
                mosym = "A"
                if uses_symmetry:
                    mosym = orca_normalizesym(info[4].split("-")[1])
                constructed_mosyms[0].append(mosym)
                line = file_handler.virtual_next()
            line = file_handler.virtual_next()
            # handle beta orbitals for UHF
            if line[17:35] == "SPIN DOWN ORBITALS":
                file_handler.skip_lines(["text"], virtual=True)
                #constructed_mooccnos.append([])
                #constructed_moenergies.append([])
                constructed_mosyms.append([])
                line = file_handler.virtual_next()
                while len(line) > 20:  # actually terminated by ------
                    info = line.split()
                    #mooccno = int(float(info[1]))
                    #moenergy = float(info[2])
                    mosym = "A"
                    if uses_symmetry:
                        mosym = orca_normalizesym(info[4].split("-")[1])
                    #constructed_mooccnos[1].append(mooccno)
                    #constructed_moenergies[1].append(moenergy)
                    constructed_mosyms[1].append(mosym)
                    line = file_handler.virtual_next()
            return {mosyms.__name__: constructed_mosyms}
        return None


    @staticmethod
    def parse(file_handler, program: str, ccdata) -> Optional[dict]:
        constructed_data = None
        if program in mosyms.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(mosyms, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
