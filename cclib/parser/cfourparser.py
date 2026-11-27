# Copyright (c) 2025-2026, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for CFOUR output files"""

import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import timedelta
from typing import TYPE_CHECKING

from cclib.parser import logfileparser, utils

import numpy as np


if TYPE_CHECKING:
    from cclib.parser.logfilewrapper import FileWrapper


class CFOUR(logfileparser.Logfile):
    """A CFOUR v2.1 log file."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(logname="CFOUR", *args, **kwargs)

    def normalisesym(self, label: str) -> str:
        # CFOUR uses A'' instead of A"
        label = label.replace("''", '"')
        # CFOUR uses SG+ and SGg+ sometimes
        label = label.replace("+", "")
        # CFOUR uses 1g, 1u, 2g, 2u,... for E1g, E1u, E2g, E2u,...
        try:
            return f"E{int(label[0])}{label[1]}"
        except ValueError:
            if len(label) >= 2:
                # CFOUR uses SG for sigma, PI for pi, DE for delta, and PH for phi
                if "u" == label:
                    label = "Eu"
                if "g" == label:
                    label = "Eg"
                if "SG" == label[:2]:
                    if len(label) == 2:
                        label = "sigma"
                    else:
                        label = f"sigma.{label[2]}"
                if "PI" == label[:2]:
                    if len(label) == 2:
                        label = "pi"
                    else:
                        label = f"pi.{label[2]}"
                if "DE" == label[:2]:
                    if len(label) == 2:
                        label = "delta"
                    else:
                        label = f"delta.{label[2]}"
                if "PH" == label[:2]:
                    if len(label) == 2:
                        label = "phi"
                    else:
                        label = f"phi.{label[2]}"
            return label

    # get coefficients and exponents of the gaussian basis set from blocks such as the following
    """   ATOM                 EXPONENT      COEFFICIENTS

 C #1  1    S
+                 1   33980.000000   0.0000910  -0.0000190   0.0000000   0.0000000   0.0000000   0.0000000
                  2    5089.000000   0.0007040  -0.0001510   0.0000000   0.0000000   0.0000000   0.0000000
                  3    1157.000000   0.0036930  -0.0007850   0.0000000   0.0000000   0.0000000   0.0000000
                  4     326.600000   0.0153600  -0.0033240   0.0000000   0.0000000   0.0000000   0.0000000
                  5     106.100000   0.0529290  -0.0115120   0.0000000   0.0000000   0.0000000   0.0000000
                  6      38.110000   0.1470430  -0.0341600   0.0000000   0.0000000   0.0000000   0.0000000
                  7      14.750000   0.3056310  -0.0771730   0.0000000   0.0000000   0.0000000   0.0000000
                  8       6.035000   0.3993450  -0.1414930   0.0000000   0.0000000   0.0000000   0.0000000
                  9       2.530000   0.2170510  -0.1180190   0.0000000   0.0000000   0.0000000   0.0000000
                 10       0.735500   0.0158940   0.2738060   1.0000000   0.0000000   0.0000000   0.0000000
                 11       0.290500  -0.0030840   0.5865100   0.0000000   1.0000000   0.0000000   0.0000000
                 12       0.111100   0.0009780   0.2854300   0.0000000   0.0000000   1.0000000   0.0000000
                 13       0.041450   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   1.0000000

 C #1  1    X
+                14      34.510000   0.0053780   0.0000000   0.0000000   0.0000000   0.0000000
                 15       7.915000   0.0361320   0.0000000   0.0000000   0.0000000   0.0000000
                 16       2.368000   0.1424930   0.0000000   0.0000000   0.0000000   0.0000000
                 17       0.813200   0.3421500   1.0000000   0.0000000   0.0000000   0.0000000
                 18       0.289000   0.4638640   0.0000000   1.0000000   0.0000000   0.0000000
                 19       0.100700   0.2500280   0.0000000   0.0000000   1.0000000   0.0000000
                 20       0.032180   0.0000000   0.0000000   0.0000000   0.0000000   1.0000000

 C #1  1    Y
+                21      34.510000   0.0053780   0.0000000   0.0000000   0.0000000   0.0000000
                 22       7.915000   0.0361320   0.0000000   0.0000000   0.0000000   0.0000000
                 23       2.368000   0.1424930   0.0000000   0.0000000   0.0000000   0.0000000
                 24       0.813200   0.3421500   1.0000000   0.0000000   0.0000000   0.0000000
                 25       0.289000   0.4638640   0.0000000   1.0000000   0.0000000   0.0000000
                 26       0.100700   0.2500280   0.0000000   0.0000000   1.0000000   0.0000000
                 27       0.032180   0.0000000   0.0000000   0.0000000   0.0000000   1.0000000

 C #1  1    Z
+                28      34.510000   0.0053780   0.0000000   0.0000000   0.0000000   0.0000000
                 29       7.915000   0.0361320   0.0000000   0.0000000   0.0000000   0.0000000
                 30       2.368000   0.1424930   0.0000000   0.0000000   0.0000000   0.0000000
                 31       0.813200   0.3421500   1.0000000   0.0000000   0.0000000   0.0000000
                 32       0.289000   0.4638640   0.0000000   1.0000000   0.0000000   0.0000000
                 33       0.100700   0.2500280   0.0000000   0.0000000   1.0000000   0.0000000
                 34       0.032180   0.0000000   0.0000000   0.0000000   0.0000000   1.0000000

 C #1  1    XX
+                35       1.848000   1.0000000   0.0000000   0.0000000   0.0000000
                 36       0.649000   0.0000000   1.0000000   0.0000000   0.0000000
                 37       0.228000   0.0000000   0.0000000   1.0000000   0.0000000
                 38       0.076600   0.0000000   0.0000000   0.0000000   1.0000000

 C #1  1    XY
+                39       1.848000   1.0000000   0.0000000   0.0000000   0.0000000
                 40       0.649000   0.0000000   1.0000000   0.0000000   0.0000000
                 41       0.228000   0.0000000   0.0000000   1.0000000   0.0000000
                 42       0.076600   0.0000000   0.0000000   0.0000000   1.0000000

 C #1  1    XZ
+                43       1.848000   1.0000000   0.0000000   0.0000000   0.0000000
                 44       0.649000   0.0000000   1.0000000   0.0000000   0.0000000
                 45       0.228000   0.0000000   0.0000000   1.0000000   0.0000000
                 46       0.076600   0.0000000   0.0000000   0.0000000   1.0000000

 C #1  1    YY
+                47       1.848000   1.0000000   0.0000000   0.0000000   0.0000000
                 48       0.649000   0.0000000   1.0000000   0.0000000   0.0000000
                 49       0.228000   0.0000000   0.0000000   1.0000000   0.0000000
                 50       0.076600   0.0000000   0.0000000   0.0000000   1.0000000

 C #1  1    YZ
+                51       1.848000   1.0000000   0.0000000   0.0000000   0.0000000
                 52       0.649000   0.0000000   1.0000000   0.0000000   0.0000000
                 53       0.228000   0.0000000   0.0000000   1.0000000   0.0000000
                 54       0.076600   0.0000000   0.0000000   0.0000000   1.0000000

 C #1  1    ZZ
+                55       1.848000   1.0000000   0.0000000   0.0000000   0.0000000
                 56       0.649000   0.0000000   1.0000000   0.0000000   0.0000000
                 57       0.228000   0.0000000   0.0000000   1.0000000   0.0000000
                 58       0.076600   0.0000000   0.0000000   0.0000000   1.0000000

 C #1  1    F300
+                59       1.419000   1.0000000   0.0000000   0.0000000
                 60       0.485000   0.0000000   1.0000000   0.0000000
                 61       0.187000   0.0000000   0.0000000   1.0000000
"""

    def parse_basis(self, inputfile: "FileWrapper", line: str):
        atom_index = {}
        line = next(inputfile)
        should_parse = True
        last_line = ""
        temp_basis_info = []
        gbasis = []
        self.ecp_labels = defaultdict(list)
        for i in range(len(self.atomic_symbols)):
            gbasis.append([])
        hashtag_in_last_line = False
        first_iter = True
        while True:
            if not line.strip():
                line = next(inputfile)
                continue
            line = line.replace("+", "")
            tokens = line.split()
            line_length = len(tokens)
            if not (
                (hashtag_in_last_line) or ("#" in line) or (len(last_line.split()) == line_length)
            ):
                break
            if "#" in line:
                if (not first_iter) and should_parse:
                    for i in temp_basis_info:
                        gbasis[atom_index[curr_atom]].append((curr_ang_mom, i))  # noqa: F821
                first_iter = False
                hashtag_in_last_line = True
                if len(tokens[0]) == 1:
                    curr_atom = tokens[0] + tokens[1] + tokens[2]
                    ecp_label = tokens[0] + tokens[1]
                else:
                    curr_atom = tokens[0] + tokens[1]
                    ecp_label = tokens[0]
                symbol_len = 0
                for i in curr_atom:
                    if not i == "#":
                        symbol_len += 1
                    else:
                        break
                if curr_atom not in atom_index:
                    for i in range(len(self.atomic_symbols)):
                        if self.atomic_symbols[i] == curr_atom[:symbol_len] and (
                            i not in atom_index.values()
                        ):
                            self.ecp_labels[ecp_label] += [i]
                            atom_index[curr_atom] = i
                            break
                if len(tokens[-1]) == 1:
                    if tokens[-1] == "S":
                        curr_ang_mom = "S"
                        should_parse = True
                    elif tokens[-1] == "X":
                        curr_ang_mom = "P"
                        should_parse = True
                    else:
                        should_parse = False
                elif len(tokens[-1]) == 2:
                    if tokens[-1] == "XX":
                        curr_ang_mom = "D"
                        should_parse = True
                    else:
                        should_parse = False
                else:
                    if (
                        (int(tokens[-1][-3]) > 2)
                        and (int(tokens[-1][-2]) == 0)
                        and (int(tokens[-1][-1]) == 0)
                    ):
                        curr_ang_mom = tokens[-1][:-3]
                        should_parse = True
                    else:
                        should_parse = False
            if should_parse and ("#" not in line):
                if hashtag_in_last_line:
                    temp_basis_info = []
                    for i in range(line_length - 2):
                        temp_basis_info.append([])
                for i in range(line_length - 2):
                    if not float(tokens[2 + i]) == 0.0:
                        temp_basis_info[i].append((float(tokens[1]), float(tokens[2 + i])))
            if hashtag_in_last_line and ("#" not in line):
                hashtag_in_last_line = False
            last_line = line
            line = next(inputfile)
        if should_parse:
            for i in temp_basis_info:
                gbasis[atom_index[curr_atom]].append((curr_ang_mom, i))
        return gbasis

    def before_parsing(self) -> None:
        # geting atomic number and symbol is different for 1 atom
        self.set_attribute("only_one_atom", False)

        # set to true so that atomic numbers are parsed on the first block of coordinates
        self.set_attribute("first_coord_block", True)

        # set to True so that alpha MOs are parsed first
        self.set_attribute("alpha_mos_to_parse", True)

        # set temp alpha mo coeffs to []
        self.set_attribute("temp_alpha_mocoeffs", [])

        # set temp beta mo coeffs to []
        self.set_attribute("temp_beta_mocoeffs", [])

        # should mo coeffs be reset the next time an mo is parsed
        self.set_attribute("mocoeffs_should_be_reset", True)

        # set to True so that the first time parsing mo coeffs it also parses ao names
        self.set_attribute("parse_aonames", True)

        # dict of ECP labels
        self.set_attribute("ecp_labels", {})

        # core electron dict
        self.set_attribute("core_electron_dict", {})

        # set sym numbering to {}
        self.set_attribute("sym_numbering", {})

        # set current symmetry to "0"
        self.set_attribute("curr_sym", "0")

        # set to True so to indicate no time was recorded
        self.set_attribute("no_time", True)

    def after_parsing(self) -> None:
        # set metadata "success" to False if no time was recorded
        if self.no_time:
            self.metadata["success"] = False

        # get optdone
        if hasattr(self, "geovalues") and hasattr(self, "geotargets"):
            for i in range(len(self.geovalues)):
                if self.geovalues[i][0] < self.geotargets[0]:
                    self.append_attribute("optdone", i)

        # get the number of atoms
        if hasattr(self, "atomcoords"):
            if len(self.atomcoords) >= 1:
                self.set_attribute("natom", len(self.atomnos))

        # get the number of atomic orbitals in the basis
        if hasattr(self, "aonames"):
            self.set_attribute("nbasis", len(self.aonames))

        # get the number of molecular orbitals
        if hasattr(self, "moenergies"):
            if len(self.moenergies) >= 1:
                self.set_attribute("nmo", len(self.moenergies[0]))

        # get core electrons
        for i in self.ecp_labels:
            if i in self.core_electron_dict:
                for j in self.ecp_labels[i]:
                    self.coreelectrons[j] = self.core_electron_dict[i]

        # sort etenergies, etoscs, etsecs, and etsyms
        if hasattr(self, "etenergies"):
            sort_inds = np.argsort(self.etenergies)
            temp_etenergies = []
            temp_etoscs = []
            temp_etsecs = []
            temp_etsyms = []
            for i in sort_inds:
                temp_etenergies.append(self.etenergies[i])
                if hasattr(self, "etoscs"):
                    if self.estate_prop_on:
                        temp_etoscs.append(self.etoscs[i])
                temp_etsecs.append(self.etsecs[i])
                temp_etsyms.append(self.etsyms[i])
            self.etenergies = temp_etenergies
            if hasattr(self, "etoscs"):
                self.etoscs = temp_etoscs
            self.etsecs = temp_etsecs
            self.etsyms = temp_etsyms

    def extract(self, inputfile: "FileWrapper", line: str) -> None:
        tokens = line.split()

        # get the version of CFOUR
        if "Version" in line:
            self.metadata["package_version"] = tokens[1]

        # parse the control parameters
        if line.strip() == "CFOUR Control Parameters":
            self.skip_lines(inputfile, ["d", "external internal value units", "name name", "d"])
            control_parameters = parse_control_parameters(inputfile)
            self.metadata["basis_set"] = control_parameters["BASIS"].value
            excite = str(control_parameters["EXCITE"].value)
            if excite != "NONE":
                self.metadata["excited_states_method"] = (
                    excite + "-" + str(control_parameters["CALCLEVEL"].value)
                )
            reference = str(control_parameters["REFERENCE"].value)
            if reference[0] == "U":
                self.metadata["unrestricted"] = True
                self.set_attribute("homos", np.array([0, 0]))
            elif reference == "ROHF":
                self.set_attribute("homos", np.array([0, 0]))
                self.metadata["unrestricted"] = False
            else:
                self.set_attribute("homos", np.array([0]))
                self.metadata["unrestricted"] = False
            self.set_attribute("scf_target_value", float(control_parameters["SCF_CONV"].value))
            self.set_attribute("geotargets", [10 ** -float(control_parameters["GEO_CONV"].value)])
            self.set_attribute("charge", int(control_parameters["CHARGE"].value))
            self.set_attribute(
                "estate_prop_on", str(control_parameters["ESTATE_PROP"].value) != "OFF"
            )
            self.set_attribute("mult", int(control_parameters["MULTIPLICTY"].value))

        if ("Minimum force:" in line) and ("RMS force:" in line):
            self.append_attribute("geovalues", [float(tokens[6])])

        # get full point group
        if "The full molecular point group is" in line:
            self.metadata["symmetry_detected"] = tokens[6].lower()

        # get used point group
        if "The computational point group is" in line:
            self.metadata["symmetry_used"] = tokens[5].lower()

        # get success, cpu time, and wall time
        if "@CHECKOUT-I, Total execution time (CPU/WALL):" in line:
            self.no_time = False
            self.metadata["success"] = True
            self.metadata["cpu_time"] = [timedelta(seconds=float(tokens[5][:-1]))]
            self.metadata["wall_time"] = [timedelta(seconds=float(tokens[6]))]

        # get coupled cluster energy
        if "A miracle has come to pass. The CC iterations have converged." in line:
            cc_lines = []
            while "@CHECKOUT-I," not in line:
                if not line.strip() == "":
                    cc_lines.append(line)
                line = next(inputfile)
            cc_tokens = cc_lines[-1].split()
            if cc_tokens[-1] == "a.u.":
                ccenergy_index = -2
            else:
                ccenergy_index = -1
            self.append_attribute("ccenergies", cc_tokens[ccenergy_index])

        # get coefficients and exponents of the gaussian basis set from blocks such as the following
        if "ATOM                 EXPONENT      COEFFICIENTS" in line:
            self.set_attribute("gbasis", self.parse_basis(inputfile, line))

        # exception for only one atom
        if "1 entries found in Z-matrix" in line:
            self.only_one_atom = True

        # The next 2 sections parse blocks like the following if there is only one atom in the calculation
        """ NUCLEAR CHARGE:                        6
    NUMBER OF SYMMETRY INDEPENDENT ATOMS:  1
    HIGHEST ORBITAL TYPE:                  G

      1 GROUPS OF CGTOS OF S TYPE
      1 GROUPS OF CGTOS OF P TYPE
      1 GROUPS OF CGTOS OF D TYPE
      1 GROUPS OF CGTOS OF F TYPE
      1 GROUPS OF CGTOS OF G TYPE

    NUCLEAR COORDINATES (IN A.U.) ARE :

    C #1        0.000000000000000        0.000000000000000        0.000000000000000

   INTERNUCLEAR DISTANCES (A) :

     FOR ATOM C #1 (COORDINATES :   0.00000   0.00000   0.00000)
"""
        if ("NUCLEAR CHARGE:" in line) and self.only_one_atom:
            self.set_attribute("atomnos", [int(tokens[2])])
        if ("NUCLEAR COORDINATES (IN A.U.) ARE :" in line) and self.only_one_atom:
            line = next(inputfile)
            line = next(inputfile)
            tokens = line.split()
            self.set_attribute("atomic_symbols", [tokens[0]])
            self.set_attribute(
                "atomcoords", [[[float(tokens[2]), float(tokens[3]), float(tokens[4])]]]
            )

        # get the coordinates at each step in a geometry optimization
        # if this is the first time parsing a block of coordinates also get the atomic numbers
        """
         ----------------------------------------------------------------
                     Coordinates used in calculation (QCOMP)
         ----------------------------------------------------------------
          Z-matrix   Atomic            Coordinates (in bohr)
           Symbol    Number           X              Y              Z
         ----------------------------------------------------------------
             C         6        -0.00000000     0.00000000     1.23064327
             O         8         0.00000000     0.00000000    -0.92327590
         ----------------------------------------------------------------
        """
        if "Coordinates used in calculation (QCOMP)" in line:
            self.skip_lines(inputfile, ["d", "header", "header", "d"])
            line = next(inputfile)
            tokens = line.split()
            atomnos = []
            atomic_symbols = []
            temp_atomcoords = []
            while "----------------------------------------------------------------" not in line:
                if self.first_coord_block:
                    atomic_number = int(tokens[1])
                    atomic_symbol = tokens[0]
                    if atomic_number > 0:
                        atomnos.append(atomic_number)
                        atomic_symbols.append(atomic_symbol)
                temp_atomcoords.append(
                    [
                        utils.convertor(float(tokens[2]), "bohr", "Angstrom"),
                        utils.convertor(float(tokens[3]), "bohr", "Angstrom"),
                        utils.convertor(float(tokens[4]), "bohr", "Angstrom"),
                    ]
                )
                line = next(inputfile)
                tokens = line.split()
            if self.first_coord_block:
                self.set_attribute("atomnos", atomnos)
                self.set_attribute("coreelectrons", np.zeros(len(atomnos)))
                self.set_attribute("atomic_symbols", atomic_symbols)
                self.first_coord_block = False
            self.append_attribute("atomcoords", temp_atomcoords)

        # get core electrons in each atoms ECP
        if "ECP PARAMETERS FOR ATOM" in line:
            num_ce_index = tokens[-1]
            line = next(inputfile)
            tokens = line.split()
            let_ce_index = line.split(":")[0]
            ce_index = let_ce_index + "#" + num_ce_index
            while "NCORE =" not in line:
                line = next(inputfile)
                tokens = line.split()
            self.core_electron_dict[ce_index] = int(tokens[2])

        # get scfenergies, scftargets, and scfvalues at each step in a geometry optimization
        if "Iteration         Total Energy            Largest Density Difference" in line:
            no_scf_energy_yet = True
            while no_scf_energy_yet:
                last_line = line
                line = next(inputfile)
                tokens = line.split()
                if "current occupation vector" in line:
                    last_tokens = last_line.split()
                    if last_tokens[0] == "0":
                        self.append_attribute("scfvalues", [])
                    self.scfvalues[-1].append(
                        [
                            float(last_tokens[2].split("D")[0])
                            * float(np.power(10.0, int(last_tokens[2].split("D")[1])))
                        ]
                    )
                if "E(SCF)=" in line:
                    self.append_attribute("scftargets", [self.scf_target_value])
                    self.scfvalues[-1].append(
                        [
                            float(tokens[2].split("D")[0])
                            * float(np.power(10.0, int(tokens[2].split("D")[1])))
                        ]
                    )
                    self.append_attribute("scfenergies", float(tokens[1]))
                    no_scf_energy_yet = False

        # get alpha mo energies of the last ran scf method
        if "ORBITAL EIGENVALUES (ALPHA)  (1H = 27.2113834 eV)" in line:
            self.skip_lines(inputfile, ["b", "header", "d"])
            line = next(inputfile)
            tokens = line.split()
            self.set_attribute("moenergies", [])
            self.set_attribute("mosyms", [])
            alpha_moenergies = []
            alpha_mosyms = []
            while not (
                ("VSCF finished." in line)
                or ("ORBITAL EIGENVALUES ( BETA)  (1H = 27.2113834 eV)" in line)
                or ("SCF failed to converge in" in line)
            ):
                if ("+++++" in line) or (line.strip() == ""):
                    if "+++++" in line:
                        self.homos[0] = int(last_line.split()[0]) - 1
                    line = next(inputfile)
                    tokens = line.split()
                    continue
                alpha_moenergies.append(utils.float(tokens[2]))
                alpha_mosyms.append(self.normalisesym(tokens[5]))
                if self.normalisesym(tokens[5]) not in self.sym_numbering.values():
                    self.sym_numbering[tokens[6][1:-1]] = self.normalisesym(tokens[5])
                last_line = line
                line = next(inputfile)
                tokens = line.split()
            self.append_attribute("moenergies", np.array(alpha_moenergies))
            self.append_attribute("mosyms", alpha_mosyms)

        # get beta mo energies of the last ran scf method if an unrestricted reference is used
        if "ORBITAL EIGENVALUES ( BETA)  (1H = 27.2113834 eV)" in line:
            self.skip_lines(inputfile, ["b", "header", "d"])
            line = next(inputfile)
            tokens = line.split()
            beta_moenergies = []
            beta_mosyms = []
            while not (("VSCF finished." in line) or ("SCF failed to converge in" in line)):
                if ("+++++" in line) or (line.strip() == ""):
                    if "+++++" in line:
                        self.homos[1] = int(last_line.split()[0]) - 1
                    line = next(inputfile)
                    tokens = line.split()
                    continue
                beta_moenergies.append(utils.float(tokens[2]))
                beta_mosyms.append(self.normalisesym(tokens[5]))
                last_line = line
                line = next(inputfile)
                tokens = line.split()
            self.append_attribute("moenergies", np.array(beta_moenergies))
            self.append_attribute("mosyms", beta_mosyms)

        # add on to the molecular orbital coefficients
        if len(tokens) >= 2:
            if "Symmetry" == tokens[0]:
                line = next(inputfile)
                tokens = line.split()
                if (
                    "--------------------------------------------------------------------------------"
                    in line
                ):
                    if self.mocoeffs_should_be_reset:
                        self.set_attribute("mocoeffs", [])
                        self.mocoeffs_should_be_reset = False
                    line = next(inputfile)
                    tokens = line.split()
                    len_split_line = len(tokens)
                    if self.alpha_mos_to_parse:
                        if self.parse_aonames:
                            aonames = []
                            atombasis = []
                            start_atom = 0
                            sub_amount = 0
                            basis_index = 0
                        for i in range(len_split_line - 2):
                            self.temp_alpha_mocoeffs.append([])
                        while "----------" not in line:
                            if self.parse_aonames:
                                if int(tokens[0]) - start_atom > 1:
                                    sub_amount += (int(tokens[0]) - start_atom) - 1
                                    start_atom = int(tokens[0])
                                    atombasis.append([])
                                    subshell_number = 0
                                elif int(tokens[0]) - start_atom == 1:
                                    start_atom = int(tokens[0])
                                    atombasis.append([])
                                    subshell_number = 0
                                if tokens[1] == "S":
                                    subshell_number += 1
                                if len(tokens[1]) == 2:
                                    if tokens[1][1] == "X":
                                        subshell_number += 1
                                if len(tokens[1]) >= 3:
                                    if tokens[1][1] == "X":
                                        try:
                                            x_power = int(tokens[1][2:])  # noqa: F841
                                            subshell_number += 1
                                        except ValueError:
                                            pass
                                aonames.append(
                                    f"{self.atomic_symbols[(int(tokens[0]) - sub_amount) - 1]}{int(tokens[0]) - sub_amount}_{subshell_number}{tokens[1]}"
                                )
                                atombasis[-1].append(basis_index)
                                basis_index += 1
                            for i in range(len_split_line - 2):
                                self.temp_alpha_mocoeffs[i - (len_split_line - 2)].append(
                                    utils.float(tokens[i + 2])
                                )
                            line = next(inputfile)
                            tokens = line.split()
                        if self.parse_aonames:
                            self.set_attribute("aonames", aonames)
                            self.set_attribute("atombasis", atombasis)
                            self.parse_aonames = False
                        line = next(inputfile)
                        tokens = line.split()
                        if ("-----------" in line) or ("" == line.strip()):
                            self.alpha_mos_to_parse = False
                            self.mocoeffs.append(self.temp_alpha_mocoeffs)
                            self.temp_alpha_mocoeffs = []
                            if "" == line.strip():
                                self.alpha_mos_to_parse = True
                                self.mocoeffs_should_be_reset = True
                    else:
                        for i in range(len_split_line - 2):
                            self.temp_beta_mocoeffs.append([])
                        while "----------------" not in line:
                            for i in range(len_split_line - 2):
                                self.temp_beta_mocoeffs[i - (len_split_line - 2)].append(
                                    utils.float(tokens[i + 2])
                                )
                            line = next(inputfile)
                            tokens = line.split()
                        line = next(inputfile)
                        tokens = line.split()
                        if "" == line.strip():
                            self.alpha_mos_to_parse = True
                            self.append_attribute("mocoeffs", self.temp_beta_mocoeffs)
                            self.temp_beta_mocoeffs = []
                            self.mocoeffs_should_be_reset = True

        # change curr_sym
        if "Beginning symmetry block" in line:
            self.curr_sym = tokens[3][:-1]

        # get excitation energies
        if "Converged eigenvalue:" in line:
            temp_etenergy = float(tokens[2])
            mult = "Singlet"
            if self.estate_prop_on:
                num_dash_lines = 0
                parse_etsecs = False
                keep_parse = True
                temp_etsecs = []
                line = next(inputfile)
                tokens = line.split()
                while "Right Transition Moment" not in line:
                    if "Converged eigenvalue:" in line:
                        temp_etenergy = float(tokens[2])
                        parse_etsecs = True
                    if parse_etsecs and (
                        "--------------------------------------------------------------------------------"
                        in line
                    ):
                        keep_parse = False
                    if parse_etsecs and keep_parse and num_dash_lines == 2:
                        i = int(tokens[0])
                        if i:
                            i -= 1
                        j = int(tokens[1])
                        if j:
                            j -= 1
                        a = int(tokens[2])
                        if a:
                            a -= 1
                        b = int(tokens[3])
                        if b:
                            b -= 1
                        temp_etsecs.append(((i, j), (a, b), float(tokens[4])))
                    if (
                        parse_etsecs
                        and keep_parse
                        and num_dash_lines == 3
                        and ("The state is a" in line)
                    ):
                        mult = tokens[4]
                        if mult == "triplet":
                            mult = "Triplet"
                        else:
                            mult = "Singlet"
                    if parse_etsecs and (
                        "--------------------------------------------------------------------------------"
                        in line
                    ):
                        num_dash_lines += 1
                        keep_parse = True
                    line = next(inputfile)
                    tokens = line.split()
                self.append_attribute("etenergies", temp_etenergy)
                self.append_attribute("etsecs", temp_etsecs)
                self.append_attribute("etsyms", f"{mult}-{self.sym_numbering[self.curr_sym]}")
            else:
                self.append_attribute("etenergies", temp_etenergy)
                self.skip_lines(inputfile, ["header", "header", "header", "d", "header", "d"])
                line = next(inputfile)
                tokens = line.split()
                temp_etsecs = []
                while (
                    "--------------------------------------------------------------------------------"
                    not in line
                ):
                    i = int(tokens[0])
                    if i:
                        i -= 1
                    j = int(tokens[1])
                    if j:
                        j -= 1
                    a = int(tokens[2])
                    if a:
                        a -= 1
                    b = int(tokens[3])
                    if b:
                        b -= 1
                    temp_etsecs.append(((i, j), (a, b), float(tokens[4])))
                    line = next(inputfile)
                    tokens = line.split()
                self.append_attribute("etsecs", temp_etsecs)
                self.append_attribute("etsyms", f"{mult}-{self.sym_numbering[self.curr_sym]}")

        # get etoscs
        if "Norm of oscillator strength :" in line:
            self.append_attribute("etoscs", float(tokens[-1]))

        # zero point vibrational energy correction
        if "Zero-point energy:" in line:
            self.set_attribute("zpve", float(tokens[5]) / 2625.5)

        # get vibrational data
        if "Normal Coordinate Analysis" in line:
            self.skip_lines(inputfile, ["d", "header", "header", "d", "header", "d"])
            line = next(inputfile)
            tokens = line.split()
            while "----------------------------------------------------------------" not in line:
                if "VIBRATION" in line:
                    self.append_attribute("vibfreqs", utils.parse_vibfreq(tokens[1]))
                    self.append_attribute("vibirs", float(tokens[2]))
                    self.append_attribute("vibsyms", self.normalisesym(tokens[0]))
                line = next(inputfile)
                tokens = line.split()

        # get vibrational displacements
        if "Normal Coordinates" in line:
            while True:
                if "VIBRATION" in line:
                    num_vibs_this_line = len(tokens)
                    for i in range(num_vibs_this_line):
                        self.append_attribute("vibdisps", [])
                    line = next(inputfile)
                    tokens = line.split()
                    while ("Gradient vector in normal coordinate representation" not in line) and (
                        not (line.strip() == "")
                    ):
                        test_tokens = tokens[1].replace("-", " ").split()
                        if len(test_tokens) == 2:
                            token1 = test_tokens[0]
                            token2 = f"-{test_tokens[1]}"
                            if tokens[1][0] == "-":
                                token1 = f"-{token1}"
                            tokens.remove(tokens[1])
                            tokens.insert(1, token2)
                            tokens.insert(1, token1)
                        for i in range(num_vibs_this_line):
                            self.vibdisps[-(i + 1)].append(
                                [
                                    float(tokens[-((3 * i) + 3)]),
                                    float(tokens[-((3 * i) + 2)]),
                                    float(tokens[-((3 * i) + 1)]),
                                ]
                            )
                        line = next(inputfile)
                        tokens = line.split()
                    if "Gradient vector in normal coordinate representation" in line:
                        break
                line = next(inputfile)
                tokens = line.split()

        # get gradients
        if "gradient from JOBARC" in line:
            line = next(inputfile)
            tokens = line.split()
            self.append_attribute("grads", [])
            while line.strip():
                self.grads[-1].append([float(tokens[0]), float(tokens[1]), float(tokens[2])])
                line = next(inputfile)
                if "Norm is" in line:
                    break
                tokens = line.split()
        if "Total MP2 energy" in line:
            self.append_attribute("mpenergies", [])
            self.mpenergies[-1].append(float(tokens[4]))
        if "D-MBPT(3)" in line:
            self.mpenergies[-1].append(float(tokens[2]))
        if "Total MBPT(4)       energy:" in line:
            self.mpenergies[-1].append(float(tokens[3]))


ControlParameterValue = str | int | float | None


@dataclass(frozen=True, kw_only=True, slots=True)
class ControlParameter:
    """A row in the "CFOUR Control Parameters" table."""

    external_name: str
    """The row's value under the "External Name" column."""

    internal_name: str
    """The row's value under the "Internal Name" column."""

    value: ControlParameterValue
    """The row's value under the "Value" column, parsed into a Python type."""

    raw_value: str
    """The row's value under the "Value" column, as-is."""

    integer: int | None

    units: str | None
    """The row's value under the "Units" column, if it exists."""


LINE_RE = re.compile(
    r"""
    ^\s*
    (?P<external>\S+) \s+
    (?P<internal>\S+) \s+
    (?P<value>
        \d+D-\s*\d+                 # 10D-  4, 10D-100 (has internal spaces)
      | [^\s\[]+(?:\ [^\s\[]+)*     # tokens joined by single spaces; stops before '['
    )
    \s*
    (?:\[\s*(?P<int>\d+)\])?        # optional [  0], with or without a space before it
    (?:\s+(?P<units>\S+(?:\ \S+)*))?  # optional units: ***, x 10-6, # of cyc., ...
    \s*$
    """,
    re.VERBOSE,
)


def _convert_control_parameter_value(raw: str) -> ControlParameterValue:
    if set(raw) == {"*"}:  # *** (unset) or ********* (Fortran overflow)
        return None
    if m := re.fullmatch(r"10D-\s*(\d+)", raw):  # CFOUR's way of printing 10**-n
        return 10.0 ** -int(m[1])
    try:
        return int(raw)
    except ValueError:
        return raw


def parse_control_parameter_line(control_parameter_line: str) -> ControlParameter | None:
    m = LINE_RE.match(control_parameter_line)
    if not m:
        return None
    units = m["units"]
    control_parameter = ControlParameter(
        external_name=m["external"],
        internal_name=m["internal"],
        value=_convert_control_parameter_value(m["value"]),
        raw_value=m["value"],
        integer=int(m["int"]) if m["int"] is not None else None,
        units=None if units is None or set(units) == {"*"} else units,
    )
    return control_parameter


def parse_control_parameters(inputfile: "FileWrapper") -> dict[str, ControlParameter]:
    control_parameters = dict()
    line = next(inputfile)
    while set(line.strip()) != {"-"}:
        control_parameter = parse_control_parameter_line(line)
        if control_parameter is not None:
            control_parameters[control_parameter.external_name] = control_parameter
        line = next(inputfile)
    return control_parameters
