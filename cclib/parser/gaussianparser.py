# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Gaussian output files"""

import datetime
import re
from typing import Dict, List, Tuple

from cclib.parser import data, logfileparser, utils
from cclib.parser.logfileparser import StopParsing

import numpy

__all__ = ("Gaussian", "parse_version")


class Gaussian(logfileparser.Logfile):
    """A Gaussian log file."""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="Gaussian", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return f"Gaussian log file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'Gaussian("{self.filename}")'

    def normalisesym(self, label):
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

    def before_parsing(self):
        # Calculations use point group symmetry by default.
        self.uses_symmetry = True

        # Extract only well-formed numbers in scientific notation.
        self.re_scinot = re.compile(r"(\w*)=\s*(-?\d\.\d{2}D[+-]\d{2})")
        # Extract only well-formed numbers in traditional
        # floating-point format.
        self.re_float = re.compile(r"(\w*-?\w*)=\s*(-?\d+\.\d{10,})")

        # Detect when a new link (Gaussian subprogram) is entered.
        self.re_link = re.compile(r"\(Enter [\w/:\\\-\.]+l(?P<link_number>\d{3,4})\.exe\)")
        # Printed in a different format, but always start here.
        self.link = 1

        # Flag for identifying Coupled Cluster runs.
        self.coupledcluster = False

        # Fragment number for counterpoise or fragment guess calculations
        # (normally zero).
        self.counterpoise = 0

        # Flag for identifying ONIOM calculations.
        self.oniom = False

        # Flag for identifying BOMD calculations.
        # These calculations have a back-integration algorithm so that not all
        # geometries should be kept.
        # We also add a "time" attribute to the parser.
        self.BOMD = False

        # Do we have high-precision polarizabilities printed from a
        # dedicated `polar` job? If so, avoid duplicate parsing.
        self.hp_polarizabilities = False

        # Used to estimate wall times from CPU times.
        self.num_cpu = 1

        # For detecting when excited states reset (because of a new level of theory,
        # or a new round in an optimisation etc).
        self.last_et = 0

        # When running an NBO analysis on a calculation with different alpha
        # and beta orbitals, the charges are printed for combined, alpha, and
        # beta spins.  This is false when in the combined section (the general
        # case) and true when in an alpha or beta section.
        self.nbo_spin_section = False

        # Define strings needed for line detection. Older Gaussian versions
        # don't always give the charge type explicitly, so we must include
        # "atomic" as a general term to catch all other atomic charge or spin
        # lines.
        self.atomprops_no_atomic = ["mulliken", "lowdin", "apt", "hirshfeld"]
        self.atomprops = self.atomprops_no_atomic + ["atomic"]
        # The (Hirshfeld) partial atomic charges and spins parsed in these
        # sections need to be swapped after parsing via the general extraction
        # code.  This is only required for version 09; Hirshfeld charges were
        # not present earlier, and 16 combines them with CM5 in a different
        # section.
        self.atomcharges_atomspins_headers_swap = [
            " spin densities, charges and dipoles using iradan="
        ]
        self.atomcharges_atomspins_headers = [
            " atomic charges:",
            " charges:",
            " charges with hydrogens summed into heavy atoms:",
            " atomic charges with hydrogens summed into heavy atoms:",
            " atomic spin densities:",
            " charges and spin densities with hydrogens summed into heavy atoms:",
            " charges and spin densities:",
        ] + self.atomcharges_atomspins_headers_swap

    def after_parsing(self):
        # atomcoords are parsed as a list of lists but it should be an array.
        # Done automatically later in arrayify, but we need it now for the
        # rest of this method.
        if hasattr(self, "atomcoords"):
            self.atomcoords = numpy.array(self.atomcoords)

        # Correct the percent values in the etsecs in the case of
        # a restricted calculation. The following has the
        # effect of including each transition twice.
        if hasattr(self, "etsecs") and len(self.homos) == 1:
            new_etsecs = [
                [(x[0], x[1], x[2] * numpy.sqrt(2)) for x in etsec] for etsec in self.etsecs
            ]
            self.etsecs = new_etsecs

        if (
            hasattr(self, "enthalpy")
            and hasattr(self, "temperature")
            and hasattr(self, "freeenergy")
        ):
            self.set_attribute("entropy", (self.enthalpy - self.freeenergy) / self.temperature)

        # This bit is needed in order to trim coordinates that are printed a second time
        # at the end of geometry optimizations. Note that we need to do this for both atomcoords
        # and inputcoords. The reason is that normally a standard orientation is printed and that
        # is what we parse into atomcoords, but inputcoords stores the input (unmodified) coordinates
        # and that is copied over to atomcoords if no standard oritentation was printed, which happens
        # for example for jobs with no symmetry. This last step, however, is now generic for all parsers.
        # Perhaps then this part should also be generic code...
        # Regression that tests this: Gaussian03/cyclopropenyl.rhf.g03.cut.log
        if hasattr(self, "optstatus") and len(self.optstatus) > 0:
            last_point = len(self.optstatus) - 1
            if hasattr(self, "atomcoords"):
                self.atomcoords = self.atomcoords[: last_point + 1]
            if hasattr(self, "inputcoords"):
                self.inputcoords = self.inputcoords[: last_point + 1]
            if hasattr(self, "rotconsts"):
                self.rotconsts = self.rotconsts[: last_point + 1]

        # If we parsed high-precision vibrational displacements, overwrite
        # lower-precision displacements in self.vibdisps
        if hasattr(self, "vibdispshp"):
            self.vibdisps = self.vibdispshp
            del self.vibdispshp
        if hasattr(self, "time"):
            self.time = [self.time[i] for i in sorted(self.time.keys())]
        if hasattr(self, "energies_BOMD"):
            self.set_attribute(
                "scfenergies", [self.energies_BOMD[i] for i in sorted(self.energies_BOMD.keys())]
            )
        if hasattr(self, "atomcoords_BOMD"):
            self.atomcoords = [self.atomcoords_BOMD[i] for i in sorted(self.atomcoords_BOMD.keys())]

        # Gaussian prints 'forces' in input orientation unlike other values such as 'moments' or 'vibdisp'.
        # Therefore, we convert 'grads' to the values in standard orientation with rotation matrix.
        if hasattr(self, "grads") and hasattr(self, "inputcoords") and hasattr(self, "atomcoords"):
            grads_std = []
            for grad, inputcoord, atomcoord in zip(self.grads, self.inputcoords, self.atomcoords):
                rotation = utils.get_rotation(numpy.array(inputcoord), numpy.array(atomcoord))
                grads_std.append(rotation.apply(grad))
            self.set_attribute("grads", numpy.array(grads_std))

        if hasattr(self, "ccenergy"):
            self.append_attribute("ccenergies", self.ccenergy)
            del self.ccenergy

        # If we have cpu times but no wall-times, we can calculate the later based on the number of CPUs used.
        if "cpu_time" in self.metadata and "wall_time" not in self.metadata:
            self.metadata["wall_time"] = []
            for cpu_time in self.metadata["cpu_time"]:
                self.metadata["wall_time"].append(cpu_time / self.num_cpu)

        # Some calculation types don't use consistent MP levels throughout.
        # For example, EOM-CCSD + opt prints MP3 and MP4 energies for the first iteration only.
        # However, mpenergies needs to be homogeneous, so we need to do something about it.
        if hasattr(self, "mpenergies"):
            max_mp = max(len(mp_e) for mp_e in self.mpenergies)
            mp_energies = []
            for index, energy in enumerate(self.mpenergies):
                if len(energy) != max_mp:
                    self.logger.warning(
                        "MP energies of order %d are incomplete and will be ignored", index + 2
                    )

                else:
                    mp_energies.append(energy)

            self.set_attribute("mpenergies", mp_energies)

        # The Hirshfeld section from version 09 prints atomspins no matter
        # what; no spin density means all printed values are zero.  Rather
        # than determine whether or not they should be present by looking at
        # (un)restricted flags and HOMO indices, check that they are the only
        # key present in atomspins, since otherwise Mulliken will be there,
        # since it isn't possible to simultaneously turn off Mulliken and
        # Hirshfeld printing.
        if hasattr(self, "atomcharges") and hasattr(self, "atomspins"):
            if (
                set(self.atomspins.keys()) == {"hirshfeld"}
                and sum(self.atomspins["hirshfeld"]) == 0.0
            ):
                delattr(self, "atomspins")

        super().after_parsing()

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        link_match = self.re_link.search(line)
        if link_match is not None:
            self.link = int(link_match.groupdict()["link_number"])

        # Extract the version number: "Gaussian 09, Revision D.01"
        # becomes "09revisionD.01".
        if line.strip() == "Cite this work as:":
            tokens = next(inputfile).split()
            self.metadata["legacy_package_version"] = "".join(
                [tokens[1][:-1], "revision", tokens[-1][:-1]]
            )

        if "Gaussian, Inc.," in line:
            self.skip_lines(inputfile, ["b", "s"])
            parsed_version = parse_version(next(inputfile))
            self.metadata.update(parsed_version)
            run_date = next(inputfile).strip()  # noqa: F841
            line = self.skip_line(inputfile, "s")[0]
            while set(line.strip()) != {"-"}:
                # Parse performance info.
                if "Will use up to" in line and "processors via shared memory." in line:
                    self.metadata["num_cpu"] = int(line.split()[4])

                elif line[1:6].lower() == "%mem=":
                    # The maximum amount of memory requested.
                    # We need to do some unit juggling.
                    mem_str = line.strip().upper()
                    # No space is allowed between units, units can probably only be uppercase but convert anyway.
                    # Supported units are: 'KB, MB, GB, TB, KW, MW, GW or TW'.
                    # TODO: Good opportunity for case here in future.
                    units = mem_str[-2:]
                    raw_mem = int(mem_str[5:-2])
                    # We are converting to bytes.
                    if units == "TB":
                        memory = raw_mem * 1e12

                    elif units == "GB":
                        memory = raw_mem * 1e9

                    elif units == "MB":
                        memory = raw_mem * 1e6

                    elif units == "KB":
                        memory = raw_mem * 1e3

                    elif units == "TW":
                        memory = raw_mem * 8e12

                    elif units == "GW":
                        memory = raw_mem * 8e9

                    elif units == "MW":
                        memory = raw_mem * 8e6

                    elif units == "KW":
                        memory = raw_mem * 8e3

                    else:
                        # No explicit units, default are single words (8-bytes)
                        memory = mem_str[5:] * 8

                    self.metadata["memory_available"] = int(memory)
                line = next(inputfile)
            line = next(inputfile)
            self.parse_keywords_route_and_comment(inputfile, line)

        if line.strip().startswith("Link1:  Proceeding to internal job step number"):
            self.new_internal_job()
            self.skip_line(inputfile, "d")
            line = next(inputfile)
            self.parse_keywords_route_and_comment(inputfile, line)

        # Parse performance info.
        if "Will use up to" in line and "processors via shared memory." in line:
            self.metadata["num_cpu"] = int(line.split()[4])

        elif "Leave Link    1" in line and "MaxMem=" in line and "num_cpu" in self.metadata:
            # Leave Link    1 at Wed Apr  4 10:49:19 2018, MaxMem=   805306368 cpu:               0.3 elap:               0.0
            # Gaussian helpfully prints the total 'available' memory for us. There are, however a few caveats here:
            # 1) This memory (in bytes) is per CPU
            # 2) The total memory here (x num_cpu) will not equal the amount requested in %mem because Gaussian (probably erroneously)
            #    interprets GB as gibibytes (1024 x 1024 x 1024 bytes) rather than gigabytes. This has the unfortunate consequence of
            #    Gaussian assigning more memory than you probably expected.
            #
            # MaxMem can appear with or without whitespace:
            # MaxMem=  805306368
            # MaxMem=174483046400
            memory_per_cpu = int(re.search(r"MaxMem=\s*(\d+)", line).group(1))
            self.metadata["memory_used"] = memory_per_cpu * self.metadata["num_cpu"]

        elif line[1:6].lower() == "%mem=":
            # The maximum amount of memory requested.
            # We need to do some unit juggling.
            mem_str = line.strip().upper()
            # No space is allowed between units, units can probably only be uppercase but convert anyway.
            # Supported units are: 'KB, MB, GB, TB, KW, MW, GW or TW'.
            # TODO: Good opportunity for case here in future.
            units = mem_str[-2:]
            raw_mem = int(mem_str[5:-2])
            # We are converting to bytes.
            if units == "TB":
                memory = raw_mem * 1e12

            elif units == "GB":
                memory = raw_mem * 1e9

            elif units == "MB":
                memory = raw_mem * 1e6

            elif units == "KB":
                memory = raw_mem * 1e3

            elif units == "TW":
                memory = raw_mem * 8e12

            elif units == "GW":
                memory = raw_mem * 8e9

            elif units == "MW":
                memory = raw_mem * 8e6

            elif units == "KW":
                memory = raw_mem * 8e3

            else:
                # No explicit units, default are single words (8-bytes)
                memory = mem_str[5:] * 8

            self.metadata["memory_available"] = int(memory)

        # This block contains some general information as well as coordinates,
        # which could be parsed in the future:
        #
        # Symbolic Z-matrix:
        # Charge =  0 Multiplicity = 1
        # C                     0.73465   0.        0.
        # C                     1.93465   0.        0.
        # C
        # ...
        #
        # It also lists fragments, if there are any, which is potentially valuable:
        #
        # Symbolic Z-matrix:
        # Charge =  0 Multiplicity = 1 in supermolecule
        # Charge =  0 Multiplicity = 1 in fragment      1.
        # Charge =  0 Multiplicity = 1 in fragment      2.
        # B(Fragment=1)         0.06457  -0.0279    0.01364
        # H(Fragment=1)         0.03117  -0.02317   1.21604
        # ...
        #
        # Note, however, that currently we only parse information for the whole system
        # or supermolecule as Gaussian calls it.
        if line.strip() == "Symbolic Z-matrix:":
            self.updateprogress(inputfile, "Symbolic Z-matrix", self.fupdate)

            line = next(inputfile)
            while line.split()[0] == "Charge":
                # For the supermolecule, we can parse the charge and multicplicity.
                regex = r".*=(.*)Mul.*=\s*-?(\d+).*"
                match = re.match(regex, line)
                assert match, f"Something unusual about the line: '{line}'"

                self.set_attribute("charge", int(match.groups()[0]))
                self.set_attribute("mult", int(match.groups()[1]))

                if line.split()[-2] == "fragment":
                    self.nfragments = int(line.split()[-1].strip("."))

                if line.strip()[-13:] == "model system.":
                    self.nmodels = getattr(self, "nmodels", 0) + 1

                line = next(inputfile)

            # The remaining part will allow us to get the atom count.
            # When coordinates are given, there is a blank line at the end, but if
            # there is a Z-matrix here, there will also be variables and we need to
            # stop at those to get the right atom count.
            # Also, in older versions there is bo blank line (G98 regressions),
            # so we need to watch out for leaving the link.
            natom = 0
            while line.split() and "Variables" not in line and "Leave Link" not in line:
                natom += 1
                line = next(inputfile)
            self.set_attribute("natom", natom)

        # Continuing from above, there is not always a symbolic matrix, for example
        # if the Z-matrix was in the input file. In such cases, try to match the
        # line and get at the charge and multiplicity.
        #
        #   Charge =  0 Multiplicity = 1 in supermolecule
        #   Charge =  0 Multiplicity = 1 in fragment  1.
        #   Charge =  0 Multiplicity = 1 in fragment  2.
        if line[1:7] == "Charge" and line.find("Multiplicity") >= 0:
            self.updateprogress(inputfile, "Charge and Multiplicity", self.fupdate)

            if line.split()[-1] == "supermolecule" or (
                "fragment" not in line and "model system" not in line
            ):
                regex = r".*=(.*)Mul.*=\s*-?(\d+).*"
                match = re.match(regex, line)
                assert match, f"Something unusual about the line: '{line}'"

                self.set_attribute("charge", int(match.groups()[0]))
                self.set_attribute("mult", int(match.groups()[1]))

            if line.split()[-2] == "fragment":
                self.nfragments = int(line.split()[-1].strip("."))

            if line.strip()[-13:] == "model system.":
                self.nmodels = getattr(self, "nmodels", 0) + 1

        # Number of atoms is also explicitely printed after the above.
        if line[1:8] == "NAtoms=":
            self.updateprogress(inputfile, "Attributes", self.fupdate)

            natom = int(re.search(r"NAtoms=\s*(\d+)", line).group(1))
            self.set_attribute("natom", natom)

            # Necessary for `if line.strip().split()[0:3] == ["Atom", "AN", "X"]:` block
            if not hasattr(self, "nqmf"):
                match = re.search(r"NQMF=\s*(\d+)", line)
                if match is not None:
                    nqmf = int(match.group(1))
                    if nqmf > 0:
                        self.set_attribute("nqmf", nqmf)

        # Basis set name
        if line[1:15] == "Standard basis":
            self.metadata["basis_set"] = line.split()[2]

        # Solvent information.
        # PCM (the default gaussian solvent method).
        if line[1:34] == "Polarizable Continuum Model (PCM)":
            # Solvent parameters are unique to each solvent model,
            # so they are packaged together to prevent clogging the namespace.
            self.metadata["solvent_params"] = {}

            # Keep looking until dashed only.
            while set(line.strip()) != set("-"):
                line = next(inputfile)

                # PCM has a few different subtypes.
                # Model                : PCM.
                if "Model" in line:
                    self.metadata["solvent_model"] = " ".join(line.split()[2:])[:-1]

                    if self.metadata["solvent_model"] == "PCM":
                        self.metadata["solvent_model"] = "IEFPCM"

                    elif self.metadata["solvent_model"] == "C-PCM":
                        self.metadata["solvent_model"] = "CPCM"

                elif "Atomic radii" in line and line.split()[-1] == "SMD-Coulomb.":
                    self.metadata["solvent_model"] = "SMD-IEFPCM"

                # Solvent by keyword.
                #  Solvent              : Toluene, Eps=   2.374100 Eps(inf)=   2.238315
                # Solvent by definition.
                #  Solvent              : Generic,
                #            Eps                           =   9.000000
                #            Eps(infinity)                 =   2.000000

                elif "Solvent" in line and "Eps=" in line and "Eps(inf)= " in line:
                    split_line = line.split()

                    # Capture the human readable name, as well as params.
                    self.metadata["solvent_name"] = split_line[2][:-1].lower()

                    self.metadata["solvent_params"]["epsilon"] = float(split_line[4])
                    self.metadata["solvent_params"]["epsilon_infinite"] = float(split_line[6])

                elif "Eps(infinity)" in line:
                    # Assume manually specified solvent.
                    self.metadata["solvent_params"]["epsilon_infinite"] = float(line.split()[-1])

                elif "Eps" in line:
                    # Assume manually specified solvent.
                    self.metadata["solvent_params"]["epsilon"] = float(line.split()[-1])

        elif "Reaction Field using a Density IsoSurface Boundary" in line:
            self.metadata["solvent_model"] = "IPCM"

        #  Epsi=   78.3000 Cont =    0.0010
        elif "Epsi=" in line and "Cont =":
            if "solvent_params" not in self.metadata:
                self.metadata["solvent_params"] = {}

            self.metadata["solvent_params"]["epsilon"] = float(line.split()[1])
            self.metadata["solvent_params"]["isovalue"] = float(line.split()[4])

        # Dipole moment
        # e.g. from G09
        #  Dipole moment (field-independent basis, Debye):
        #    X=              0.0000    Y=              0.0000    Z=              0.0930
        # e.g. from G03
        #     X=     0.0000    Y=     0.0000    Z=    -1.6735  Tot=     1.6735
        # need the "field independent" part - ONIOM and other calc use diff formats
        if line[1:39] == "Dipole moment (field-independent basis":
            self.updateprogress(inputfile, "Dipole and Higher Moments", self.fupdate)

            self.reference = [0.0, 0.0, 0.0]
            self.moments = [self.reference]

            tokens = next(inputfile).split()
            # split - dipole would need to be *huge* to fail a split
            # and G03 and G09 use different spacing
            if len(tokens) >= 6:
                dipole = (float(tokens[1]), float(tokens[3]), float(tokens[5]))

            if not hasattr(self, "moments"):
                self.moments = [self.reference, dipole]
            else:
                self.moments.append(dipole)

        if line[1:43] == "Quadrupole moment (field-independent basis":
            # e.g. (g09)
            # Quadrupole moment (field-independent basis, Debye-Ang):
            #   XX=             -6.1213   YY=             -4.2950   ZZ=             -5.4175
            #   XY=              0.0000   XZ=              0.0000   YZ=              0.0000
            # or from g03
            #   XX=    -6.1213   YY=    -4.2950   ZZ=    -5.4175
            quadrupole = {}
            for j in range(2):  # two rows
                line = next(inputfile)
                if line[22] == "=":  # g03 file
                    for i in (1, 18, 35):
                        quadrupole[line[i : i + 4]] = float(line[i + 5 : i + 16])
                else:
                    for i in (1, 27, 53):
                        quadrupole[line[i : i + 4]] = float(line[i + 5 : i + 25])

            lex = sorted(quadrupole.keys())
            quadrupole = [quadrupole[key] for key in lex]

            if not hasattr(self, "moments") or len(self.moments) < 2:
                self.logger.warning("Found quadrupole moments but no previous dipole")
                self.reference = [0.0, 0.0, 0.0]
                self.moments = [self.reference, None, quadrupole]
            else:
                if len(self.moments) == 2:
                    self.moments.append(quadrupole)
                else:
                    assert self.moments[2] == quadrupole

        if line[1:41] == "Octapole moment (field-independent basis":
            # e.g.
            # Octapole moment (field-independent basis, Debye-Ang**2):
            #  XXX=              0.0000  YYY=              0.0000  ZZZ=             -0.1457  XYY=              0.0000
            #  XXY=              0.0000  XXZ=              0.0136  XZZ=              0.0000  YZZ=              0.0000
            #  YYZ=             -0.5848  XYZ=              0.0000
            octapole = {}
            for j in range(2):  # two rows
                line = next(inputfile)
                if line[22] == "=":  # g03 file
                    for i in (1, 18, 35, 52):
                        octapole[line[i : i + 4]] = float(line[i + 5 : i + 16])
                else:
                    for i in (1, 27, 53, 79):
                        octapole[line[i : i + 4]] = float(line[i + 5 : i + 25])

            # last line only 2 moments
            line = next(inputfile)
            if line[22] == "=":  # g03 file
                for i in (1, 18):
                    octapole[line[i : i + 4]] = float(line[i + 5 : i + 16])
            else:
                for i in (1, 27):
                    octapole[line[i : i + 4]] = float(line[i + 5 : i + 25])

            lex = sorted(octapole.keys())
            octapole = [octapole[key] for key in lex]

            if not hasattr(self, "moments") or len(self.moments) < 3:
                self.logger.warning("Found octapole moments but no previous dipole or quadrupole")
                self.reference = [0.0, 0.0, 0.0]
                self.moments = [self.reference, None, None, octapole]
            else:
                if len(self.moments) == 3:
                    self.moments.append(octapole)
                else:
                    assert self.moments[3] == octapole

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
                line = next(inputfile)
                if line[22] == "=":  # g03 file
                    for i in (1, 18, 35, 52):
                        hexadecapole[line[i : i + 4]] = float(line[i + 5 : i + 16])
                else:
                    for i in (1, 27, 53, 79):
                        hexadecapole[line[i : i + 4]] = float(line[i + 5 : i + 25])

            # last line only 3 moments
            line = next(inputfile)
            if line[22] == "=":  # g03 file
                for i in (1, 18, 35):
                    hexadecapole[line[i : i + 4]] = float(line[i + 5 : i + 16])
            else:
                for i in (1, 27, 53):
                    hexadecapole[line[i : i + 4]] = float(line[i + 5 : i + 25])

            lex = sorted(hexadecapole.keys())
            hexadecapole = [hexadecapole[key] for key in lex]

            if not hasattr(self, "moments") or len(self.moments) < 4:
                self.reference = [0.0, 0.0, 0.0]
                self.moments = [self.reference, None, None, None, hexadecapole]
            else:
                if len(self.moments) == 4:
                    self.append_attribute("moments", hexadecapole)
                else:
                    try:
                        numpy.testing.assert_equal(self.moments[4], hexadecapole)
                    except AssertionError:
                        self.logger.warning(
                            "Attribute hexadecapole changed value (%s -> %s)",
                            self.moments[4],
                            hexadecapole,
                        )
                    self.append_attribute("moments", hexadecapole)

        # Catch message about completed optimization.
        if line[1:23] == "Optimization completed":
            self.append_attribute("optdone", len(self.geovalues) - 1)

            assert hasattr(self, "optstatus") and len(self.optstatus) > 0
            self.optstatus[-1] += data.ccData.OPT_DONE

        # Catch message about stopped optimization (not converged).
        if line[1:21] == "Optimization stopped":
            if not hasattr(self, "optdone"):
                self.optdone = []

            assert hasattr(self, "optstatus") and len(self.optstatus) > 0
            self.optstatus[-1] += data.ccData.OPT_UNCONVERGED

        # Extract the atomic numbers and coordinates from the input orientation,
        #   in the event the standard orientation isn't available.
        # Don't extract from Input or Z-matrix orientation in a BOMD run, as only
        #   the final geometry should be kept but extract inputatoms.
        # We also use "inputcoords" to convert "grads" from input orientation
        #   to standard orientation
        if line.find("Input orientation") > -1 or line.find("Z-Matrix orientation") > -1:
            # If this is a counterpoise calculation, this output means that
            #   the supermolecule is now being considered, so we can set:
            self.counterpoise = 0

            self.updateprogress(inputfile, "Attributes", self.cupdate)

            if not self.BOMD and not hasattr(self, "inputcoords"):
                self.inputcoords = []
            self.inputatoms = []

            self.skip_lines(inputfile, ["d", "cols", "cols", "d"])

            atomcoords = []
            line = next(inputfile)
            while list(set(line.strip())) != ["-"]:
                broken = line.split()
                atomno = int(broken[1])
                # Atom with atomno -1 only appears on "Z-Matrix orientation", and excluded on
                #   "Input orientation" or "Standard orientation".
                # We remove this line to keep shape consistency of "atomcoords" and "inputcoords",
                #   so that we can convert "grads" from input orientation to standard orientaion
                #   with rotation matrix calculated from "atomcoords" and "inputcoords"
                if atomno != -1:
                    self.inputatoms.append(atomno)
                    atomcoords.append(list(map(float, broken[3:6])))
                line = next(inputfile)

            if not self.BOMD:
                self.inputcoords.append(atomcoords)

            self.set_attribute("atomnos", numpy.array(self.inputatoms))
            self.set_attribute("natom", len(self.inputatoms))

        if self.BOMD and line.startswith(" Summary information for step"):
            # We keep time and energies_BOMD and coordinates in a dictionary
            #  because steps can be recalculated, and we need to overwite the
            #  previous data

            broken = line.split()
            step = int(broken[-1])
            line = next(inputfile)
            broken = line.split()
            if not hasattr(self, "time"):
                self.set_attribute("time", {step: float(broken[-1])})
            else:
                self.time[step] = float(broken[-1])

            line = next(inputfile)
            broken = line.split(";")[1].split()
            ene = utils.float(broken[-1])
            if not hasattr(self, "energies_BOMD"):
                self.set_attribute("energies_BOMD", {step: ene})
            else:
                self.energies_BOMD[step] = ene

            self.updateprogress(inputfile, "Attributes", self.cupdate)

            if not hasattr(self, "atomcoords_BOMD"):
                self.atomcoords_BOMD = {}
            # self.inputatoms = []

            self.skip_lines(inputfile, ["EKin", "Angular", "JX", "Total", "Total", "Cartesian"])

            atomcoords = []
            line = next(inputfile)
            while "MW cartesian" not in line:
                broken = line.split()
                atomcoords.append(list(map(utils.float, (broken[3], broken[5], broken[7]))))
                #    self.inputatoms.append(int(broken[1]))
                line = next(inputfile)

            self.atomcoords_BOMD[step] = atomcoords

            # self.set_attribute('atomnos', self.inputatoms)
            # self.set_attribute('natom', len(self.inputatoms))

        # Extract the atomic masses.
        # Typical section:
        #                    Isotopes and Nuclear Properties:
        # (Nuclear quadrupole moments (NQMom) in fm**2, nuclear magnetic moments (NMagM)
        # in nuclear magnetons)
        #
        #  Atom         1           2           3           4           5           6           7           8           9          10
        # IAtWgt=          12          12          12          12          12           1           1           1          12          12
        # AtmWgt=  12.0000000  12.0000000  12.0000000  12.0000000  12.0000000   1.0078250   1.0078250   1.0078250  12.0000000  12.0000000
        # NucSpn=           0           0           0           0           0           1           1           1           0           0
        # AtZEff=  -3.6000000  -3.6000000  -3.6000000  -3.6000000  -3.6000000  -1.0000000  -1.0000000  -1.0000000  -3.6000000  -3.6000000
        # NQMom=    0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
        # NMagM=    0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   2.7928460   2.7928460   2.7928460   0.0000000   0.0000000
        # ... with blank lines dividing blocks of ten, and Leave Link 101 at the end.
        # This is generally parsed before coordinates, so atomnos is not defined.
        # Note that in Gaussian03 the comments are not there yet and the labels are different.
        if line.strip() == "Isotopes and Nuclear Properties:":
            line = next(inputfile)
            while line[1:16] != "Leave Link  101":
                if line[1:8] == "AtmWgt=":
                    self.extend_attribute("atommasses", list(map(float, line.split()[1:])))
                line = next(inputfile)

        # Symmetry: point group
        if "Symmetry turned off" in line:
            self.set_attribute("uses_symmetry", False)
        if "Full point group" in line:
            point_group_detected = line.split()[3].lower()
            if self.uses_symmetry:
                while "Largest Abelian subgroup" not in line:
                    line = next(inputfile)
                    if "Leave Link" in line:
                        # TODO To handle this properly, it needs to be
                        # calculated from the full point group.
                        point_group_used = point_group_detected
                        break
                if "Leave Link" not in line:
                    point_group_used = line.split()[3].lower()
            else:
                point_group_used = "c1"
            self.metadata["symmetry_detected"] = point_group_detected
            self.metadata["symmetry_used"] = point_group_used

        # Symmetry: ordering of irreducible representations
        if "symmetry adapted cartesian basis functions" in line:
            irrep = self.normalisesym(line.split()[-2])
            self.append_attribute("symlabels", irrep)

        # Extract the atomic numbers and coordinates of the atoms.
        if line.strip() == "Standard orientation:":
            self.updateprogress(inputfile, "Attributes", self.cupdate)

            # If this is a counterpoise calculation, this output means that
            #   the supermolecule is now being considered, so we can set:
            self.counterpoise = 0

            if not hasattr(self, "atomcoords"):
                self.atomcoords = []

            self.skip_lines(inputfile, ["d", "cols", "cols", "d"])

            atomnos = []
            atomcoords = []
            line = next(inputfile)
            while list(set(line.strip())) != ["-"]:
                broken = line.split()
                atomnos.append(int(broken[1]))
                atomcoords.append(list(map(float, broken[-3:])))
                line = next(inputfile)
            self.atomcoords.append(atomcoords)

            self.set_attribute("natom", len(atomnos))
            self.set_attribute("atomnos", atomnos)

        # This is a bit of a hack for regression Gaussian09/BH3_fragment_guess.pop_minimal.log
        # to skip output for all fragments, assuming the supermolecule is always printed first.
        # Eventually we want to make this more general, or even better parse the output for
        # all fragment, but that will happen in a newer version of cclib.
        if line[1:16] == "Fragment guess:" and getattr(self, "nfragments", 0) > 1:
            if "full" not in line:
                raise StopParsing()

        # Another hack for regression Gaussian03/ortho_prod_freq.log, which is an ONIOM job.
        # Basically for now we stop parsing after the output for the real system, because
        # currently we don't support changes in system size or fragments in cclib. When we do,
        # we will want to parse the model systems, too, and that is what nmodels could track.
        if (
            "ONIOM: generating point" in line
            and line.strip()[-13:] == "model system."
            and getattr(self, "nmodels", 0) > 0
        ):
            while not line[1:30] == "ONIOM: Integrating ONIOM file":
                line = next(inputfile)

        # With the gfinput keyword, the atomic basis set functions are:
        #
        # AO basis set in the form of general basis input (Overlap normalization):
        #  1 0
        # S   3 1.00       0.000000000000
        #      0.7161683735D+02  0.1543289673D+00
        #      0.1304509632D+02  0.5353281423D+00
        #      0.3530512160D+01  0.4446345422D+00
        # SP   3 1.00       0.000000000000
        #      0.2941249355D+01 -0.9996722919D-01  0.1559162750D+00
        #      0.6834830964D+00  0.3995128261D+00  0.6076837186D+00
        #      0.2222899159D+00  0.7001154689D+00  0.3919573931D+00
        # ****
        #  2 0
        # S   3 1.00       0.000000000000
        #      0.7161683735D+02  0.1543289673D+00
        # ...
        #
        # The same is also printed when the gfprint keyword is used, but the
        # interstitial lines differ and there are no stars between atoms:
        #
        # AO basis set (Overlap normalization):
        # Atom C1       Shell     1 S   3     bf    1 -     1          0.509245180608         -2.664678875191          0.000000000000
        #       0.7161683735D+02  0.1543289673D+00
        #       0.1304509632D+02  0.5353281423D+00
        #       0.3530512160D+01  0.4446345422D+00
        # Atom C1       Shell     2 SP   3    bf    2 -     5          0.509245180608         -2.664678875191          0.000000000000
        #       0.2941249355D+01 -0.9996722919D-01  0.1559162750D+00
        # ...

        # ONIOM calculations result basis sets reported for atoms that are not in order of atom number which breaks this code (line 390 relies on atoms coming in order)
        if line[1:13] == "AO basis set" and not self.oniom:
            self.gbasis = []

            # For counterpoise fragment calcualtions, skip these lines.
            if self.counterpoise != 0:
                return

            atom_line = next(inputfile)
            self.gfprint = atom_line.split()[0] == "Atom"
            self.gfinput = not self.gfprint

            # Note how the shell information is on a separate line for gfinput,
            # whereas for gfprint it is on the same line as atom information.
            if self.gfinput:
                shell_line = next(inputfile)

            shell = []
            while len(self.gbasis) < self.natom:
                if self.gfprint:
                    cols = atom_line.split()
                    subshells = cols[4]
                    ngauss = int(cols[5])
                else:
                    cols = shell_line.split()
                    subshells = cols[0]
                    ngauss = int(cols[1])

                parameters = []
                for ig in range(ngauss):
                    line = next(inputfile)
                    parameters.append(list(map(utils.float, line.split())))
                for iss, ss in enumerate(subshells):
                    contractions = []
                    for param in parameters:
                        exponent = param[0]
                        coefficient = param[iss + 1]
                        contractions.append((exponent, coefficient))
                    subshell = (ss, contractions)
                    shell.append(subshell)

                if self.gfprint:
                    line = next(inputfile)
                    if line.split()[0] == "Atom":
                        atomnum = int(re.sub(r"\D", "", line.split()[1]))
                        if atomnum == len(self.gbasis) + 2:
                            self.gbasis.append(shell)
                            shell = []
                        atom_line = line
                    else:
                        self.gbasis.append(shell)
                else:
                    line = next(inputfile)
                    if line.strip() == "****":
                        self.gbasis.append(shell)
                        shell = []
                        atom_line = next(inputfile)
                        shell_line = next(inputfile)
                    else:
                        shell_line = line

        if "Dispersion energy=" in line:
            self.append_attribute("dispersionenergies", float(line.split()[-2]))

        # Find the targets for SCF convergence (QM calcs).
        # Not for BOMD as targets are not available in the summary
        if not self.BOMD and line[1:44] == "Requested convergence on RMS density matrix":
            if not hasattr(self, "scftargets"):
                self.scftargets = []
            # The following can happen with ONIOM which are mixed SCF
            # and semi-empirical
            if isinstance(self.scftargets, numpy.ndarray):
                self.scftargets = []

            scftargets = []
            # The RMS density matrix.
            scftargets.append(utils.float(line.split("=")[1].split()[0]))
            line = next(inputfile)
            # The MAX density matrix.
            scftargets.append(utils.float(line.strip().split("=")[1][:-1]))
            line = next(inputfile)
            # For G03, there's also the energy (not for G98).
            if line[1:10] == "Requested":
                scftargets.append(utils.float(line.strip().split("=")[1][:-1]))

            self.scftargets.append(scftargets)

        # Extract SCF convergence information (QM calcs).
        if line[1:10] == "Cycle   1":
            if not hasattr(self, "scfvalues"):
                self.scfvalues = []

            scfvalues = []
            line = next(inputfile)
            while line.find("SCF Done") == -1:
                self.updateprogress(inputfile, "QM convergence", self.fupdate)

                # SCI-PCM solvent info appears in each SCF section...
                #  Compute SCI-PCM surface.
                if "Compute SCI-PCM surface" in line:
                    self.metadata["solvent_model"] = "SCIPCM"

                # For SCI-PCM.
                # Dielectric constant of solvent =     2.374100"
                if line[1:33] == "Dielectric constant of solvent =":
                    if "solvent_params" not in self.metadata:
                        self.metadata["solvent_params"] = {}

                    self.metadata["solvent_params"]["epsilon"] = float(line.split()[-1])

                if line.find(" E=") == 0:
                    self.logger.debug(line)

                #  RMSDP=3.74D-06 MaxDP=7.27D-05 DE=-1.73D-07 OVMax= 3.67D-05
                # or
                #  RMSDP=1.13D-05 MaxDP=1.08D-04              OVMax= 1.66D-04
                if line.find(" RMSDP") == 0:
                    # Fields of interest:
                    # RMSDP
                    # MaxDP
                    # (DE) -> Only add the energy if it's a target criteria

                    matches = self.re_scinot.findall(line)
                    matches = {match[0]: utils.float(match[1]) for match in matches}
                    scfvalues_step = [
                        matches.get("RMSDP", numpy.nan),
                        matches.get("MaxDP", numpy.nan),
                    ]
                    if hasattr(self, "scftargets") and len(self.scftargets[0]) == 3:
                        scfvalues_step.append(matches.get("DE", numpy.nan))
                    scfvalues.append(scfvalues_step)

                try:
                    line = next(inputfile)
                # May be interupted by EOF.
                except StopIteration:
                    self.logger.warning("File terminated before end of last SCF!")
                    break

            self.scfvalues.append(scfvalues)

        # Extract SCF convergence information (AM1, INDO and other semi-empirical calcs).
        # The output (for AM1) looks like this:
        # Ext34=T Pulay=F Camp-King=F BShift= 0.00D+00
        # It=  1 PL= 0.103D+01 DiagD=T ESCF=     31.564733 Diff= 0.272D+02 RMSDP= 0.152D+00.
        # It=  2 PL= 0.114D+00 DiagD=T ESCF=      7.265370 Diff=-0.243D+02 RMSDP= 0.589D-02.
        # ...
        # It= 11 PL= 0.184D-04 DiagD=F ESCF=      4.687669 Diff= 0.260D-05 RMSDP= 0.134D-05.
        # It= 12 PL= 0.105D-04 DiagD=F ESCF=      4.687669 Diff=-0.686D-07 RMSDP= 0.215D-05.
        # 4-point extrapolation.
        # It= 13 PL= 0.110D-05 DiagD=F ESCF=      4.687669 Diff=-0.111D-06 RMSDP= 0.653D-07.
        # Energy=    0.172272018655 NIter=  14.
        if line[1:4] == "It=":
            scftargets = numpy.array([1e-7], "d")  # This is the target value for the rms
            scfvalues = [[]]

            while line.find(" Energy") == -1:
                self.updateprogress(inputfile, "AM1 Convergence")

                if line[1:4] == "It=":
                    parts = line.strip().split()
                    scfvalues[0].append(utils.float(parts[-1][:-1]))

                line = next(inputfile)

                # If an AM1 or INDO guess is used (Guess=INDO in the input, for example),
                # this will be printed after a single iteration, so that is the line
                # that should trigger a break from this loop. At least that's what we see
                # for regression Gaussian/Gaussian09/guessIndo_modified_ALT.out
                if line[:14] == " Initial guess":
                    break

            # Attach the attributes to the object Only after the energy is found .
            if line.find(" Energy") == 0:
                self.scftargets = scftargets
                self.scfvalues = scfvalues

        # Note: this needs to follow the section where 'SCF Done' is used
        #   to terminate a loop when extracting SCF convergence information.
        if not self.BOMD and line[1:9] == "SCF Done":
            t1 = line.split()[2]
            if t1 == "E(RHF)":
                self.metadata["methods"].append("HF")
            else:
                self.metadata["methods"].append("DFT")
                self.metadata["functional"] = t1[t1.index("(") + 2 : t1.rindex(")")]

            self.append_attribute("scfenergies", utils.float(line.split()[4]))
        # gmagoon 5/27/09: added scfenergies reading for PM3 case
        # Example line: " Energy=   -0.077520562724 NIter=  14."
        # See regression Gaussian03/QVGXLLKOCUKJST-UHFFFAOYAJmult3Fixed.out
        if line[1:8] == "Energy=":
            self.append_attribute("scfenergies", utils.float(line.split()[1]))

        # Total energies after Moller-Plesset corrections.
        # Second order correction is always first, so its first occurance
        #   triggers creation of mpenergies (list of lists of energies).
        # Further MP2 corrections are appended as found.
        #
        # Example MP2 output line:
        #  E2 =    -0.9505918144D+00 EUMP2 =    -0.28670924198852D+03
        # Warning! this output line is subtly different for MP3/4/5 runs
        # Newer versions of gausian introduced a space between 'EUMP2' and '='...
        if "EUMP2 =" in line[27:36] or "EUMP2=" in line[27:35]:
            self.metadata["methods"].append("MP2")
            self.append_attribute("mpenergies", [utils.float(line.split("=")[2])])

        # Example MP3 output line:
        #  E3=       -0.10518801D-01     EUMP3=      -0.75012800924D+02
        if line[34:40] == "EUMP3=":
            self.metadata["methods"].append("MP3")
            self.mpenergies[-1].append(utils.float(line.split("=")[2]))

        # Example MP4 output lines:
        #  E4(DQ)=   -0.31002157D-02        UMP4(DQ)=   -0.75015901139D+02
        #  E4(SDQ)=  -0.32127241D-02        UMP4(SDQ)=  -0.75016013648D+02
        #  E4(SDTQ)= -0.32671209D-02        UMP4(SDTQ)= -0.75016068045D+02
        # Energy for most substitutions is used only (SDTQ by default)
        if line[34:43] == "UMP4(DQ)=":
            self.metadata["methods"].append("MP4")

            mp4energy = utils.float(line.split("=")[2])
            line = next(inputfile)
            if line[34:44] == "UMP4(SDQ)=":
                mp4energy = utils.float(line.split("=")[2])
                line = next(inputfile)
                if line[34:45] == "UMP4(SDTQ)=":
                    mp4energy = utils.float(line.split("=")[2])
            self.mpenergies[-1].append(mp4energy)

        # Example MP5 output line:
        #  DEMP5 =  -0.11048812312D-02 MP5 =  -0.75017172926D+02
        if line[29:34] == "MP5 =":
            self.metadata["methods"].append("MP5")
            self.mpenergies[-1].append(utils.float(line.split("=")[2]))

        # Total energies after Coupled Cluster corrections.
        # Second order MBPT energies (MP2) are also calculated for these runs,
        # but the output is the same as when parsing for mpenergies.
        # Read the consecutive correlated energies
        # but append only the last one to ccenergies.
        # Only the highest level energy is appended - ex. CCSD(T), not CCSD.
        if line[1:10] == "DE(Corr)=" and line[27:35] == "E(CORR)=":
            self.metadata["methods"].append("CCSD")
            self.ccenergy = utils.float(line.split()[3])
        if line[1:14] == "T1 Diagnostic":
            self.metadata["t1_diagnostic"] = utils.float(line.split()[-1])
        if line[1:10] == "T5(CCSD)=":
            line = next(inputfile)
            if line[1:9] == "CCSD(T)=":
                self.metadata["methods"].append("CCSD-T")
                self.ccenergy = utils.float(line.split()[1])
        # Find step number for current optimization/IRC
        # Matches "Step number  123", "Pt XX Step number 123" and "PtXXX Step number 123"
        if " Step number" in line:
            step = int(line.split()[line.split().index("Step") + 2])
            self.append_attribute("optstatus", data.ccData.OPT_UNKNOWN)
            if step == 1:
                self.optstatus[-1] += data.ccData.OPT_NEW

        # Geometry convergence information.
        if line[49:59] == "Converged?":
            if not hasattr(self, "geotargets"):
                self.geovalues = []
                self.geotargets = numpy.array([0.0, 0.0, 0.0, 0.0], "d")
            allconverged = True
            newlist = [0] * 4
            for i in range(4):
                line = next(inputfile)
                self.logger.debug(line)
                parts = line.split()
                if "NO" in parts[-1]:
                    allconverged = False  # noqa: F841
                try:
                    value = utils.float(parts[2])
                except ValueError:
                    self.logger.error(
                        "Problem parsing the value for geometry optimisation: %s is not a number.",
                        parts[2],
                    )
                else:
                    newlist[i] = value
                self.geotargets[i] = utils.float(parts[3])

            self.geovalues.append(newlist)

        # Gradients.
        # Read in the cartesian energy gradients (forces) from a block like this:
        # -------------------------------------------------------------------
        # Center     Atomic                   Forces (Hartrees/Bohr)
        # Number     Number              X              Y              Z
        # -------------------------------------------------------------------
        # 1          1          -0.012534744   -0.021754635   -0.008346094
        # 2          6           0.018984731    0.032948887   -0.038003451
        # 3          1          -0.002133484   -0.006226040    0.023174772
        # 4          1          -0.004316502   -0.004968213    0.023174772
        #           -2          -0.001830728   -0.000743108   -0.000196625
        # ------------------------------------------------------------------
        #
        # The "-2" line is for a dummy atom
        #
        # Then optimization is done in internal coordinates, Gaussian also
        # print the forces in internal coordinates, which can be produced from
        # the above. This block looks like this:
        # Variable       Old X    -DE/DX   Delta X   Delta X   Delta X     New X
        #                                 (Linear)    (Quad)   (Total)
        #   ch        2.05980   0.01260   0.00000   0.01134   0.01134   2.07114
        #   hch        1.75406   0.09547   0.00000   0.24861   0.24861   2.00267
        #   hchh       2.09614   0.01261   0.00000   0.16875   0.16875   2.26489
        #         Item               Value     Threshold  Converged?

        # We could get the gradients in BOMD, but it is more complex because
        #   they are not in the summary, and they are not as relevant as for
        #   an optimization
        if not self.BOMD and line[37:43] == "Forces":
            if not hasattr(self, "grads"):
                self.grads = []

            self.skip_lines(inputfile, ["header", "d"])

            forces = []
            line = next(inputfile)
            while list(set(line.strip())) != ["-"]:
                tmpforces = []
                for N in range(3):  # Fx, Fy, Fz
                    force = line[23 + N * 15 : 38 + N * 15]
                    if force.startswith("*"):
                        force = "NaN"
                    tmpforces.append(float(force))
                forces.append(tmpforces)
                line = next(inputfile)
            self.grads.append(forces)

        if "Number of optimizations in scan" in line:
            self.scan_length = int(line.split()[-1])

        # All PES scans have a list of initial parameters from which we
        # can get the names and more.
        #
        #                            ----------------------------
        #                           !    Initial Parameters    !
        #                           ! (Angstroms and Degrees)  !
        # --------------------------                            --------------------------
        # ! Name  Definition              Value          Derivative Info.                !
        # --------------------------------------------------------------------------------
        # ! R1    R(1,2)                  1.4212         estimate D2E/DX2                !
        # ! R2    R(1,14)                 1.4976         estimate D2E/DX2                !
        # ...
        if "Initial Parameters" in line:
            self.scannames_all = []
            self.scannames_scanned = []
            self.skip_lines(inputfile, ["units", "d", "header", "d"])
            line = next(inputfile)
            while line.strip()[0] == "!":
                name = line.split()[1]
                definition = line.split()[2]
                self.scannames_all.append(name)
                if line.split()[4] == "Scan":
                    self.scannames_scanned.append(name)
                    self.append_attribute("scannames", definition)
                line = next(inputfile)

        # Extract unrelaxed PES scan data, which looks something like:
        #
        # Summary of the potential surface scan:
        #   N       A          SCF
        # ----  ---------  -----------
        #    1   109.0000    -76.43373
        #    2   119.0000    -76.43011
        #    3   129.0000    -76.42311
        #    4   139.0000    -76.41398
        #    5   149.0000    -76.40420
        #    6   159.0000    -76.39541
        #    7   169.0000    -76.38916
        #    8   179.0000    -76.38664
        #    9   189.0000    -76.38833
        #   10   199.0000    -76.39391
        #   11   209.0000    -76.40231
        # ----  ---------  -----------
        if "Summary of the potential surface scan:" in line:
            colmnames = next(inputfile)
            if not hasattr(self, "scannames"):
                self.set_attribute("scannames", colmnames.split()[1:-1])

            hyphens = next(inputfile)
            line = next(inputfile)
            scanparm = [[] for _ in range(len(self.scannames))]
            while line != hyphens:
                broken = line.split()
                self.append_attribute("scanenergies", float(broken[-1]))
                for idx, p in enumerate(broken[1:-1]):
                    scanparm[idx].append(float(p))
                # self.append_attribute('scanparm', [float(p) for p in broken[1:-1]])
                line = next(inputfile)
            self.set_attribute("scanparm", scanparm)

        # Extract relaxed (optimized) PES scan data, for which the form
        # of the output is transposed:
        #
        #  Summary of Optimized Potential Surface Scan (add -382.0 to energies):
        #                            1         2         3         4         5
        #      Eigenvalues --    -0.30827  -0.30695  -0.30265  -0.29955  -0.30260
        #            R1           1.42115   1.42152   1.42162   1.42070   1.42071
        #            R2           1.49761   1.49787   1.49855   1.49901   1.49858
        #            R3           1.42245   1.42185   1.42062   1.42048   1.42147
        #            R4           1.40217   1.40236   1.40306   1.40441   1.40412
        # ...
        if "Summary of Optimized Potential Surface Scan" in line:
            # The base energy separation is version dependent, and first
            # appears in Gaussian16.
            base_energy = 0.0
            if "add" in line and "to energies" in line:
                base_energy = float(line.split()[-3])

            scanenergies = []
            scanparm = [[] for _ in range(len(self.scannames))]
            while len(scanenergies) != self.scan_length:
                line = next(inputfile)
                indices = [int(i) for i in line.split()]
                widths = [10] * len(indices)
                splitter = utils.WidthSplitter(widths)

                line = next(inputfile)
                eigenvalues_in_line = line[21:].rstrip()
                assert len(eigenvalues_in_line) == sum(widths)
                cols = list(splitter.split(eigenvalues_in_line))
                try:
                    eigenvalues = [float(e) for e in cols]
                    eigenvalues = [base_energy + e for e in eigenvalues]
                except ValueError:
                    eigenvalues = [numpy.nan for _ in cols]
                assert len(eigenvalues) == len(indices)
                scanenergies.extend(eigenvalues)

                for name in self.scannames_all:
                    line = next(inputfile)
                    assert line.split()[0] == name
                    if name in self.scannames_scanned:
                        iname = self.scannames_scanned.index(name)
                        params_in_line = line[21:].rstrip()
                        assert len(params_in_line) == sum(widths)
                        params = [float(v) for v in splitter.split(params_in_line)]
                        scanparm[iname].extend(params)

            self.set_attribute("scanenergies", scanenergies)
            self.set_attribute("scanparm", scanparm)

        # Orbital symmetries.
        if line[1:20] == "Orbital symmetries:" and not hasattr(self, "mosyms"):
            # For counterpoise fragments, skip these lines.
            if self.counterpoise != 0:
                return

            self.updateprogress(inputfile, "MO Symmetries", self.fupdate)

            self.mosyms = [[]]
            line = next(inputfile)
            unres = False
            if line.find("Alpha Orbitals") == 1:
                unres = True
                line = next(inputfile)
            i = 0
            while len(line) > 18 and line[17] == "(":
                if line.find("Virtual") >= 0:
                    self.homos = [i - 1]
                parts = line[17:].split()
                for x in parts:
                    self.mosyms[0].append(self.normalisesym(x.strip("()")))
                    i += 1
                line = next(inputfile)
            if unres:
                line = next(inputfile)
                # Repeat with beta orbital information
                i = 0
                self.mosyms.append([])
                while len(line) > 18 and line[17] == "(":
                    if line.find("Virtual") >= 0:
                        # Here we consider beta
                        # If there was also an alpha virtual orbital,
                        #  we will store two indices in the array
                        # Otherwise there is no alpha virtual orbital,
                        #  only beta virtual orbitals, and we initialize
                        #  the array with one element. See the regression
                        #  QVGXLLKOCUKJST-UHFFFAOYAJmult3Fixed.out
                        #  donated by Gregory Magoon (gmagoon).
                        if hasattr(self, "homos"):
                            # Extend the array to two elements
                            # 'HOMO' indexes the HOMO in the arrays
                            self.homos.append(i - 1)
                        else:
                            # 'HOMO' indexes the HOMO in the arrays
                            self.homos = [i - 1]
                    parts = line[17:].split()
                    for x in parts:
                        self.mosyms[1].append(self.normalisesym(x.strip("()")))
                        i += 1
                    line = next(inputfile)

            # Some calculations won't explicitely print the number of basis sets used,
            # and will occasionally drop some without warning. We can infer the number,
            # however, from the MO symmetries printed here. Specifically, this fixes
            # regression Gaussian/Gaussian09/dvb_sp_terse.log (#23 on github).
            self.set_attribute("nmo", len(self.mosyms[-1]))

        # Alpha/Beta electron eigenvalues.
        if line[1:6] == "Alpha" and line.find("eigenvalues") >= 0:
            # For counterpoise fragments, skip these lines.
            if self.counterpoise != 0:
                return

            # For ONIOM calcs, ignore this section in order to bypass assertion failure.
            if self.oniom:
                return

            self.updateprogress(inputfile, "Eigenvalues", self.fupdate)
            self.moenergies = [[]]
            HOMO = -2

            while line.find("Alpha") == 1:
                if line.split()[1] == "virt." and HOMO == -2:
                    # If there aren't any symmetries, this is a good way to find the HOMO.
                    HOMO = len(self.moenergies[0]) - 1
                    self.homos = [HOMO]

                # Convert to floats and append to moenergies, but sometimes Gaussian
                #  doesn't print correctly so test for ValueError (bug 1756789).
                part = line[28:]
                i = 0
                while i * 10 + 4 < len(part):
                    s = part[i * 10 : (i + 1) * 10]
                    try:
                        x = utils.float(s)
                    except ValueError:
                        x = numpy.nan
                    self.moenergies[0].append(x)
                    i += 1
                line = next(inputfile)

            # If, at this point, self.homos is unset, then there were not
            # any alpha virtual orbitals
            if not hasattr(self, "homos"):
                HOMO = len(self.moenergies[0]) - 1
                self.homos = [HOMO]

            if line.find("Beta") == 2:
                self.moenergies.append([])

            HOMO = -2
            while line.find("Beta") == 2:
                if line.split()[1] == "virt." and HOMO == -2:
                    # If there aren't any symmetries, this is a good way to find the HOMO.
                    # Also, check for consistency if homos was already parsed.
                    HOMO = len(self.moenergies[1]) - 1
                    self.homos.append(HOMO)

                part = line[28:]
                i = 0
                while i * 10 + 4 < len(part):
                    x = part[i * 10 : (i + 1) * 10]
                    self.moenergies[1].append(utils.float(x))
                    i += 1
                line = next(inputfile)

            self.moenergies = [numpy.array(x, "d") for x in self.moenergies]

        # Start of the IR/Raman frequency section.
        # Caution is advised here, as additional frequency blocks
        #   can be printed by Gaussian (with slightly different formats),
        #   often doubling the information printed.
        # See, for a non-standard exmaple, regression Gaussian98/test_H2.log
        # If either the Gaussian freq=hpmodes keyword or IOP(7/33=1) is used,
        # an extra frequency block with higher-precision vibdisps is
        # printed before the normal frequency block.
        # Note that the code parses only the vibsyms and vibdisps
        # from the high-precision block, but parses vibsyms, vibfreqs, vibfconsts,
        # vibramans, vibrmasses and vibirs from the normal block. vibsyms parsed
        # from the high-precision block are discarded and replaced by those
        # from the normal block while the high-precision vibdisps, if present,
        # are used to overwrite default-precision vibdisps at the end of the parse.
        if line[1:14] == "Harmonic freq":  # This matches in both freq block types
            self.updateprogress(inputfile, "Frequency Information", self.fupdate)

            # The whole block should not have any blank lines.
            while line.strip() != "":
                # The line with indices
                if line[1:15].strip() == "" and line[15:60].split()[0].isdigit():
                    freqbase = int(line[15:60].split()[0])
                    if freqbase == 1 and hasattr(self, "vibsyms"):
                        # we are coming accross duplicated information.
                        # We might be be parsing a default-precision block having
                        # already parsed (only) vibsyms and displacements from
                        # the high-precision block, or might be encountering
                        # a second low-precision block (see e.g. 25DMF_HRANH.log
                        # regression).
                        self.vibsyms = []
                        if hasattr(self, "vibirs"):
                            self.vibirs = []
                        if hasattr(self, "vibfreqs"):
                            self.vibfreqs = []
                        if hasattr(self, "vibramans"):
                            self.vibramans = []
                        if hasattr(self, "vibrmasses"):
                            self.vibrmasses = []
                        if hasattr(self, "vibfconsts"):
                            self.vibfconsts = []
                        if hasattr(self, "vibdisps"):
                            self.vibdisps = []

                # Lines with symmetries and symm. indices begin with whitespace.
                if line[1:15].strip() == "" and not line[15:60].split()[0].isdigit():
                    if not hasattr(self, "vibsyms"):
                        self.vibsyms = []
                    syms = line.split()
                    self.vibsyms.extend(syms)

                if line[1:15] == "Frequencies --":  # note: matches low-precision block, and
                    if not hasattr(self, "vibfreqs"):
                        self.vibfreqs = []

                    freqs = [utils.float(f) for f in line[15:].split()]
                    self.vibfreqs.extend(freqs)

                if line[1:15] == "Red. masses --":  # note: matches only low-precision block
                    if not hasattr(self, "vibrmasses"):
                        self.vibrmasses = []

                    rmasses = [utils.float(f) for f in line[15:].split()]
                    self.vibrmasses.extend(rmasses)

                if line[1:15] == "Frc consts  --":  # note: matches only low-precision block
                    if not hasattr(self, "vibfconsts"):
                        self.vibfconsts = []

                    fconsts = [utils.float(f) for f in line[15:].split()]
                    self.vibfconsts.extend(fconsts)

                if line[1:15] == "IR Inten    --":  # note: matches only low-precision block
                    if not hasattr(self, "vibirs"):
                        self.vibirs = []

                    irs = []
                    for ir in line[15:].split():
                        try:
                            irs.append(utils.float(ir))
                        except ValueError:
                            irs.append(utils.float("nan"))
                    self.vibirs.extend(irs)

                if line[1:15] == "Raman Activ --":  # note: matches only low-precision block
                    if not hasattr(self, "vibramans"):
                        self.vibramans = []

                    ramans = []
                    for raman in line[15:].split():
                        try:
                            ramans.append(utils.float(raman))
                        except ValueError:
                            ramans.append(utils.float("nan"))

                    self.vibramans.extend(ramans)

                # Block with (default-precision) displacements should start with this.
                #                     1                      2                      3
                #                     A                      A                      A
                # Frequencies --   370.7936               370.7987               618.0103
                # Red. masses --     2.3022                 2.3023                 1.9355
                # Frc consts  --     0.1865                 0.1865                 0.4355
                # IR Inten    --     0.0000                 0.0000                 0.0000
                #  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z
                #     1   6     0.00   0.00  -0.04     0.00   0.00   0.19     0.00   0.00   0.12
                #     2   6     0.00   0.00   0.19     0.00   0.00  -0.06     0.00   0.00  -0.12

                if line.strip().split()[0:3] == ["Atom", "AN", "X"]:
                    if not hasattr(self, "vibdisps"):
                        self.vibdisps = []
                    disps = []
                    if not hasattr(self, "nqmf"):
                        self.set_attribute("nqmf", self.natom)
                    for n in range(self.nqmf):
                        line = next(inputfile)
                        numbers = [float(s) for s in line[10:].split()]
                        N = len(numbers) // 3
                        if not disps:
                            for n in range(N):
                                disps.append([])
                        for n in range(N):
                            disps[n].append(numbers[3 * n : 3 * n + 3])
                    self.vibdisps.extend(disps)

                # Block with high-precision (freq=hpmodes) displacements should start with this.
                #                           1         2         3         4         5
                #                           A         A         A         A         A
                #       Frequencies ---   370.7936  370.7987  618.0103  647.7864  647.7895
                #    Reduced masses ---     2.3022    2.3023    1.9355    6.4600    6.4600
                #   Force constants ---     0.1865    0.1865    0.4355    1.5971    1.5972
                #    IR Intensities ---     0.0000    0.0000    0.0000    0.0000    0.0000
                # Coord Atom Element:
                #   1     1     6          0.00000   0.00000   0.00000  -0.18677   0.05592
                #   2     1     6          0.00000   0.00000   0.00000   0.28440   0.21550
                #   3     1     6         -0.04497   0.19296   0.11859   0.00000   0.00000
                #   1     2     6          0.00000   0.00000   0.00000   0.03243   0.37351
                #   2     2     6          0.00000   0.00000   0.00000   0.14503  -0.06117
                #   3     2     6          0.18959  -0.05753  -0.11859   0.00000   0.00000
                if line.strip().split()[0:3] == ["Coord", "Atom", "Element:"]:
                    # Wait until very end of parsing to assign vibdispshp to self.vibdisps
                    # as otherwise the higher precision displacements will be overwritten
                    # by low precision displacements which are printed further down file
                    if not hasattr(self, "vibdispshp"):
                        self.vibdispshp = []

                    disps = []
                    for n in range(3 * self.natom):
                        line = next(inputfile)
                        numbers = [float(s) for s in line[16:].split()]
                        atomindex = int(line[4:10]) - 1  # atom index, starting at zero
                        numbermodes = len(numbers)

                        if not disps:
                            for mode in range(numbermodes):
                                # For each mode, make list of list [atom][coord_index]
                                disps.append([[] for x in range(0, self.natom)])
                        for mode in range(numbermodes):
                            disps[mode][atomindex].append(numbers[mode])
                    self.vibdispshp.extend(disps)

                line = next(inputfile)

        # Metadata for excited states methods
        #
        # For HF/DFT level ES, this is our trigger line:
        #  MDV=  1342177280 DFT=T DoStab=F Mixed=T DoRPA=F DoScal=F NonHer=F
        # Comparing DFT=T/F and DoRPA=T/F allows us to distinguish CIS, RPA, TD-DFT and TDA.
        if "MDV=" in line and "DFT=" in line and "DoRPA=" in line:
            if "DFT=T" in line:
                if "RPA=T" in line:
                    method = "TD-DFT"

                else:
                    method = "TDA"

            else:
                if "RPA=T" in line:
                    method = "RPA"

                else:
                    method = "CIS"

            self.metadata["excited_states_method"] = method

        if line.strip() == "EOM-CCSD":
            self.metadata["excited_states_method"] = "EOM-CCSD"

        # Electronic transitions.
        if line[1:14] == "Excited State":
            # Excited State 1:
            et_index = float(line.split()[2][:-1])

            if not hasattr(self, "etenergies") or et_index <= self.last_et:
                self.etenergies = []
                self.etoscs = []
                self.etsyms = []
                self.etsecs = []

            # Keep track of the highest excited state, so we can detect when we enter a new
            # section (the 'highest' excited state will be the same or lower as the last one).
            self.last_et = et_index

            # Need to deal with lines like:
            # (restricted calc)
            # Excited State   1:   Singlet-BU     5.3351 eV  232.39 nm  f=0.1695
            # (unrestricted calc) (first excited state is 2!)
            # Excited State   2:   ?Spin  -A      0.1222 eV 10148.75 nm  f=0.0000
            # (Gaussian 09 ZINDO)
            # Excited State   1:      Singlet-?Sym    2.5938 eV  478.01 nm  f=0.0000  <S**2>=0.000
            p = re.compile(r":(?P<sym>.*?)(?P<energy>-?\d*\.\d*) eV")
            groups = p.search(line).groups()
            self.append_attribute(
                "etenergies", utils.convertor(utils.float(groups[1]), "eV", "hartree")
            )
            self.append_attribute("etoscs", utils.float(line.split("f=")[-1].split()[0]))
            # Fix Gaussian's weird capitalisation.
            mult, symm = groups[0].strip().split("-")
            self.append_attribute("etsyms", f"{mult}-{self.normalisesym(symm)}")

            line = next(inputfile)

            p = re.compile(r"(\d+)")
            CIScontrib = []
            while line.find(" ->") >= 0:  # This is a contribution to the transition
                parts = line.split("->")
                self.logger.debug(parts)
                # Has to deal with lines like:
                #       32 -> 38         0.04990
                #      35A -> 45A        0.01921
                frommoindex = 0  # For restricted or alpha unrestricted
                fromMO = parts[0].strip()
                if fromMO[-1] == "B":
                    frommoindex = 1  # For beta unrestricted
                fromMO = (
                    int(p.match(fromMO).group()) - 1
                )  # subtract 1 so that it is an index into moenergies

                t = parts[1].split()
                tomoindex = 0
                toMO = t[0]
                if toMO[-1] == "B":
                    tomoindex = 1
                toMO = (
                    int(p.match(toMO).group()) - 1
                )  # subtract 1 so that it is an index into moenergies

                percent = utils.float(t[1])
                # For restricted calculations, the percentage will be corrected
                # after parsing (see after_parsing() above).
                CIScontrib.append([(fromMO, frommoindex), (toMO, tomoindex), percent])
                line = next(inputfile)
            self.append_attribute("etsecs", CIScontrib)

            # Skip over 'de-excitation' contributions (these are typically hidden but can be revealed
            # by iop(9/40=2)).
            while line.find(" <-") >= 0:
                # These are not processed atm.
                line = next(inputfile)

            # Check if this state is our 'state of interest' (for optimisations etc).
            if "This state for optimization and/or second-order correction" in line:
                # Index to the current excited state.
                self.metadata["opt_state"] = len(self.etenergies) - 1

        # Electronic transition transition-dipole data
        #
        # Ground to excited state transition electric dipole moments (Au):
        #       state          X           Y           Z        Dip. S.      Osc.
        #         1         0.0008     -0.0963      0.0000      0.0093      0.0005
        #         2         0.0553      0.0002      0.0000      0.0031      0.0002
        #         3         0.0019      2.3193      0.0000      5.3790      0.4456
        # Ground to excited state transition velocity dipole moments (Au):
        #       state          X           Y           Z        Dip. S.      Osc.
        #         1        -0.0001      0.0032      0.0000      0.0000      0.0001
        #         2        -0.0081      0.0000      0.0000      0.0001      0.0005
        #         3        -0.0002     -0.2692      0.0000      0.0725      0.3887
        # Ground to excited state transition magnetic dipole moments (Au):
        #       state          X           Y           Z
        #         1         0.0000      0.0000     -0.0003
        #         2         0.0000      0.0000      0.0000
        #         3         0.0000      0.0000      0.0035
        #
        # NOTE: In Gaussian03, there were some inconsitancies in the use of
        # small / capital letters: e.g.
        # Ground to excited state Transition electric dipole moments (Au)
        # Ground to excited state transition velocity dipole Moments (Au)
        # so to look for a match, we will lower() everything.
        #
        # EOM-CCSD looks similar, but has more data.
        #  ==============================================
        #
        #          EOM-CCSD transition properties
        #
        #  ==============================================
        #  Ground to excited state transition electric dipole moments (Au):
        #        state          X           Y           Z        Dip. S.      Osc.
        #          1         0.1021     -0.0000     -0.0000      0.0107      0.0030
        #          2         0.0000      0.0000      0.0000      0.0000      0.0000
        #  Excited to ground state transition electric dipole moments (Au):
        #        state          X           Y           Z        Dip. S.      Osc.
        #          1         0.1046     -0.0000     -0.0000      0.0107      0.0030
        #          2         0.0000      0.0000      0.0000      0.0000      0.0000
        #  Ground to excited state transition velocity dipole moments (Au):
        #        state          X           Y           Z        Dip. S.      Osc.
        #          1        -0.1351     -0.0000     -0.0000      0.0186      0.0296
        #          2         0.0000      0.0000      0.0000      0.0000      0.0000
        #  Excited to ground state transition velocity dipole moments (Au):
        #        state          X           Y           Z        Dip. S.      Osc.
        #          1        -0.1378     -0.0000     -0.0000      0.0186      0.0296
        #          2        -0.0000      0.0000      0.0000      0.0000      0.0000
        #  Ground to excited state transition magnetic dipole moments (Au):
        #        state          X           Y           Z
        #          1         0.0000      0.5361      0.0000
        #          2        -0.0000      0.0000      0.6910
        #  Excited to ground state transition magnetic dipole moments (Au):
        #        state          X           Y           Z
        #          1         0.0000      0.5473      0.0000
        #          2        -0.0000     -0.0000      0.7057

        if line[1:51].lower() == "ground to excited state transition electric dipole":
            # In EOM-CCSD we have multiple levels of theory, so we always want to reset.
            self.etdips = []
            self.etveldips = []
            self.etmagdips = []

            # now loop over lines reading eteltrdips until we find eteltrdipvel
            line = next(inputfile)  # state          X ...
            line = next(inputfile)  # 1        -0.0001 ...

            # Older versions have fewer fields.
            while len(line.split()) in [5, 6]:
                self.etdips.append(list(map(float, line.split()[1:4])))
                line = next(inputfile)

        if line[1:51].lower() == "ground to excited state transition velocity dipole":
            # now loop over lines reading etveleltrdips until we find
            # etmagtrdip
            line = next(inputfile)  # state          X ...
            line = next(inputfile)  # 1        -0.0001 ...

            # Older versions have fewer fields.
            while len(line.split()) in [5, 6]:
                self.etveldips.append(list(map(float, line.split()[1:4])))
                line = next(inputfile)

        if line[1:51].lower() == "ground to excited state transition magnetic dipole":
            # now loop over lines while the line starts with at least 3 spaces
            line = next(inputfile)  # state          X ...
            line = next(inputfile)  # 1        -0.0001 ...
            while len(line.split()) == 4:
                self.etmagdips.append(list(map(float, line.split()[1:4])))
                line = next(inputfile)

        # Circular dichroism data (different for G03 vs G09)
        #
        # G03
        #
        # ## <0|r|b> * <b|rxdel|0>  (Au), Rotatory Strengths (R) in
        # ## cgs (10**-40 erg-esu-cm/Gauss)
        # ##       state          X           Y           Z     R(length)
        # ##         1         0.0006      0.0096     -0.0082     -0.4568
        # ##         2         0.0251     -0.0025      0.0002     -5.3846
        # ##         3         0.0168      0.4204     -0.3707    -15.6580
        # ##         4         0.0721      0.9196     -0.9775     -3.3553
        #
        # G09
        #
        # ## 1/2[<0|r|b>*<b|rxdel|0> + (<0|rxdel|b>*<b|r|0>)*]
        # ## Rotatory Strengths (R) in cgs (10**-40 erg-esu-cm/Gauss)
        # ##       state          XX          YY          ZZ     R(length)     R(au)
        # ##         1        -0.3893     -6.7546      5.7736     -0.4568     -0.0010
        # ##         2       -17.7437      1.7335     -0.1435     -5.3845     -0.0114
        # ##         3       -11.8655   -297.2604    262.1519    -15.6580     -0.0332
        if (
            line[1:52] == "<0|r|b> * <b|rxdel|0>  (Au), Rotatory Strengths (R)"
            or line[1:50] == "1/2[<0|r|b>*<b|rxdel|0> + (<0|rxdel|b>*<b|r|0>)*]"
        ):
            self.etrotats = []

            self.skip_lines(inputfile, ["units"])

            headers = next(inputfile)
            Ncolms = len(headers.split())
            line = next(inputfile)
            parts = line.strip().split()
            while len(parts) == Ncolms:
                try:
                    R = utils.float(parts[4])
                except ValueError:
                    # nan or -nan if there is no first excited state
                    # (for unrestricted calculations)
                    pass
                else:
                    self.etrotats.append(R)
                line = next(inputfile)
                temp = line.strip().split()
                parts = line.strip().split()
            self.etrotats = numpy.array(self.etrotats, "d")

        # Number of basis sets functions.
        # Has to deal with lines like:
        #  NBasis =   434 NAE=    97 NBE=    97 NFC=    34 NFV=     0
        # and...
        #  NBasis = 148  MinDer = 0  MaxDer = 0
        # Although the former is in every file, it doesn't occur before
        #   the overlap matrix is printed.
        if line[1:7] == "NBasis" or line[4:10] == "NBasis":
            # For counterpoise fragment, skip these lines.
            if self.counterpoise != 0:
                return

            # For ONIOM calcs, ignore this section in order to bypass assertion failure.
            if self.oniom:
                return

            # If nbasis was already parsed, check if it changed. If it did, issue a warning.
            # In the future, we will probably want to have nbasis, as well as nmo below,
            # as a list so that we don't need to pick one value when it changes.
            nbasis = int(line.split("=")[1].split()[0])
            if hasattr(self, "nbasis"):
                try:
                    assert nbasis == self.nbasis
                except AssertionError:
                    self.logger.warning(
                        "Number of basis functions (nbasis) has changed from %s to %s",
                        self.nbasis,
                        nbasis,
                    )
            self.nbasis = nbasis

        # Number of linearly-independent basis sets.
        if line[1:7] == "NBsUse":
            # For counterpoise fragment, skip these lines.
            if self.counterpoise != 0:
                return

            # For ONIOM calcs, ignore this section in order to bypass assertion failure.
            if self.oniom:
                return

            nmo = int(line.split("=")[1].split()[0])
            self.set_attribute("nmo", nmo)

        # For AM1 calculations, set nbasis by a second method,
        #   as nmo may not always be explicitly stated.
        if line[7:22] == "basis functions, ":
            nbasis = int(line.split()[0])
            self.set_attribute("nbasis", nbasis)

        # Molecular orbital overlap matrix.
        # Has to deal with lines such as:
        #   *** Overlap ***
        #   ****** Overlap ******
        # Note that Gaussian sometimes drops basis functions,
        #  causing the overlap matrix as parsed below to not be
        #  symmetric (which is a problem for population analyses, etc.)
        if line[1:4] == "***" and (line[5:12] == "Overlap" or line[8:15] == "Overlap"):
            # Ensure that this is the main calc and not a fragment
            if self.counterpoise != 0:
                return

            self.aooverlaps = numpy.zeros((self.nbasis, self.nbasis), "d")
            # Overlap integrals for basis fn#1 are in aooverlaps[0]
            base = 0
            colmNames = next(inputfile)
            while base < self.nbasis:
                self.updateprogress(inputfile, "Overlap", self.fupdate)

                for i in range(self.nbasis - base):  # Fewer lines this time
                    line = next(inputfile)
                    parts = line.split()
                    for j in range(len(parts) - 1):  # Some lines are longer than others
                        k = float(parts[j + 1].replace("D", "E"))
                        self.aooverlaps[base + j, i + base] = k
                        self.aooverlaps[i + base, base + j] = k
                base += 5
                colmNames = next(inputfile)
            self.aooverlaps = numpy.array(self.aooverlaps, "d")

        # Molecular orbital coefficients (mocoeffs).
        # Essentially only produced for SCF calculations.
        # This is also the place where aonames and atombasis are parsed.
        if (
            line[5:35] == "Molecular Orbital Coefficients"
            or line[5:41] == "Alpha Molecular Orbital Coefficients"
            or line[5:40] == "Beta Molecular Orbital Coefficients"
        ):
            # If counterpoise fragment, return without parsing orbital info
            if self.counterpoise != 0:
                return
            # Skip this for ONIOM calcs
            if self.oniom:
                return

            if line[5:40] == "Beta Molecular Orbital Coefficients":
                beta = True
                if self.popregular:
                    return
                    # This was continue before refactoring the parsers.
                    # continue # Not going to extract mocoeffs
                # Need to add an extra array to self.mocoeffs
                self.mocoeffs.append(numpy.zeros((self.nmo, self.nbasis), "d"))
            else:
                beta = False
                self.aonames = []
                self.atombasis = []
                mocoeffs = [numpy.zeros((self.nmo, self.nbasis), "d")]

            base = 0
            self.popregular = False
            for base in range(0, self.nmo, 5):
                self.updateprogress(inputfile, "Coefficients", self.fupdate)

                colmNames = next(inputfile)

                if any(
                    name in colmNames
                    for name in [
                        "Density Matrix:",
                        "DENSITY MATRIX.",
                        "Beta Molecular Orbital Coefficients",
                    ]
                ):
                    # Reached end of mocoeff section early, this implies pop was not full.
                    self.popregular = True
                    # We can stop processing.
                    break

                if not colmNames.split():
                    self.logger.warning("Molecular coefficients header found but no coefficients.")
                    break

                symmetries = next(inputfile)  # noqa: F841
                eigenvalues = next(inputfile)
                for i in range(self.nbasis):
                    line = next(inputfile)
                    if i == 0:
                        # Find location of the start of the basis function name
                        start_of_basis_fn_name = line.find(line.split()[3]) - 1
                    if base == 0 and not beta:  # Just do this the first time 'round
                        parts = line[:start_of_basis_fn_name].split()
                        if len(parts) > 1:  # New atom
                            if i > 0:
                                self.atombasis.append(atombasis)  # noqa: F821
                            atombasis = []
                            atomname = f"{parts[2]}{parts[1]}"
                        orbital = line[start_of_basis_fn_name:20].strip()
                        self.aonames.append(f"{atomname}_{orbital}")
                        atombasis.append(i)

                    part = line[21:].replace("D", "E").rstrip()
                    temp = []
                    for j in range(0, len(part), 10):
                        temp.append(float(part[j : j + 10]))
                    if beta:
                        self.mocoeffs[1][base : base + len(part) // 10, i] = temp
                    else:
                        mocoeffs[0][base : base + len(part) // 10, i] = temp

                if base == 0 and not beta:  # Do the last update of atombasis
                    self.atombasis.append(atombasis)
            if not self.popregular and not beta:
                self.mocoeffs = mocoeffs

        # Natural orbital coefficients (nocoeffs) and occupation numbers (nooccnos),
        # which are respectively define the eigenvectors and eigenvalues of the
        # diagonalized one-electron density matrix. These orbitals are formed after
        # configuration interaction (CI) calculations, but not only. Similarly to mocoeffs,
        # we can parse and check aonames and atombasis here.
        #
        #     Natural Orbital Coefficients:
        #                           1         2         3         4         5
        #     Eigenvalues --     2.01580   2.00363   2.00000   2.00000   1.00000
        #   1 1   O  1S          0.00000  -0.15731  -0.28062   0.97330   0.00000
        #   2        2S          0.00000   0.75440   0.57746   0.07245   0.00000
        # ...
        #
        def natural_orbital_single_spin_parsing(inputfile, updateprogress_title):
            coeffs = numpy.zeros((self.nmo, self.nbasis), "d")
            occnos = []
            aonames = []
            atombasis = []
            for base in range(0, self.nmo, 5):
                self.updateprogress(inputfile, updateprogress_title, self.fupdate)
                self.skip_line(inputfile, "column numbers")
                eigenvalues = next(inputfile)
                occnos.extend(map(float, eigenvalues.split()[2:]))
                for i in range(self.nbasis):
                    line = next(inputfile)
                    # Just do this the first time 'round.
                    if base == 0:
                        parts = line[:12].split()
                        # New atom.
                        if len(parts) > 1:
                            if i > 0:
                                atombasis.append(basisonatom)  # noqa: F821
                            basisonatom = []
                            atomname = f"{parts[2]}{parts[1]}"
                        orbital = line[11:20].strip()
                        aonames.append(f"{atomname}_{orbital}")
                        basisonatom.append(i)
                    part = line[21:].replace("D", "E").rstrip()
                    temp = []
                    for j in range(0, len(part), 10):
                        temp.append(float(part[j : j + 10]))
                    coeffs[base : base + len(part) // 10, i] = temp
                # Do the last update of atombasis.
                if base == 0:
                    atombasis.append(basisonatom)
            return occnos, coeffs, aonames, atombasis

        if line[5:33] == "Natural Orbital Coefficients":
            updateprogress_title = "Natural orbitals"
            nooccnos, nocoeffs, aonames, atombasis = natural_orbital_single_spin_parsing(
                inputfile, updateprogress_title
            )
            self.set_attribute("nocoeffs", nocoeffs)
            self.set_attribute("nooccnos", nooccnos)
            self.set_attribute("atombasis", atombasis)
            self.set_attribute("aonames", aonames)

        # Natural spin orbital coefficients (nsocoeffs) and occupation numbers (nsooccnos)
        # Parsed attributes are similar to the natural orbitals above except
        # the natural spin orbitals and occupation numbers are the eigenvalues
        # and eigenvectors of the one particles spin density matrices
        #     Alpha Natural Orbital Coefficients:
        #                           1         2         3         4         5
        #     Eigenvalues --     1.00000   1.00000   0.99615   0.99320   0.99107
        #   1 1   O  1S          0.70425   0.70600  -0.16844  -0.14996  -0.00000
        #   2        2S          0.01499   0.01209   0.36089   0.34940  -0.00000
        # ...
        #     Beta Natural Orbital Coefficients:
        #                           1         2         3         4         5
        #     Eigenvalues --     1.00000   1.00000   0.99429   0.98790   0.98506
        #   1 1   O  1S          0.70822   0.70798  -0.15316  -0.13458   0.00465
        #   2        2S          0.00521   0.00532   0.33837   0.33189  -0.01301
        #   3        3S         -0.02542  -0.00841   0.28649   0.53224   0.18902
        # ...

        if line[5:39] == "Alpha Natural Orbital Coefficients":
            updateprogress_title = "Natural Spin orbitals (alpha)"
            nsooccnos, nsocoeffs, aonames, atombasis = natural_orbital_single_spin_parsing(
                inputfile, updateprogress_title
            )
            if self.unified_no_nso:
                self.append_attribute("nocoeffs", nsocoeffs)
                self.append_attribute("nooccnos", nsooccnos)
            else:
                self.append_attribute("nsocoeffs", nsocoeffs)
                self.append_attribute("nsooccnos", nsooccnos)
            self.set_attribute("atombasis", atombasis)
            self.set_attribute("aonames", aonames)
        if line[5:38] == "Beta Natural Orbital Coefficients":
            updateprogress_title = "Natural Spin orbitals (beta)"
            nsooccnos, nsocoeffs, aonames, atombasis = natural_orbital_single_spin_parsing(
                inputfile, updateprogress_title
            )
            if self.unified_no_nso:
                self.append_attribute("nocoeffs", nsocoeffs)
                self.append_attribute("nooccnos", nsooccnos)
            else:
                self.append_attribute("nsocoeffs", nsocoeffs)
                self.append_attribute("nsooccnos", nsooccnos)
            self.set_attribute("atombasis", atombasis)
            self.set_attribute("aonames", aonames)

        # For FREQ=Anharm, extract anharmonicity constants
        if line[1:40] == "X matrix of Anharmonic Constants (cm-1)":
            Nvibs = len(self.vibfreqs)
            self.vibanharms = numpy.zeros((Nvibs, Nvibs), "d")

            base = 0
            colmNames = next(inputfile)
            while base < Nvibs:
                for i in range(Nvibs - base):  # Fewer lines this time
                    line = next(inputfile)
                    parts = line.split()
                    for j in range(len(parts) - 1):  # Some lines are longer than others
                        k = float(parts[j + 1].replace("D", "E"))
                        self.vibanharms[base + j, i + base] = k
                        self.vibanharms[i + base, base + j] = k
                base += 5
                colmNames = next(inputfile)

        # Pseudopotential charges.
        if line.find("Pseudopotential Parameters") > -1:
            self.skip_lines(inputfile, ["e", "label1", "label2", "e"])

            line = next(inputfile)
            if line.find("Centers:") < 0:
                return
                # This was continue before parser refactoring.
                # continue

            # Needs to handle code like the following:
            #
            #  Center     Atomic      Valence      Angular      Power                                                       Coordinates
            #  Number     Number     Electrons     Momentum     of R      Exponent        Coefficient                X           Y           Z
            # ===================================================================================================================================
            # Centers:   1
            # Centers:  16
            # Centers:  21 24
            # Centers:  99100101102
            #    1         44           16                                                                      -4.012684 -0.696698  0.006750
            #                                      F and up
            #                                                     0      554.3796303       -0.05152700
            centers = []
            while line.find("Centers:") >= 0:
                temp = line[10:]
                for i in range(0, len(temp) - 3, 3):
                    centers.append(int(temp[i : i + 3]))
                line = next(inputfile)
            centers.sort()  # Not always in increasing order

            self.coreelectrons = numpy.zeros(self.natom, "i")

            for center in centers:
                front = line[:10].strip()
                while not (front and int(front) == center):
                    line = next(inputfile)
                    front = line[:10].strip()
                info = line.split()
                self.coreelectrons[center - 1] = int(info[1]) - int(info[2])
                line = next(inputfile)

        # This will be printed for counterpoise calcualtions only.
        # To prevent crashing, we need to know which fragment is being considered.
        # Other information is also printed in lines that start like this.
        if line[1:14] == "Counterpoise:":
            if line[42:50] == "fragment":
                self.counterpoise = int(line[51:54])

        # This will be printed only during ONIOM calcs; use it to set a flag
        # that will allow assertion failures to be bypassed in the code.
        if line[1:7] == "ONIOM:":
            self.oniom = True

        # This will be printed only during BOMD calcs;
        if line.startswith(" INPUT DATA FOR L118"):
            self.BOMD = True

        # Atomic charges are straightforward to parse, although the header
        # has changed over time somewhat.
        #
        # Mulliken charges:
        #                1
        #     1  C   -0.004513
        #     2  C   -0.077156
        # ...
        # Sum of Mulliken charges =   0.00000
        # Mulliken charges with hydrogens summed into heavy atoms:
        #               1
        #     1  C   -0.004513
        #     2  C    0.002063
        # ...
        #
        # Spins may be included in the same section.
        #
        # Mulliken charges and spin densities:
        #               1          2
        #     1  O   -0.596188   0.019133
        #     2  C    0.320624   0.000869
        #
        # APT and Lowdin charges are also displayed in this way.
        def extract_charges_spins(line: str, prop: str) -> None:
            """Extracts atomic charges and spin de strnsities into
               self.atomcharges and self.atomspins dictionaries.

            Inputs:
                line - line header marking the beginning of a particular set of charges or spins.
                prop - property type to be extracted as a string (e.g. Mulliken, Lowdin, APT).
            """
            has_spin = "spin" in line.lower()
            has_charges = "charges" in line.lower()
            if has_charges and not hasattr(self, "atomcharges"):
                self.atomcharges = {}
            if has_spin and not hasattr(self, "atomspins"):
                self.atomspins = {}
            _ = next(inputfile)
            charges = []
            spins = []
            # It's possible that this could match on a header other than how
            # the match happen in the code before calling this function, but
            # the effect should be the same.
            swap_indices = any(
                header in line.lower() for header in self.atomcharges_atomspins_headers_swap
            )
            is_sum = "summed" in line
            # Iterate over each line and append values to a list
            # based on whether they are charges or spins.
            if is_sum:
                for i in self.atomnos:
                    # currently bug exists where files with translation vectors report
                    # an extra atom with atomnumber -2 in self.atomnos, so must ignore
                    # this by passing whenever i in self.atomnos == -2.
                    if i == -2:
                        pass
                    # For lists of summed charges or spins, a value
                    # of 0 is added if the atom is a hydrogen.
                    elif i == 1:
                        if has_charges:
                            charges.append(float(0))
                            if has_spin:
                                spins.append(float(0))
                        elif has_spin:
                            spins.append(float(0))
                    else:
                        nline = next(inputfile)
                        # Some older versions of Gaussian already include
                        # hydrogens with value 0 for summed charges or
                        # spins, so these should be ignored.
                        while nline.split()[1] == "H":
                            nline = next(inputfile)
                        _append_charges_and_spins(
                            nline, has_charges, charges, has_spin, spins, swap_indices
                        )
            else:
                for i in self.atomnos:
                    # Ignore translation vectors.
                    if i == -2:
                        pass
                    else:
                        nline = next(inputfile)
                        _append_charges_and_spins(
                            nline, has_charges, charges, has_spin, spins, swap_indices
                        )
            # When the charge type is not given explicitly we
            # must find it from the bottom line, which always
            # has the format: "Sum of Mulliken charges=   0.00000"
            # so we can extract the type by splitting each
            # line until we get a valid charge type.
            while prop.lower() not in self.atomprops_no_atomic:
                nline = next(inputfile)
                prop = nline.split()[2].lower()
            # Input extracted values into self.atomcharges.
            if has_charges:
                if is_sum:
                    self.atomcharges[f"{prop}_sum"] = charges
                else:
                    self.atomcharges[f"{prop}"] = charges
            if has_spin:
                if is_sum:
                    self.atomspins[f"{prop}_sum"] = spins
                else:
                    self.atomspins[f"{prop}"] = spins

            return None

        if hasattr(self, "atomnos"):
            # Combine props and headers to find lines heading lists
            # of atom charges or spins.
            for prop in self.atomprops:
                for header in self.atomcharges_atomspins_headers:
                    if f"{prop}{header}".lower() in line.lower():
                        # When we use "atomic" as the property, only
                        # extract if the charge type isn't given explicity.
                        # This prevents us from reading some lines twice.
                        # e.g. "Mulliken atomic charges:" is caught by
                        # "mulliken atomic charges:" and " atomic charges:"
                        if prop == "atomic":
                            if all(
                                atomprop not in line.lower()
                                for atomprop in self.atomprops_no_atomic
                            ):
                                extract_charges_spins(line, prop)
                        else:
                            extract_charges_spins(line, prop)

        if "N A T U R A L   A T O M I C   O R B I T A L" in line:
            self.nbo_spin_section = False
        elif "Alpha spin orbitals" in line:
            self.nbo_spin_section = True
        elif "Beta spin orbitals" in line:
            self.nbo_spin_section = True

        if line.strip() == "Natural Population":
            if not hasattr(self, "atomcharges"):
                self.atomcharges = {}
            if not self.nbo_spin_section:
                line1 = next(inputfile)
                line2 = next(inputfile)
                if line1.split()[0] == "Natural" and line2.split()[2] == "Charge":
                    self.skip_line(inputfile, "d")
                    charges = []
                    for _ in range(self.natom):
                        nline = next(inputfile)
                        charges.append(float(nline.split()[2]))
                    self.atomcharges["natural"] = charges

        # Combined Hirshfeld/CM5 is different enough that we don't try and
        # reuse extract_charges_spins, at least for now.
        if "Hirshfeld charges, spin densities, dipoles, and CM5 charges" in line:
            if not hasattr(self, "atomcharges"):
                self.atomcharges = {}
            has_spins = len(self.homos) == 2
            if has_spins and not hasattr(self, "atomspins"):
                self.atomspins = {}
            self.skip_line(inputfile, "Q-H        S-H        Dx")
            line = next(inputfile)
            atomcharges_hirshfeld = []
            atomcharges_cm5 = []
            atomspins_hirshfeld = []
            while "Tot" not in line:
                tokens = line.split()
                atomcharges_hirshfeld.append(float(tokens[2]))
                atomcharges_cm5.append(float(tokens[7]))
                if has_spins:
                    atomspins_hirshfeld.append(float(tokens[3]))
                line = next(inputfile)
            self.atomcharges["hirshfeld"] = atomcharges_hirshfeld
            self.atomcharges["cm5"] = atomcharges_cm5
            if has_spins:
                self.atomspins["hirshfeld"] = atomspins_hirshfeld
            self.skip_lines(
                inputfile,
                ["Hirshfeld charges with hydrogens summed into heavy atoms:", "Q-H       Q-CM5"],
            )
            line = next(inputfile)
            atomcharges_hirshfeld_summed = []
            atomcharges_cm5_summed = []
            for i in self.atomnos:
                if i == 1:
                    atomcharges_hirshfeld_summed.append(0.0)
                    atomcharges_cm5_summed.append(0.0)
                else:
                    tokens = line.split()
                    atomcharges_hirshfeld_summed.append(float(tokens[2]))
                    atomcharges_cm5_summed.append(float(tokens[3]))
                    line = next(inputfile)
            self.atomcharges["hirshfeld_sum"] = atomcharges_hirshfeld_summed
            self.atomcharges["cm5_sum"] = atomcharges_cm5_summed

        # Determine information about ESP-derived charges.
        if self.link == 602:
            # Some are not prefixed or postfixed for consistency with the
            # Q-Chem parser.
            radii_to_esp_model = {
                "Breneman": "chelpg",
                "Chirlian-Francl": "chelp",
                "Francl": "chelp",
                # Ideally would be "esp_hly".
                "Hu-Lu-Yang": "esp",
                # Ideally would be "esp_mk".
                "Merz-Kollman": "esp",
                # Ideally would be "esp_mkuff".
                "UFF": "esp",
            }
            esp_iop = self.metadata["routes"][-1][6][20]
            if line.startswith(" Read replacement radii"):
                # This section still generally depends on a parent set of
                # radii being defined.  This parent set determines the
                # atomcharges key.
                re_radii_line = re.compile(
                    r"Rad\((?P<atnum>\d+)\)\s+=\s+(?P<radius>\d\.\d{12}D\+\d{2})"
                )
                line = next(inputfile)
                while line.startswith(" Rad("):
                    mtch = re_radii_line.search(line)
                    assert mtch is not None
                    _atnum = int(mtch.group("atnum"))
                    _radius = utils.float(mtch.group("radius"))
                    line = next(inputfile)
            elif "radii" in line:
                self.esp_model = radii_to_esp_model[line.split()[0]]
            elif line.startswith(" Generate Potential Derived Charges using the"):
                self.esp_model = radii_to_esp_model[line.split()[6]]
            elif line.strip() == "Electrostatic Properties Using The SCF Density":
                key = self.esp_model
                _ = self.skip_lines(inputfile, ["b", "s", "b"])
                # atomic positions
                for _ in range(self.natom):
                    line = next(inputfile)
                line = next(inputfile)
                assert "points will be used for fitting atomic charges" in line
                _npoints = int(line.split()[0])
                line = next(inputfile)
                assert "Fitting point charges" in line
                fitting_point_dipoles = False
                if "Fitting point charges and point dipoles" in line:
                    fitting_point_dipoles = True
                line = next(inputfile)
                if line.strip() == "The dipole moment will be constrained to the correct value":
                    # This doesn't make sense for ChElP(G), so avoid setting
                    # it there.
                    if key == "esp":
                        key = "resp"
                    if fitting_point_dipoles:
                        _key = "esp_mk_dipole_atom"
                    line = next(inputfile)
                assert line.startswith(" Charges from ESP fit")
                line = next(inputfile)
                assert line.startswith(" Charge=") or line.strip() in (
                    "ESP charges:",
                    "ESP charges         Point Dipoles (au):",
                )
                line = next(inputfile)
                line = next(inputfile)
                charges = []
                for _ in range(self.natom):
                    charges.append(float(line.split()[2]))
                    line = next(inputfile)
                # Gaussian atomic densities instead of those from the HLY scheme
                if esp_iop == 15:
                    _key = "esp_hlygat"
                self.atomcharges[key] = charges

        # Extract Thermochemistry
        # Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
        # Zero-point correction=                           0.342233 (Hartree/
        # Thermal correction to Energy=                    0.
        # Thermal correction to Enthalpy=                  0.
        # Thermal correction to Gibbs Free Energy=         0.302940
        # Sum of electronic and zero-point Energies=           -563.649744
        # Sum of electronic and thermal Energies=              -563.636699
        # Sum of electronic and thermal Enthalpies=            -563.635755
        # Sum of electronic and thermal Free Energies=         -563.689037
        if "Zero-point correction" in line:
            self.set_attribute("zpve", float(line.split()[2]))
        if "Sum of electronic and thermal Enthalpies" in line:
            self.set_attribute("enthalpy", float(line.split()[6]))
        if "Sum of electronic and thermal Free Energies=" in line:
            self.set_attribute("freeenergy", float(line.split()[7]))
        if line[1:13] == "Temperature ":
            self.set_attribute("temperature", float(line.split()[1]))
            self.set_attribute("pressure", float(line.split()[4]))

        # Static polarizability (from `polar`), lower triangular
        # matrix.
        if line[1:26] == "SCF Polarizability for W=":
            self.hp_polarizabilities = True
            polarizability = numpy.zeros(shape=(3, 3))
            self.skip_line(inputfile, "directions")
            for i in range(3):
                line = next(inputfile)
                polarizability[i, : i + 1] = [utils.float(x) for x in line.split()[1:]]
            self.append_attribute(
                "polarizabilities", utils.symmetrize(polarizability, use_triangle="lower")
            )

        # Static polarizability (from `freq`), lower triangular matrix.
        if line[1:16] == "Polarizability=":
            self.hp_polarizabilities = True
            polarizability = numpy.zeros(shape=(3, 3))
            polarizability_list = []
            polarizability_list.extend([line[16:31], line[31:46], line[46:61]])
            line = next(inputfile)
            polarizability_list.extend([line[16:31], line[31:46], line[46:61]])
            indices = numpy.tril_indices(3)
            polarizability[indices] = [utils.float(x) for x in polarizability_list]
            self.append_attribute(
                "polarizabilities", utils.symmetrize(polarizability, use_triangle="lower")
            )

        # Static polarizability, compressed into a single line from
        # terse printing.
        # Order is XX, YX, YY, ZX, ZY, ZZ (lower triangle).
        if line[2:23] == "Exact polarizability:":
            if not self.hp_polarizabilities:
                polarizability = numpy.zeros(shape=(3, 3))
                indices = numpy.tril_indices(3)
                try:
                    # G16 C01 changes polarizability printing
                    # Sample:
                    #       Exact polarizability:      68.238       6.777     143.018       0.000       0.000      11.343
                    polarizability[indices] = [
                        utils.float(x)
                        for x in [
                            line[23:35],
                            line[35:47],
                            line[47:59],
                            line[59:71],
                            line[71:83],
                            line[83:95],
                        ]
                    ]
                except:  # noqa: E722
                    # G16A03 and older
                    # Sample:
                    #       Exact polarizability:  68.238  -6.777 143.018   0.000   0.000  11.343
                    polarizability[indices] = [
                        utils.float(x)
                        for x in [
                            line[23:31],
                            line[31:39],
                            line[39:47],
                            line[47:55],
                            line[55:63],
                            line[63:71],
                        ]
                    ]
                self.append_attribute(
                    "polarizabilities", utils.symmetrize(polarizability, use_triangle="lower")
                )

        # IRC Computation convergence checks.
        #
        # -------- Sample extract for IRC step --------
        #
        # IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC
        # ------------------------------------------------------------------------
        # INPUT DATA FOR L123
        # ------------------------------------------------------------------------
        # GENERAL PARAMETERS:
        # Follow reaction path in both directions.
        # Maximum points per path      = 200
        # Step size                    =   0.100 bohr
        # Integration scheme           = HPC
        #    Redo corrector integration= Yes
        #    DWI Weight Power          =  2
        #    DWI will use Hessian update vectors when possible.
        #    Max correction cycles     =  50
        # Initial Hessian              = CalcFC
        # Hessian evaluation           = Analytic every   5 predictor steps
        #                              = Analytic every   5 corrector steps
        # Hessian updating method      = Bofill
        # ------------------------------------------------------------------------
        # IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC
        #
        # -------- Sample extract for converged step --------
        # IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC
        # Pt  1 Step number   1 out of a maximum of  50
        # Modified Bulirsch-Stoer Extrapolation Cycles:
        #    EPS =    0.000010000000000
        # Maximum DWI energy   std dev =  0.000000595 at pt     1
        # Maximum DWI gradient std dev =  0.135684493 at pt     2
        # CORRECTOR INTEGRATION CONVERGENCE:
        #   Recorrection delta-x convergence threshold:    0.010000
        #   Delta-x Convergence Met
        # Point Number:   1          Path Number:   1
        #   CHANGE IN THE REACTION COORDINATE =    0.16730
        #   NET REACTION COORDINATE UP TO THIS POINT =    0.16730
        #  # OF POINTS ALONG THE PATH =   1
        #  # OF STEPS =   1
        #
        # Calculating another point on the path.
        # IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC
        #
        # -------- Sample extract for unconverged intermediate step --------
        # IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC
        # Error in corrector energy =          -0.0000457166
        # Magnitude of corrector gradient =     0.0140183779
        # Magnitude of analytic gradient =      0.0144021969
        # Magnitude of difference =             0.0078709968
        # Angle between gradients (degrees)=   32.1199
        # Pt 40 Step number   2 out of a maximum of  20
        # Modified Bulirsch-Stoer Extrapolation Cycles:
        #    EPS =    0.000010000000000
        # Maximum DWI energy   std dev =  0.000007300 at pt    31
        # Maximum DWI gradient std dev =  0.085197906 at pt    59
        # CORRECTOR INTEGRATION CONVERGENCE:
        #   Recorrection delta-x convergence threshold:    0.010000
        #   Delta-x Convergence NOT Met
        # IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC
        if line[1:20] == "INPUT DATA FOR L123":  # First IRC step
            self.append_attribute("optstatus", data.ccData.OPT_NEW)
        if line[3:22] == "Delta-x Convergence":
            if line[23:30] == "NOT Met":
                self.optstatus[-1] += data.ccData.OPT_UNCONVERGED
            elif line[23:26] == "Met":
                self.optstatus[-1] += data.ccData.OPT_DONE
                self.append_attribute("optdone", len(self.optstatus) - 1)

        # Save num CPUs incase we have an old version of Gaussian which doesn't print wall times.
        if "Will use up to" in line:
            self.num_cpu = int(line.split()[4])

        # Extract total elapsed (wall) and CPU job times
        if line[:14] == " Elapsed time:" or line[:14] == " Job cpu time:":
            # create empty list for the times to be stored in
            if line[:14] == " Elapsed time:" and "wall_time" not in self.metadata:
                self.metadata["wall_time"] = []
            if line[:14] == " Job cpu time:" and "cpu_time" not in self.metadata:
                self.metadata["cpu_time"] = []
            # the line format is " Elapsed time:       0 days  0 hours  0 minutes 47.5 seconds." at the end of each job ran.
            # the line format is " Job cpu time:       0 days  0 hours  8 minutes 45.7 seconds." at the end of each job ran.
            try:
                n = 2
                key = "wall_time"
                # if parsing a cpu time change key and shift n
                if line[:14] == " Job cpu time:":
                    n += 1
                    key = "cpu_time"
                # split the line by white space
                split_line = line.split()
                # cast the time elements as floats for use in timedelta data structure
                time = datetime.timedelta(
                    days=float(split_line[n]),
                    hours=float(split_line[n + 2]),
                    minutes=float(split_line[n + 4]),
                    seconds=float(split_line[n + 6]),
                )
                self.metadata[key].append(time)
            except:  # noqa: E722
                pass

        # Extract Rotational Constants
        # Example:
        # Rotational constants (GHZ):           3.13081     1.24272     0.88960
        # OR for linear molecules:
        # Rotational constants (GHZ): ************ 12.73690 12.73690
        # Note: rotational constant will be converted to wavenumber units (1/cm) to standardize across parsers
        if line[:28] == " Rotational constants (GHZ):":
            splits = line[28:].split()
            self.append_attribute("rotconsts", [float(splits[i]) for i in (-3, -2, -1)])

        # Extract Molecular Mass (in amu)
        # Example:
        # Molecular mass:   128.06260 amu.
        if line[:16] == " Molecular mass:":
            splits = line.split()
            self.molmass = utils.float(splits[2])

        if line[:31] == " Normal termination of Gaussian":
            self.metadata["success"] = True

        if line.startswith(" Error termination via"):
            self.metadata["success"] = False

    def parse_keywords_route_and_comment(self, inputfile, line: str) -> None:
        """Parse the keyword line and comment from the echoed input file into metadata."""
        # example:
        #
        # --------------------------------
        # #p hf/sto-3g polar symmetry=none
        # --------------------------------
        # 1/38=1,172=1/1;
        # 2/12=2,15=3,17=6,18=5,40=1/2;
        # 3/6=3,11=9,25=1,30=1/1,2,3;
        # 4//1;
        # 5/5=2,38=5,98=1/2;
        # 8/6=4,10=90,11=11/1;
        # 10/6=1,13=10,31=1/2;
        # 6/7=2,8=2,9=2,10=2,28=1/1;
        # 99/5=1,9=1/99;
        # Leave Link    1 at Wed Apr  4 10:15:34 2018, MaxMem=           0 cpu:               0.0 elap:               0.1
        # (Enter /software/Gaussian16/g16_sse4/g16/l101.exe)
        # --------------------------------
        # tryptophan static polarizability
        # --------------------------------
        keywords = []
        while set(line.strip()) != {"-"}:
            keywords.append(line.strip())
            line = next(inputfile)
        if "keywords" not in self.metadata:
            self.metadata["keywords"] = []
        self.metadata["keywords"].append("".join(keywords))
        line = next(inputfile)
        route_lines = []
        while set(line.strip()) != {"-"}:
            if "Leave Link    1" in line and "MaxMem=" in line and "num_cpu" in self.metadata:
                # Leave Link    1 at Wed Apr  4 10:49:19 2018, MaxMem=   805306368 cpu:               0.3 elap:               0.0
                # Gaussian helpfully prints the total 'available' memory for us. There are, however a few caveats here:
                # 1) This memory (in bytes) is per CPU
                # 2) The total memory here (x num_cpu) will not equal the amount requested in %mem because Gaussian (probably erroneously)
                #    interprets GB as gibibytes (1024 x 1024 x 1024 bytes) rather than gigabytes. This has the unfortunate consequence of
                #    Gaussian assigning more memory than you probably expected.
                #
                # MaxMem can appear with or without whitespace:
                # MaxMem=  805306368
                # MaxMem=174483046400
                maxmem = re.search(r"MaxMem=\s*(\d+)", line)
                assert maxmem is not None
                memory_per_cpu = int(maxmem.group(1))
                self.metadata["memory_used"] = memory_per_cpu * self.metadata["num_cpu"]
            route_lines.append(line)
            line = next(inputfile)
        route_lines = [line for line in route_lines if re.match(r"\s\d", line) is not None]
        route = parse_route_lines(route_lines)
        if "routes" not in self.metadata:
            self.metadata["routes"] = []
        self.metadata["routes"].append(route)
        line = next(inputfile)
        comments = []
        while set(line.strip()) != {"-"}:
            comments.append(line.strip())
            line = next(inputfile)
        if "comments" not in self.metadata:
            self.metadata["comments"] = []
        self.metadata["comments"].append("".join(comments))


# Examples:
#  Gaussian 16:  Apple M1-G16RevC.02  7-Dec-2021
#  Gaussian 16:  ES64L-G16RevC.01  3-Jul-2019
#  Gaussian 98:  IBM-RS6000-G98RevA.7 11-Apr-1999
#  Gaussian 98:  SGI64-G98RevA.7 11-Apr-1999
RE_PLATFORM_AND_VERSION = re.compile(
    r"""
    (?:\ Gaussian\ (?P<year>\d{2}):\ {2})?  # not present in fchk files
    (?P<platform>\w*\ ?\w*(?:-\w*)?)-G(?P<year_suffix>\d{2})Rev(?P<revision>[A-Z]\.\d{1,2}(?:\.\d)?)
    (?:\ +(?P<compile_date>\d{1,2}-[A-Z][a-z]{2}-\d{4}))?  # not present in fchk files
    """,
    re.VERBOSE,
)

# Use to map from the usual year suffixed to full years so package
# versions can be sorted properly after parsing with
# `packaging.parse.version`.
YEAR_SUFFIXES_TO_YEARS = {
    "70": "1970",
    "76": "1976",
    "80": "1980",
    "82": "1982",
    "86": "1986",
    "88": "1988",
    "90": "1990",
    "92": "1992",
    "94": "1994",
    "98": "1998",
    "03": "2003",
    "09": "2009",
    "16": "2016",
}


def parse_version(line: str) -> Dict[str, str]:
    """Extract the version number from Gaussian log and formatted checkpoint
    files.

    From log files, "Gaussian 98: x86-Linux-G98RevA.11.3 5-Feb-2002" becomes
    "1998+A.11.3", and "Gaussian 16: ES64L-G16RevA.03 25-Dec-2016" becomes
    "2016+A.03".

    In formatted checkpoint files, only the middle part of the string is
    present.

    The compile date, which is the last part of the string in log files, is
    not used.
    """
    mtch = RE_PLATFORM_AND_VERSION.match(line)
    assert mtch is not None
    groupdict = mtch.groupdict()
    year_suffix = groupdict["year_suffix"]
    revision = groupdict["revision"]
    return {
        "package_version": f"{YEAR_SUFFIXES_TO_YEARS[year_suffix]}+{revision}",
        "platform": groupdict["platform"],
    }


def _append_charges_and_spins(
    nline: str,
    has_charges: bool,
    charges: List[float],
    has_spin: bool,
    spins: List[float],
    swap_indices: bool,
) -> None:
    """Append partial atomic charges and spins to their respective lists, if
    they were determined beforehand to exist.
    """
    split_line = nline.split()
    if has_charges:
        charge_index = 2 if not swap_indices else 3
        charges.append(float(split_line[charge_index]))
        if has_spin:
            spin_index = 3 if not swap_indices else 2
            spins.append(float(split_line[spin_index]))
    elif has_spin:
        # This will never need swapping, since the particular section where
        # swapping is relevant 1. will always have both charges and spins and
        # 2. has spins first.
        spins.append(float(split_line[2]))


Overlay = int
Option = int
Value = int
Route = Dict[Overlay, Dict[Option, Value]]


def parse_route_line(route_line: str) -> Tuple[Overlay, Dict[Option, Value]]:
    """Parse a single line of the route section.

    This line should be preserved as present in the input, with a single
    leading space, a trailing semicolon, and a possible trailing newline.
    """
    cleaned = route_line.rstrip()
    # remove single leading space and single trailing semicolon
    assert cleaned[-1] == ";"
    cleaned = cleaned[1:-1]
    # it's not clear the stuff at the end means anything
    route, options, _ = cleaned.split("/")
    if not options:
        return int(route), dict()
    return int(route), dict([[int(x) for x in pair.split("=")] for pair in options.split(",")])


def parse_route_lines(route_lines: List[str]) -> Route:
    """Parse the route lines appearing after the input keyword line.

    The route is the section of global internal options (IOps) derived from
    the input keywords and dictates what calculation will be performed
    (https://gaussian.com/iops/).
    """
    route = dict()
    for route_line in route_lines:
        overlay, mapping = parse_route_line(route_line)
        route[overlay] = mapping
    return route
