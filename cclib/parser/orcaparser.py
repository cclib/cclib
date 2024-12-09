# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for ORCA output files"""

import datetime
import re
from itertools import zip_longest
from typing import Callable, Optional, Tuple

from cclib.parser import logfileparser, utils

import numpy
from packaging.version import parse as parse_version


class ORCA(logfileparser.Logfile):
    """An ORCA log file."""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="ORCA", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return f"ORCA log file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'ORCA("{self.filename}")'

    def normalisesym(self, label):
        """ORCA does not require normalizing symmetry labels."""
        return label

    def before_parsing(self):
        self.uses_symmetry = False

        # A geometry optimization is started only when
        # we parse a cycle (so it will be larger than zero().
        self.gopt_cycle = 0

        # Keep track of whether this is a relaxed scan calculation
        self.is_relaxed_scan = False

        # Flag for whether this calc is DFT.
        self.is_DFT = False

        # Used to estimate CPU time from wall time.
        self.metadata["num_cpu"] = 1

        # The excited state multiplicity for post-HF excited states
        self.mdci_et_mult = None

        # needs to be here so regression tests pass
        self.reference = [0.0, 0.0, 0.0]
        
    def sort_et(self):
        # ORCA prints singlet and triplet excited states separately, so the energies are out of order.
        if hasattr(self, "etenergies"):
            prop_names = ("etenergies", "etsyms", "etoscs", "etsecs", "etrotats")

            # First, set energies properly, keeping track of each energy's old index.
            energy_index = sorted(
                [(energy, index) for index, energy in enumerate(self.etenergies)],
                key=lambda energy_index: energy_index[0],
            )

            props = {}
            for prop_name in prop_names:
                if hasattr(self, prop_name):
                    # Check this property and etenergies are the same length (otherwise we can accidentally and silently truncate a list that's too long).
                    if len(getattr(self, prop_name)) != len(self.etenergies):
                        raise Exception(
                            f"Parsed different number of {prop_name} "
                            f"({len(getattr(self, prop_name))}) than "
                            f"etenergies ({len(self.etenergies)})"
                        )

                    # Reorder based on our mapping.
                    props[prop_name] = [
                        getattr(self, prop_name)[old_index] for energy, old_index in energy_index
                    ]

            # Assign back again
            for prop_name in props:
                setattr(self, prop_name, props[prop_name])

    def after_parsing(self):
        super().after_parsing()
        # ORCA doesn't add the dispersion energy to the "Total energy" (which
        # we parse), only to the "FINAL SINGLE POINT ENERGY" (which we don't
        # parse).
        if hasattr(self, "scfenergies") and hasattr(self, "dispersionenergies"):
            for i, (scfenergy, dispersionenergy) in enumerate(
                zip_longest(self.scfenergies, self.dispersionenergies)
            ):
                # It isn't as problematic if there are more dispersion than
                # SCF energies, since all dispersion energies can still be
                # added to the SCF energies, hence the difference in log level.
                if dispersionenergy is None:
                    self.logger.error(
                        "The number of SCF and dispersion energies are not equal: %d vs. %d, "
                        "can't add dispersion energy to all SCF energies",
                        len(self.scfenergies),
                        len(self.dispersionenergies),
                    )
                    break
                if scfenergy is None:
                    self.logger.warning(
                        "The number of SCF and dispersion energies are not equal: %d vs. %d, "
                        "can't add dispersion energy to all SCF energies",
                        len(self.scfenergies),
                        len(self.dispersionenergies),
                    )
                    break
                self.scfenergies[i] += dispersionenergy

        self.sort_et()

        # If we previously stored the mem per cpu, add the total mem now.
        if hasattr(self, "mem_per_cpu"):
            self.metadata["memory_available"] = int(self.mem_per_cpu * self.metadata["num_cpu"])

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # Extract the version number.
        if "Program Version" == line.strip()[:15]:
            # Handle development versions.
            self.metadata["legacy_package_version"] = line.split()[2]
            self.metadata["package_version"] = self.metadata["legacy_package_version"].replace(
                ".x", "dev"
            )
            possible_revision_line = next(inputfile)
            if "SVN: $Rev" in possible_revision_line:
                version = re.search(r"\d+", possible_revision_line).group()
                self.metadata["package_version"] += f"+{version}"
            
            self.version = parse_version(self.metadata["package_version"]).release

        # Extract basis-set info.
        # ----- Orbital basis set information -----
        # Your calculation utilizes the basis: cc-pVDZ
        if "Your calculation utilizes the basis:" == line[:36]:
            self.metadata["basis_set"] = line[37:].strip()

        # ================================================================================
        #                                         WARNINGS
        #                        Please study these warnings very carefully!
        # ================================================================================
        #
        # Warning: TCutStore was < 0. Adjusted to Thresh (uncritical)
        #
        # WARNING: your system is open-shell and RHF/RKS was chosen
        #   ===> : WILL SWITCH to UHF/UKS
        #
        #
        # INFO   : the flag for use of LIBINT has been found!
        #
        # ================================================================================
        if "WARNINGS" == line.strip():
            self.skip_lines(inputfile, ["text", "=", "blank"])
            if "warnings" not in self.metadata:
                self.metadata["warnings"] = []
            if "info" not in self.metadata:
                self.metadata["info"] = []

            line = next(inputfile)
            while line[0] != "=":
                if line.lower()[:7] == "warning":
                    self.metadata["warnings"].append("")
                    while len(line) > 1 and set(line.strip()) != {"="}:
                        self.metadata["warnings"][-1] += line[9:].strip()
                        line = next(inputfile)
                elif line.lower()[:4] == "info":
                    self.metadata["info"].append("")
                    while len(line) > 1 and set(line.strip()) != {"="}:
                        self.metadata["info"][-1] += line[9:].strip()
                        line = next(inputfile)
                else:
                    line = next(inputfile)

        # ================================================================================
        #                                        INPUT FILE
        # ================================================================================
        # NAME = input.dat
        # |  1> %pal nprocs 4 end
        # |  2> ! B3LYP def2-svp
        # |  3> ! Grid4
        # |  4>
        # |  5> *xyz 0 3
        # |  6>     O   0   0   0
        # |  7>     O   0   0   1.5
        # |  8> *
        # |  9>
        # | 10>                          ****END OF INPUT****
        # ================================================================================
        if "INPUT FILE" == line.strip():
            self.skip_line(inputfile, "=")
            self.metadata["input_file_name"] = next(inputfile).split()[-1]

            # First, collect all the lines...
            lines = []
            for line in inputfile:
                if line[0] != "|":
                    break
                lines.append(line[line.find("> ") + 2 :])

            self.metadata["input_file_contents"] = "".join(lines[:-1])
            lines_iter = iter(lines[:-1])

            keywords = []
            coords = []
            # ...then parse them separately.
            for line in lines_iter:
                line = line.strip()
                if not line:
                    continue

                # Keywords block
                if line[0] == "!":
                    keywords += line[1:].split()

                elif line[0:8] == "%MaxCore":
                    self.mem_per_cpu = int(float(line.split()[1]) * 1e6)

                # Impossible to parse without knowing whether a keyword opens a new block
                elif line[0] == "%":
                    pass
                # Geometry block
                elif line[0] == "*":
                    coord_type, charge, multiplicity = line[1:].split()[:3]
                    self.set_attribute("charge", int(charge))
                    self.set_attribute("multiplicity", int(multiplicity))
                    coord_type = coord_type.lower()
                    self.metadata["coord_type"] = coord_type
                    if coord_type == "xyz":

                        def splitter(line):
                            atom, x, y, z = line.split()[:4]
                            return [atom, float(x), float(y), float(z)]
                    elif coord_type in ["int", "internal"]:

                        def splitter(line):
                            atom, a1, a2, a3, bond, angle, dihedral = line.split()[:7]
                            # This could be some combination of floats and variables
                            # C                  0    0    0       0.0                  0.0                   0.0
                            # C                  3    2    1       {B3}                 {A2}                  {D1}
                            try:
                                return [
                                    atom,
                                    int(a1),
                                    int(a2),
                                    int(a3),
                                    float(bond),
                                    float(angle),
                                    float(dihedral),
                                ]
                            except:  # noqa: E722
                                return [
                                    atom,
                                    int(a1),
                                    int(a2),
                                    int(a3),
                                    str(bond),
                                    str(angle),
                                    str(dihedral),
                                ]
                    elif coord_type == "gzmt":

                        def splitter(line):
                            vals = line.split()[:7]
                            if len(vals) == 7:
                                atom, a1, bond, a2, angle, a3, dihedral = vals
                                return [
                                    atom,
                                    int(a1),
                                    float(bond),
                                    int(a2),
                                    float(angle),
                                    int(a3),
                                    float(dihedral),
                                ]
                            elif len(vals) == 5:
                                return [
                                    vals[0],
                                    int(vals[1]),
                                    float(vals[2]),
                                    int(vals[3]),
                                    float(vals[4]),
                                ]
                            elif len(vals) == 3:
                                return [vals[0], int(vals[1]), float(vals[2])]
                            elif len(vals) == 1:
                                return [vals[0]]
                            self.logger.warning("Incorrect number of atoms in input geometry.")
                    elif "file" in coord_type:
                        pass
                    else:
                        self.logger.warning("Invalid coordinate type.")

                    if "file" not in coord_type:
                        for line in lines_iter:
                            if not line:
                                continue
                            if line[0] == "#" or line.strip(" ") == "\n":
                                continue
                            if line.strip()[0] == "*" or line.strip() == "end":
                                break
                            # Strip basis specification that can appear after coordinates
                            line = line.split("newGTO")[0].strip()
                            coords.append(splitter(line))
            self.metadata["keywords"] = keywords
            self.metadata["coords"] = coords

        # Semiempirical methods use a minimal basis fit to Slater functions,
        # not def2-SVP or whatever default is given before the input file is
        # echoed.
        if "FIT TO SLATER BASIS" in line:
            self.metadata["basis_set"] = line.split()[0][4:]

        # If the calculations is a unrelaxed parameter scan then immediately following the
        # input file block is the following section:

        # | 12> **                         ****END OF INPUT****
        # ================================================================================
        #
        #                        ******************************
        #                        * Parameter Scan Calculation *
        #                        ******************************
        #
        # Trajectory settings:
        #     -> SCF surface will be mapped
        #
        # There are 1 parameter(s) to be scanned
        #              R: range=   0.58220000 ..   5.08220000  steps=   46
        # There will be   46 energy evaluations
        #
        #
        # following this each calculation has the following block at the start
        #
        #         *************************************************************
        #                                TRAJECTORY STEP   1
        #                  R  :   0.58220000
        #         *************************************************************

        if "Parameter Scan Calculation" in line:
            self.skip_lines(
                inputfile, ["s", "b", "Trajectory settings", "Surface information", "b"]
            )
            line = next(inputfile)
            num_params = int(line.strip().split()[2])
            for i in range(num_params):
                line = next(inputfile).strip()
                self.append_attribute("scannames", line.split(":")[0])
        if "TRAJECTORY STEP" in line:
            current_params = []
            for i in range(len(self.scannames)):
                line = next(inputfile)
                current_params.append(float(line.split(":")[-1].strip()))
            self.append_attribute("scanparm", tuple(current_params))

        # If the calculations is a relaxed parameter scan then immediately following the
        # input file block is the following section:

        #                        ******************************
        #                        *    Relaxed Surface Scan    *
        #                        ******************************
        #
        #         Dihedral (  9,   8,   3,   2):   range=   0.00000000 .. 360.00000000  steps =   12
        #
        # There is 1 parameter to be scanned.
        # There will be   12 constrained geometry optimizations.
        #
        #
        #          *************************************************************
        #          *               RELAXED SURFACE SCAN STEP   1               *
        #          *                                                           *
        #          *   Dihedral (  9,   8,   3,   2)  :   0.00000000           *
        #          *************************************************************

        if "Relaxed Surface Scan" in line:
            self.skip_lines(inputfile, ["s", "b"])
            line = next(inputfile)
            while not line.isspace():
                line = line.strip()
                self.append_attribute("scannames", line.split(":")[0])
                line = next(inputfile)
            line = next(inputfile)
            num_params = int(line.strip().split()[2])

        if line[0:15] == "Number of atoms":
            natom = int(line.split()[-1])
            self.set_attribute("natom", natom)

        if line[1:13] == "Total Charge":
            charge = int(line.split()[-1])
            self.set_attribute("charge", charge)

            line = next(inputfile)

            mult = int(line.split()[-1])
            self.set_attribute("mult", mult)

        if line[1:18] == "Symmetry handling":
            self.parse_symmetry_section(inputfile)

        if "Density Functional" == line[1:19]:
            self.is_DFT = True
            # In theory we could also parse the functional from this section,
            # but sadly ORCA doesn't print simple functional names.

        # --------------------
        # CPCM SOLVATION MODEL
        # --------------------
        # CPCM parameters:
        #   Epsilon                                         ...       2.3741
        #   Refrac                                          ...       1.4970
        #   Rsolv                                           ...       1.3000
        #   Surface type                                    ... GAUSSIAN VDW
        #   Epsilon function type                           ...         CPCM
        # Radii:
        #  Radius for O  used is    3.4469 Bohr (=   1.8240 Ang.)
        #  Radius for H  used is    2.4944 Bohr (=   1.3200 Ang.)
        # Calculating surface                               ...        done! (  0.0s)
        # GEPOL surface points                              ...          244
        # GEPOL Volume                                      ...     194.1477
        # GEPOL Surface-area                                ...     165.6341
        # Calculating surface distance matrix               ...        done! (  0.0s)
        # Performing Cholesky decomposition & store         ...        done! (  0.0s)
        # Overall time for CPCM initialization              ...                 0.0s
        if line.strip() == "CPCM SOLVATION MODEL":
            # We can assume we're using CPCM if we see this line.
            # SMD also uses this line, but we can update later.
            self.metadata["solvent_model"] = "CPCM"
            self.metadata["solvent_params"] = {}

            line = next(inputfile)
            line = next(inputfile)

            while set(line.strip()) != set("-"):
                line = next(inputfile)

                if "Epsilon function type" in line:
                    if line.split()[-1] == "COSMO":
                        self.metadata["solvent_model"] = "CPCM-COSMO"

                elif "Epsilon" in line:
                    self.metadata["solvent_params"]["epsilon"] = float(line.split()[-1])

                elif "Refrac" in line:
                    self.metadata["solvent_params"]["refractive_index"] = float(line.split()[-1])

                elif "SMD-CDS solvent descriptors" in line:
                    self.metadata["solvent_model"] = "SMD-CPCM"

                elif "Solvent:" in line:
                    # Only get this for SMD.
                    self.metadata["solvent_name"] = line.split()[-1].lower()

        # In Orca <= 5.x, SCF convergence output begins with:
        #
        # --------------
        # SCF ITERATIONS
        # --------------
        #
        # However, there are two common formats which need to be handled, implemented as separate functions.
        # In Orca 6.x, the SCF ITERATIONS header has been removed, and there is no warning before the
        # convergence data is printed.
        if line.strip() == "SCF ITERATIONS" or \
            ("Iteration" in line and "Energy" in line and "Delta-E" in line):
            
            # Skip for old Orca.
            if line.strip() == "SCF ITERATIONS":
                self.skip_line(inputfile, "dashes")
                line = next(inputfile)
            
            columns = line.split()
            # "Starting incremental Fock matrix formation" doesn't
            # necessarily appear before the extended format.
            if not columns:
                self.parse_scf_expanded_format(inputfile, columns)
            # A header with distinct columns indicates the condensed
            # format.
            elif columns[1] == "Energy":
                self.parse_scf_condensed_format(inputfile, columns)
            # Assume the extended format.
            else:
                self.parse_scf_expanded_format(inputfile, columns)

        # Information about the final iteration, which also includes the convergence
        # targets and the convergence values, is printed separately, in a section like this:
        #
        #       *****************************************************
        #       *                     SUCCESS                       *
        #       *           SCF CONVERGED AFTER   9 CYCLES          *
        #       *****************************************************
        #
        # ...
        #
        # Total Energy       :         -382.04963064 Eh          -10396.09898 eV
        #
        # ...
        #
        # -------------------------   ----------------
        # FINAL SINGLE POINT ENERGY     -382.049630637
        # -------------------------   ----------------
        #
        # We cannot use this last message as a stop condition in general, because
        # often there is vibrational output before it. So we use the 'Total Energy'
        # line. However, what comes after that is different for single point calculations
        # and in the inner steps of geometry optimizations.
        if "SCF CONVERGED AFTER" in line:
            if not hasattr(self, "scfvalues"):
                self.scfvalues = []
            if not hasattr(self, "scftargets"):
                self.scftargets = []

            while "Total Energy       :" not in line:
                line = next(inputfile)
            self.append_attribute("scfenergies", utils.float(line.split()[3]))
            if self.is_DFT:
                method = "DFT"
            else:
                semiempirical_methods = _METHODS_SEMIEMPIRICAL & {
                    keyword.upper() for keyword in self.metadata["keywords"]
                }
                assert len(semiempirical_methods) in (0, 1)
                if semiempirical_methods:
                    method = semiempirical_methods.pop()
                else:
                    method = "HF"
            self.metadata["methods"].append(method)

            self._append_scfvalues_scftargets(inputfile, line)

        # Sometimes the SCF does not converge, but does not halt the
        # the run (like in bug 3184890). In this this case, we should
        # remain consistent and use the energy from the last reported
        # SCF cycle. In this case, ORCA print a banner like this:
        #
        #       *****************************************************
        #       *                     ERROR                         *
        #       *           SCF NOT CONVERGED AFTER   8 CYCLES      *
        #       *****************************************************
        if "SCF NOT CONVERGED AFTER" in line:
            if not hasattr(self, "scfvalues"):
                self.scfvalues = []
            if not hasattr(self, "scftargets"):
                self.scftargets = []

            self.append_attribute("scfenergies", self.scfvalues[-1][-1][0])
            self.metadata["methods"].append("HF" if not self.is_DFT else "DFT")

            self._append_scfvalues_scftargets(inputfile, line)

        """
-------------------------------------------------------------------------------
                          DFT DISPERSION CORRECTION

                              DFTD3 V3.1  Rev 1
                              USING zero damping
-------------------------------------------------------------------------------
The omegaB97X-D3 functional is recognized. Fit by Chai et al.
Active option DFTDOPT                   ...         3

molecular C6(AA) [au] = 9563.878941


            DFT-D V3
 parameters
 s6 scaling factor         :     1.0000
 rs6 scaling factor        :     1.2810
 s8 scaling factor         :     1.0000
 rs8 scaling factor        :     1.0940
 Damping factor alpha6     :    14.0000
 Damping factor alpha8     :    16.0000
 ad hoc parameters k1-k3   :    16.0000     1.3333    -4.0000

 Edisp/kcal,au: -10.165629059768  -0.016199959356
 E6   /kcal   :  -4.994512983
 E8   /kcal   :  -5.171116077
 % E8         :  50.868628459

-------------------------   ----------------
Dispersion correction           -0.016199959
-------------------------   ----------------
"""
        if "DFT DISPERSION CORRECTION" in line:
            # A bunch of parameters are printed the first time dispersion is called
            # However, they vary wildly in form and number, making parsing problematic
            line = next(inputfile)
            while "Dispersion correction" not in line:
                line = next(inputfile)
            self.append_attribute("dispersionenergies", utils.float(line.split()[-1]))

        # The convergence targets for geometry optimizations are printed at the
        # beginning of the output, although the order and their description is
        # different than later on. So, try to standardize the names of the criteria
        # and save them for later so that we can get the order right.
        #
        #                        *****************************
        #                        * Geometry Optimization Run *
        #                        *****************************
        #
        # Geometry optimization settings:
        # Update method            Update   .... BFGS
        # Choice of coordinates    CoordSys .... Redundant Internals
        # Initial Hessian          InHess   .... Almoef's Model
        #
        # Convergence Tolerances:
        # Energy Change            TolE     ....  5.0000e-06 Eh
        # Max. Gradient            TolMAXG  ....  3.0000e-04 Eh/bohr
        # RMS Gradient             TolRMSG  ....  1.0000e-04 Eh/bohr
        # Max. Displacement        TolMAXD  ....  4.0000e-03 bohr
        # RMS Displacement         TolRMSD  ....  2.0000e-03 bohr
        #
        if line[25:50] == "Geometry Optimization Run":
            self.skip_lines(inputfile, ["s", "b"])
            line = next(inputfile)
            while line[0:23] != "Convergence Tolerances:":
                line = next(inputfile)

            if hasattr(self, "geotargets"):
                self.logger.warning(
                    "The geotargets attribute should not exist yet. There is a problem in the parser."
                )
            self.geotargets = []
            self.geotargets_names = []

            # There should always be five tolerance values printed here.
            for i in range(5):
                line = next(inputfile)
                name = line[:25].strip().lower().replace(".", "").replace("displacement", "step")
                target = float(line.split()[-2])
                self.geotargets_names.append(name)
                self.geotargets.append(target)

        # The convergence targets for relaxed surface scan steps are printed at the
        # beginning of the output, although the order and their description is
        # different than later on. So, try to standardize the names of the criteria
        # and save them for later so that we can get the order right.
        #
        #         *************************************************************
        #         *               RELAXED SURFACE SCAN STEP  12               *
        #         *                                                           *
        #         *   Dihedral ( 11,  10,   3,   4)  : 180.00000000           *
        #         *************************************************************
        #
        # Geometry optimization settings:
        # Update method            Update   .... BFGS
        # Choice of coordinates    CoordSys .... Redundant Internals
        # Initial Hessian          InHess   .... Almoef's Model
        #
        # Convergence Tolerances:
        # Energy Change            TolE     ....  5.0000e-06 Eh
        # Max. Gradient            TolMAXG  ....  3.0000e-04 Eh/bohr
        # RMS Gradient             TolRMSG  ....  1.0000e-04 Eh/bohr
        # Max. Displacement        TolMAXD  ....  4.0000e-03 bohr
        # RMS Displacement         TolRMSD  ....  2.0000e-03 bohr
        if "RELAXED SURFACE SCAN STEP" in line:
            self.skip_lines(inputfile, ["b"])
            current_params = []
            for i in range(len(self.scannames)):
                line = next(inputfile)
                line = line.replace("*", "")
                current_params.append(float(line.split(":")[-1].strip()))
            self.append_attribute("scanparm", tuple(current_params))

            self.is_relaxed_scan = True
            while "Convergence Tolerances:" not in line:
                line = next(inputfile)

            self.geotargets = []
            self.geotargets_names = []

            # There should always be five tolerance values printed here.
            for i in range(5):
                line = next(inputfile)
                name = line[:25].strip().lower().replace(".", "").replace("displacement", "step")
                target = float(line.split()[-2])
                self.geotargets_names.append(name)
                self.geotargets.append(target)

        # Moller-Plesset energies.
        #
        # ---------------------------------------
        # MP2 TOTAL ENERGY:      -76.112119693 Eh
        # ---------------------------------------
        if "MP2 TOTAL ENERGY" in line[:16]:
            if not hasattr(self, "mpenergies"):
                self.metadata["methods"].append("MP2")

            self.append_attribute("mpenergies", [utils.float(line.split()[-2])])

        # MP2 energy output line is different for MP3, since it uses the MDCI
        # code, which is also in charge of coupled cluster.
        #
        # MP3 calculation:
        # E(MP2)  =    -76.112119775   EC(MP2)=    -0.128216417
        # E(MP3)  =    -76.113783480   EC(MP3)=    -0.129880122  E3=    -0.001663705
        #
        # CCSD calculation:
        # E(MP2)                                     ...     -0.393722942
        # Initial E(tot)                             ...  -1639.631576169
        # <T|T>                                      ...      0.087231847
        # Number of pairs included                   ... 55
        # Total number of pairs                      ... 55
        if "E(MP2)" in line:
            self.append_attribute("mpenergies", [utils.float(line.split()[-1])])

            line = next(inputfile)
            if line[:6] == "E(MP3)":
                self.metadata["methods"].append("MP3")
                self.mpenergies[-1].append(utils.float(line.split()[2]))
            else:
                assert line[:14] == "Initial E(tot)"

        # ----------------------
        # COUPLED CLUSTER ENERGY
        # ----------------------
        #
        # E(0)                                       ...  -1639.237853227
        # E(CORR)                                    ...     -0.360153516
        # E(TOT)                                     ...  -1639.598006742
        # Singles Norm <S|S>**1/2                    ...      0.176406354
        # T1 diagnostic                              ...      0.039445660
        if line[:22] == "COUPLED CLUSTER ENERGY":
            self.skip_lines(inputfile, ["d", "b"])
            line = next(inputfile)
            assert line[:4] == "E(0)"
            scfenergy = float(line.split()[-1])  # noqa: F841
            line = next(inputfile)
            assert line[:7] == "E(CORR)"
            while "E(TOT)" not in line:
                line = next(inputfile)
            self.append_attribute("ccenergies", float(line.split()[-1]))
            self.metadata["methods"].append("CCSD")
            line = next(inputfile)
            assert line[:23] == "Singles Norm <S|S>**1/2"
            line = next(inputfile)
            self.metadata["t1_diagnostic"] = float(line.split()[-1])

        # Most of the "TRIPLES CORRECTION" correction block can be ignored.
        if line[:10] == "E(CCSD(T))":
            self.ccenergies[-1] = float(line.split()[-1])
            assert self.metadata["methods"][-1] == "CCSD"
            self.metadata["methods"].append("CCSD(T)")

        # ------------------
        # CARTESIAN GRADIENT
        # ------------------
        #
        # 1   H   :    0.000000004    0.019501450   -0.021537091
        # 2   O   :    0.000000054   -0.042431648    0.042431420
        # 3   H   :    0.000000004    0.021537179   -0.019501388
        #
        # ORCA MP2 module has different signal than 'CARTESIAN GRADIENT'.
        #
        # The final MP2 gradient
        # 0:   0.01527469  -0.00292883   0.01125000
        # 1:   0.00098782  -0.00040549   0.00196825
        # 2:  -0.01626251   0.00333431  -0.01321825
        if line[:18] == "CARTESIAN GRADIENT" or line[:22] == "The final MP2 gradient":
            grads = []
            if line[:18] == "CARTESIAN GRADIENT":
                self.skip_lines(inputfile, ["dashes", "blank"])

            line = next(inputfile).strip()
            if "CONSTRAINED CARTESIAN COORDINATES" in line:
                self.skip_line(inputfile, "constrained Cartesian coordinate warning")
                line = next(inputfile).strip()

            while line:
                tokens = line.split()
                x, y, z = float(tokens[-3]), float(tokens[-2]), float(tokens[-1])
                grads.append((x, y, z))
                line = next(inputfile).strip()

            if not hasattr(self, "grads"):
                self.grads = []
            self.grads.append(grads)

        # After each geometry optimization step, ORCA prints the current convergence
        # parameters and the targets (again), so it is a good idea to check that they
        # have not changed. Note that the order of these criteria here are different
        # than at the beginning of the output, so make use of the geotargets_names created
        # before and save the new geovalues in correct order.
        #
        #          ----------------------|Geometry convergence|---------------------
        #          Item                value                 Tolerance   Converged
        #          -----------------------------------------------------------------
        #          Energy change       0.00006021            0.00000500      NO
        #          RMS gradient        0.00031313            0.00010000      NO
        #          RMS step            0.01596159            0.00200000      NO
        #          MAX step            0.04324586            0.00400000      NO
        #          ....................................................
        #          Max(Bonds)      0.0218      Max(Angles)    2.48
        #          Max(Dihed)        0.00      Max(Improp)    0.00
        #          -----------------------------------------------------------------
        #
        if line[33:53] == "Geometry convergence":
            self.skip_lines(inputfile, ["header", "d"])
            names = []
            values = []
            targets = []
            line = next(inputfile)
            # Handle both the dots only and dashes only cases
            while len(list(set(line.strip()))) != 1:
                name = line[10:28].strip().lower()
                tokens = line.split()
                value = float(tokens[2])
                target = float(tokens[3])
                names.append(name)
                values.append(value)
                targets.append(target)
                line = next(inputfile)

            # The energy change is normally not printed in the first iteration, because
            # there was no previous energy -- in that case assume zero. There are also some
            # edge cases where the energy change is not printed, for example when internal
            # angles become improper and internal coordinates are rebuilt as in regression
            # CuI-MePY2-CH3CN_optxes, and in such cases use NaN.
            newvalues = []
            for i, n in enumerate(self.geotargets_names):
                if (n == "energy change") and (n not in names):
                    if self.is_relaxed_scan:
                        newvalues.append(0.0)
                    else:
                        newvalues.append(numpy.nan)
                else:
                    newvalues.append(values[names.index(n)])
                    assert targets[names.index(n)] == self.geotargets[i]

            self.append_attribute("geovalues", newvalues)

        """ Grab cartesian coordinates
        ---------------------------------
        CARTESIAN COORDINATES (ANGSTROEM)
        ---------------------------------
        H      0.000000    0.000000    0.000000
        O      0.000000    0.000000    1.000000
        H      0.000000    1.000000    1.000000
        """
        if line[0:33] == "CARTESIAN COORDINATES (ANGSTROEM)":
            next(inputfile)

            atomnos = []
            atomcoords = []
            line = next(inputfile)
            while len(line) > 1:
                atom, x, y, z = line.split()
                if atom[-1] != ">":
                    atomnos.append(self.table.number[atom])
                    atomcoords.append([float(x), float(y), float(z)])
                line = next(inputfile)

            self.set_attribute("natom", len(atomnos))
            self.set_attribute("atomnos", atomnos)
            self.append_attribute("atomcoords", atomcoords)

        """ Grab atom masses
        ----------------------------
        CARTESIAN COORDINATES (A.U.)
        ----------------------------
        NO LB      ZA    FRAG     MASS         X           Y           Z
        0 H     1.0000    0     1.008    0.000000    0.000000    0.000000
        1 O     8.0000    0    15.999    0.000000    0.000000    1.889726
        2 H     1.0000    0     1.008    0.000000    1.889726    1.889726
        """
        if line[0:28] == "CARTESIAN COORDINATES (A.U.)" and not hasattr(self, "atommasses"):
            next(inputfile)
            next(inputfile)

            line = next(inputfile)
            self.atommasses = []
            while len(line) > 1:
                if line[:32] == "* core charge reduced due to ECP":
                    break
                if line.strip() == "> coreless ECP center with (optional) point charge":
                    break
                no, lb, za, frag, mass, x, y, z = line.split()
                if lb[-1] != ">":
                    self.atommasses.append(float(mass))
                line = next(inputfile)

        if line[21:68] == "FINAL ENERGY EVALUATION AT THE STATIONARY POINT":
            if not hasattr(self, "optdone"):
                self.optdone = []
            self.optdone.append(len(self.atomcoords))

        if "The optimization did not converge" in line:
            if not hasattr(self, "optdone"):
                self.optdone = []

        if line[0:16] == "ORBITAL ENERGIES":
            self.skip_lines(inputfile, ["d", "text", "text"])

            self.mooccnos = [[]]
            self.moenergies = [[]]
            self.mosyms = [[]]

            line = next(inputfile)
            # restricted calcs are terminated by ------
            # OR has *Only the first 10 virtual orbitals were printed.
            while len(line) > 20 and line[:5] != '*Only':
                info = line.split()
                mooccno = int(float(info[1]))
                moenergy = float(info[2])
                mosym = "A"
                if self.uses_symmetry:
                    mosym = self.normalisesym(info[4].split("-")[1])
                self.mooccnos[0].append(mooccno)
                self.moenergies[0].append(moenergy)
                self.mosyms[0].append(mosym)
                line = next(inputfile)

            line = next(inputfile)

            # handle beta orbitals for UHF
            if line[17:35] == "SPIN DOWN ORBITALS":
                self.skip_line(inputfile, "text")
                self.mooccnos.append([])
                self.moenergies.append([])
                self.mosyms.append([])

                line = next(inputfile)
                # actually terminated by ------
                # OR has *Only the first 10 virtual orbitals were printed.
                while len(line) > 20 and line[:5] != '*Only':
                    info = line.split()
                    mooccno = int(float(info[1]))
                    moenergy = float(info[2])
                    mosym = "A"
                    if self.uses_symmetry:
                        mosym = self.normalisesym(info[4].split("-")[1])
                    self.mooccnos[1].append(mooccno)
                    self.moenergies[1].append(moenergy)
                    self.mosyms[1].append(mosym)
                    line = next(inputfile)

            if not hasattr(self, "homos"):
                doubly_occupied = self.mooccnos[0].count(2)
                singly_occupied = self.mooccnos[0].count(1)
                # Restricted closed-shell.
                if doubly_occupied > 0 and singly_occupied == 0:
                    self.set_attribute("homos", [doubly_occupied - 1])
                # Restricted open-shell.
                elif doubly_occupied > 0 and singly_occupied > 0:
                    self.set_attribute(
                        "homos", [doubly_occupied + singly_occupied - 1, doubly_occupied - 1]
                    )
                # Unrestricted.
                else:
                    assert len(self.moenergies) == 2
                    assert doubly_occupied == 0
                    assert self.mooccnos[1].count(2) == 0
                    nbeta = self.mooccnos[1].count(1)
                    self.set_attribute("homos", [singly_occupied - 1, nbeta - 1])

        # So nbasis was parsed at first with the first pattern, but it turns out that
        # semiempirical methods (at least AM1 as reported by Julien Idé) do not use this.
        # For this reason, also check for the second patterns, and use it as an assert
        # if nbasis was already parsed. Regression PCB_1_122.out covers this test case.
        if line[1:32] == "# of contracted basis functions":
            self.set_attribute("nbasis", int(line.split()[-1]))
        if line[1:27] == "Basis Dimension        Dim":
            self.set_attribute("nbasis", int(line.split()[-1]))

        if line[0:14] == "OVERLAP MATRIX":
            self.skip_line(inputfile, "dashes")

            self.aooverlaps = numpy.zeros((self.nbasis, self.nbasis), "d")
            for i in range(0, self.nbasis, 6):
                self.updateprogress(inputfile, "Overlap")

                header = next(inputfile)
                size = len(header.split())

                for j in range(self.nbasis):
                    line = next(inputfile)
                    broken = line.split()
                    self.aooverlaps[j, i : i + size] = list(map(float, broken[1 : size + 1]))

        # Molecular orbital coefficients are parsed here, but also related things
        # like atombasis and aonames if possible.
        #
        # Normally the output is easy to parse like this:
        # ------------------
        # MOLECULAR ORBITALS
        # ------------------
        #                       0         1         2         3         4         5
        #                  -19.28527 -19.26828 -19.26356 -19.25801 -19.25765 -19.21471
        #                    2.00000   2.00000   2.00000   2.00000   2.00000   2.00000
        #                   --------  --------  --------  --------  --------  --------
        #   0C   1s         0.000002 -0.000001  0.000000  0.000000 -0.000000  0.000001
        #   0C   2s        -0.000007  0.000006 -0.000002 -0.000000  0.000001 -0.000003
        #   0C   3s        -0.000086 -0.000061  0.000058 -0.000033 -0.000027 -0.000058
        # ...
        #
        # But when the numbers get big, things get yucky since ORCA does not use
        # fixed width formatting for the floats, and does not insert extra spaces
        # when the numbers get wider. So things get stuck together overflowing columns,
        # like this:
        #   12C   6s       -11.608845-53.775398161.302640-76.633779 29.914985 22.083999
        #
        # One assumption that seems to hold is that there are always six significant
        # digits in the coefficients, so we can try to use that to delineate numbers
        # when the parsing gets rough. This is what we do below with a regex, and a case
        # like this is tested in regression ORCA/ORCA4.0/invalid-literal-for-float.out
        # which was reported in https://github.com/cclib/cclib/issues/629
        if line[0:18] == "MOLECULAR ORBITALS":
            self.skip_line(inputfile, "dashes")

            aonames = []
            atombasis = [[] for i in range(self.natom)]
            mocoeffs = [numpy.zeros((self.nbasis, self.nbasis), "d")]
            
            moenergies = []

            for spin in range(len(self.moenergies)):
                moenergies.append([]);
                if spin == 1:
                    self.skip_line(inputfile, "blank")
                    mocoeffs.append(numpy.zeros((self.nbasis, self.nbasis), "d"))

                for i in range(0, self.nbasis, 6):
                    self.updateprogress(inputfile, "Coefficients")

                    line = self.skip_lines(inputfile, ["numbers", "energies"])[-1]
                    moenergies[-1].extend([float(energy) for energy in line.split()])
                    self.skip_lines(inputfile, ["occs", "d"])
#                    self.skip_lines(inputfile, ["numbers", "energies", "occs", "d"])
                    

                    for j in range(self.nbasis):
                        line = next(inputfile)

                        # Only need this in the first iteration.
                        if spin == 0 and i == 0:
                            atomname = line[3:5].split()[0]
                            num = int(line[0:3])
                            orbital = line.split()[1].upper()

                            aonames.append(f"{atomname}{int(num + 1)}_{orbital}")
                            atombasis[num].append(j)

                        # This regex will tease out all number with exactly
                        # six digits after the decimal point.
                        coeffs = re.findall(r"-?\d+\.\d{6}", line)

                        # Something is very wrong if this does not hold.
                        assert len(coeffs) <= 6

                        mocoeffs[spin][i : i + len(coeffs), j] = [float(c) for c in coeffs]

            self.set_attribute("aonames", aonames)
            self.set_attribute("atombasis", atombasis)
            self.set_attribute("mocoeffs", mocoeffs)
            if hasattr(self, 'moenergies') and len(self.moenergies[0]) != self.nbasis:
                self.logger.warning(
                    "Only {} of {} orbital energies parsed from orbital table; switching to lower precision coefficients table".format(len(self.moenergies[0]), self.nbasis)
                )
                # The previously parsed orbital energies can be cut off in Orca 6
                # (only the first 10 virtual orbitals are printed)
                self.set_attribute("moenergies", moenergies, False)

        # Basis set information
        # ORCA prints this out in a somewhat indirect fashion.
        # Therefore, parsing occurs in several steps:
        # 1. read which atom belongs to which basis set group
        if line[0:21] == "BASIS SET INFORMATION":
            line = next(inputfile)
            line = next(inputfile)

            self.tmp_atnames = []  # temporary attribute, needed later
            while not line[0:5] == "-----":
                if line[0:4] == "Atom":
                    self.tmp_atnames.append(line[8:12].strip())
                line = next(inputfile)

        # 2. Read information for the basis set groups
        if line[0:25] == "BASIS SET IN INPUT FORMAT":
            line = next(inputfile)
            line = next(inputfile)

            # loop over basis set groups
            gbasis_tmp = {}
            while not line[0:5] == "-----":
                if line[1:7] == "NewGTO":
                    bas_atname = line.split()[1]
                    gbasis_tmp[bas_atname] = []

                    line = next(inputfile)
                    # loop over contracted GTOs
                    while not line[0:6] == "  end;":
                        words = line.split()
                        ang = words[0]
                        nprim = int(words[1])

                        # loop over primitives
                        coeff = []
                        for iprim in range(nprim):
                            words = next(inputfile).split()
                            coeff.append((float(words[1]), float(words[2])))
                        gbasis_tmp[bas_atname].append((ang, coeff))

                        line = next(inputfile)
                line = next(inputfile)

            # 3. Assign the basis sets to gbasis
            self.gbasis = []
            for bas_atname in self.tmp_atnames:
                self.gbasis.append(gbasis_tmp[bas_atname])
            del self.tmp_atnames

        """
        --------------------------
        THERMOCHEMISTRY AT 298.15K
        --------------------------

        Temperature         ... 298.15 K
        Pressure            ... 1.00 atm
        Total Mass          ... 130.19 AMU

        Throughout the following assumptions are being made:
          (1) The electronic state is orbitally nondegenerate
          ...

        freq.      45.75  E(vib)   ...       0.53
        freq.      78.40  E(vib)   ...       0.49
        ...


        ------------
        INNER ENERGY
        ------------

        The inner energy is: U= E(el) + E(ZPE) + E(vib) + E(rot) + E(trans)
             E(el)   - is the total energy from the electronic structure calc
             ...

        Summary of contributions to the inner energy U:
        Electronic energy                ...   -382.05075804 Eh
        ...
        """
        if line.strip().startswith("THERMOCHEMISTRY AT"):
            self.skip_lines(inputfile, ["dashes", "blank"])
            self.set_attribute("temperature", float(next(inputfile).split()[2]))
            self.set_attribute("pressure", float(next(inputfile).split()[2]))
            total_mass = float(next(inputfile).split()[3])  # noqa: F841

            # Vibrations, rotations, and translations
            line = next(inputfile)
            while line[:17] != "Electronic energy":
                line = next(inputfile)
            self.electronic_energy = float(line.split()[3])
            self.set_attribute("zpve", float(next(inputfile).split()[4]))
            thermal_vibrational_correction = float(next(inputfile).split()[4])  # noqa: F841
            thermal_rotional_correction = float(next(inputfile).split()[4])  # noqa: F841
            thermal_translational_correction = float(next(inputfile).split()[4])
            self.skip_lines(inputfile, ["dashes"])
            total_thermal_energy = float(next(inputfile).split()[3])  # noqa: F841

            # Enthalpy
            # In Orca 6.x, this line gets renamed to Total thermal energy
            while line[:20].strip() not in ["Total free energy", "Total thermal energy"]:
                line = next(inputfile)
            thermal_enthalpy_correction = float(next(inputfile).split()[4])  # noqa: F841
            next(inputfile)

            # For a single atom, ORCA provides the total free energy or inner energy
            # which includes a spurious vibrational correction (see #817 for details).
            if self.natom > 1:
                enthalpy = float(next(inputfile).split()[3])
            else:
                enthalpy = self.electronic_energy + thermal_translational_correction
            self.set_attribute("enthalpy", enthalpy)

            # Entropy
            while line[:18] != "Electronic entropy":
                line = next(inputfile)
            electronic_entropy = float(line.split()[3])
            vibrational_entropy = float(next(inputfile).split()[3])  # noqa: F841
            rotational_entropy = float(next(inputfile).split()[3])  # noqa: F841
            translational_entropy = float(next(inputfile).split()[3])
            self.skip_lines(inputfile, ["dashes"])

            # ORCA prints -inf for single atom entropy.
            if self.natom > 1:
                entropy = float(next(inputfile).split()[4]) / self.temperature
            else:
                entropy = (electronic_entropy + translational_entropy) / self.temperature
            self.set_attribute("entropy", entropy)

            while (line[:25] != "Final Gibbs free enthalpy") and (
                line[:23] != "Final Gibbs free energy"
            ):
                line = next(inputfile)
            self.skip_lines(inputfile, ["dashes"])

            # ORCA prints -inf for single atom free energy, in which case it
            # will be computed after parsing.
            if self.natom > 1:
                self.set_attribute("freeenergy", float(line.split()[5]))

        if line.strip() in (
            "ORCA TD-DFT/TDA CALCULATION",
            "ORCA TD-DFT CALCULATION",
            "ORCA CIS CALCULATION",
            "ORCA ROCIS CALCULATION",
        ):
            # Start of excited states, reset our attributes in case this is an optimised excited state calc
            # (or another type of calc where excited states are calculated multiple times).
            for attr in ("etenergies", "etsyms", "etoscs", "etsecs", "etrotats"):
                if hasattr(self, attr):
                    delattr(self, attr)

            # Excited state metadata.
            if line.strip() == "ORCA ROCIS CALCULATION":
                # Here we consider ROCIS the same as CIS (?)
                self.metadata["excited_states_method"] = "CIS"

            else:
                if "TD-DFT" in line:
                    method = "TD-DFT"

                else:
                    method = "RPA"

                while "Tamm-Dancoff approximation" not in line:
                    line = next(inputfile)

                if line.split()[-1] == "operative":
                    if method == "TD-DFT":
                        method = "TDA"

                    else:
                        method = "CIS"

                self.metadata["excited_states_method"] = method

        # Read TDDFT information
        if any(
            x in line
            for x in ("TD-DFT/TDA EXCITED", "TD-DFT EXCITED", "CIS-EXCITED", "CIS EXCITED")
        ):
            # Could be singlets or triplets
            if line.find("SINGLETS") >= 0:
                mult = "Singlet"
            elif line.find("TRIPLETS") >= 0:
                mult = "Triplet"
            else:
                # This behaviour matches the output Gaussian produces when it encounters an unfamiliar multiplicity.
                mult = "???"

            etsecs = []
            etenergies = []
            etsyms = []

            lookup = {"a": 0, "b": 1}
            line = next(inputfile)
            while line.find("STATE") < 0:
                line = next(inputfile)
            # Contains STATE or is blank
            while line.find("STATE") >= 0:
                broken = line.split()
                etenergies.append(utils.float(broken[3]))
                # In Orca 6, symmetry is printed at the end of the line.
                if len(broken) >= 14 and broken[12] == "Sym:":
                    symm = broken[13]
                    
                else:
                    symm = ""
                
                line = next(inputfile)
                sec = []
                # Contains SEC or is blank
                while line.strip():
                    start = line[0:8].strip()
                    start = (int(start[:-1]), lookup[start[-1]])
                    end = line[10:17].strip()
                    end = (int(end[:-1]), lookup[end[-1]])
                    # Coeffients are not printed for RPA, only
                    # TDA/CIS.
                    contrib = line[35:47].strip()
                    try:
                        contrib = float(contrib)
                    except ValueError:
                        contrib = numpy.nan
                    sec.append([start, end, contrib])
                    line = next(inputfile)
                    # ORCA 5.0 seems to print symmetry at end of block listing transitions
                    if "Symmetry" in line:
                        symm = line.split()[-1]
                        line = next(inputfile)
                    
                etsecs.append(sec)
                if mult != "" and symm != "":
                    etsyms.append(mult + "-" + symm)
                elif mult != "" or symm != "":
                    etsyms.append(mult + symm)
                line = next(inputfile)

            self.extend_attribute("etenergies", etenergies)
            self.extend_attribute("etsecs", etsecs)
            if len(etsyms) > 0:
                self.extend_attribute("etsyms", etsyms)

        # Parse the various absorption spectra for TDDFT and ROCIS.
        if "CD SPECTRUM" not in line and ( "ABSORPTION SPECTRUM" in line or "ELECTRIC DIPOLE" in line ):
            # CASSCF has an anomalous printing of ABSORPTION SPECTRUM.
            if line[:-1] == "ABSORPTION SPECTRUM":
                return

            line = line.strip()

            # Standard header, occasionally changes
            header = ["d", "header", "header", "d"]
            energy_intensity: Optional[Callable[[str], Tuple[float, float]]] = None

            
            if line == "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" and self.version < (6, 0):

                def energy_intensity(line: str) -> Tuple[float, float]:
                    """TDDFT and related methods standard method of output
                    -----------------------------------------------------------------------------
                             ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
                    -----------------------------------------------------------------------------
                    State   Energy  Wavelength   fosc         T2         TX        TY        TZ
                            (cm-1)    (nm)                  (au**2)     (au)      (au)      (au)
                    -----------------------------------------------------------------------------
                       1 5184116.7      1.9   0.040578220   0.00258  -0.05076  -0.00000  -0.00000"""
                    try:
                        state, energy, wavelength, intensity, t2, tx, ty, tz = (
                            utils.float(x) for x in line.split()
                        )
                    except ValueError:
                        # Must be spin forbidden and thus no intensity
                        energy = utils.float(line.split()[1])
                        intensity = 0
                    energy = utils.convertor(energy, "wavenumber", "hartree")
                    return energy, intensity

            elif line == "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" and self.version > (6, 0):
                def energy_intensity(line: str) -> Tuple[float, float]:
                    """TDDFT and related methods standard method of output
                    ----------------------------------------------------------------------------------------------------
                     ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
                    ----------------------------------------------------------------------------------------------------
                         Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ
                                          (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)
                    ----------------------------------------------------------------------------------------------------
                      0-1Ag ->  1-3Bu   3.129277   25239.3   396.2   0.000000000   0.00000   0.00000   0.00000   0.00000"""
                    state1, arrow, state2, energy_ev, energy_wavenumber, wavelength, intensity, t2, tx, ty, tz = line.split()

                    energy = utils.convertor(utils.float(energy_wavenumber), "wavenumber", "hartree")
                    return energy, utils.float(intensity)

            # Check for variations
            elif (
                line == "COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM"
                or line
                == "COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (origin adjusted)"
            ):

                def energy_intensity(line: str) -> Tuple[float, float]:
                    """TDDFT with DoQuad == True
                    ------------------------------------------------------------------------------------------------------
                                    COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM
                    ------------------------------------------------------------------------------------------------------
                    State   Energy Wavelength    D2        m2        Q2         D2+m2+Q2       D2/TOT    m2/TOT    Q2/TOT
                            (cm-1)   (nm)                (*1e6)    (*1e6)
                    ------------------------------------------------------------------------------------------------------
                       1 61784150.6      0.2   0.00000   0.00000   3.23572   0.00000323571519   0.00000   0.00000   1.00000"""
                    (
                        state,
                        energy,
                        wavelength,
                        d2,
                        m2,
                        q2,
                        intensity,
                        d2_contrib,
                        m2_contrib,
                        q2_contrib,
                    ) = (utils.float(x) for x in line.split())
                    energy = utils.convertor(energy, "wavenumber", "hartree")
                    return energy, intensity

            elif (
                line
                == "COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (Origin Independent, Length Representation)"
            ):

                def energy_intensity(line: str) -> Tuple[float, float]:
                    """TDDFT with doQuad == True (Origin Independent Length Representation)
                    -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                                        COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM (Origin Independent, Length Representation)
                    -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    State  Energy   Wavelength       D2            m2              Q2               DM             DO               D2+m2+Q2+DM+DO          D2/TOT          m2/TOT          Q2/TOT         DM/TOT          DO/TOT
                           (cm-1)      (nm)                      (*1e6)          (*1e6)           (*1e6)         (*1e6)
                    -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                       1 61784150.6      0.2      0.00000         0.00000         3.23572         0.00000         0.00000         0.00000323571519         0.00000         0.00000         1.00000         0.00000          0.00000
                       2 61793079.3      0.2      0.00000         0.00000         2.85949         0.00000        -0.00000         0.00000285948800         0.00000         0.00000         1.00000         0.00000         -0.00000"""
                    vals = [utils.float(x) for x in line.split()]
                    energy = utils.convertor(vals[1], "wavenumber", "hartree")
                    if len(vals) < 14:
                        return energy, 0
                    return energy, vals[8]

            elif line[:5] == "X-RAY" and (
                line[6:23] == "EMISSION SPECTRUM" or line[6:25] == "ABSORPTION SPECTRUM"
            ):

                def energy_intensity(line: str) -> Tuple[float, float]:
                    """X-Ray from XES (emission or absorption, electric or velocity dipole moments)
                    -------------------------------------------------------------------------------------
                              X-RAY ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
                    -------------------------------------------------------------------------------------
                           Transition          Energy           INT             TX        TY        TZ
                                                (eV)        (normalized)       (au)      (au)      (au)
                    -------------------------------------------------------------------------------------
                        1   90a ->    0a      8748.824     0.000002678629     0.00004  -0.00001   0.00003"""
                    state, start, arrow, end, energy, intensity, tx, ty, tz = line.split()
                    energy = utils.convertor(utils.float(energy), "eV", "hartree")
                    return energy, intensity

            elif (
                line[:70]
                == "COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE X-RAY"
            ):
                header = ["header", "d", "header", "d", "header", "header", "d"]

                def energy_intensity(line: str) -> Tuple[float, float]:
                    """XAS with quadrupole (origin adjusted)
                    -------------------------------------------------------------------------------------------------------------------------------
                              COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE X-RAY ABSORPTION SPECTRUM
                                                          (origin adjusted)
                    -------------------------------------------------------------------------------------------------------------------------------
                                                                            INT (normalized)
                                                         ---------------------------------------------------------
                           Transition         Energy        D2             M2             Q2           D2+M2+Q2       D2/TOT     M2/TOT     Q2/TOT
                                               (eV)                      (*1e6)         (*1e6)
                    -------------------------------------------------------------------------------------------------------------------------------
                        1   90a ->    0a     8748.824    0.000000       0.000292       0.003615     0.000000027512   0.858012   0.010602   0.131386"""
                    (
                        state,
                        start,
                        arrow,
                        end,
                        energy,
                        d2,
                        m2,
                        q2,
                        intensity,
                        d2_contrib,
                        m2_contrib,
                        q2_contrib,
                    ) = line.split()
                    energy = utils.convertor(utils.float(energy), "eV", "hartree")
                    return energy, intensity

            elif line[:55] == "SPIN ORBIT CORRECTED ABSORPTION SPECTRUM VIA TRANSITION":

                def energy_intensity(line: str) -> Tuple[float, float]:
                    """ROCIS dipole approximation with SOC == True (electric or velocity dipole moments)
                    -------------------------------------------------------------------------------
                    SPIN ORBIT CORRECTED ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
                    -------------------------------------------------------------------------------
                    States    Energy  Wavelength   fosc         T2         TX        TY        TZ
                              (cm-1)    (nm)                  (au**2)     (au)      (au)      (au)
                    -------------------------------------------------------------------------------
                     0  1       0.0      0.0   0.000000000   0.00000   0.00000   0.00000   0.00000
                     0  2 5184116.4      1.9   0.020288451   0.00258   0.05076   0.00003   0.00000"""
                    state, state2, energy, wavelength, intensity, t2, tx, ty, tz = (
                        utils.float(x) for x in line.split()
                    )
                    energy = utils.convertor(energy, "wavenumber", "hartree")
                    return energy, intensity

            elif (
                line[:79]
                == "ROCIS COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM"
                or line[:87]
                == "SOC CORRECTED COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM"
            ):

                def energy_intensity(line: str) -> Tuple[float, float]:
                    """ROCIS with DoQuad = True and SOC = True (also does origin adjusted)
                    ------------------------------------------------------------------------------------------------------
                              ROCIS COMBINED ELECTRIC DIPOLE + MAGNETIC DIPOLE + ELECTRIC QUADRUPOLE SPECTRUM
                    ------------------------------------------------------------------------------------------------------
                    States  Energy Wavelength    D2        m2        Q2         D2+m2+Q2       D2/TOT    m2/TOT    Q2/TOT
                            (cm-1)   (nm)                (*1e6)    (*1e6)     (*population)
                    ------------------------------------------------------------------------------------------------------
                     0  1       0.0      0.0   0.00000   0.00000   0.00000   0.00000000000000   0.00000   0.00000   0.00000
                     0  2 669388066.6      0.0   0.00000   0.00000   0.00876   0.00000000437784   0.00000   0.00000   1.00000"""
                    (
                        state,
                        state2,
                        energy,
                        wavelength,
                        d2,
                        m2,
                        q2,
                        intensity,
                        d2_contrib,
                        m2_contrib,
                        q2_contrib,
                    ) = (utils.float(x) for x in line.split())
                    energy = utils.convertor(energy, "wavenumber", "hartree")
                    return energy, intensity

            # Clashes with Orca 2.6 (and presumably before) TDDFT absorption spectrum printing
            elif line == "ABSORPTION SPECTRUM" and parse_version(
                self.metadata["package_version"]
            ).release > (2, 6):

                def energy_intensity(line: str) -> Tuple[float, float]:
                    """CASSCF absorption spectrum
                    ------------------------------------------------------------------------------------------
                                                    ABSORPTION SPECTRUM
                    ------------------------------------------------------------------------------------------
                      States           Energy   Wavelength   fosc          T2        TX         TY        TZ
                                       (cm-1)     (nm)                   (D**2)      (D)        (D)       (D)
                    ------------------------------------------------------------------------------------------
                      0( 0)-> 1( 0) 1   83163.2    120.2   0.088250385   2.25340   0.00000   0.00000   1.50113"""
                    reg = (
                        r"(\d+)\( ?(\d+)\)-> ?(\d+)\( ?(\d+)\) (\d+)"
                        + r"\s+(\d+\.\d+)" * 4
                        + r"\s+(-?\d+\.\d+)" * 3
                    )
                    res = re.search(reg, line)
                    (
                        jstate,
                        jblock,
                        istate,
                        iblock,
                        mult,
                        energy,
                        wavelength,
                        intensity,
                        t2,
                        tx,
                        ty,
                        tz,
                    ) = res.groups()
                    energy = utils.convertor(utils.float(energy), "wavenumber", "hartree")
                    return energy, intensity

            name = line
            self.skip_lines(inputfile, header)

            if not hasattr(self, "transprop"):
                self.transprop = {}

            # A spectrum section was found, so a function for parsing the energy and intensity is available.
            if energy_intensity is not None:
                etenergies = []
                etoscs = []
                line = next(inputfile)
                # The sections are occasionally ended with dashed lines
                # other times they are blank (other than a new line)
                while len(line.strip("-")) > 2:
                    energy, intensity = energy_intensity(line)
                    etenergies.append(energy)
                    etoscs.append(intensity)

                    line = next(inputfile)

                # Some of these sections contain data that we probably do not want to be populating etenergies
                # and/or etoscs with.  For example, the SOC corrected spectra are for mixed singlet/triplet states,
                # so they do not correspond to the symmetries given in etsyms, and the energy values given are
                # probably not what the user would expect to find in etenergies anyway?
                # Also, there are twice as many SOC states as true spin states, so half of the etenergies wouldn't
                # have a symmetry in etsyms at all...
                #
                # Don't parse from SOC sections.
                # ROCIS COMBINED is combination of SOC and ROCIS (we still parse the normal ROCIS section).
                if not any(
                    [
                        soc_header in name
                        for soc_header in [
                            "SPIN ORBIT CORRECTED",
                            "SOC CORRECTED",
                            "ROCIS COMBINED",
                        ]
                    ]
                ):
                    # We need to be careful about how we parse etenergies from these spectrum sections.
                    # First, and in most cases, energies printed here will be the same as those printed in
                    # previous sections. The energies in cm-1 aught to match exactly to those parsed previously,
                    # but other units may have rounding errors. Occasionally even cm-1 does not match exactly.
                    # Secondly, some methods (ROCIS, CASSCF, SOC to name a few) may only print their final excited state
                    # energies in this spectrum section, in which case the energies will not match those previously parsed
                    # (which will be from lower levels of theory that we're not interested in). This means we cannot simply
                    # ignore the energies printed. Also, in this case we must decide whether to discard other previously
                    # parsed etdata (etsyms, etsecs etc).
                    # Thirdly, SOC prints spin-mixed excited state spectra. This is interesting, but does not match the
                    # number of states or symmetry of data parsed in previous sections, so is not used to overwrite etenergies.

                    # If we have no previously parsed etenergies, there's nothing to worry about.
                    if not hasattr(self, "etenergies"):
                        self.set_attribute("etenergies", etenergies)
                        
                    elif self.version >= (6, 0):
                        # Uniquely (so far), Orca 6 reorders the spectrum states in terms of energy.
                        # Fix our internal states to match.
                        self.sort_et()

                    # Determine if these energies are same as those previously parsed.
                    if len(etenergies) == len(self.etenergies) and numpy.allclose(
                        etenergies, self.etenergies
                    ):
                        pass

                    # New energies.
                    else:
                        # Because these energies are new, we do not know if they correspond to the same level of theory
                        # as the previously parsed etsyms etc.
                        self.logger.warning(
                            "New excited state energies encountered in spectrum section, resetting excited state attributes"
                        )

                        for attr in ("etenergies", "etsyms", "etoscs", "etsecs", "etrotats"):
                            if hasattr(self, attr):
                                delattr(self, attr)

                        self.set_attribute("etenergies", etenergies)

                    self.set_attribute("etoscs", etoscs)

                # Save everything to transprop.
                self.transprop[name] = (numpy.asarray(etenergies), numpy.asarray(etoscs))

        if line.strip() in ["CD SPECTRUM", "CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"]: 
            # -------------------------------------------------------------------
            #                              CD SPECTRUM
            # -------------------------------------------------------------------
            # State   Energy Wavelength       R         MX        MY        MZ
            #         (cm-1)   (nm)       (1e40*cgs)   (au)      (au)      (au)
            # -------------------------------------------------------------------
            #    1   43167.6    231.7      0.00000   0.00000  -0.00000   0.00000
            # ...
            #    6   25291.1    395.4 spin forbidden
            #
            # OR (from 4.2.0 onwards)
            # ------------------------------------------------------------------------------
            #                             CD SPECTRUM
            # ------------------------------------------------------------------------------
            #      States        Energy   Wavelength   R*T        RX        RY        RZ
            #                    (cm-1)      (nm)   (1e40*sgs)   (au)      (au)      (au)
            # ------------------------------------------------------------------------------
            #  0( 1)-> 1( 1) 1   37192.8    268.9     0.00000  -0.00000  -0.34085   0.00000
            # ...
            # ------------------------------------------------------------------------------
            etenergies = []
            etrotats = []
            self.skip_lines(inputfile, ["d", "State   Energy Wavelength", "(cm-1)   (nm)", "d"])
            line = next(inputfile)
            while line.strip() and not utils.str_contains_only(line.strip(), ["-"]):
                tokens = line.split()
                if "spin forbidden" in line:
                    etrotat, mx, my, mz = 0.0, 0.0, 0.0, 0.0
                    etenergy_wavenumber = utils.float(tokens[-4])
                else:
                    etrotat, mx, my, mz = (utils.float(t) for t in tokens[-4:])
                    etenergy_wavenumber = utils.float(tokens[-6])
                etenergies.append(utils.convertor(etenergy_wavenumber, "wavenumber", "hartree"))
                etrotats.append(etrotat)
                line = next(inputfile)
            self.set_attribute("etrotats", etrotats)
            if not hasattr(self, "etenergies"):
                self.logger.warning(
                    "etenergies not parsed before ECD section, the output file may be malformed"
                )
                self.set_attribute("etenergies", etenergies)

        # Read higher-level excited states (EOM-CCSD etc).
        # Multiplicity is in a different section to energies.
        # We can only calculate one type of mult at a time.
        if "Multiplicity                               ..." in line:
            self.mdci_et_mult = line.split()[-1].capitalize()

        if any(
            x in line
            for x in ("CIS RESULTS", "ADC(2) RESULTS", "EOM-CCSD RESULTS", "STEOM-CCSD RESULTS")
        ):
            if "ADC(2)" in line:
                self.metadata["excited_states_method"] = "ADC(2)"

            elif "STEOM-CCSD RESULTS" in line:
                self.metadata["excited_states_method"] = "STEOM-CCSD"

            elif "EOM-CCSD RESULTS" in line:
                self.metadata["excited_states_method"] = "EOM-CCSD"

            if self.mdci_et_mult is None and ("EOM-CCSD" in line or "ADC(2)" in line):
                # These methods can only do singlets.
                # Think this is safe?.
                self.mdci_et_mult = "Singlet"

            # CIS prints orbital contributions different to everone else.
            cis = "CIS RESULTS" in line

            etsecs = []
            etenergies = []
            etsyms = []

            self.skip_lines(inputfile, ["dashes", "blank"])
            line = next(inputfile)
            while line.find("IROOT=") >= 0:
                # ------------------
                # STEOM-CCSD RESULTS
                # ------------------
                #
                # IROOT=  1:  0.120159 au     3.270 eV   26371.8 cm**-1
                #   Amplitude    Excitation
                #    0.104201    61 ->  76
                #    ...
                #   -0.362075    69 ->  73
                #   Ground state amplitude:  0.000000
                #
                # Percentage Active Character     99.28
                #
                #   Amplitude    Excitation in Canonical Basis
                #   -0.158598    64 ->  72
                #   ...
                #   -0.111587    69 ->  75
                #
                # IROOT=  2:  0.123788 au     3.368 eV   27168.4 cm**-1
                #
                # or:
                #
                #  ----------------------
                #  EOM-CCSD RESULTS (RHS)
                #  ----------------------
                #
                #   IROOT=  1: -0.001688 au    -0.046 eV    -370.5 cm**-1
                #     Amplitude    Excitation
                #     -0.693399     x ->  70
                #   Percentage singles character=     93.46
                #
                #   IROOT=  2:  0.061276 au     1.667 eV   13448.5 cm**-1
                etenergies.append(float(line.split()[2]))
                if self.mdci_et_mult is not None:
                    etsyms.append(self.mdci_et_mult)
                sec = []
                # Header line.
                # There is no header for CIS.
                if not cis:
                    line = next(inputfile)
                # First orbital contribution line.
                line = next(inputfile)
                while "->" in line:
                    coeff_split = line.split()
                    # CIS prints coefficients after the orbitals, other modules are reversed.
                    if not cis:
                        contrib = coeff_split[0]
                        start = coeff_split[1]
                        end = coeff_split[3]
                    else:
                        # CIS looks like this:
                        # 32 ->  37    0.037502 (-0.193653)
                        # 34 ->  35    0.935796 ( 0.967365)
                        # 34 ->  36    0.014167 (-0.119024)
                        # Note that sometimes the coefficient has whitespace between the brackets, sometimes not.
                        start = coeff_split[0]
                        end = coeff_split[2]
                        contrib = "".join(coeff_split[4:])[1:-1]

                    try:
                        # TODO: Unrestricted?
                        sec.append([(int(start), 0), (int(end), 0), float(contrib)])

                    except ValueError:
                        # Sometimes we come across lines like:
                        # 0.690372    69 -> x
                        # Ignore these for now.
                        pass

                    line = next(inputfile)

                # Sort contributions so largest is first.
                etsecs.append(sorted(sec, key=lambda sec_item: sec_item[2] ** 2, reverse=True))

                if "Ground state amplitude" in line:
                    # Data currently not parsed. Just skip.
                    line = self.next_filled_line(inputfile)

                if "Percentage singles character" in line:
                    # Data currently not parsed. Just skip.
                    line = self.next_filled_line(inputfile)

                # (Possibly) blank line
                if "IROOT=" in line:
                    continue
                elif line.strip() == "":
                    line = self.next_filled_line(inputfile)

                if "Percentage Active Character " in line:
                    # Data currently not parsed. Just skip.
                    line = self.next_filled_line(inputfile)

                if (
                    "Warning:: the state may have not converged with respect to active space"
                    in line
                ):
                    # Skip this line and the next (which both contain warnings).
                    self.logger.warning(line)
                    line = next(inputfile)
                    line = self.next_filled_line(inputfile)

                if "Amplitude    Excitation in Canonical Basis" in line:
                    # Data currently not parsed. Just skip.
                    # Header line.
                    line = next(inputfile)
                    while "->" in line:
                        line = next(inputfile)

                    line = self.next_filled_line(inputfile)

            # High level excited states will calculate excited states at a number of levels iteratively.
            # We only care about the highest, so overwrite anything from before.
            self.set_attribute("etenergies", etenergies)
            if sum(len(item) for item in etsecs) != 0:
                self.set_attribute("etsecs", etsecs)
            else:
                self.del_attribute("etsecs")

            if len(etsyms) > 0:
                self.set_attribute("etsyms", etsyms)

            else:
                self.del_attribute("etsyms")

        # ---------------
        # CHEMICAL SHIFTS
        # ---------------
        #
        # Note: using conversion factor for au to ppm alpha^2/2 =   26.625677252
        # GIAO: Doing para- and diamagnetic shielding integrals analytically     ...done
        #  --------------
        #  Nucleus   0C :
        #  --------------
        #
        # Diamagnetic contribution to the shielding tensor (ppm) :
        #            261.007         -0.294        -0.000
        #             -0.093        266.836         0.000
        #             -0.000         -0.000       244.502
        #
        # Paramagnetic contribution to the shielding tensor (ppm):
        #           -179.969          6.352        -0.000
        #              3.571       -210.368        -0.000
        #              0.000          0.000       -31.944
        #
        # Total shielding tensor (ppm):
        #             81.038          6.058        -0.000
        #              3.478         56.468        -0.000
        #              0.000          0.000       212.558
        #
        #
        #  Diagonalized sT*s matrix:
        #
        #  sDSO           266.692          261.151          244.502  iso=     257.448
        #  sPSO          -211.114         -179.223          -31.944  iso=    -140.760
        #         ---------------  ---------------  ---------------
        #  Total           55.577           81.929          212.558  iso=     116.688
        # ...
        #
        # --------------------------
        # CHEMICAL SHIELDING SUMMARY (ppm)
        # --------------------------
        #
        #
        #  Nucleus  Element    Isotropic     Anisotropy
        #  -------  -------  ------------   ------------
        #      0       C          116.686        143.809
        #      1       C          122.158        130.692
        # ...
        if line.strip() == "CHEMICAL SHIFTS":
            nmrtensors = dict()
            while line.strip() != "CHEMICAL SHIELDING SUMMARY (ppm)":
                if line[:8] == " Nucleus":
                    atom = int(re.search(r"Nucleus\s+(\d+)\w", line).groups()[0])
                    atomtensors = dict()

                    while "Diagonalized sT*s matrix:" not in line:
                        if "contribution" in line or "Total shielding tensor" in line:
                            # Tensor section.
                            t_type = line.split()[0].lower()

                            # Read the tensor.
                            tensor = numpy.zeros((3, 3))
                            for j, row in zip(range(3), inputfile):
                                tensor[j] = list(map(float, row.split()))

                            atomtensors[t_type] = tensor

                        line = next(inputfile)

                    while "Total" not in line:
                        line = next(inputfile)

                    atomtensors["isotropic"] = float(line.split()[-1])
                    nmrtensors[atom] = atomtensors

                line = next(inputfile)

            self.set_attribute("nmrtensors", nmrtensors)

        # -----------------------------------------------------------
        #  NUCLEUS A = C    0 NUCLEUS B = C    1
        #  ( 13C  gnA =  1.405  13C  gnB =  1.405) r(AB) =     2.8677
        # -----------------------------------------------------------
        #
        # Diamagnetic contribution (Hz):
        #         0.4891        -0.1270       -0.0000
        #        -0.1270        -0.2550        0.0000
        #        -0.0000         0.0000       -0.1388
        # Paramagnetic contribution (Hz):
        #         1.1869         0.2802        0.0000
        #         0.2802        -0.5515        0.0000
        #        -0.0000         0.0000       -0.0408
        # Fermi-contact contribution (Hz):
        #         7.4196         0.0000        0.0000
        #         0.0000         7.4196        0.0000
        #         0.0000         0.0000        7.4196
        # Spin-dipolar contribution (Hz):
        #         0.7215         0.0394       -0.0000
        #         0.0394         1.0985       -0.0000
        #         0.0000         0.0000        3.6092
        # Spin-dipolar/Fermi contact cross term contribution (Hz):
        #         1.9743         0.0164        0.0000
        #         0.0164         2.2237       -0.0000
        #         0.0000        -0.0000       -4.1983
        #
        # Total spin-spin coupling tensor  (Hz):
        #        11.7914         0.2090       -0.0000
        #         0.2090         9.9353       -0.0000
        #         0.0000         0.0000        6.6509
        #
        #  Diagonalized sT*s matrix:
        #
        #  ssDSO           -0.139           -0.218            0.452  iso=       0.032
        #  ssPSO           -0.041           -0.592            1.227  iso=       0.198
        #  ssFC             7.420            7.420            7.420  iso=       7.420
        #  ssSD             3.609            1.085            0.735  iso=       1.810
        #  ssSD/FC         -4.198            2.217            1.981  iso=      -0.000
        #         ---------------  ---------------  ---------------  ----------------
        #  Total            6.651            9.912           11.815  iso=       9.459
        #
        # Sections for NMR spin-spin couplings.
        if "NMR SPIN-SPIN COUPLING CONSTANTS" in line:
            # Reset attributes for upcoming section.
            setattr(self, "nmrcouplingtensors", dict())

        if "NUCLEUS A =" in line and "NUCLEUS B =" in line:
            line_split = line.split()
            # Here we're relying on whitespace between the element symbol and index.
            # For two character elements (eg Cu) and big molecules (>1000 atoms) this space may disappear...
            atoms = (int(line_split[4]), int(line_split[9]))

            # Even though our atom indices reference back to atomnos/atommasses etc, we also need to record
            # the NMR isotope (this isn't recorded anywhere else, and multiple isotopes might get printed).
            line = next(inputfile)
            line_split = line.split()
            # We might have similar whitespace problems here.
            isotopes = (
                int(re.search(r"\d+", line_split[1])[0]),
                int(re.search(r"\d+", line_split[5])[0]),
            )

            # Look for tensor sections.
            # The order and number of tensors is not guaranteed (because different tensors can be
            # explicitly requested).
            tensors = dict()
            while line.strip() not in ["Diagonalized sT*s matrix:", "Diagonalized JT*J matrix:"]:
                if "contribution" in line or "Total spin-spin coupling tensor" in line:
                    # Tensor section.
                    t_type = line.split()[0].lower()

                    # Do some name-nudging.
                    if t_type == "fermi-contact":
                        t_type = "fermi"

                    elif t_type == "spin-dipolar/fermi":
                        t_type = "spin-dipolar-fermi"

                    # Read the tensor.
                    tensor = numpy.zeros((3, 3))
                    for j, row in zip(range(3), inputfile):
                        tensor[j] = list(map(float, row.split()))

                    tensors[t_type] = tensor

                line = next(inputfile)

            while "Total" not in line:
                line = next(inputfile)

            tensors["isotropic"] = float(line.split()[-1])

            if atoms not in self.nmrcouplingtensors:
                self.nmrcouplingtensors[atoms] = {}

            self.nmrcouplingtensors[atoms][isotopes] = tensors

        if line[:23] == "VIBRATIONAL FREQUENCIES":
            self.skip_lines(inputfile, ["d", "b"])

            # Starting with 6.0, the point group is printed.
            if float(self.metadata["package_version"][:3]) >= 6.0:
                self.skip_lines(inputfile, [
                    "Scaling factor for frequencies",
                    "Point group:",
                    "Irrep"
                ])
                
            # Starting with 4.1, a scaling factor for frequencies is printed
            elif float(self.metadata["package_version"][:3]) > 4.0:
                self.skip_lines(inputfile, ["Scaling factor for frequencies", "b"])
            
            

            if self.natom > 1:
                vibfreqs = numpy.zeros(3 * self.natom)
                for i, line in zip(range(3 * self.natom), inputfile):
                    vibfreqs[i] = float(line.split()[1])

                nonzero = numpy.nonzero(vibfreqs)[0]
                self.first_mode = nonzero[0]
                # Take all modes after first
                # Mode between imaginary and real modes could be 0
                self.num_modes = 3 * self.natom - self.first_mode
                if self.num_modes > 3 * self.natom - 6:
                    msg = "Modes corresponding to rotations/translations may be non-zero."
                    if self.num_modes == 3 * self.natom - 5:
                        msg += "\n You can ignore this if the molecule is linear."
                self.set_attribute("vibfreqs", vibfreqs[self.first_mode :])
            else:
                # we have a single atom
                self.set_attribute("vibfreqs", numpy.array([]))

        # NORMAL MODES
        # ------------
        #
        # These modes are the cartesian displacements weighted by the diagonal matrix
        # M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom
        # Thus, these vectors are normalized but *not* orthogonal
        #
        #                   0          1          2          3          4          5
        #       0       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
        #       1       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
        #       2       0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
        # ...
        if line[:12] == "NORMAL MODES":
            if self.natom > 1:
                all_vibdisps = numpy.zeros((3 * self.natom, self.natom, 3), "d")
                
                if self.version >= (6, 0):
                    # Orca 6 once again prints the point group.
                    self.skip_lines(inputfile, ["d", "b", "text", "text", "text", "b", "Point group:", "b"])
                    # And has a wider matrix.
                    matrix_columns = 10
                
                else:
                    self.skip_lines(inputfile, ["d", "b", "text", "text", "text", "b"])
                    matrix_columns = 6
                    

                for mode in range(0, 3 * self.natom, matrix_columns):
                    header = next(inputfile)
                    if self.version >= (6, 0):
                        irreps = next(inputfile)
                    
                    for atom in range(self.natom):
                        all_vibdisps[mode : mode + matrix_columns, atom, 0] = next(inputfile).split()[1:]
                        all_vibdisps[mode : mode + matrix_columns, atom, 1] = next(inputfile).split()[1:]
                        all_vibdisps[mode : mode + matrix_columns, atom, 2] = next(inputfile).split()[1:]
                        
                    if self.version >= (6, 0):
                        self.skip_lines(inputfile, ["b"])

                self.set_attribute("vibdisps", all_vibdisps[self.first_mode :])
            else:
                # we have a single atom
                self.set_attribute("vibdisps", numpy.array([]))

        # ORCA 4 example
        # -----------
        # IR SPECTRUM
        # -----------
        #
        #  Mode    freq (cm**-1)   T**2         TX         TY         TZ
        # -------------------------------------------------------------------
        #    6:      2069.36    1.674341  ( -0.000000   0.914970  -0.914970)
        #    7:      3978.81   76.856228  (  0.000000   6.199041  -6.199042)
        #    8:      4113.34   61.077784  ( -0.000000   5.526201   5.526200)
        # ...

        # ORCA 5 example
        # -----------
        # IR SPECTRUM
        # -----------
        #
        # Mode   freq       eps      Int      T**2         TX        TY        TZ
        #       cm**-1   L/(mol*cm) km/mol    a.u.
        # ----------------------------------------------------------------------------
        #  6:     45.66   0.000006    0.03  0.000039  ( 0.000000  0.000000  0.006256)
        #  7:     78.63   0.000000    0.00  0.000000  ( 0.000000  0.000000  0.000000)
        # ...
        if line[:11] == "IR SPECTRUM":
            package_version = self.metadata.get("package_version", None)
            if package_version is None:
                package_version = "5.x.x"
                self.logger.warning(
                    "package_version has not been set, assuming %s", package_version
                )
            major_version = int(package_version[0])
            if major_version >= 5:
                self.skip_lines(inputfile, ["d", "b", "header", "units", "d"])
                regex = r"\s+(?P<num>\d+):\s+(?P<frequency>\d+\.\d+)\s+(?P<eps>\d+\.\d+)\s+(?P<intensity>\d+\.\d+)"
                
            else:
                self.skip_lines(inputfile, ["d", "b", "header", "d"])
                regex = r"\s+(?P<num>\d+):\s+(?P<frequency>\d+\.\d+)\s+(?P<intensity>\d+\.\d+)"
                

            if self.natom > 1:
                all_vibirs = numpy.zeros((3 * self.natom,), "d")

                line = next(inputfile)
                matches = re.match(regex, line)
                while matches:
                    num = int(matches.group("num"))
                    intensity = float(matches.group("intensity"))
                    all_vibirs[num] = intensity
                    line = next(inputfile)
                    matches = re.match(regex, line)

                self.set_attribute("vibirs", all_vibirs[self.first_mode :])
            else:
                # we have a single atom
                self.set_attribute("vibirs", numpy.array([]))

        # --------------
        # RAMAN SPECTRUM
        # --------------
        #
        #  Mode    freq (cm**-1)   Activity   Depolarization
        # -------------------------------------------------------------------
        #    6:       296.23      5.291229      0.399982
        #    7:       356.70      0.000000      0.749764
        #    8:       368.27      0.000000      0.202068
        if line[:14] == "RAMAN SPECTRUM":
            self.skip_lines(inputfile, ["d", "b", "header", "d"])

            if self.natom > 1:
                all_vibramans = numpy.zeros(3 * self.natom)

                line = next(inputfile)
                while len(line) > 2:
                    num = int(line[0:4])
                    all_vibramans[num] = float(line.split()[2])
                    line = next(inputfile)

                self.set_attribute("vibramans", all_vibramans[self.first_mode :])
            else:
                # we have a single atom
                self.set_attribute("vibramans", numpy.array([]))

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
            self.parse_charge_section(line, inputfile, "mulliken")
        # Things are the same for Lowdin populations, except that the sums
        #   are not printed (there is a blank line at the end).
        if line[:22] == "LOEWDIN ATOMIC CHARGES":
            self.parse_charge_section(line, inputfile, "lowdin")
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
            self.parse_charge_section(line, inputfile, "hirshfeld")
        # CHELPG Charges
        # --------------------------------
        #  0   C   :       0.363939
        #  1   H   :       0.025695
        # ...
        # --------------------------------
        # Total charge:    -0.000000
        # --------------------------------
        if line.startswith("CHELPG Charges"):
            self.parse_charge_section(line, inputfile, "chelpg")

        # The center of mass is used as the origin
        # It is not stated explicitely, but the dipole moment components printed by ORCA
        # seem to be in atomic units, so they will need to be converted.

        # example:
        # The origin for moment calculation is the CENTER OF MASS  = (-1.651256, -1.258772 -1.572312)

        # -------------
        # DIPOLE MOMENT
        # -------------
        #                                 X             Y             Z
        # Electronic contribution:      0.00000      -0.00000      -0.00000
        # Nuclear contribution   :      0.00000       0.00000       0.00000
        #                         -----------------------------------------
        # Total Dipole Moment    :      0.00000      -0.00000      -0.00000
        #                         -----------------------------------------
        # Magnitude (a.u.)       :      0.00000
        # Magnitude (Debye)      :      0.00000
        #
        # TODO: add quadrupole moment parsing, which can be optionally calculated with ORCA

        # the origin/reference might be printed in multiple places in the output file
        # depending on the calculation type
        if line.startswith("The origin for moment calculation is"):
            tmp_reference = line.split()[-3:]
            reference_x = float(tmp_reference[0].replace("(", "").replace(",", ""))
            reference_y = float(tmp_reference[1])
            reference_z = float(tmp_reference[2].replace(")", ""))
            self.reference = numpy.array([reference_x, reference_y, reference_z])

        if line.startswith("DIPOLE MOMENT"):
            self.skip_lines(inputfile, 'd')
            line = next(inputfile) # blank or XYZ
            if line.strip() == '':
                self.skip_lines(inputfile, ['Method', 'Type', 'Multiplicity', 'Irrep', 'Energy', 'Relativity', 'Basis', 'XYZ'])
            self.skip_lines(inputfile, ["electronic", "nuclear", "d"])
            total = next(inputfile)
            assert "Total Dipole Moment" in total

            dipole = numpy.array([float(d) for d in total.split()[-3:]])
            dipole = utils.convertor(dipole, "ebohr", "Debye")

            if not hasattr(self, "moments"):
                self.set_attribute("moments", [self.reference, dipole])
            else:
                try:
                    assert numpy.all(self.moments[1] == dipole)
                except AssertionError:
                    self.logger.warning("Overwriting previous multipole moments with new values")
                    self.set_attribute("moments", [self.reference, dipole])

        if "Molecular Dynamics Iteration" in line:
            self.skip_lines(inputfile, ["d", "ORCA MD", "d", "New Coordinates"])
            line = next(inputfile)
            tokens = line.split()
            assert tokens[0] == "time"
            time = utils.convertor(float(tokens[2]), "time_au", "fs")
            self.append_attribute("time", time)

        # Static polarizability.
        if line.strip() == "THE POLARIZABILITY TENSOR":
            if not hasattr(self, "polarizabilities"):
                self.polarizabilities = []
            self.skip_lines(inputfile, ["d", "b"])
            line = next(inputfile)
            assert line.strip() == "The raw cartesian tensor (atomic units):"
            polarizability = []
            for _ in range(3):
                line = next(inputfile)
                polarizability.append(line.split())
            self.polarizabilities.append(numpy.array(polarizability))

        if line.strip() == "ORCA-CASSCF":
            # -------------------------------------------------------------------------------
            #                               ORCA-CASSCF
            # -------------------------------------------------------------------------------
            #
            # Symmetry handling      UseSym     ... ON
            # Point group                       ... C2
            # Used point group                  ... C2
            # Number of irreps                  ... 2
            #    Irrep    A has   10 SALCs (ofs=   0) #(closed)=   0 #(active)=   2
            #    Irrep    B has   10 SALCs (ofs=  10) #(closed)=   0 #(active)=   2
            #  Symmetries of active orbitals:
            #    MO =    0  IRREP= 0 (A)
            #    MO =    1  IRREP= 1 (B)
            self.skip_lines(inputfile, ["d", "b"])
            vals = next(inputfile).split()
            # Symmetry section is only printed if symmetry is used.
            if vals[0] == "Symmetry":
                assert vals[-1] == "ON"
                point_group = next(inputfile).split()[-1]  # noqa: F841
                used_point_group = next(inputfile).split()[-1]  # noqa: F841
                num_irreps = int(next(inputfile).split()[-1])
                num_active = 0
                # Parse the irreps.
                for i, line in zip(range(num_irreps), inputfile):
                    reg = r"Irrep\s+(\w+) has\s+(\d+) SALCs \(ofs=\s*(\d+)\) #\(closed\)=\s*(\d+) #\(active\)=\s*(\d+)"
                    groups = re.search(reg, line).groups()
                    irrep = groups[0]
                    salcs, ofs, closed, active = map(int, groups[1:])
                    num_active += active
                self.skip_line(inputfile, "Symmetries")
                # Parse the symmetries of the active orbitals.
                for i, line in zip(range(num_active), inputfile):
                    reg = r"(\d+)  IRREP= (\d+) \((\w+)\)"
                    groups = re.search(reg, line).groups()
                    mo, irrep_idx, irrep = groups

            # Skip until the system specific settings.
            # This will align the cases of symmetry on and off.
            line = next(inputfile)
            while line[:25] != "SYSTEM-SPECIFIC SETTINGS:":
                line = next(inputfile)

            # SYSTEM-SPECIFIC SETTINGS:
            # Number of active electrons          ...    4
            # Number of active orbitals           ...    4
            # Total number of electrons           ...    4
            # Total number of orbitals            ...   20
            num_el = int(next(inputfile).split()[-1])  # noqa: F841
            num_orbs = int(next(inputfile).split()[-1])
            total_el = int(next(inputfile).split()[-1])  # noqa: F841
            total_orbs = int(next(inputfile).split()[-1])  # noqa: F841

            line = utils.skip_until_no_match(
                inputfile, r"^\s*$|^Total number aux.*$|^Determined.*$"
            )

            # Determined orbital ranges:
            #    Internal       0 -   -1 (   0 orbitals)
            #    Active         0 -    3 (   4 orbitals)
            #    External       4 -   19 (  16 orbitals)
            orbital_ranges = []
            # Parse the orbital ranges for: Internal, Active, and External orbitals.
            for i in range(3):
                vals = line.split()
                start, end, num = int(vals[1]), int(vals[3]), int(vals[5])
                # Change from inclusive to exclusive in order to match python.
                end = end + 1
                assert end - start == num
                orbital_ranges.append((start, end, num))

            line = next(inputfile)
            while line[:8] != "CI-STEP:":
                line = next(inputfile)

            # CI-STEP:
            # CI strategy                         ... General CI
            # Number of symmetry/multplity blocks ...    1
            # BLOCK  1 WEIGHT=   1.0000
            #   Multiplicity                      ...    1
            #   Irrep                             ...    0 (A)
            #   #(Configurations)                 ...   11
            #   #(CSFs)                           ...   12
            #   #(Roots)                          ...    1
            #     ROOT=0 WEIGHT=    1.000000
            self.skip_line(inputfile, "CI strategy")
            num_blocks = int(next(inputfile).split()[-1])
            for b in range(1, num_blocks + 1):
                line = utils.skip_until_no_match(inputfile, r"^\s*$")
                vals = line.split()
                block = int(vals[1])
                weight = float(vals[3])
                assert b == block
                mult = int(next(inputfile).split()[-1])
                vals = next(inputfile).split()
                # The irrep will only be printed if using symmetry.
                if vals[0] == "Irrep":
                    irrep_idx = int(vals[-2])  # noqa: F841
                    irrep = vals[-1].strip("()")
                    vals = next(inputfile).split()
                num_confs = int(vals[-1])  # noqa: F841
                num_csfs = int(next(inputfile).split()[-1])  # noqa: F841
                num_roots = int(next(inputfile).split()[-1])
                # Parse the roots.
                for r, line in zip(range(num_roots), inputfile):
                    reg = r"=(\d+) WEIGHT=\s*(\d\.\d+)"
                    groups = re.search(reg, line).groups()
                    root = int(groups[0])
                    weight = float(groups[1])  # noqa: F841
                    assert r == root

            # Skip additional setup printing and CASSCF iterations.
            line = next(inputfile).strip()
            while line != "CASSCF RESULTS":
                line = next(inputfile).strip()

            # --------------
            # CASSCF RESULTS
            # --------------
            #
            # Final CASSCF energy       : -14.597120777 Eh    -397.2078 eV
            self.skip_lines(inputfile, ["d", "b"])
            casscf_energy = float(next(inputfile).split()[4])  # noqa: F841

            # This is only printed for first and last step of geometry optimization.
            # ----------------
            # ORBITAL ENERGIES
            # ----------------
            #
            #   NO   OCC          E(Eh)            E(eV)    Irrep
            #    0   0.0868       0.257841         7.0162    1-A
            self.skip_lines(inputfile, ["b", "d"])
            if next(inputfile).strip() == "ORBITAL ENERGIES":
                self.skip_lines(inputfile, ["d", "b", "NO"])
                orbitals = []
                vals = next(inputfile).split()
                while vals:
                    occ, eh, ev = map(float, vals[1:4])
                    # The irrep will only be printed if using symmetry.
                    if len(vals) == 5:
                        idx, irrep = vals[4].split("-")
                        orbitals.append((occ, ev, int(idx), irrep))
                    else:
                        orbitals.append((occ, ev))
                    vals = next(inputfile).split()
                self.skip_lines(inputfile, ["b", "d"])

            # Orbital Compositions
            # ---------------------------------------------
            # CAS-SCF STATES FOR BLOCK  1 MULT= 1 IRREP= Ag NROOTS= 2
            # ---------------------------------------------
            #
            # ROOT   0:  E=     -14.5950507665 Eh
            #       0.89724 [     0]: 2000
            for b in range(num_blocks):
                line = utils.skip_until_no_match(inputfile, r"^\s*$|^-*$")
                # Parse the block data.
                reg = r"BLOCK\s+(\d+) MULT=\s*(\d+) (IRREP=\s*\w+ )?(NROOTS=\s*(\d+))?"
                groups = re.search(reg, line).groups()
                block = int(groups[0])
                mult = int(groups[1])
                # The irrep will only be printed if using symmetry.
                if groups[2] is not None:
                    irrep = groups[2].split("=")[1].strip()
                nroots = int(groups[3].split("=")[1])  # noqa: F841

                self.skip_lines(inputfile, ["d", "b"])

                line = next(inputfile).strip()
                while line:
                    if line[:4] == "ROOT":
                        # Parse the root section.
                        reg = r"(\d+):\s*E=\s*(-?\d+.\d+) Eh(\s+\d+\.\d+ eV)?(\s+\d+\.\d+)?"
                        groups = re.search(reg, line).groups()
                        root = int(groups[0])
                        energy = float(groups[1])
                        # Excitation energies are only printed for excited state roots.
                        if groups[2] is not None:
                            excitation_energy_ev = float(groups[2].split()[0])  # noqa: F841
                            excitation_energy_cm = float(groups[3])  # noqa: F841
                    else:
                        # Parse the occupations section.
                        reg = r"(\d+\.\d+) \[\s*(\d+)\]: (\d+)"
                        groups = re.search(reg, line).groups()
                        coeff = float(groups[0])
                        number = float(groups[1])  # noqa: F841
                        occupations = list(map(int, groups[2]))  # noqa: F841

                    line = next(inputfile).strip()

            # Skip any extended wavefunction printing.
            while line != "DENSITY MATRIX":
                line = next(inputfile).strip()

            self.skip_lines(inputfile, ["d", "b"])
            # --------------
            # DENSITY MATRIX
            # --------------
            #
            #                    0          1          2          3
            #       0       0.897244   0.000000   0.000000   0.000000
            #       1       0.000000   0.533964   0.000000   0.000000
            density = numpy.zeros((num_orbs, num_orbs))
            for i in range(0, num_orbs, 6):
                next(inputfile)
                for j, line in zip(range(num_orbs), inputfile):
                    density[j][i : i + 6] = list(map(float, line.split()[1:]))

            line = utils.skip_until_no_match(inputfile, r"^\s*$|^-*$|^Trace.*$|^Extracting.*$")

            # This is only printed for open-shells.
            # -------------------
            # SPIN-DENSITY MATRIX
            # -------------------
            #
            #                   0          1          2          3          4          5
            #       0      -0.003709   0.001410   0.000074  -0.000564  -0.007978   0.000735
            #       1       0.001410  -0.001750  -0.000544  -0.003815   0.008462  -0.004529
            if line.strip() == "SPIN-DENSITY MATRIX":
                self.skip_lines(inputfile, ["d", "b"])
                spin_density = numpy.zeros((num_orbs, num_orbs))
                for i in range(0, num_orbs, 6):
                    next(inputfile)
                    for j, line in zip(range(num_orbs), inputfile):
                        spin_density[j][i : i + 6] = list(map(float, line.split()[1:]))
                self.skip_lines(inputfile, ["Trace", "b", "d", "ENERGY"])
            self.skip_lines(inputfile, ["d", "b"])

            # -----------------
            # ENERGY COMPONENTS
            # -----------------
            #
            # One electron energy          :    -18.811767801 Eh        -511.8942 eV
            # Two electron energy          :      4.367616074 Eh         118.8489 eV
            # Nuclear repuslion energy     :      0.000000000 Eh           0.0000 eV
            #                                ----------------
            #                                   -14.444151727
            #
            # Kinetic energy               :     14.371970266 Eh         391.0812 eV
            # Potential energy             :    -28.816121993 Eh        -784.1265 eV
            # Virial ratio                 :     -2.005022378
            #                                ----------------
            #                                   -14.444151727
            #
            # Core energy                  :    -13.604678408 Eh     -370.2021 eV
            one_el_energy = float(next(inputfile).split()[4])  # noqa: F841
            two_el_energy = float(next(inputfile).split()[4])  # noqa: F841
            nuclear_repulsion_energy = float(next(inputfile).split()[4])  # noqa: F841
            self.skip_line(inputfile, "dashes")
            energy = float(next(inputfile).strip())
            self.skip_line(inputfile, "blank")
            kinetic_energy = float(next(inputfile).split()[3])  # noqa: F841
            potential_energy = float(next(inputfile).split()[3])  # noqa: F841
            virial_ratio = float(next(inputfile).split()[3])  # noqa: F841
            self.skip_line(inputfile, "dashes")
            energy = float(next(inputfile).strip())
            self.skip_line(inputfile, "blank")
            core_energy = float(next(inputfile).split()[3])  # noqa: F841

        if "Program running with" in line and "parallel MPI-processes" in line:
            # ************************************************************
            # *        Program running with 4 parallel MPI-processes     *
            # *              working on a common directory               *
            # ************************************************************
            self.metadata["num_cpu"] = int(line.split()[4])

        elif "Memory available" in line:
            split_line = line.split()
            if len(split_line) == 5:
                # This is the amount of memory, per cpu. The units are printed afterwards, although
                # it always seems to be in MB...
                #
                # ORCA has a strange relationship with memory management. Some modules appear to be
                # able to use the full allocated amount, some slightly less, some only half (?):
                # Memory available                           ...   2500.00 MB
                # Memory available                           ...   1250.00 MB
                # Memory available                       ... 2291 MB
                # To counter this, we'll always try and store the largest amount available.
                if split_line[4] == "MB":
                    memory = int(float(split_line[3]) * 1e6) * self.metadata["num_cpu"]

                if memory > self.metadata.get("memory_available", 0):
                    self.metadata["memory_available"] = memory

        elif "Maximum memory used throughout the entire" in line:
            # Memory used, making an educated guess that this is per CPU.
            # This is probably also always in MB
            mem_split = line.split()
            memory = float(mem_split[-2])
            mem_units = mem_split[-1]

            if mem_units == "MB":
                memory *= 1e6

            memory *= self.metadata["num_cpu"]

            if memory > self.metadata.get("memory_used", 0):
                self.metadata["memory_used"] = int(memory)

        if line[:15] == "TOTAL RUN TIME:":
            # TOTAL RUN TIME: 0 days 0 hours 0 minutes 11 seconds 901 msec
            self.metadata["success"] = True

            # Parse timings.
            # We also have timings for individual modules (SCF, MDCI etc) which we could use instead?
            time_split = line.split()
            days = int(time_split[3])
            hours = int(time_split[5])
            minutes = int(time_split[7])
            seconds = int(time_split[9])
            milliseconds = int(time_split[11])

            if "wall_time" not in self.metadata:
                self.metadata["wall_time"] = []
            if "cpu_time" not in self.metadata:
                self.metadata["cpu_time"] = []

            self.metadata["wall_time"].append(
                datetime.timedelta(
                    days=days,
                    hours=hours,
                    minutes=minutes,
                    seconds=seconds,
                    milliseconds=milliseconds,
                )
            )

            self.metadata["cpu_time"].append(
                datetime.timedelta(
                    days=days,
                    hours=hours,
                    minutes=minutes,
                    seconds=seconds,
                    milliseconds=milliseconds,
                )
                * self.metadata["num_cpu"]
            )

    def parse_symmetry_section(self, inputfile):
        self.uses_symmetry = True

        line = next(inputfile)
        assert "Point group" in line
        point_group_full = line.split()[3].lower()
        line = next(inputfile)
        # ORCA < 6
        if "Used point group" in line:
            point_group_abelian = line.split()[4].lower()
            line = next(inputfile)
        # ORCA >= 6
        elif "Symmetry-adapted orbitals":
            point_group_abelian = line.split()[3].lower()
            next(inputfile)
            line = next(inputfile)
        assert "Number of irreps" in line
        nirrep = int(line.split()[4])
        for n in range(nirrep):
            line = next(inputfile)
            assert "symmetry adapted basis functions" in line
            irrep = line[8:13]
            if not hasattr(self, "symlabels"):
                self.symlabels = []
            self.symlabels.append(self.normalisesym(irrep))

        self.metadata["symmetry_detected"] = point_group_full
        self.metadata["symmetry_used"] = point_group_abelian
    def parse_charge_section(self, line, inputfile, chargestype):
        """Parse a charge section, modifies class in place

        Parameters
        ----------
        line : str
          the line which triggered entry here
        inputfile : file
          handle to file object
        chargestype : str
          what type of charge we're dealing with, must be one of
          'mulliken', 'lowdin', 'chelpg' or 'hirshfeld'
        """
        has_spins = "AND SPIN POPULATIONS" in line

        if not hasattr(self, "atomcharges"):
            self.atomcharges = {}
        if has_spins and not hasattr(self, "atomspins"):
            self.atomspins = {}

        self.skip_line(inputfile, "dashes")

        # depending on chargestype, decide when to stop parsing lines
        # start, stop - indices for slicing lines and grabbing values
        # should_stop: when to stop parsing
        if chargestype == "mulliken":

            def should_stop(x: str) -> bool:
                return x.startswith("Sum of atomic charges")

            start, stop = 8, 20
        elif chargestype == "lowdin":

            def should_stop(x: str) -> bool:
                return not bool(x.strip())

            start, stop = 8, 20
        elif chargestype == "chelpg":

            def should_stop(x: str) -> bool:
                return x.startswith("---")

            start, stop = 11, 26
        elif chargestype == "hirshfeld":

            def should_stop(x: str) -> bool:
                return not bool(x.strip())

            start, stop = 9, 18
            self.skip_lines(
                inputfile,
                [
                    "d",
                    "b",
                    "Total integrated alpha density",
                    "Total integrated beta density",
                    "header",
                ],
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
            
    def parse_scf_condensed_format(self, inputfile, splitline):
        """"""
        # Possible formats:
        # Orca 5:
        # ITER       Energy         Delta-E        Max-DP      RMS-DP      [F,P]     Damp
        #                ***  Starting incremental Fock matrix formation  ***
        #   0   -384.5203638934   0.000000000000 0.03375012  0.00223249  0.1351565 0.7000
        #   1   -384.5792776162  -0.058913722842 0.02841696  0.00175952  0.0734529 0.7000
        #                                ***Turning on DIIS***
        #   2   -384.6074211837  -0.028143567475 0.04968025  0.00326114  0.0310435 0.0000
        #   3   -384.6479682063  -0.040547022616 0.02097477  0.00121132  0.0361982 0.0000
        #   4   -384.6571124353  -0.009144228947 0.00576471  0.00035160  0.0061205 0.0000
        #   5   -384.6574659959  -0.000353560584 0.00191156  0.00010160  0.0025838 0.0000
        #   6   -384.6574990782  -0.000033082375 0.00052492  0.00003800  0.0002061 0.0000
        #   7   -384.6575005762  -0.000001497987 0.00020257  0.00001146  0.0001652 0.0000
        #   8   -384.6575007321  -0.000000155848 0.00008572  0.00000435  0.0000745 0.0000
        #          **** Energy Check signals convergence ****
        # Orca 6:
        # ----------------------------------------D-I-I-S--------------------------------------------
        # Iteration    Energy (Eh)           Delta-E    RMSDP     MaxDP     DIISErr   Damp  Time(sec)
        # -------------------------------------------------------------------------------------------
        #                ***  Starting incremental Fock matrix formation  ***
        #     1    -382.0551222939103582     0.00e+00  1.60e-04  1.99e-03  3.38e-03  0.700   1.0
        #                               *** Initializing SOSCF ***
        # ---------------------------------------S-O-S-C-F--------------------------------------
        # Iteration    Energy (Eh)           Delta-E    RMSDP     MaxDP     MaxGrad    Time(sec)
        # --------------------------------------------------------------------------------------
        #     2    -382.0551463039915916    -2.40e-05  3.97e-04  5.10e-03  1.52e-03     0.4
        #                *** Restarting incremental Fock matrix formation ***
        #     3    -382.0551558700066153    -9.57e-06  5.54e-04  9.87e-03  1.10e-03     1.0
        #     4    -382.0551154272437770     4.04e-05  3.63e-04  6.82e-03  2.09e-03     0.8
        #     5    -382.0551748409881156    -5.94e-05  5.18e-05  4.86e-04  8.23e-05     0.8
        #     6    -382.0551748903446878    -4.94e-08  2.54e-05  3.18e-04  9.86e-05     0.7
        #                  **** Energy Check signals convergence ****
        #
        
        # Decide which version this is
        headers = splitline[0][0:5] == "-----"
        # Based on that, decide on our column positions:
        if headers:
            # New, Orca 6 format:
            index = {
                'energy':1,
                'deltaE': 2,
                'maxDP': 4,
                'rmsDP': 3
            }
        else:
            # Old, Orca 5 format.
            index = {
                'energy':1,
                'deltaE': 2,
                'maxDP': 3,
                'rmsDP': 4
            }
        
        self.append_attribute("scfvalues", [])
        #diis_active = True
        
        # Stop on newline.
        while splitline:
            
            maxDP = None
            # Taken from old code. I don't think these statements can ever get triggered
            # (looks like they want to check the whole line contents rather than any of
            # the split elements). Furthermore we don't have any tests to check their
            # function or intent...
#             if "Newton-Raphson" in splitline:
#                 diis_active = False
#                 
#             elif "SOSCF" in splitline:
#                 diis_active = False
#                 
#             el
            if splitline[0].isdigit():
                # In some (Orca 4 at least) versions of Orca there appears to be a printing error, and columns can overlap each other
                # looks like:
                # %3i %17.10f%12.12f%11.8f %11.8f
                if splitline[1].count(".") == 2:
                    integer1, decimal1_integer2, decimal2 = splitline[1].split(".")
                    decimal1, integer2 = decimal1_integer2[:10], decimal1_integer2[10:]
                    splitline = [
                        splitline[0],
                        integer1 + "." + decimal1,
                        integer2 + "." + decimal2
                    ] + splitline[2:]
                    
                elif splitline[1].count(".") == 3:
                    integer1, decimal1_integer2, decimal2_integer3, decimal3 = splitline[1].split(
                        "."
                    )
                    decimal1, integer2 = decimal1_integer2[:10], decimal1_integer2[10:]
                    decimal2, integer3 = decimal2_integer3[:12], decimal2_integer3[12:]
                    splitline = [
                        splitline[0],
                        integer1 + "." + decimal1,
                        integer2 + "." + decimal2,
                        integer3 + "." + decimal3
                    ] + splitline[2:]
                    
                elif splitline[2].count(".") == 2:
                    integer1, decimal1_integer2, decimal2 = splitline[2].split(".")
                    decimal1, integer2 = decimal1_integer2[:12], decimal1_integer2[12:]
                    
                    splitline = [
                        splitline[0],
                        splitline[1],
                        integer1 + "." + decimal1,
                        integer2 + "." + decimal2
                    ] + splitline[3:]
                
                
                deltaE = float(splitline[index['deltaE']])
                maxDP = float(splitline[index['maxDP']]) #+ int(not diis_active)])
                rmsDP = float(splitline[index['rmsDP']]) #+ int(not diis_active)])

                self.scfvalues[-1].append([deltaE, maxDP, rmsDP])
            
            try:
                splitline = next(inputfile).split()
            
            except StopIteration:
                self.logger.warning(f"File terminated before end of last SCF! Last Max-DP: {maxDP}")
                break

    def parse_scf_expanded_format(self, inputfile, line):
        """Parse SCF convergence when in expanded format."""

        # The following is an example of the format
        # -----------------------------------------
        #
        #               ***  Starting incremental Fock matrix formation  ***
        #
        #                         ----------------------------
        #                         !        ITERATION     0   !
        #                         ----------------------------
        #   Total Energy        :    -377.960836651297 Eh
        #   Energy Change       :    -377.960836651297 Eh
        #   MAX-DP              :       0.100175793695
        #   RMS-DP              :       0.004437973661
        #   Actual Damping      :       0.7000
        #   Actual Level Shift  :       0.2500 Eh
        #   Int. Num. El.       :    43.99982197 (UP=   21.99991099 DN=   21.99991099)
        #   Exchange            :   -34.27550826
        #   Correlation         :    -2.02540957
        #
        #
        #                         ----------------------------
        #                         !        ITERATION     1   !
        #                         ----------------------------
        #   Total Energy        :    -378.118458080109 Eh
        #   Energy Change       :      -0.157621428812 Eh
        #   MAX-DP              :       0.053240648588
        #   RMS-DP              :       0.002375092508
        #   Actual Damping      :       0.7000
        #   Actual Level Shift  :       0.2500 Eh
        #   Int. Num. El.       :    43.99994143 (UP=   21.99997071 DN=   21.99997071)
        #   Exchange            :   -34.00291075
        #   Correlation         :    -2.01607243
        #
        #                               ***Turning on DIIS***
        #
        #                         ----------------------------
        #                         !        ITERATION     2   !
        #                         ----------------------------
        # ....
        #
        if not hasattr(self, "scfvalues"):
            self.scfvalues = []

        self.scfvalues.append([])

        line = "Foo"  # dummy argument to enter loop
        while line.find("******") < 0:
            try:
                line = next(inputfile)
            except StopIteration:
                self.logger.warning("File terminated before end of last SCF!")
                break
            info = line.split()
            if len(info) > 1 and info[1] == "ITERATION":
                self.skip_line(inputfile, "d")
                energy_line = next(inputfile).split()
                energy = float(energy_line[3])
                deltaE_line = next(inputfile).split()
                deltaE = float(deltaE_line[3])
                if energy == deltaE:
                    deltaE = 0
                maxDP_line = next(inputfile).split()
                maxDP = float(maxDP_line[2])
                rmsDP_line = next(inputfile).split()
                rmsDP = float(rmsDP_line[2])
                self.scfvalues[-1].append([deltaE, maxDP, rmsDP])

        return

    # end of parse_scf_expanded_format

    def _append_scfvalues_scftargets(self, inputfile, line):
        # The SCF convergence targets are always printed in this next section
        # but which targets are available depends on the SCF method in use,
        # among other things.
        while "Last Energy change" not in line:
            line = next(inputfile)

        deltaE_value = float(line.split()[4])
        deltaE_target = float(line.split()[7])
        maxDP_value = None
        rmsDP_value = None
        maxDP_target = None
        rmsDP_target = None

        line = next(inputfile)
        if "Last MAX-Density change" in line:
            maxDP_value = float(line.split()[4])
            maxDP_target = float(line.split()[7])
            line = next(inputfile)
            if "Last RMS-Density change" in line:
                rmsDP_value = float(line.split()[4])
                rmsDP_target = float(line.split()[7])
            else:
                rmsDP_value = None
                rmsDP_target = None

        if len(self.scfvalues) > 0:
            self.scfvalues[-1].append([deltaE_value, maxDP_value, rmsDP_value])
            self.scftargets.append([deltaE_target, maxDP_target, rmsDP_target])
        else:
            self.logger.warning('No SCF values when parsing final changes')


_METHODS_SEMIEMPIRICAL = {"AM1", "MNDO", "PM3", "ZINDO/1", "ZINDO/S"}
