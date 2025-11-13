# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Serenity output files"""

import datetime
from pathlib import Path

from cclib.parser import data, logfileparser, utils
from cclib.parser.logfileparser import StopParsing

import numpy


class Serenity(logfileparser.Logfile):
    """A Serenity output file"""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="Serenity", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return f"Serenity output file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'Serenity("{self.filename}")'

    def normalisesym(self, label):
        """Serenity does not require normalizing symmetry labels."""
        return label

    def before_parsing(self):
        self.skipSystem = False
        self.skipSCF = False
        self.populationdict = {}
        # TODO Serenity also has Becke, IAO, and CHELPG population analyses that are currently not supported.
        self.populationtypes = ["Mulliken", "CM5", "Hirshfeld"]
        self.metadata["unrestricted"] = False
        self.beta_parsing = False
        self.path = Path(self.inputfile.filenames[0]).resolve()

    def after_parsing(self):
        self.populationdict = {k.lower(): v for k, v in self.populationdict.items()}
        self.set_attribute("atomcharges", self.populationdict)
        # Get molecular orbital information
        assert utils.find_package("h5py"), (
            "h5py is needed to read in molecular orbital info from Serenity."
        )
        import h5py

        if self.metadata["unrestricted"]:
            orbpath = self.path.parent / self.systemname / f"{self.systemname}.orbs.unres.h5"
            if orbpath.is_file():
                with h5py.File(orbpath, "r") as orbfile:
                    eigenvalues = []
                    coeffs = []
                    eigenvalues.append(orbfile["eigenvalues_alpha"][:].flatten())
                    eigenvalues.append(orbfile["eigenvalues_beta"][:].flatten())
                    coeffs.append(orbfile["coefficients_alpha"][:])
                    coeffs.append(orbfile["coefficients_alpha"][:])
                    self.set_attribute("moenergies", eigenvalues)
                    self.set_attribute("mocoeffs", coeffs)
                    self.set_attribute("nmo", len(eigenvalues[0]))
        else:
            orbpath = self.path.parent / self.systemname / f"{self.systemname}.orbs.res.h5"
            if orbpath.is_file():
                with h5py.File(orbpath, "r") as orbfile:
                    coeffs = [orbfile["coefficients"][:]]
                    eigenvalues = [orbfile["eigenvalues"][:].flatten()]
                    self.set_attribute("moenergies", eigenvalues)
                    self.set_attribute("mocoeffs", coeffs)
                    self.set_attribute("nmo", len(eigenvalues[0]))

        super().after_parsing()

    def convert_to_spin(self, symbol):
        if symbol == "a":
            return 0
        elif symbol == "b":
            return 1
        else:
            self.logger.error(
                "Unexpected symbol encountered while parsing for spins of transitions."
            )

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        ### SYSTEM specific data
        # Extract system name
        if line.strip().startswith("------------------------------------------------------------"):
            line = next(inputfile)
            if line.strip().startswith("System"):
                name = line.split()[1]
                if not hasattr(self, "systemname"):
                    self.systemname = name
                elif name[-5:] == "_FREE":
                    self.logger.warning(
                        "Skipping the parsing of the atom SCFs done by some population analysis methods in Serenity."
                    )
                    self.skipSystem = True
                    self.skipSCF = True
                elif name != self.systemname:
                    self.logger.error("Multiple systems detected which is not currently supported.")
                    raise StopParsing()

        if not self.skipSystem:
            # Extract charge and multiplicity
            if line[5:11] == "Charge":
                self.set_attribute("charge", int(line.split()[1]))

            # Extract multiplicity
            if line[5:9] == "Spin":
                self.set_attribute("mult", int(line.split()[1]) + 1)

            # Extract from atoms: number of atoms, elements, and coordinates
            if line.strip().startswith("Current Geometry (Angstrom):"):
                line = next(inputfile)
                line = next(inputfile)
                atomnos = []
                coords = []
                while line.strip():
                    atominfo = line.split()
                    element = atominfo[1]
                    x, y, z = map(float, atominfo[2:5])
                    atomnos.append(self.table.number[element])
                    coords.append([x, y, z])
                    line = next(inputfile)

                self.set_attribute("atomnos", atomnos)
                self.set_attribute("natom", len(atomnos))
                if hasattr(self, "optstatus"):
                    if (
                        not self.optstatus[-1] & data.ccData.OPT_DONE
                        and not len(self.optstatus) == 1
                    ):
                        self.append_attribute("atomcoords", coords)
                else:
                    self.append_attribute("atomcoords", coords)

        if not self.skipSCF:
            ### SCF data
            if line[5:21] == "Basis Functions:":
                self.set_attribute("nbasis", int(line.split()[2]))

            if line[5:15] == "Basis Set:":
                self.metadata["basis_set"] = line.split()[2]
            if line[5:16] == "Functional:":
                self.metadata["functional"] = line.split()[1]

            # Extract SCF thresholds
            if line.strip().startswith("Energy Threshold:"):
                ethresh = float(line.split()[2])
                line = next(inputfile)
                if "RMSD[D]" in line:
                    rmsd = float(line.split()[2])
                    line = next(inputfile)
                    if "DIIS" in line:
                        diis = float(line.split()[2])
                        scftargets = [ethresh, rmsd, diis]
                        self.append_attribute("scftargets", scftargets)

            if line[5:13] == "SCF Mode":
                if line.split()[2] == "RESTRICTED":
                    self.metadata["unrestricted"] = False
                elif line.split()[2] == "UNRESTRICTED":
                    self.metadata["unrestricted"] = True
                line = next(inputfile)
                if line.startswith("Method"):
                    self.metadata["methods"].append(line.split()[1])

            if "Cycle" in line and "Mode" in line:
                line = next(inputfile)
                values = []
                while not line.strip().startswith("Converged after"):
                    linedata = line.split()
                    c1, c2, c3 = map(float, linedata[2:5])
                    values.append([c1, c2, c3])
                    line = next(inputfile)
                self.append_attribute("scfvalues", numpy.vstack(numpy.array(values)))

            if line.startswith("Total Energy ("):
                self.append_attribute("scfenergies", float(line.split()[3]))
                if not hasattr(self, "optstatus") and hasattr(self, "scfenergies"):
                    self.logger.warning(
                        "Multiple instances of scfenergies despite no geometry optimization being done. This Serenity calculation possibly has several systems."
                    )

            if "Total Supersystem Energy" in line:
                self.logger.warning(
                    "Supersystem energy encountered. Serenity calculations involving subsystems (and by extension, supersystems) are currently not supported."
                )
                raise StopParsing()

            if "Dispersion Correction (" in line:
                self.append_attribute("dispersionenergies", float(line.split()[3]))

            # Extract index of HOMO(s)
            if line.strip().startswith("Orbital Energies:"):
                self.skip_line(inputfile, ["Orbital"])
                self.skip_line(inputfile, ["dashes"])
                if self.metadata["unrestricted"] and not self.beta_parsing:
                    self.skip_line(inputfile, ["Alpha:"])
                    self.beta_parsing = True
                self.skip_line(inputfile, ["#   Occ."])
                # self.skip_lines(inputfile, ["Orbital","dashes","#   Occ."]) # TODO test results in warnings
                line = next(inputfile)
                homos = None
                occ_number = None
                if self.metadata["unrestricted"]:
                    occ_number = "1.00"
                else:
                    occ_number = "2.00"
                while line.split()[1] == occ_number:
                    homos = int(line.split()[0])
                    line = next(inputfile)
                self.append_attribute("homos", homos - 1)  # Serenity starts at 1, python at 0

                if self.beta_parsing:
                    while line.split()[0] != "Beta:":
                        line = next(inputfile)
                    self.skip_line(inputfile, ["Beta:"])
                    self.skip_line(inputfile, ["#   Occ."])
                    line = next(inputfile)
                    while line.split()[1] == occ_number:
                        homos = int(line.split()[0])
                        line = next(inputfile)
                    self.append_attribute("homos", homos - 1)

        ### Multipole moments
        if line.strip().startswith("Origin chosen as:"):
            line = self.skip_line(inputfile, "Origin chosen as:")[0]
            origin_data = line.replace("(", "").replace(")", "").replace(",", "").split()
            x, y, z = map(float, origin_data)
            origin = [x, y, z]
            self.append_attribute("moments", origin)

        if line.strip().startswith("Dipole Moment:"):
            self.skip_line(inputfile, ["Dipole Moment"])
            self.skip_line(inputfile, ["dashes"])
            # self.skip_lines(inputfile, ["Dipole Moment","dashes"]) # TODO test results in warnings
            line = self.skip_line(inputfile, "x")[0]
            dipole_data = line.split()
            x, y, z = map(float, dipole_data[:3])
            dipole_raw = [x, y, z]
            dipole = [utils.convertor(value, "ebohr", "Debye") for value in dipole_raw]
            self.append_attribute("moments", dipole)

        if line.strip().startswith("Quadrupole Moment:"):
            self.skip_line(inputfile, "Quadrupole Moment")
            self.skip_line(inputfile, ["dashes"])
            line = self.skip_line(inputfile, "x")[0]
            q_data = line.split()
            xx, xy, xz = map(float, q_data[1:4])
            line = self.skip_line(inputfile, "x")[0]
            q_data = line.split()
            yy, yz = map(float, q_data[2:4])
            line = self.skip_line(inputfile, "y")[0]
            q_data = line.split()
            zz = float(q_data[3])
            quadrupole_raw = [xx, xy, xz, yy, yz, zz]
            quadrupole = [
                utils.convertor(value, "ebohr2", "Buckingham") for value in quadrupole_raw
            ]
            self.append_attribute("moments", quadrupole)

        # Extract charges from population analysis
        if line.strip().startswith(tuple(self.populationtypes)) and line.split()[1] == "Population":
            key = line.split()[0]
            line = next(inputfile)
            self.skip_line(inputfile, "dashes")
            self.skip_line(inputfile, "blank")
            line = next(inputfile)
            line = next(inputfile)
            chargelist = []
            while line.strip() and not line.strip().startswith("---"):
                charge = float(line.split()[3])
                chargelist.append(charge)
                line = next(inputfile)
            value = numpy.array(chargelist)
            self.populationdict[key] = value

        if "Total Local-CCSD Energy" in line:
            self.set_attribute("ccenergies", float(line.split()[3]))
            self.metadata["methods"].append("Local CCSD")  # TODO might not work
        if "Total Local-CCSD(T0) Energy" in line:
            self.set_attribute("ccenergies", float(line.split()[3]))
            self.metadata["methods"].append("Local CCSD(T0)")  # TODO might not work
        if "Total CCSD Energy" in line:
            self.set_attribute("ccenergies", float(line.split()[3]))
            self.metadata["methods"].append("CCSD")
        if "Total CCSD(T) Energy" in line:
            self.set_attribute("ccenergies", float(line.split()[3]))
            self.metadata["methods"].append("CCSD(T)")

        # for additional robustness.
        # this should already be stopped by the presence of several systems in the output
        if line.strip().startswith("Freeze-and-Thaw Cycle"):
            self.logger.error(
                "Freeze-and-Thaw subsystem calculation detected which is not currently supported."
            )
            raise StopParsing()

        ### geometry optimization
        if line.strip().startswith("Cycle:"):
            self.append_attribute("optstatus", data.ccData.OPT_UNKNOWN)

        if hasattr(self, "optstatus"):
            if line.split() == ["Cycle:", "1"]:
                self.optstatus[-1] += data.ccData.OPT_NEW

            if line.strip().startswith("WARNING: Geometry Optimization not yet converged!"):
                self.set_attribute("optdone", [])
                self.optstatus[-1] += data.ccData.OPT_UNCONVERGED

            if line.strip().startswith("Convergence reached after"):
                if not self.optstatus[-1] & data.ccData.OPT_DONE:
                    self.append_attribute("optdone", len(self.atomcoords) - 1)
                    self.optstatus[-1] += data.ccData.OPT_DONE
                    # removing unnecessary entries of HOMOs in case of geom. opt. etc.
                    if self.metadata["unrestricted"]:
                        self.homos = self.homos[0:2]
                    else:
                        self.homos = self.homos[0]

            if line.strip().startswith("Current Geometry Gradients (a.u.):"):
                if not self.optstatus[-1] & data.ccData.OPT_DONE:
                    self.skip_line(inputfile, ["Current"])
                    line = next(inputfile)
                    grad = []
                    for i in range(self.natom):
                        grad_data_raw = line.split()
                        x, y, z = map(float, grad_data_raw[2:5])
                        grad.append([x, y, z])
                        line = next(inputfile)
                    self.append_attribute("grads", grad)

            # The 5 convergence criteria in Serenity, of which 3 must be met, are (in order):
            # Energy Change, RMS Gradient, Max Gradient, RMS Step, Max Step
            if line.strip().startswith("Geometry Relaxation:"):
                criteria = []
                self.skip_line(inputfile, ["Geometry Relaxation:"])
                self.skip_line(inputfile, ["dashes"])
                line = next(inputfile)
                criteria.append(line.split()[2])
                line = next(inputfile)
                line = next(inputfile)
                criteria.append(line.split()[2])
                line = next(inputfile)
                criteria.append(line.split()[2])
                line = next(inputfile)
                criteria.append(line.split()[2])
                line = next(inputfile)
                criteria.append(line.split()[2])
                self.append_attribute("geovalues", criteria)

        if line.split()[1:3] == ["MP2", "Results"] or line.split()[1:3] == [
            "(Local-)MP2",
            "Results",
        ]:
            # Serenity has no higher order than MP2 and cannot do geometry optimization with it,
            # but still may contain several MP2 energies in one file.
            if hasattr(self, "mpenergies"):
                self.logger.warning("Warning: Multiple MP2 energies in Serenity!")
            line = next(inputfile)
            # skip forward to string "Total Energy", but only for max 20 lines
            i = 0
            while not line.strip().startswith("Total Energy") and i < 20:
                line = next(inputfile)
                i += 1
            self.append_attribute("mpenergies", [line.split()[2]])
            self.metadata["methods"].append("MP2")

        if line.strip().startswith("Polarizability Tensor / a.u.:"):
            self.skip_line(inputfile, "Polarizability")
            self.skip_line(inputfile, ["dashes"])
            line = next(inputfile)
            polarizability = []
            polarizability.append([float(i) for i in line.split()[-3:]])
            line = next(inputfile)
            polarizability.append([float(i) for i in line.split()[-3:]])
            line = next(inputfile)
            polarizability.append([float(i) for i in line.split()[-3:]])
            self.append_attribute("polarizabilities", numpy.array(polarizability))

        ### metadata
        if line[4:34] == "Time taken for the entire run:":
            # TODO this condition is probably too straigthforward and will take some more testing.
            self.metadata["success"] = True

        # note: cpu time is not printed straightforwardly in Serenity.
        if line[4:23] == "Time taken for task":
            if "wall_time" not in self.metadata:
                self.metadata["wall_time"] = []
            if line.split()[6] == "s.":
                walltime = datetime.timedelta(seconds=float(line.split()[5]))
            elif line.split()[6] == "min.":
                walltime = datetime.timedelta(
                    minutes=float(line.split()[5].split(":")[0]),
                    seconds=float(line.split()[5].split(":")[1]),
                )
            self.metadata["wall_time"].append(walltime)

        if line.strip().startswith("Warning") or line.strip().startswith("WARNING"):
            if "warnings" not in self.metadata:
                self.metadata["warnings"] = []
            # TODO just adding the entire warning line for now. Warnings may be longer than this line.
            self.metadata["warning"].append(line)

        if line[4:23] == "Version           :":
            self.metadata["package_version"] = line.split()[2]

        ### EXCITED STATE
        # TODO add CC2 etc
        if line.strip().startswith("TDDFT Summary"):
            self.metadata["excited_states_method"] = "TD-DFT"

        if line.strip().startswith("TDA Summary"):
            # TODO need to differentiate between TDA and CIS depending on DFT or HF for ground state
            self.metadata["excited_states_method"] = "TDA"

        if line.strip().startswith("CC2 Summary"):
            self.metadata["excited_states_method"] = "CC2"

        if line.strip().startswith("ADC(2) Summary"):
            self.metadata["excited_states_method"] = "ADC2"

        if line.strip().startswith("CIS(D) Summary"):
            self.metadata["excited_states_method"] = "CIS(D)"

        # excitation energies and singly-excited configuration data
        if line.strip().startswith("Dominant Contributions"):
            self.skip_line(inputfile, ["Dominant"])
            self.skip_line(inputfile, ["dashes"])
            self.skip_line(inputfile, ["state"])
            self.skip_line(inputfile, ["(a.u.)"])
            line = next(inputfile)

            # TODO maybe modify the method of warning. having multiple excited state calculations  obscures the metadata
            if hasattr(self, "etenergies"):
                self.logger.warning("Warning: Multiple Excited state calculations in Serenity!")

            exc_iterator = 1
            transition_data = []
            while not line.strip().startswith("--"):
                if line.strip().startswith(str(exc_iterator)):  # new exc
                    if exc_iterator > 1:
                        self.append_attribute("etsecs", transition_data)
                        transition_data = []
                    self.append_attribute("etenergies", line.split()[1])
                    line_data = line.split()
                    i = int(line_data[6])
                    a = int(line_data[8])
                    i_spin = self.convert_to_spin(line_data[7])
                    a_spin = self.convert_to_spin(line_data[9])
                    # Serenity prints coef as coef^2 * 100
                    coef = numpy.sqrt(float(line_data[10]) / 100.00)
                    transition_data.append([(i, i_spin), (a, a_spin), coef])
                    exc_iterator += 1
                else:
                    line_data = line.split()
                    i = int(line_data[1])
                    a = int(line_data[3])
                    i_spin = self.convert_to_spin(line_data[2])
                    a_spin = self.convert_to_spin(line_data[4])
                    coef = numpy.sqrt(float(line_data[5]) / 100.00)
                    transition_data.append([(i, i_spin), (a, a_spin), coef])
                line = next(inputfile)
            self.append_attribute("etsecs", transition_data)

        # oscillator strengths and transition dipoles (length gauge)
        if line.strip().startswith("Absorption Spectrum (dipole-length)"):
            self.skip_line(inputfile, ["Absorption"])
            self.skip_line(inputfile, ["dashes"])
            self.skip_line(inputfile, ["state"])
            self.skip_line(inputfile, ["(eV)"])
            line = next(inputfile)
            while not line.strip().startswith("--"):
                line_data = line.split()
                x, y, z = map(float, line_data[4:])
                self.append_attribute("etdips", [x, y, z])
                self.append_attribute("etoscs", line_data[3])
                line = next(inputfile)

        # transition dipoles (velocity gauge)
        if line.strip().startswith("Absorption Spectrum (dipole-velocity)"):
            self.skip_line(inputfile, ["Absorption"])
            self.skip_line(inputfile, ["dashes"])
            self.skip_line(inputfile, ["state"])
            self.skip_line(inputfile, ["(eV)"])
            line = next(inputfile)
            while not line.strip().startswith("--"):
                line_data = line.split()
                x, y, z = map(float, line_data[4:])
                self.append_attribute("etveldips", [x, y, z])
                line = next(inputfile)

        # rotatory strengths and magnetic transition dipoles (length gauge)
        if line.strip().startswith("CD Spectrum (dipole-length)"):
            self.skip_line(inputfile, ["CD"])
            self.skip_line(inputfile, ["dashes"])
            self.skip_line(inputfile, ["state"])
            self.skip_line(inputfile, ["(eV)"])
            line = next(inputfile)
            while not line.strip().startswith("--"):
                line_data = line.split()
                x, y, z = map(float, line_data[4:])
                self.append_attribute("etmagdips", [x, y, z])
                self.append_attribute("etrotats", line_data[3])
                line = next(inputfile)
