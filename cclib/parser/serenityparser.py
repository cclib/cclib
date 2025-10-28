# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Serenity output files"""

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
        self.unrestricted = False
        self.path = Path(self.inputfile.filenames[0]).resolve()

    def after_parsing(self):
        # Get molecular orbital information
        orbpath = self.path.parent / self.systemname / f"{self.systemname}.orbs.res.h5"
        if orbpath.is_file():
            assert utils.find_package("h5py"), (
                "h5py is needed to read in molecular orbital info from Serenity."
            )
            import h5py

            with h5py.File(orbpath, "r") as orbfile:
                coeffs = [orbfile["coefficients"][:]]
                eigenvalues = [orbfile["eigenvalues"][:].flatten()]
                self.set_attribute("moenergies", eigenvalues)
                self.set_attribute("mocoeffs", coeffs)
                self.set_attribute("nmo", len(eigenvalues[0]))

        super().after_parsing()

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # Extract system name
        if line.strip().startswith("------------------------------------------------------------"):
            line = next(inputfile)
            if line.strip().startswith("System"):
                name = line.split()[1]
                if not hasattr(self, "systemname"):
                    self.systemname = name
                elif name != self.systemname:
                    self.logger.error("Multiple systems detected which is not currently supported.")
                    raise StopParsing()

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
                if not self.optstatus[-1] & data.ccData.OPT_DONE and not len(self.optstatus) == 1:
                    self.append_attribute("atomcoords", coords)
            else:
                self.append_attribute("atomcoords", coords)

        if line[5:21] == "Basis Functions:":
            self.set_attribute("nbasis", int(line.split()[2]))

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

        if "Total Energy (" in line:
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
            dipoleRaw = [x, y, z]
            dipole = [utils.convertor(value, "ebohr", "Debye") for value in dipoleRaw]
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
            quadrupoleRaw = [xx, xy, xz, yy, yz, zz]
            quadrupole = [utils.convertor(value, "ebohr2", "Buckingham") for value in quadrupoleRaw]
            self.append_attribute("moments", quadrupole)

        if "Cycle" in line and "Mode" in line:
            line = next(inputfile)
            values = []
            while not line.strip().startswith("Converged after"):
                linedata = line.split()
                c1, c2, c3 = map(float, linedata[2:5])
                values.append([c1, c2, c3])
                line = next(inputfile)
            self.append_attribute("scfvalues", numpy.vstack(numpy.array(values)))

        if line.strip().startswith("Dispersion Correction ("):
            self.append_attribute("dispersionenergies", float(line.split()[3]))

        if "Total Local-CCSD Energy" in line:
            self.set_attribute("ccenergies", float(line.split()[3]))
            self.metadata["methods"].append("Local CCSD")
        if "Total Local-CCSD(T0) Energy" in line:
            self.set_attribute("ccenergies", float(line.split()[3]))
            self.metadata["methods"].append("Local CCSD(T0)")
        if "Total CCSD Energy" in line:
            self.set_attribute("ccenergies", float(line.split()[3]))
            self.metadata["methods"].append("CCSD")
        if "Total CCSD(T) Energy" in line:
            self.set_attribute("ccenergies", float(line.split()[3]))
            self.metadata["methods"].append("CCSD(T)")

        # Extract index of HOMO
        if line.strip().startswith("Orbital Energies:"):
            self.skip_line(inputfile, ["Orbital"])
            self.skip_line(inputfile, ["dashes"])
            self.skip_line(inputfile, ["#   Occ."])
            # self.skip_lines(inputfile, ["Orbital","dashes","#   Occ."]) # TODO test results in warnings
            homos = None
            line = next(inputfile)
            while line.split()[1] == "2.00":
                homos = int(line.split()[0])
                line = next(inputfile)
            self.set_attribute("homos", [homos - 1])  # Serenity starts at 1, python at 0

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
