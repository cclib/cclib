# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Serenity output files"""

from cclib.parser import logfileparser, utils


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

    def after_parsing(self):
        super().after_parsing()

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # Extract charge and multiplicity
        if line[5:11] == "Charge":
            self.set_attribute("charge", int(line.split()[1]))

        # Extract multiplicity
        if line[5:9] == "Spin":
            self.set_attribute("mult", int(line.split()[1]) + 1)

        if line[5:21] == "Basis Functions:":
            self.set_attribute("nbasis", int(line.split()[2]))

        # Extract from atoms: number of atoms, elements, and coordinates
        if line.strip().startswith("Origin chosen as:"):
            line = next(inputfile)
            origin_data = line.replace("(", "").replace(")", "").replace(",", "").split()
            x, y, z = map(float, origin_data)
            origin = [x, y, z]
            self.append_attribute("moments", origin)

        if line.strip().startswith("Dipole Moment:"):
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            dipole_data = line.split()
            x, y, z = map(float, dipole_data[:3])
            dipoleRaw = [x, y, z]
            dipole = [utils.convertor(value, "ebohr", "Debye") for value in dipoleRaw]
            self.append_attribute("moments", dipole)

        if line.strip().startswith("Quadrupole Moment:"):
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            q_data = line.split()
            xx, xy, xz = map(float, q_data[1:4])
            line = next(inputfile)
            q_data = line.split()
            yy, yz = map(float, q_data[2:4])
            line = next(inputfile)
            q_data = line.split()
            zz = float(q_data[3])
            quadrupoleRaw = [xx, xy, xz, yy, yz, zz]
            quadrupole = [utils.convertor(value, "ebohr2", "Buckingham") for value in quadrupoleRaw]
            self.append_attribute("moments", quadrupole)
