# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Serenity output files"""

from cclib.parser import logfileparser

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

        # Extract from atoms: number of atoms, elements, and coordinates
        if line.strip().startswith("Current Geometry (Angstrom):"):
            line = next(inputfile)
            line = next(inputfile)
            atomnos = []
            coords = []
            while line.strip():
                print(line)
                atominfo = line.split()
                element = atominfo[1]
                x, y, z = map(float, atominfo[2:5])
                atomnos.append(self.table.number[element])
                coords.append([x, y, z])
                line = next(inputfile)

            self.set_attribute("atomnos", atomnos)
            self.set_attribute("natom", len(atomnos))
            self.append_attribute("atomcoords", coords)

        if line[5:21] == "Basis Functions:":
            self.set_attribute("nbasis", int(line.split()[2]))

        # Extract SCF thresholds
        if line.strip().startswith("Energy Threshold:"):
            scftargets = []
            ethresh = float(line.split()[2])
            line = next(inputfile)
            if "RMSD[D]" in line:
                rmsd = float(line.split()[2])
                line = next(inputfile)
                if "DIIS" in line:
                    diis = float(line.split()[2])
                    scftargets.append(numpy.array([ethresh, rmsd, diis]))
                    self.set_attribute("scftargets", scftargets)
