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
        if line.strip().startswith("Current Geometry (Angstrom):") and not getattr(
            self, "atomcoords", []
        ):
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
            self.append_attribute("atomcoords", coords)

        if line[5:21] == "Basis Functions:":
            self.set_attribute("nbasis", int(line.split()[2]))

        if "Cycle" in line and "Mode" in line:
            line = next(inputfile)
            values = []
            while not line.strip().startswith("Converged after"):
                linedata = line.split()
                c1, c2, c3 = map(float, linedata[2:5])
                values.append([c1, c2, c3])
                line = next(inputfile)
            self.append_attribute("scfvalues", numpy.vstack(numpy.array(values)))

        if "Dispersion Correction (" in line:
            self.append_attribute("dispersionenergies", float(line.split()[3]))

        # create dict

        populationtypes = ["CM5", "Mulliken", "Hirshfeld"]
        thisdict = {}
        if line.strip().startswith(tuple(populationtypes)) and line.split()[1] == "Population":
            key = line.split()[0]
            line = next(inputfile)
            self.skip_line(inputfile, "dashes")
            self.skip_line(inputfile, "blank")
            line = next(inputfile)
            line = next(inputfile)
            chargeList = []

            while line.strip():
                charge = float(line.split()[3])
                chargeList.append(charge)
                line = next(inputfile)

            value = numpy.array(chargeList)
            thisdict[key] = value
            print(thisdict)
        self.set_attribute("atomcharges", thisdict)
