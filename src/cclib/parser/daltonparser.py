# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2015, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Parser for DALTON output files"""


from __future__ import print_function

import numpy


from . import logfileparser
from . import utils


class DALTON(logfileparser.Logfile):
    """A DALTON log file."""

    # Used to index self.scftargets[].
    SCFRMS, SCFMAX, SCFENERGY = list(range(3))

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(DALTON, self).__init__(logname="DALTON", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "DALTON log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'DALTON("%s")' % (self.filename)

    def normalisesym(self, label):
        """Normalise the symmetries used by DALTON."""

        # it appears that DALTON is using the correct labels
        return label

    def before_parsing(self):

        self.firststdorient = True # Used to decide whether to wipe the atomcoords clean
        self.scftype = "none" # Type of SCF calculation: BLYP, RHF, ROHF, etc.

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # Since DALTON sometimes uses symmetry labels (Ag, Au, etc.) and
        # sometimes the symmetry group index, we need to parse and keep
        # a mapping between these two for later.
        #
        #  Symmetry Orbitals
        #  -----------------
        #
        #  Number of orbitals in each symmetry:          25    5   25    5
        #
        #
        #  Symmetry  Ag ( 1)
        #
        #    1     C        1s         1 +    2
        #    2     C        1s         3 +    4
        # ...
        #
        if line.strip() == "Symmetry Orbitals":

            self.skip_lines(inputfile, ['d', 'b'])

            line = inputfile.next()
            self.symcounts = [int(c) for c in line.split()[-4:]]

            self.symlabels = []
            for sc in self.symcounts:

                self.skip_lines(inputfile, ['b', 'b'])

                line = inputfile.next()
                assert line.split()[0] == "Symmetry"
                self.symlabels.append(line.split()[1])

                self.skip_line(inputfile, 'blank')
                for i in range(sc):
                    orbital = inputfile.next()

        # -------------------------------------------------
        if "Final DFT energy" in line or "Final HF energy" in line:
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            temp = line.split()
            self.scfenergies.append(utils.convertor(float(temp[-1]), "hartree", "eV"))

        # -------------------------------------------------
        if not hasattr(self, "natom") and "Total number of atoms:" in line:
            temp = line.split()
            self.set_attribute("natom", int(temp[-1]))

        # -------------------------------------------------
        if not hasattr(self, "charge") and "@    Total charge of the molecule" in line:
            temp = line.split()
            self.set_attribute("charge", int(temp[-1]))

        # -------------------------------------------------
        if not hasattr(self, "mult") and "@    Spin multiplicity and 2 M_S                1         0" in line:
            temp = line.split()
            self.set_attribute("mult", int(temp[-2]))

        # -------------------------------------------------
        if not hasattr(self, "nbasis") and "Total number of orbitals" in line:
            temp = line.split()
            self.set_attribute("nbasis", int(temp[4]))

        # -------------------------------------------------
        # SCF target threshold
        if "Threshold for SCF convergence" in line:
            if not hasattr(self, "scftargets"):
                self.scftargets = []
            temp = line.split()
            scftarget = self.float(temp[-1])
            self.scftargets.append([scftarget])

        # -------------------------------------------------
        # starting the SCF iterations
        #
        # with and without symmetry, the "Total energy" line is shifted a little.
        if "Iter" in line and "Total energy" in line:
            iteration = 0
            converged = False
            values = []
            if not hasattr(self, "scfvalues"):
                self.scfvalues = []

            while not converged:
                line = next(inputfile)
                # each iteration is bracketed by "-------------"
                if "-------------------" in line:
                    iteration += 1
                    continue

                # the first hit of @ n where n is the current iteration
                strcompare = "@{0:>3d}".format(iteration)
                if strcompare in line:
                    temp = line.split()
                    values.append([float(temp[2])])
                    #print(line.split())
                if line[0] == "@" and "converged in" in line:
                    converged = True

            self.scfvalues.append(values)

        if "@    Occupied SCF orbitals" in line and not hasattr(self, 'homos'):
            temp = line.split()
            homos = int(temp[4])
            self.set_attribute('homos', [homos-1]) # it is the index (python counting, so -1)

        # -------------------------------------------------
        # the molecular geometry requires the use of of
        # .RUN PROPERTIES
        if "Molecular geometry (au)" in line:
            if not hasattr(self, "atomcoords"):
                self.atomcoords = []

            if self.firststdorient:
                self.firststdorient = False

            line = next(inputfile)
            line = next(inputfile)
            #line = next(inputfile)

            atomcoords = []
            atomnos = []
            for i in range(self.natom):
                line = next(inputfile)
                temp = line.split()
                atomnos.append(self.table.number[temp[0]])

                # if symmetry has been enabled, extra labels are printed. if not, the list is one shorter
                coords = [1, 2, 3]
                try:
                    float(temp[1])
                except ValueError:
                    coords = [2, 3, 4]


                atomcoords.append([utils.convertor(float(temp[i]), "bohr", "Angstrom") for i in coords])
            self.atomcoords.append(atomcoords)
            self.set_attribute('atomnos', atomnos)

        # DALTON organizes the energies by symmetry, so we need to parser first,
        # and then sort the energies (and labels) before we store them.
        #
        # The formatting varies depending on RHF/DFT and/or version. Here is
        # an example from a DFT job:
        #
        #  *** SCF orbital energy analysis ***
        #
        #  Only the five lowest virtual orbital energies printed in each symmetry.
        #
        #  Number of electrons :   70
        #  Orbital occupations :   15    2   15    3
        #
        #  Sym       Kohn-Sham orbital energies
        #
        # 1 Ag    -10.01616533   -10.00394288   -10.00288640   -10.00209612    -9.98818062
        #          -0.80583154    -0.71422407    -0.58487249    -0.55551093    -0.50630125
        # ...
        #
        # Here is an example from an RHF job that only has symmetry group indices:
        #
        #  *** SCF orbital energy analysis ***
        #
        #  Only the five lowest virtual orbital energies printed in each symmetry.
        #
        #  Number of electrons :   70
        #  Orbital occupations :   15    2   15    3
        #
        #  Sym       Hartree-Fock orbital energies
        #
        #   1    -11.04052518   -11.03158921   -11.02882211   -11.02858563   -11.01747921
        #         -1.09029777    -0.97492511    -0.79988247    -0.76282547    -0.69677619
        # ...
        #
        if "*** SCF orbital energy analysis ***" in line:

            # to get ALL orbital energies, the .PRINTLEVELS keyword needs
            # to be at least 0,10 (up from 0,5). I know, obvious, right?
            # this, however, will conflict with the scfvalues output that
            # changes into some weird form of DIIS debug output.

            mosyms = []
            moenergies = []

            line = next(inputfile)
            # if the second line has "Only the five" present in it, then we have a reduced output in
            # the number of orbital energies.
            line = next(inputfile)  # has number of electrons if we need it
            if "Only the five" in line:
                line = next(inputfile)
                line = next(inputfile)
            
            line = next(inputfile)  # has # of orbital symmetries
            temp = line.split()
            nsym = len(temp[3:])

            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)

            # now parse nsym symmetries
            while nsym > 0:

                line = next(inputfile)
                temp = line.split()

                # there seems to be an errornous printout of a newline character in DALTON
                # if the number of orbital energies printed for a given symmetry is 5
                if len(temp) == 0:
                    continue

                # The first line has the orbital symmetry information, but sometimes
                # it's the label and sometimes it's the index. There are always five
                # energies per line, though, so we can deduce if we have the labels or
                # not just the index. In the latter case, we depend on the labels
                # being read earlier into the list `symlabels`.
                if len(temp) == 7:
                    sym = self.normalisesym(temp[1])
                    temp = [float(t) for t in temp[2:]]
                if len(temp) == 6:
                    sym = self.normalisesym(self.symlabels[int(temp[0])-1])
                    temp = [float(t) for t in temp[1:]]

                moenergies.extend(temp)
                mosyms.extend(len(temp)*[sym])
                while len(temp) > 0:
                    line = next(inputfile)
                    temp = [float(col) for col in line.split()]

                    if len(temp) == 0:
                        continue

                    moenergies.extend(temp)
                    mosyms.extend(len(temp)*[sym])

                nsym -= 1

            # now sort the data about energies and symmetries. see the following post for the magic
            # http://stackoverflow.com/questions/19339/a-transpose-unzip-function-in-python-inverse-of-zip
            sdata = sorted(zip(moenergies, mosyms), key=lambda x: x[0])
            moenergies, mosyms = zip(*sdata)

            self.moenergies = [[]]
            self.moenergies[0] = moenergies
            self.mosyms = [[]]
            self.mosyms[0] = mosyms

            if not hasattr(self, "nmo"):
                self.nmo = self.nbasis
                if len(self.moenergies[0]) != self.nmo:
                    self.set_attribute('nmo', len(self.moenergies[0]))

            #self.mocoeffs = [numpy.zeros((self.nmo, self.nbasis), "d")]

        # -------------------------------------------------
        # extract the center of mass line
        if "Center-of-mass coordinates (a.u.):" in line:
            temp = line.split()
            reference = [utils.convertor(float(temp[i]), "bohr", "Angstrom") for i in [3, 4, 5]]
            if not hasattr(self, 'moments'):
                self.moments = [reference]

        # -------------------------------------------------
        # Extract the dipole moment
        if "Dipole moment components" in line:
            dipole = numpy.zeros(3)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            if not "zero by symmetry" in line:
                line = next(inputfile)

                line = next(inputfile)
                temp = line.split()
                for i in range(3):
                    dipole[i] = float(temp[2])  # store the Debye value
            if hasattr(self, 'moments'):
                self.moments.append(dipole)


if __name__ == "__main__":
    import doctest, daltonparser, sys
    if len(sys.argv) == 1:
        doctest.testmod(daltonparser, verbose=False)

    if len(sys.argv) >= 2:
        parser = daltonparser.DALTON(sys.argv[1])
        data = parser.parse()

    if len(sys.argv) > 2:
        for i in range(len(sys.argv[2:])):
            if hasattr(data, sys.argv[2 + i]):
                print(getattr(data, sys.argv[2 + i]))
