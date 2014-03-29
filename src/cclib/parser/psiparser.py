# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2008-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.


import re

import numpy

from . import logfileparser
from . import utils


class Psi(logfileparser.Logfile):
    """A Psi log file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(Psi, self).__init__(logname="Psi", *args, **kwargs)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Psi log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Psi("%s")' % (self.filename)
    
    def normalisesym(self, label):
        """Use standard symmetry labels instead of NWChem labels.

        To normalise:
        (1) If label is one of [SG, PI, PHI, DLTA], replace by [sigma, pi, phi, delta]
        (2) replace any G or U by their lowercase equivalent

        >>> sym = NWChem("dummyfile").normalisesym
        >>> labels = ['A1', 'AG', 'A1G', "SG", "PI", "PHI", "DLTA", 'DLTU', 'SGG']
        >>> map(sym, labels)
        ['A1', 'Ag', 'A1g', 'sigma', 'pi', 'phi', 'delta', 'delta.u', 'sigma.g']
        """
        # FIXME if necessary
        return label

    def before_parsing(self):

        # Set any global variables for the parser here
        pass

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        #  ==> Geometry <==

        if "Geometry (in Angstrom)" in line:
            assert line.split()[3] == "charge"
            self.charge = int(line.split()[5].strip(','))
            assert line.split()[6] == "multiplicity"
            self.mult = int(line.split()[8].strip(':'))
            blank = next(inputfile)
            header = next(inputfile)
            dashes = next(inputfile)
            line = next(inputfile)
            if not hasattr(self, 'atomcoords'):
                self.atomcoords = []
            coords = []
            while line.strip():
                el, x, y, z = line.split()
                coords.append([float(x), float(y), float(z)])
                line = next(inputfile)
            self.atomcoords.append(coords)

        # At the beginning, in the geometry section, both Psi3 and Psi4 print
        # charge and multiplicity, differing in some lower/upper characters.
        if line[2:16].lower() == "charge       =":
            self.charge = int(line.split()[-1])
        if line[2:16].lower() == "multiplicity =":
            self.mult = int(line.split()[-1])

        #  ==> Algorithm <==

        # This is for Psi3, and might not be the only possible trigger for SCF parameters.
        if line.strip() == "CSCF3.0: An SCF program written in C":
            blank = next(inputfile)
            authors = next(inputfile)
            blank = next(inputfile)
            dashes = next(inputfile)
            blank = next(inputfile)
            mult = next(inputfile)
            mult_comment = next(inputfile)
            blank = next(inputfile)
            line = next(inputfile)
            while line.strip():
                if line.split()[0] == "convergence":
                    conv = float(line.split()[-1])
                line = next(inputfile)
            if not hasattr(self, 'scftargets'):
                self.scftargets = []
            self.scftargets.append([conv])

        # The printout for Psi4 has a more obvious trigger in this case.
        if line.strip() == "==> Algorithm <==":
            blank = next(inputfile)
            line = next(inputfile)
            if not hasattr(self, "scftargets"):
                self.scftargets = []
            while line.strip():
                if "Energy threshold" in line:
                    etarget = float(line.split()[-1])
                if "Density threshold" in line:
                    dtarget = float(line.split()[-1])
                line = next(inputfile)
            self.scftargets.append([etarget, dtarget])

        #  ==> Pre-Iterations <==

        # A block called 'Calculation Information' prints these before starting SCF.
        if "Number of atoms" in line:
            self.natom = int(line.split()[-1])
        if "Number of atomic orbitals" in line:
            self.nbasis = int(line.split()[-1])

        #  ==> Iterations <==

        # Psi3 converges just the density elements, although it reports in the iterations
        # changes in the energy as well as the DIIS error.
        if line.strip() == "iter       total energy        delta E         delta P          diiser":
            line = next(inputfile)
            if not hasattr(self, 'scfvalues'):
                self.scfvalues = []
            self.scfvalues.append([])
            while line.strip():
                ddensity = float(line.split()[-3])
                self.scfvalues[-1].append([ddensity])
                line = next(inputfile)

        # Psi4 converges both the SCF energy and density elements and reports both in the
        # iterations printout. However, the default convergence scheme involves a density-fitted
        # algorithm for efficiency, and this is often followed by a something with exact electron
        # repulsion integrals. In that case, there are actually two convergence cycles performed,
        # one for the density-fitted algorithm and one for the exact one, and the iterations are
        # printed in two blocks separated by some set-up information.
        if line.strip() == "==> Iterations <==":
            if not hasattr(self, 'scfvalues'):
                self.scfvalues = []
            blank = next(inputfile)
            header = next(inputfile)
            assert header.strip() == "Total Energy        Delta E     RMS |[F,P]|"
            blank = next(inputfile)
            line = next(inputfile)
            scfvalues = []
            while line.strip() != "==> Post-Iterations <==":
                if line.strip() and line.split()[0] in ["@DF-RHF", "@RHF", "@DF-RKS", "@RKS"]:
                    denergy = float(line.split()[4])
                    ddensity = float(line.split()[5])
                    scfvalues.append([denergy, ddensity])
                line = next(inputfile)
            self.scfvalues.append(scfvalues)

        #  ==> Post-Iterations <==

        # This section, from which we parse molecular orbital symmetries and
        # orbital energies, is quite similar for both Psi3 and Psi4, and in fact
        # the format for orbtials is the same, although the headers and spacers
        # are a bit different. Let's try to get both parsed with one code block.
        #
        # Here is how the block looks like for Psi4:
        #
        #	Orbital Energies (a.u.)
        #	-----------------------
        #
        #	Doubly Occupied:                                                      
        #
        #	   1Bu   -11.040586     1Ag   -11.040524     2Bu   -11.031589  
        #	   2Ag   -11.031589     3Bu   -11.028950     3Ag   -11.028820 
        # (...)
        #	  15Ag    -0.415620     1Bg    -0.376962     2Au    -0.315126  
        #	   2Bg    -0.278361     3Bg    -0.222189  
        #
        #	Virtual:                                                              
        #
        #	   3Au     0.198995     4Au     0.268517     4Bg     0.308826  
        #	   5Au     0.397078     5Bg     0.521759    16Ag     0.565017 
        # (...)
        #	  24Ag     0.990287    24Bu     1.027266    25Ag     1.107702  
        #	  25Bu     1.124938
        #
        # The case is different in the trigger string.
        if "orbital energies (a.u.)" in line.lower():

            self.moenergies = [[]]
            self.mosyms = [[]]

            # Psi4 has dashes under the trigger line, but Psi3 did not.
            blank_or_dashes = next(inputfile)
            if list(set(blank_or_dashes.strip())) == ["-"]:
                blank = next(inputfile)

            # Both versions have this case insensisitive substring.
            doubly = next(inputfile)
            assert "doubly occupied" in doubly.lower()

            # Psi4 now has a blank line, Psi3 does not.
            line = next(inputfile)
            if not line.strip():
                line = next(inputfile)

            while line.strip():
                for i in range(len(line.split())//2):
                    self.mosyms[0].append(line.split()[i*2][-2:])
                    self.moenergies[0].append(line.split()[i*2+1])
                line = next(inputfile)

            # The last orbital energy represented the HOMO.
            self.homos = [len(self.moenergies[0])-1]

            # Different numbers of blank lines in Psi3 and Psi4.
            line = next(inputfile)
            while not line.strip():
                line = next(inputfile)

            # The header for virtual orbitals is different for the two versions.
            assert (("unoccupied orbitals" in line.lower()) or ("virtual" in line.lower()))

            # Psi4 now has a blank line, Psi3 does not.
            line = next(inputfile)
            if not line.strip():
                line = next(inputfile)

            while line.strip():
                for i in range(len(line.split())//2):
                    self.mosyms[0].append(line.split()[i*2][-2:])
                    self.moenergies[0].append(line.split()[i*2+1])
                line = next(inputfile)

        # Both Psi3 and Psi4 print the final SCF energy right after the orbital energies,
        # but the label is different. Psi4 also does DFT, and here the label is also different.
        if "* SCF total energy" in line or "@RHF Final Energy:" in line or "@RKS Final Energy" in line:
            e = float(line.split()[-1])
            self.scfenergies = [utils.convertor(e, 'hartree', 'eV')]


if __name__ == "__main__":
    import doctest, psiparser
    doctest.testmod(psiparser, verbose=False)
