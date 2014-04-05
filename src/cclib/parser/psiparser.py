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

    def before_parsing(self):

        # There are some major differences between the output of Psi3 and Psi4,
        # so it will be useful to register which one we are dealing with.
        self.version = None

        # This is just used to track which part of the output we are in for Psi4,
        # with changes triggered by ==> things like this <== (Psi3 does not have this)
        self.section = None

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

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # The version should always be detected.
        if "PSI3: An Open-Source Ab Initio" in line:
            self.version = 3
        if "PSI4: An Open-Source Ab Initio" in line:
            self.version = 4

        # This will automatically change the section attribute, when encountering
        # a line that <== looks like this ==>, to whatever is in between.
        if (line.strip()[:3] == "==>") and (line.strip()[-3:] == "<=="):
            self.section = line.strip()[4:-4]

        #  ==> Geometry <==
        #
        #    Molecular point group: c2h
        #    Full point group: C2h
        #
        #    Geometry (in Angstrom), charge = 0, multiplicity = 1:
        #
        #       Center              X                  Y                   Z       
        #    ------------   -----------------  -----------------  -----------------
        #           C         -1.415253322400     0.230221785400     0.000000000000
        #           C          1.415253322400    -0.230221785400     0.000000000000
        # ...
        #
        if (self.section == "Geometry") and ("Geometry (in Angstrom), charge" in line):

            assert line.split()[3] == "charge"
            charge = int(line.split()[5].strip(','))
            self.set_scalar('charge', charge)

            assert line.split()[6] == "multiplicity"
            mult = int(line.split()[8].strip(':'))
            self.set_scalar('mult', mult)

            self.skip_line(inputfile, "blank")
            line = next(inputfile)

            # Usually there is the header and dashes, but, for example, the coordinates
            # printed when a geometry optimization finishes do not have it.
            if line.split()[0] == "Center":
                self.skip_line(inputfile, "dashes")
                line = next(inputfile)

            coords = []
            while line.strip():
                el, x, y, z = line.split()
                coords.append([float(x), float(y), float(z)])
                line = next(inputfile)

            if not hasattr(self, 'atomcoords'):
                self.atomcoords = []
            self.atomcoords.append(coords)

        # Psi4 repeats the charge and multiplicity after the geometry.
        if (self.section == "Geometry") and (line[2:16].lower() == "charge       ="):
            charge = int(line.split()[-1])
            self.set_scalar('charge', charge)
        if (self.section == "Geometry") and (line[2:16].lower() == "multiplicity ="):
            mult = int(line.split()[-1])
            self.set_scalar('mult', mult)

        # In Psi3, the integrals program prints useful information when invoked.
        if (self.version == 3) and (line.strip() == "CINTS: An integrals program written in C"):

            self.skip_lines(inputfile, ['authors', 'd', 'b', 'b'])

            line = next(inputfile)
            assert line.strip() == "-OPTIONS:"
            while line.strip():
                line = next(inputfile)

            line = next(inputfile)
            assert line.strip() == "-CALCULATION CONSTANTS:"
            while line.strip():
                if "Number of atoms" in line:
                    natom = int(line.split()[-1])
                    self.set_scalar('natom', natom)
                if "Number of atomic orbitals" in line:
                    nbasis = int(line.split()[-1])
                    self.set_scalar('nbasis', nbasis)
                line = next(inputfile)

        # In Psi3, this part contains alot of important data pertaining to the SCF, but not only:
        if (self.version == 3) and (line.strip() == "CSCF3.0: An SCF program written in C"):

            self.skip_lines(inputfile, ['b', 'authors', 'b', 'd', 'b', 'mult', 'mult_comment', 'b'])

            line = next(inputfile)
            while line.strip():
                if line.split()[0] == "multiplicity":
                    mult = int(line.split()[-1])
                    self.set_scalar('mult', mult)
                if line.split()[0] == "charge":
                    charge = int(line.split()[-1])
                    self.set_scalar('charge', charge)
                if line.split()[0] == "convergence":
                    conv = float(line.split()[-1])
                line = next(inputfile)

            if not hasattr(self, 'scftargets'):
                self.scftargets = []
            self.scftargets.append([conv])

        # The printout for Psi4 has a more obvious trigger for the SCF parameter printout.
        if (self.section == "Algorithm") and (line.strip() == "==> Algorithm <=="):

            self.skip_line(inputfile, 'blank')

            line = next(inputfile)
            while line.strip():
                if "Energy threshold" in line:
                    etarget = float(line.split()[-1])
                if "Density threshold" in line:
                    dtarget = float(line.split()[-1])
                line = next(inputfile)

            if not hasattr(self, "scftargets"):
                self.scftargets = []
            self.scftargets.append([etarget, dtarget])

        # A block called 'Calculation Information' prints these before starting the SCF.
        if (self.section == "Pre-Iterations") and ("Number of atoms" in line):
            natom = int(line.split()[-1])
            self.set_scalar('natom', natom)
        if (self.section == "Pre-Iterations") and ("Number of atomic orbitals" in line):
            nbasis = int(line.split()[-1])
            self.set_scalar('nbasis', nbasis)

        #  ==> Iterations <==

        # Psi3 converges just the density elements, although it reports in the iterations
        # changes in the energy as well as the DIIS error.
        psi3_iterations_header = "iter       total energy        delta E         delta P          diiser"
        if (self.version == 3) and (line.strip() == psi3_iterations_header):

            if not hasattr(self, 'scfvalues'):
                self.scfvalues = []
            self.scfvalues.append([])

            line = next(inputfile)
            while line.strip():
                ddensity = float(line.split()[-2])
                self.scfvalues[-1].append([ddensity])
                line = next(inputfile)

        # Psi4 converges both the SCF energy and density elements and reports both in the
        # iterations printout. However, the default convergence scheme involves a density-fitted
        # algorithm for efficiency, and this is often followed by a something with exact electron
        # repulsion integrals. In that case, there are actually two convergence cycles performed,
        # one for the density-fitted algorithm and one for the exact one, and the iterations are
        # printed in two blocks separated by some set-up information.
        if (self.section == "Iterations") and (line.strip() == "==> Iterations <=="):

            if not hasattr(self, 'scfvalues'):
                self.scfvalues = []

            self.skip_line(inputfile, 'blank')
            header = next(inputfile)
            assert header.strip() == "Total Energy        Delta E     RMS |[F,P]|"

            scfvals = []
            self.skip_line(inputfile, 'blank')
            line = next(inputfile)
            while line.strip() != "==> Post-Iterations <==":
                if line.strip() and line.split()[0] in ["@DF-RHF", "@RHF", "@DF-RKS", "@RKS"]:
                    denergy = float(line.split()[4])
                    ddensity = float(line.split()[5])
                    scfvals.append([denergy, ddensity])
                line = next(inputfile)
            self.section = "Post-Iterations"
            self.scfvalues.append(scfvals)

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

            # If this is Psi4, we will be in the appropriate section.
            assert (self.version == 3) or (self.section == "Post-Iterations")

            self.moenergies = [[]]
            self.mosyms = [[]]

            # Psi4 has dashes under the trigger line, but Psi3 did not.
            if self.version == 4:
                self.skip_line(inputfile, 'dashes')
            self.skip_line(inputfile, 'blank')

            # Both versions have this case insensisitive substring.
            doubly = next(inputfile)
            assert "doubly occupied" in doubly.lower()

            # Psi4 now has a blank line, Psi3 does not.
            if self.version == 4:
                self.skip_line(inputfile, 'blank')

            line = next(inputfile)
            while line.strip():
                for i in range(len(line.split())//2):
                    self.mosyms[0].append(line.split()[i*2][-2:])
                    self.moenergies[0].append(line.split()[i*2+1])
                line = next(inputfile)

            # The last orbital energy here represented the HOMO.
            self.homos = [len(self.moenergies[0])-1]

            # Different numbers of blank lines in Psi3 and Psi4.
            if self.version == 3:
                self.skip_line(inputfile, 'blank')

            # The header for virtual orbitals is different for the two versions.
            unoccupied = next(inputfile)
            if self.version == 3:
                assert unoccupied.strip() == "Unoccupied orbitals"
            else:
                assert unoccupied.strip() == "Virtual:"

            # Psi4 now has a blank line, Psi3 does not.
            if self.version == 4:
                self.skip_line(inputfile, 'blank')

            line = next(inputfile)
            while line.strip():
                for i in range(len(line.split())//2):
                    self.mosyms[0].append(line.split()[i*2][-2:])
                    self.moenergies[0].append(line.split()[i*2+1])
                line = next(inputfile)

        # Both Psi3 and Psi4 print the final SCF energy right after the orbital energies,
        # but the label is different. Psi4 also does DFT, and the label is also different in that case.
        if (self.version == 3 and "* SCF total energy" in line) or \
           (self.section == "Post-Iterations" and ("@RHF Final Energy:" in line or "@RKS Final Energy" in line)):
            e = float(line.split()[-1])
            self.scfenergies = [utils.convertor(e, 'hartree', 'eV')]

        #  ==> Molecular Orbitals <==
        #
        #                 1            2            3            4            5
        #
        #    1    0.7014827    0.7015412    0.0096801    0.0100168    0.0016438
        #    2    0.0252630    0.0251793   -0.0037890   -0.0037346    0.0016447
        # ...
        #   59    0.0000133   -0.0000067    0.0000005   -0.0047455   -0.0047455
        #   60    0.0000133    0.0000067    0.0000005    0.0047455   -0.0047455
        #
        # Ene   -11.0288198  -11.0286067  -11.0285837  -11.0174766  -11.0174764
        # Sym            Ag           Bu           Ag           Bu           Ag
        # Occ             2            2            2            2            2
        #
        #
        #                11           12           13           14           15
        #
        #    1    0.1066946    0.1012709    0.0029709    0.0120562    0.1002765
        #    2   -0.2753689   -0.2708037   -0.0102079   -0.0329973   -0.2790813
        # ...
        #
        if (self.section == "Molecular Orbitals") and (line.strip() == "==> Molecular Orbitals <=="):

            self.skip_line(inputfile, 'blank')

            mocoeffs = []
            indices = next(inputfile)
            while indices.strip():

                indices = [int(i) for i in indices.split()]

                if len(mocoeffs) < indices[-1]:
                    for i in range(len(indices)):
                        mocoeffs.append([])
                else:
                    assert len(mocoeffs) == indices[-1]

                self.skip_line(inputfile, 'blank')

                line = next(inputfile)
                while line.strip():
                    iao = int(line.split()[0])
                    coeffs = [float(c) for c in line.split()[1:]]
                    for i,c in enumerate(coeffs):
                        mocoeffs[indices[i]-1].append(c)
                    line = next(inputfile)

                energies = next(inputfile)
                symmetries = next(inputfile)
                occupancies = next(inputfile)

                self.skip_lines(inputfile, ['b', 'b'])
                indices = next(inputfile)

            if not hasattr(self, 'mocoeffs'):
                self.mocoeffs = []
            self.mocoeffs.append(mocoeffs)

if __name__ == "__main__":
    import doctest, psiparser
    doctest.testmod(psiparser, verbose=False)
