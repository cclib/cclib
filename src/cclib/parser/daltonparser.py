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

        # It appears that DALTON is using the correct labels.
        return label

    def before_parsing(self):

        # Used to decide whether to wipe the atomcoords clean.
        self.firststdorient = True

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # This section is close to the beginning of the file, and can be used
        # to parse natom, nbasis and atomnos. Note that DALTON operates on the
        # idea of atom type, which are not necessarily unique.
        #
        #  Atoms and basis sets
        #  --------------------
        #
        #  Number of atom types :    6
        #  Total number of atoms:   20
        #
        #  Basis set used is "STO-3G" from the basis set library.
        #
        #  label    atoms   charge   prim   cont     basis
        #  ----------------------------------------------------------------------
        #  C           6    6.0000    15     5      [6s3p|2s1p]                                        
        #  H           4    1.0000     3     1      [3s|1s]                                            
        #  C           2    6.0000    15     5      [6s3p|2s1p]                                        
        #  H           2    1.0000     3     1      [3s|1s]                                            
        #  C           2    6.0000    15     5      [6s3p|2s1p]                                        
        #  H           4    1.0000     3     1      [3s|1s]                                            
        #  ----------------------------------------------------------------------
        #  total:     20   70.0000   180    60
        #  ----------------------------------------------------------------------
        #
        #  Threshold for neglecting AO integrals:  1.00D-12
        #
        if line.strip() == "Atoms and basis sets":

            self.skip_lines(inputfile, ['d', 'b'])

            line = next(inputfile)
            assert "Number of atom types" in line
            self.ntypes = int(line.split()[-1])

            line = next(inputfile)
            assert "Total number of atoms:" in line
            self.set_attribute("natom", int(line.split()[-1]))

            self.skip_lines(inputfile, ['b', 'basisname', 'b'])

            line = next(inputfile)
            cols = line.split()
            iatoms = cols.index('atoms')
            icharge = cols.index('charge')
            icont = cols.index('cont')

            self.skip_line(inputfile, 'dashes')

            atomnos = []
            atombasis = []
            nbasis = 0
            for itype in range(self.ntypes):
                line = next(inputfile)
                cols = line.split()
                atoms = int(cols[iatoms])
                charge = float(cols[icharge])
                assert int(charge) == charge
                charge = int(charge)
                cont = int(cols[icont])
                for at in range(atoms):
                    atomnos.append(charge)
                    atombasis.append(list(range(nbasis, nbasis + cont)))
                    nbasis += cont
            self.set_attribute('atomnos', atomnos)
            self.set_attribute('atombasis', atombasis)

        # Since DALTON sometimes uses symmetry labels (Ag, Au, etc.) and sometimes
        # just the symmetry group index, we need to parse and keep a mapping between
        # these two for later.
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

                # If the number of orbitals for a symmetry is zero, the printout
                # is different (see MP2 unittest logfile for an example).
                line = inputfile.next()
                if sc == 0:
                    assert "No orbitals in symmetry" in line
                else:
                    assert line.split()[0] == "Symmetry"
                    self.symlabels.append(line.split()[1])

                    self.skip_line(inputfile, 'blank')
                    for i in range(sc):
                        orbital = inputfile.next()

        #      Wave function specification
        #      ============================
        # @    Wave function type        >>> KS-DFT <<<
        # @    Number of closed shell electrons          70
        # @    Number of electrons in active shells       0
        # @    Total charge of the molecule               0
        #
        # @    Spin multiplicity and 2 M_S                1         0
        # @    Total number of symmetries                 4 (point group: C2h)
        # @    Reference state symmetry                   1 (irrep name : Ag )
        # 
        #     This is a DFT calculation of type: B3LYP
        # ...
        #
        if "@    Total charge of the molecule" in line:
            self.set_attribute("charge", int(line.split()[-1]))
        if "@    Spin multiplicity and 2 M_S                1         0" in line:
            self.set_attribute("mult", int(line.split()[-2]))

        #     Orbital specifications
        #     ======================
        #     Abelian symmetry species          All |    1    2    3    4
        #                                           |  Ag   Au   Bu   Bg 
        #                                       --- |  ---  ---  ---  ---
        #     Total number of orbitals           60 |   25    5   25    5
        #     Number of basis functions          60 |   25    5   25    5
        #
        #      ** Automatic occupation of RKS orbitals **
        #
        #      -- Initial occupation of symmetries is determined from extended Huckel guess.           
        #      -- Initial occupation of symmetries is :
        # @    Occupied SCF orbitals              35 |   15    2   15    3
        #
        #     Maximum number of Fock   iterations      0
        #     Maximum number of DIIS   iterations     60
        #     Maximum number of QC-SCF iterations     60
        #     Threshold for SCF convergence     1.00D-05
        #     This is a DFT calculation of type: B3LYP
        # ...
        #
        if "Total number of orbitals" in line:
            # DALTON 2015 adds a @ in front of number of orbitals
            index = 4
            if "@" in line:
                index = 5
            self.set_attribute("nbasis", int(line.split()[index]))
        if "@    Occupied SCF orbitals" in line and not hasattr(self, 'homos'):
            temp = line.split()
            homos = int(temp[4])
            self.set_attribute('homos', [homos-1]) # it is the index (python counting, so -1)
        if "Threshold for SCF convergence" in line:
            if not hasattr(self, "scftargets"):
                self.scftargets = []
            scftarget = self.float(line.split()[-1])
            self.scftargets.append([scftarget])

        #  *********************************************
        #  ***** DIIS optimization of Hartree-Fock *****
        #  *********************************************
        # 
        #  C1-DIIS algorithm; max error vectors =    8
        #
        #  Automatic occupation of symmetries with  70 electrons.
        #
        #  Iter     Total energy    Error norm  Delta(E)    SCF occupation
        #  -----------------------------------------------------------------------------
        #       K-S energy, electrons, error :    -46.547567739269  69.9999799123   -2.01D-05
        # @  1  -381.645762476       4.00D+00  -3.82D+02    15   2  15   3
        #       Virial theorem: -V/T =      2.008993
        # @      MULPOP C   _1  0.15; C   _2  0.15; C   _1  0.12; C   _2  0.12; C   _1  0.11; C   _2  0.11; H   _1 -0.15; H   _2 -0.15; H   _1 -0.14; H   _2 -0.14; 
        # @             C   _1  0.23; C   _2  0.23; H   _1 -0.15; H   _2 -0.15; C   _1  0.08; C   _2  0.08; H   _1 -0.12; H   _2 -0.12; H   _1 -0.13; H   _2 -0.13; 
        #  -----------------------------------------------------------------------------
        #       K-S energy, electrons, error :    -46.647668038900  69.9999810430   -1.90D-05
        # @  2  -381.949410128       1.05D+00  -3.04D-01    15   2  15   3
        #       Virial theorem: -V/T =      2.013393
        # ...
        #
        # Wwith and without symmetry, the "Total energy" line is shifted a little.
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

            self.skip_line(inputfile, 'blank')

            # if the second line has "Only the five" present in it, then we have a reduced
            # number of virtual orbital energies printed. If it does not, then it contains
            # the number of electrons.
            # DALTON 2015 increases the output to 20 virtuals, so only check
            # for a slightly smaller part of the string
            line = next(inputfile)
            if "Only the" in line:
                self.skip_line(inputfile, 'blank')
                line = next(inputfile)
            nelectrons = int(line.split()[-1])
            
            line = next(inputfile)
            occupations = [int(o) for o in line.split()[3:]]
            nsym = len(occupations)

            self.skip_lines(inputfile, ['b', 'header', 'b'])

            # now parse nsym symmetries
            for isym in range(nsym):

                # For unoccupied symmetries, nothing is printed here.
                if occupations[isym] == 0:
                    continue

                # When there are exactly five energies printed (on just one line), it seems
                # an extra blank line is printed after a block.
                line = next(inputfile)
                if not line.strip():
                    line = next(inputfile)
                cols = line.split()

                # The first line has the orbital symmetry information, but sometimes
                # it's the label and sometimes it's the index. There are always five
                # energies per line, though, so we can deduce if we have the labels or
                # not just the index. In the latter case, we depend on the labels
                # being read earlier into the list `symlabels`.
                if 'A' in cols[1] or 'B' in cols[1]:
                    sym = self.normalisesym(cols[1])
                    energies = [float(t) for t in cols[2:]]
                else:
                    sym = self.normalisesym(self.symlabels[int(cols[0]) - 1])
                    energies = [float(t) for t in cols[1:]]

                while len(energies) > 0:
                    moenergies.extend(energies)
                    mosyms.extend(len(energies)*[sym])
                    line = next(inputfile)
                    energies = [float(col) for col in line.split()]

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

        #                       .-----------------------------------.
        #                       | >>> Final results from SIRIUS <<< |
        #                       `-----------------------------------'
        #
        #
        # @    Spin multiplicity:           1
        # @    Spatial symmetry:            1 ( irrep  Ag  in C2h )
        # @    Total charge of molecule:    0
        #
        # @    Final DFT energy:           -382.050716652387                 
        # @    Nuclear repulsion:           445.936979976608
        # @    Electronic energy:          -827.987696628995
        #
        # @    Final gradient norm:           0.000003746706
        # ...
        #
        if "Final DFT energy" in line or "Final HF energy" in line:
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            temp = line.split()
            self.scfenergies.append(utils.convertor(float(temp[-1]), "hartree", "eV"))

        if "@   = MP2 second order energy" in line:
            energ = utils.convertor(float(line.split()[-1]), 'hartree', 'eV')
            if not hasattr(self, "mpenergies"):
                self.mpenergies = []
            self.mpenergies.append([])
            self.mpenergies[-1].append(energ)

        if "Total energy CCSD(T)" in line:
            energ = utils.convertor(float(line.split()[-1]), 'hartree', 'eV')
            if not hasattr(self, "ccenergies"):
                self.ccenergies = []
            self.ccenergies.append(energ)

        # The molecular geometry requires the use of .RUN PROPERTIES in the input.
        # Note that the second column is not the nuclear charge, but the atom type
        # index used internally by DALTON.
        #
        #                             Molecular geometry (au)
        #                             -----------------------
        #
        # C   _1     1.3498778652            2.3494125195            0.0000000000
        # C   _2    -1.3498778652           -2.3494125195            0.0000000000
        # C   _1     2.6543517307            0.0000000000            0.0000000000
        # ...
        #
        if "Molecular geometry (au)" in line:
            if not hasattr(self, "atomcoords"):
                self.atomcoords = []

            if self.firststdorient:
                self.firststdorient = False

            line = next(inputfile)
            line = next(inputfile)
            #line = next(inputfile)

            atomcoords = []
            for i in range(self.natom):
                line = next(inputfile)
                temp = line.split()

                # if symmetry has been enabled, extra labels are printed. if not, the list is one shorter
                coords = [1, 2, 3]
                try:
                    float(temp[1])
                except ValueError:
                    coords = [2, 3, 4]


                atomcoords.append([utils.convertor(float(temp[i]), "bohr", "Angstrom") for i in coords])
            self.atomcoords.append(atomcoords)

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
