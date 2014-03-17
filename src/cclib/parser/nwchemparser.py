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


class NWChem(logfileparser.Logfile):
    """An NWChem log file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(NWChem, self).__init__(logname="NWChem", *args, **kwargs)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "NWChem log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'NWChem("%s")' % (self.filename)
    
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

    def set_scalar(self, name, value, check=True):
        if hasattr(self, name):
            if check:
                assert getattr(self, name) == value
        else:
            setattr(self, name, value)

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # This is printed in the input module, so should always be the first coordinates,
        # and contains some basic information we want to parse as well. However, this is not
        # the only place where the coordinates are printed during geometry optimization,
        # since the gradients module has a separate coordinate printout, which happens
        # alongside the coordinate gradients. This geometry printout happens at the
        # beginning of each optimization step only.
        if line.strip() == 'Geometry "geometry" -> ""' or line.strip() == 'Geometry "geometry" -> "geometry"':

            dashes = next(inputfile)
            blank = next(inputfile)
            concerning_units = next(inputfile)
            blank = next(inputfile)
            header = next(inputfile)
            dashes = next(inputfile)

            if not hasattr(self, 'atomcoords'):
                self.atomcoords = []

            line = next(inputfile)
            coords = []
            atomnos = []
            while line.strip():
                # The column labeled 'tag' is usually empty, but I'm not sure whether it can have spaces,
                # so for now assume that it can and that there will be seven columns in that case.
                if len(line.split()) == 6:
                    index, atomname, nuclear, x, y, z = line.split()
                else:
                    index, atomname, tag, nuclear, x, y, z = line.split()
                coords.append(list(map(float, [x,y,z])))
                atomnos.append(int(float(nuclear)))
                line = next(inputfile)
            self.atomcoords.append(coords)
            if hasattr(self, 'atomnos'):
                assert atomnos == self.atomnos
            else:
                self.atomnos = atomnos

        # If the geometry is printed in XYZ format, it will have the number of atoms.
        if line[12:31] == "XYZ format geometry":

            dashes = next(inputfile)
            natom = int(next(inputfile).strip())
            if hasattr(self, "natom"):
                assert self.natom == natom
            else:
                self.natom = natom

        # NWChem does not normally print the basis set for each atom, but rather
        # chooses the concise option of printing Gaussian coefficients for each
        # atom type/element only once. Therefore, we need to first parse those
        # coefficients and afterwards build the appropriate gbasis attribute based
        # on that and atom types/elements already parsed (atomnos).
        if line.strip() == """Basis "ao basis" -> "ao basis" (cartesian)""":
            dashes = next(inputfile)
            gbasis_dict = {}
            line = next(inputfile)
            while line.strip():
                atomtype = line.split()[0]
                gbasis_dict[atomtype] = []
                dashes = next(inputfile)
                labels = next(inputfile)
                dashes = next(inputfile)
                shells = []
                line = next(inputfile)
                while line.strip() and line.split()[0].isdigit():
                    shell = None
                    while line.strip():
                        nshell, type, exp, coeff = line.split()
                        nshell = int(nshell)
                        assert len(shells) == nshell - 1
                        if not shell:
                            shell = (type, [])
                        else:
                            assert shell[0] == type
                        exp = float(exp)
                        coeff = float(coeff)
                        shell[1].append((exp,coeff))
                        line = next(inputfile)
                    shells.append(shell)
                    line = next(inputfile)
                gbasis_dict[atomtype].extend(shells)
            gbasis = []
            for i in range(self.natom):
                atomtype = utils.PeriodicTable().element[self.atomnos[i]]
                gbasis.append(gbasis_dict[atomtype])
            if not hasattr(self, 'gbasis'):
                self.gbasis = gbasis
            else:
                assert self.gbasis == gbasis

        # Normally the indexes of AOs assigned to specific atoms are also not printed,
        # so we need to infer that. We could do that from the previous section,
        # it might be worthwhile to take numbers from two different places, hence
        # the code below, which builds atombasis based on the number of functions
        # listed in this summary of the AO basis.
        if line.strip() == """Summary of "ao basis" -> "ao basis" (cartesian)""":
            dashes = next(inputfile)
            headers = next(inputfile)
            dashes = next(inputfile)
            atombasis_dict = {}
            line = next(inputfile)
            while line.strip():
                atomtype, desc, shells, funcs, types = line.split()
                atombasis_dict[atomtype] = int(funcs)
                line = next(inputfile)
            atombasis = []
            last = 0
            for i in range(self.natom):
                atomtype = utils.PeriodicTable().element[self.atomnos[i]]
                nfuncs = atombasis_dict[atomtype]
                atombasis.append(list(range(last,last+nfuncs)))
                last = atombasis[-1][-1] + 1
            if not hasattr(self, 'atombasis'):
                self.atombasis = atombasis
            else:
                assert self.atombasis == atombasis

        # This section contains general parameters for Hartree-Fock calculations.
        if line.strip() == "NWChem SCF Module":
            dashes = next(inputfile)
            blank = next(inputfile)
            blank = next(inputfile)
            title = next(inputfile)
            blank = next(inputfile)
            blank = next(inputfile)
            blank = next(inputfile)
            line = next(inputfile)
            while line.strip():
                if line[2:8] == "charge":
                    self.charge = int(float(line.split()[-1]))
                if line[2:13] == "open shells":
                    unpaired = int(line.split()[-1])
                    self.mult = 2*unpaired + 1
                if line[2:7] == "atoms":
                    natom = int(line.split()[-1])
                    if hasattr(self, 'natom'):
                        assert self.natom == natom
                    else:
                        self.natom = natom
                if line[2:11] == "functions":
                    functions = int(line.split()[-1])
                    self.nbasis = functions
                line = next(inputfile)

        # This section contains general parameters for DFT calculations, including
        # SCF convergence targets.
        if line.strip() == "General Information":
            while line.strip():

                if "No. of atoms" in line:
                    self.set_scalar('natom', int(line.split()[-1]))
                if "Charge" in line:
                    self.set_scalar('charge', int(line.split()[-1]))
                if "Spin multiplicity" in line:
                    self.set_scalar('mult', int(line.split()[-1]))

                if "Convergence on energy requested" in line:
                    target_energy = float(line.split()[-1].replace('D', 'E'))
                if "Convergence on density requested" in line:
                    target_density = float(line.split()[-1].replace('D', 'E'))
                if "Convergence on gradient requested" in line:
                    target_gradient = float(line.split()[-1].replace('D', 'E'))

                line = next(inputfile)

            if not hasattr(self, 'scftargets'):
                self.scftargets = []
            self.scftargets.append([target_energy, target_density, target_gradient])

        # The default (only?) SCF algorithm for Hartree-Fock is a preconditioned conjugate
        # gradient method that apparently "always" converges, so this header should reliably
        # signal a start of the SCF cycle. The convergence targets are also printed here.
        if line.strip() == "Quadratically convergent ROHF":

            while not "Final" in line:

                # Only the norm of the orbital gradient is used to test convergence.
                if line[:22] == " Convergence threshold":
                    target = float(line.split()[-1])
                    if not hasattr(self, "scftargets"):
                        self.scftargets = []
                    self.scftargets.append([target])

                    # This is critical for the stop condition of the section,
                    # because the 'Final Fock-matrix accuracy' is along the way.
                    # It would be prudent to find a more robust stop condition.
                    while list(set(line.strip())) != ["-"]:
                        line = next(inputfile)

                if line.split() == ['iter', 'energy', 'gnorm', 'gmax', 'time']:
                    values = []
                    dashes = next(inputfile)
                    line = next(inputfile)
                    while line.strip():
                        iter,energy,gnorm,gmax,time = line.split()
                        gnorm = float(gnorm.replace('D','E'))
                        values.append([gnorm])
                        line = next(inputfile)
                    if not hasattr(self, 'scfvalues'):
                        self.scfvalues = []
                    self.scfvalues.append(values)

                line = next(inputfile)

        # The SCF for DFT does not use the same algorithm as Hartree-Fock, but always
        # seems to use the following format to report SCF convergence:
        #   convergence    iter        energy       DeltaE   RMS-Dens  Diis-err    time
        # ---------------- ----- ----------------- --------- --------- ---------  ------
        # d= 0,ls=0.0,diis     1   -382.2544324446 -8.28D+02  1.42D-02  3.78D-01    23.2
        # d= 0,ls=0.0,diis     2   -382.3017298534 -4.73D-02  6.99D-03  3.82D-02    39.3
        # d= 0,ls=0.0,diis     3   -382.2954343173  6.30D-03  4.21D-03  7.95D-02    55.3
        # ...
        if line.split() == ['convergence', 'iter', 'energy', 'DeltaE', 'RMS-Dens', 'Diis-err', 'time']:
            dashes = next(inputfile)
            line = next(inputfile)
            values = []
            while line.strip():

                # Sometimes there are things in between iterations with fewer columns,
                # and we want to skip those lines, most probably. An exception might
                # unrestricted calcualtions, which show extra RMS density and DIIS
                # errors, although it is not clear yet whether these are for the
                # beta orbitals or somethine else. The iterations look like this in that case:
                #   convergence    iter        energy       DeltaE   RMS-Dens  Diis-err    time
                # ---------------- ----- ----------------- --------- --------- ---------  ------
                # d= 0,ls=0.0,diis     1   -382.0243202601 -8.28D+02  7.77D-03  1.04D-01    30.0
                #                                                     7.68D-03  1.02D-01
                # d= 0,ls=0.0,diis     2   -382.0647539758 -4.04D-02  4.64D-03  1.95D-02    59.2
                #                                                     5.39D-03  2.36D-02
                # ...
                if len(line[17:].split()) == 6:
                    iter, energy, deltaE, dens, diis, time = line[17:].split()
                    val_energy = float(deltaE.replace('D', 'E'))
                    val_density = float(dens.replace('D', 'E'))
                    val_gradient = float(diis.replace('D', 'E'))
                    values.append([val_energy, val_density, val_gradient])

                line = next(inputfile)

            if not hasattr(self, 'scfvalues'):
                self.scfvalues = []

            self.scfvalues.append(values)

        # The line containing the final SCF energy seems to be always identifiable like this.
        if "Total SCF energy" in line or "Total DFT energy" in line:
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            energy = float(line.split()[-1])
            energy = utils.convertor(energy, "hartree", "eV")
            self.scfenergies.append(energy)

        # The same format is used for HF and DFT molecular orbital analysis. We want to parse
        # the MO energies from this section, although it is printed already before this with
        # less precision (might be useful to parse that if this is not available). Also, this
        # section contains coefficients for the leading AO contributions, so it would also
        # be good to parse and use those values if the full vectors are not printed.
        if "Final Molecular Orbital Analysis" in line:
            if not hasattr(self, "moenergies"):
                self.moenergies = []
            dashes = next(inputfile)
            blank = next(inputfile)

            energies = []
            line = next(inputfile)
            homo = 0
            while line[:7] == " Vector":

                nvector = int(line[7:12])
                if "Occ=2.0" in line:
                    homo = nvector-1

                # If the printout does not start from the first MO, assume None for all previous orbitals.
                if len(energies) == 0 and nvector > 1:
                    for i in range(1,nvector):
                        energies.append(None)

                energy = float(line[34:].replace('D','E'))
                energy = utils.convertor(energy, "hartree", "eV")
                energies.append(energy)
                line = next(inputfile)
                if "MO Center" in line:
                    line = next(inputfile)
                if "Bfn." in line:
                    line = next(inputfile)
                if "-----" in line:
                    line = next(inputfile)
                while line.strip():
                    line = next(inputfile)
                line = next(inputfile)

            self.moenergies.append(energies)
            if hasattr(self, 'nmo'):
                assert self.nmo == nvector
            else:
                self.nmo = nvector

            if not hasattr(self, 'homos'):
                self.homos = []
            self.homos.append(homo)

        # This is where the full MO vectors are printed, but a special directive is needed for it.
        if line.strip() == "Final MO vectors":

            dashes = next(inputfile)
            blank = next(inputfile)
            blank = next(inputfile)

            # The columns are MOs, rows AOs, but I'm guessing the order of this array since no
            # atom information is printed alongside the indices.
            array_info = next(inputfile)
            size = array_info.split('[')[1].split(']')[0]
            nbasis = int(size.split(',')[0].split(':')[1])
            nmo = int(size.split(',')[1].split(':')[1])
            if hasattr(self, 'nbasis'):
                assert self.nbasis == nbasis
            else:
                self.nbasis = nbasis
            if hasattr(self, 'nmo'):
                assert self.nmo == nmo
            else:
                self.nmo = nmo
            
            blank = next(inputfile)
            mocoeffs = []
            while len(mocoeffs) < self.nmo:
                nmos = list(map(int,next(inputfile).split()))
                assert len(mocoeffs) == nmos[0]-1
                for n in nmos:
                    mocoeffs.append([])
                dashes = next(inputfile)
                for nb in range(nbasis):                
                    line = next(inputfile)
                    index = int(line.split()[0])
                    assert index == nb+1
                    coefficients = list(map(float,line.split()[1:]))
                    assert len(coefficients) == len(nmos)
                    for i,c in enumerate(coefficients):
                        mocoeffs[nmos[i]-1].append(c)
                blank = next(inputfile)
            self.mocoeffs = [mocoeffs]

        # The atomic Mulliken charges are typically printed like this:
        #  Mulliken analysis of the total density
        #  --------------------------------------
        #
        #    Atom       Charge   Shell Charges
        # -----------   ------   -------------------------------------------------------
        #    1 C    6     6.00   1.99  1.14  2.87
        #    2 C    6     6.00   1.99  1.14  2.87
        # ...
        if line.strip() == "Mulliken analysis of the total density":

            if not hasattr(self, "atomcharges"):
                self.atomcharges = {}

            dashes = next(inputfile)
            blank = next(inputfile)
            header = next(inputfile)
            dashes = next(inputfile)

            charges = []
            line = next(inputfile)
            while line.strip():
                index, atomname, nuclear, atom = line.split()[:4]
                shells = line.split()[4:]
                charges.append(float(atom)-float(nuclear))
                line = next(inputfile)
            self.atomcharges['mulliken'] = charges

        # If the full overlap matrix is printed, it looks like this:
        #          ----------------------------
        #          Mulliken population analysis
        #          ----------------------------
        #
        #          ----- Total      overlap population -----
        #
        #                               1              2              3              4              5              6              7
        #
        #    1   1 C  s            2.0694818227  -0.0535883400  -0.0000000000  -0.0000000000  -0.0000000000  -0.0000000000   0.0000039991
        #    2   1 C  s           -0.0535883400   0.8281341291   0.0000000000  -0.0000000000   0.0000000000   0.0000039991  -0.0009906747
        # ...
        if line.strip() == "Mulliken population analysis":

            dashes = next(inputfile)
            blank = next(inputfile)
            total_overlap_population = next(inputfile)
            blank = next(inputfile)

            overlaps = []
            line= next(inputfile)
            while all([c.isdigit() for c in line.split()]):

                # There is always a line with the MO indices printed in thie block.
                indices = [int(i)-1 for i in line.split()]
                for i in indices:
                    overlaps.append([])

                # There is usually a blank line after the MO indices, but
                # there are exceptions, so check if line is blank first.
                line = next(inputfile)
                if not line.strip():
                    line = next(inputfile)

                # Now we can iterate or atomic orbitals.
                for nao in range(self.nbasis):
                    data = list(map(float, line.split()[4:]))
                    for i,d in enumerate(data):
                        overlaps[indices[i]].append(d)
                    line = next(inputfile)

                line = next(inputfile)

            self.aooverlaps = overlaps


if __name__ == "__main__":
    import doctest, nwchemparser
    doctest.testmod(nwchemparser, verbose=False)
