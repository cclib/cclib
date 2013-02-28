# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

__revision__ = "$Revision$"

import re

import numpy

from . import logfileparser
from . import utils


class Gaussian(logfileparser.Logfile):
    """A Gaussian 98/03 log file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(Gaussian, self).__init__(logname="Gaussian", *args, **kwargs)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Gaussian log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Gaussian("%s")' % (self.filename)
    
    def normalisesym(self, label):
        """Use standard symmetry labels instead of Gaussian labels.

        To normalise:
        (1) If label is one of [SG, PI, PHI, DLTA], replace by [sigma, pi, phi, delta]
        (2) replace any G or U by their lowercase equivalent

        >>> sym = Gaussian("dummyfile").normalisesym
        >>> labels = ['A1', 'AG', 'A1G', "SG", "PI", "PHI", "DLTA", 'DLTU', 'SGG']
        >>> map(sym, labels)
        ['A1', 'Ag', 'A1g', 'sigma', 'pi', 'phi', 'delta', 'delta.u', 'sigma.g']
        """
        # note: DLT must come after DLTA
        greek = [('SG', 'sigma'), ('PI', 'pi'), ('PHI', 'phi'),
                 ('DLTA', 'delta'), ('DLT', 'delta')]
        for k, v in greek:
            if label.startswith(k):
                tmp = label[len(k):]
                label = v
                if tmp:
                    label = v + "." + tmp
        
        ans = label.replace("U", "u").replace("G", "g") 
        return ans

    def before_parsing(self):

        # Used to index self.scftargets[].
        SCFRMS, SCFMAX, SCFENERGY = list(range(3))

        # Flag that indicates whether it has reached the end of a geoopt.
        self.optfinished = False
        
        # Flag for identifying Coupled Cluster runs.
        self.coupledcluster = False

        # Fragment number for counterpoise or fragment guess calculations 
        # (normally zero).
        self.counterpoise = 0

        # Flag for identifying ONIOM calculations.
        self.oniom = False

    def after_parsing(self):

        # Correct the percent values in the etsecs in the case of
        # a restricted calculation. The following has the
        # effect of including each transition twice.
        if hasattr(self, "etsecs") and len(self.homos) == 1:
            new_etsecs = [[(x[0], x[1], x[2] * numpy.sqrt(2)) for x in etsec]
                          for etsec in self.etsecs]
            self.etsecs = new_etsecs
        if hasattr(self, "scanenergies"):
            self.scancoords = []
            self.scancoords = self.atomcoords
        if (hasattr(self, 'enthaply') and hasattr(self, 'temperature') 
                and hasattr(self, 'freeenergy')):
            self.entropy = (self.enthaply - self.freeenergy)/self.temperature
            
    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        #Extract PES scan data
        #Summary of the potential surface scan:
        #  N       A          SCF
        #----  ---------  -----------
        #   1   109.0000    -76.43373
        #   2   119.0000    -76.43011
        #   3   129.0000    -76.42311
        #   4   139.0000    -76.41398
        #   5   149.0000    -76.40420
        #   6   159.0000    -76.39541
        #   7   169.0000    -76.38916
        #   8   179.0000    -76.38664
        #   9   189.0000    -76.38833
        #  10   199.0000    -76.39391
        #  11   209.0000    -76.40231
        #----  ---------  -----------
        if "Summary of the potential surface scan:" in line:
            scanenergies = []
            scanparm = []
            colmnames = next(inputfile)
            hyphens = next(inputfile)
            line = next(inputfile)
            while line != hyphens:
                broken = line.split()
                scanenergies.append(float(broken[-1]))
                scanparm.append(map(float, broken[1:-1]))
                line = next(inputfile)
            if not hasattr(self, "scanenergies"):
                self.scanenergies = []
                self.scanenergies = scanenergies
            if not hasattr(self, "scanparm"):
                self.scanparm = []
                self.scanparm = scanparm
            if not hasattr(self, "scannames"):
                self.scannames = colmnames.split()[1:-1]

        #Extract Thermochemistry
        #Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
        #Zero-point correction=                           0.342233 (Hartree/
        #Thermal correction to Energy=                    0.
        #Thermal correction to Enthalpy=                  0.
        #Thermal correction to Gibbs Free Energy=         0.302940
        #Sum of electronic and zero-point Energies=           -563.649744
        #Sum of electronic and thermal Energies=              -563.636699
        #Sum of electronic and thermal Enthalpies=            -563.635755
        #Sum of electronic and thermal Free Energies=         -563.689037
        if "Sum of electronic and thermal Enthalpies" in line:
            if not hasattr(self, 'enthaply'):
                self.enthaply = float(line.split()[6])
        if "Sum of electronic and thermal Free Energies=" in line:
            if not hasattr(self, 'freeenergy'):
                self.freeenergy = float(line.split()[7])
        if line[1:12] == "Temperature":
            if not hasattr(self, 'temperature'):
                self.temperature = float(line.split()[1])
        
        # Number of atoms.
        if line[1:8] == "NAtoms=":

            self.updateprogress(inputfile, "Attributes", self.fupdate)
                    
            natom = int(line.split()[1])
            if not hasattr(self, "natom"):
                self.natom = natom

        # Catch message about completed optimization.
        if line[1:23] == "Optimization completed":
            self.optfinished = True
            self.optdone = True
        
        # Extract the atomic numbers and coordinates from the input orientation,
        #   in the event the standard orientation isn't available.
        if not self.optfinished and line.find("Input orientation") > -1 or line.find("Z-Matrix orientation") > -1:

            # If this is a counterpoise calculation, this output means that
            #   the supermolecule is now being considered, so we can set:
            self.counterpoise = 0

            self.updateprogress(inputfile, "Attributes", self.cupdate)
            
            if not hasattr(self, "inputcoords"):
                self.inputcoords = []
            self.inputatoms = []
            
            hyphens = next(inputfile)
            colmNames = next(inputfile)
            colmNames = next(inputfile)
            hyphens = next(inputfile)
            
            atomcoords = []
            line = next(inputfile)
            while line != hyphens:
                broken = line.split()
                self.inputatoms.append(int(broken[1]))
                atomcoords.append(list(map(float, broken[3:6])))
                line = next(inputfile)

            self.inputcoords.append(atomcoords)

            if not hasattr(self, "atomnos"):
                self.atomnos = numpy.array(self.inputatoms, 'i')
                self.natom = len(self.atomnos)

        # Extract the atomic masses.
        # Typical section:
        #                    Isotopes and Nuclear Properties:
        #(Nuclear quadrupole moments (NQMom) in fm**2, nuclear magnetic moments (NMagM)
        # in nuclear magnetons)
        #
        #  Atom         1           2           3           4           5           6           7           8           9          10
        # IAtWgt=          12          12          12          12          12           1           1           1          12          12
        # AtmWgt=  12.0000000  12.0000000  12.0000000  12.0000000  12.0000000   1.0078250   1.0078250   1.0078250  12.0000000  12.0000000
        # NucSpn=           0           0           0           0           0           1           1           1           0           0
        # AtZEff=  -3.6000000  -3.6000000  -3.6000000  -3.6000000  -3.6000000  -1.0000000  -1.0000000  -1.0000000  -3.6000000  -3.6000000
        # NQMom=    0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
        # NMagM=    0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   2.7928460   2.7928460   2.7928460   0.0000000   0.0000000
        # ... with blank lines dividing blocks of ten, and Leave Link 101 at the end.
        # This is generally parsed before coordinates, so atomnos is not defined.
        # Note that in Gaussian03 the comments are not there yet and the labels are different.
        if line.strip() == "Isotopes and Nuclear Properties:":

            if not hasattr(self, "atommasses"):
                self.atommasses = []

            line = next(inputfile)
            while line[1:16] != "Leave Link  101":
                if line[1:8] == "AtmWgt=":
                    self.atommasses.extend(list(map(float,line.split()[1:])))
                line = next(inputfile)

        # Extract the atomic numbers and coordinates of the atoms.
        if not self.optfinished and line.strip() == "Standard orientation:":

            self.updateprogress(inputfile, "Attributes", self.cupdate)

            # If this is a counterpoise calculation, this output means that
            #   the supermolecule is now being considered, so we can set:
            self.counterpoise = 0

            if not hasattr(self, "atomcoords"):
                self.atomcoords = []
            
            hyphens = next(inputfile)
            colmNames = next(inputfile)
            colmNames = next(inputfile)
            hyphens = next(inputfile)
            
            atomnos = []
            atomcoords = []
            line = next(inputfile)
            while line != hyphens:
                broken = line.split()
                atomnos.append(int(broken[1]))
                atomcoords.append(list(map(float, broken[-3:])))
                line = next(inputfile)
            self.atomcoords.append(atomcoords)

            if not hasattr(self, "natom"):
                self.atomnos = numpy.array(atomnos, 'i')
                self.natom = len(self.atomnos)

            # make sure atomnos is added for the case where natom has already been set
            elif not hasattr(self, "atomnos"):
                self.atomnos = numpy.array(atomnos, 'i')

        # Find the targets for SCF convergence (QM calcs).
        if line[1:44] == 'Requested convergence on RMS density matrix':

            if not hasattr(self, "scftargets"):
                self.scftargets = []
            # The following can happen with ONIOM which are mixed SCF
            # and semi-empirical
            if type(self.scftargets) == type(numpy.array([])):
                self.scftargets = []

            scftargets = []
            # The RMS density matrix.
            scftargets.append(self.float(line.split('=')[1].split()[0]))
            line = next(inputfile)
            # The MAX density matrix.
            scftargets.append(self.float(line.strip().split('=')[1][:-1]))
            line = next(inputfile)
            # For G03, there's also the energy (not for G98).
            if line[1:10] == "Requested":
                scftargets.append(self.float(line.strip().split('=')[1][:-1]))

            self.scftargets.append(scftargets)

        # Extract SCF convergence information (QM calcs).
        if line[1:10] == 'Cycle   1':
                    
            if not hasattr(self, "scfvalues"):
                self.scfvalues = []

            scfvalues = []
            line = next(inputfile)
            while line.find("SCF Done") == -1:
            
                self.updateprogress(inputfile, "QM convergence", self.fupdate)
                      
                if line.find(' E=') == 0:
                    self.logger.debug(line)

                #  RMSDP=3.74D-06 MaxDP=7.27D-05 DE=-1.73D-07 OVMax= 3.67D-05
                # or
                #  RMSDP=1.13D-05 MaxDP=1.08D-04              OVMax= 1.66D-04
                if line.find(" RMSDP") == 0:

                    parts = line.split()
                    newlist = [self.float(x.split('=')[1]) for x in parts[0:2]]
                    energy = 1.0
                    if len(parts) > 4:
                        energy = parts[2].split('=')[1]
                        if energy == "":
                            energy = self.float(parts[3])
                        else:
                            energy = self.float(energy)
                    if len(self.scftargets[0]) == 3: # Only add the energy if it's a target criteria
                        newlist.append(energy)
                    scfvalues.append(newlist)

                try:
                    line = next(inputfile)
                # May be interupted by EOF.
                except StopIteration:
                    break

            self.scfvalues.append(scfvalues)

        # Extract SCF convergence information (AM1 calcs).
        if line[1:4] == 'It=':
                    
            self.scftargets = numpy.array([1E-7], "d") # This is the target value for the rms
            self.scfvalues = [[]]

            line = next(inputfile)
            while line.find(" Energy") == -1:
            
                if self.progress:
                    step = inputfile.tell()
                    if step != oldstep:
                        self.progress.update(step, "AM1 Convergence")
                        oldstep = step
                        
                if line[1:4] == "It=":
                    parts = line.strip().split()
                    self.scfvalues[0].append(self.float(parts[-1][:-1]))
                line = next(inputfile)

        # Note: this needs to follow the section where 'SCF Done' is used
        #   to terminate a loop when extracting SCF convergence information.
        if line[1:9] == 'SCF Done':

            if not hasattr(self, "scfenergies"):
                self.scfenergies = []

            self.scfenergies.append(utils.convertor(self.float(line.split()[4]), "hartree", "eV"))
        # gmagoon 5/27/09: added scfenergies reading for PM3 case
        # Example line: " Energy=   -0.077520562724 NIter=  14."
        # See regression Gaussian03/QVGXLLKOCUKJST-UHFFFAOYAJmult3Fixed.out
        if line[1:8] == 'Energy=':
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            self.scfenergies.append(utils.convertor(self.float(line.split()[1]), "hartree", "eV"))
        
        # Total energies after Moller-Plesset corrections.
        # Second order correction is always first, so its first occurance
        #   triggers creation of mpenergies (list of lists of energies).
        # Further MP2 corrections are appended as found.
        #
        # Example MP2 output line:
        #  E2 =    -0.9505918144D+00 EUMP2 =    -0.28670924198852D+03
        # Warning! this output line is subtly different for MP3/4/5 runs
        if "EUMP2" in line[27:34]:

            if not hasattr(self, "mpenergies"):
                self.mpenergies = []
            self.mpenergies.append([])
            mp2energy = self.float(line.split("=")[2])
            self.mpenergies[-1].append(utils.convertor(mp2energy, "hartree", "eV"))

        # Example MP3 output line:
        #  E3=       -0.10518801D-01     EUMP3=      -0.75012800924D+02
        if line[34:39] == "EUMP3":

            mp3energy = self.float(line.split("=")[2])
            self.mpenergies[-1].append(utils.convertor(mp3energy, "hartree", "eV"))

        # Example MP4 output lines:
        #  E4(DQ)=   -0.31002157D-02        UMP4(DQ)=   -0.75015901139D+02
        #  E4(SDQ)=  -0.32127241D-02        UMP4(SDQ)=  -0.75016013648D+02
        #  E4(SDTQ)= -0.32671209D-02        UMP4(SDTQ)= -0.75016068045D+02
        # Energy for most substitutions is used only (SDTQ by default)
        if line[34:42] == "UMP4(DQ)":

            mp4energy = self.float(line.split("=")[2])
            line = next(inputfile)
            if line[34:43] == "UMP4(SDQ)":
              mp4energy = self.float(line.split("=")[2])
              line = next(inputfile)
              if line[34:44] == "UMP4(SDTQ)":
                mp4energy = self.float(line.split("=")[2])
            self.mpenergies[-1].append(utils.convertor(mp4energy, "hartree", "eV"))

        # Example MP5 output line:
        #  DEMP5 =  -0.11048812312D-02 MP5 =  -0.75017172926D+02
        if line[29:32] == "MP5":
            mp5energy = self.float(line.split("=")[2])
            self.mpenergies[-1].append(utils.convertor(mp5energy, "hartree", "eV"))

        # Total energies after Coupled Cluster corrections.
        # Second order MBPT energies (MP2) are also calculated for these runs,
        # but the output is the same as when parsing for mpenergies.
        # Read the consecutive correlated energies
        # but append only the last one to ccenergies.
        # Only the highest level energy is appended - ex. CCSD(T), not CCSD.
        if line[1:10] == "DE(Corr)=" and line[27:35] == "E(CORR)=":
            self.ccenergy = self.float(line.split()[3])
        if line[1:10] == "T5(CCSD)=":
            line = next(inputfile)
            if line[1:9] == "CCSD(T)=":
                self.ccenergy = self.float(line.split()[1])
        if line[12:53] == "Population analysis using the SCF density":
            if hasattr(self, "ccenergy"):
                if not hasattr(self, "ccenergies"):
                    self.ccenergies = []
                self.ccenergies.append(utils.convertor(self.ccenergy, "hartree", "eV"))
                del self.ccenergy

        # Geometry convergence information.
        if line[49:59] == 'Converged?':

            if not hasattr(self, "geotargets"):
                self.geovalues = []
                self.geotargets = numpy.array([0.0, 0.0, 0.0, 0.0], "d")

            newlist = [0]*4
            for i in range(4):
                line = next(inputfile)
                self.logger.debug(line)
                parts = line.split()
                try:
                    value = self.float(parts[2])
                except ValueError:
                    self.logger.error("Problem parsing the value for geometry optimisation: %s is not a number." % parts[2])
                else:
                    newlist[i] = value
                self.geotargets[i] = self.float(parts[3])

            self.geovalues.append(newlist)

        # Gradients.
        # Read in the cartesian energy gradients (forces) from a block like this:
        # -------------------------------------------------------------------
        # Center     Atomic                   Forces (Hartrees/Bohr)
        # Number     Number              X              Y              Z
        # -------------------------------------------------------------------
        # 1          1          -0.012534744   -0.021754635   -0.008346094
        # 2          6           0.018984731    0.032948887   -0.038003451
        # 3          1          -0.002133484   -0.006226040    0.023174772
        # 4          1          -0.004316502   -0.004968213    0.023174772
        #           -2          -0.001830728   -0.000743108   -0.000196625
        # ------------------------------------------------------------------
        #
        # The "-2" line is for a dummy atom
        #
        # Then optimization is done in internal coordinates, Gaussian also
        # print the forces in internal coordinates, which can be produced from 
        # the above. This block looks like this:
        # Variable       Old X    -DE/DX   Delta X   Delta X   Delta X     New X
        #                                 (Linear)    (Quad)   (Total)
        #   ch        2.05980   0.01260   0.00000   0.01134   0.01134   2.07114
        #   hch        1.75406   0.09547   0.00000   0.24861   0.24861   2.00267
        #   hchh       2.09614   0.01261   0.00000   0.16875   0.16875   2.26489
        #         Item               Value     Threshold  Converged?
        if line[37:43] == "Forces":

            if not hasattr(self, "grads"):
                self.grads = []

            header = next(inputfile)
            dashes = next(inputfile)
            line = next(inputfile)
            forces = []
            while line != dashes:
                broken = line.split()
                Fx, Fy, Fz = broken[-3:]
                forces.append([float(Fx), float(Fy), float(Fz)])
                line = next(inputfile)
            self.grads.append(forces)                

        # Charge and multiplicity.
        # If counterpoise correction is used, multiple lines match.
        # The first one contains charge/multiplicity of the whole molecule.:
        #   Charge =  0 Multiplicity = 1 in supermolecule
        #   Charge =  0 Multiplicity = 1 in fragment  1.
        #   Charge =  0 Multiplicity = 1 in fragment  2.
        if line[1:7] == 'Charge' and line.find("Multiplicity")>=0:

            regex = ".*=(.*)Mul.*=\s*-?(\d+).*"
            match = re.match(regex, line)
            assert match, "Something unusual about the line: '%s'" % line
            
            if not hasattr(self, "charge"):
                self.charge = int(match.groups()[0])
            if not hasattr(self, "mult"):
                self.mult = int(match.groups()[1])

        # Orbital symmetries.
        if line[1:20] == 'Orbital symmetries:' and not hasattr(self, "mosyms"):

            # For counterpoise fragments, skip these lines.
            if self.counterpoise != 0: return

            self.updateprogress(inputfile, "MO Symmetries", self.fupdate)
                    
            self.mosyms = [[]]
            line = next(inputfile)
            unres = False
            if line.find("Alpha Orbitals") == 1:
                unres = True
                line = next(inputfile)
            i = 0
            while len(line) > 18 and line[17] == '(':
                if line.find('Virtual') >= 0:
                    self.homos = numpy.array([i-1], "i") # 'HOMO' indexes the HOMO in the arrays
                parts = line[17:].split()
                for x in parts:
                    self.mosyms[0].append(self.normalisesym(x.strip('()')))
                    i += 1 
                line = next(inputfile)
            if unres:
                line = next(inputfile)
                # Repeat with beta orbital information
                i = 0
                self.mosyms.append([])
                while len(line) > 18 and line[17] == '(':
                    if line.find('Virtual')>=0:
                        # Here we consider beta
                        # If there was also an alpha virtual orbital,
                        #  we will store two indices in the array
                        # Otherwise there is no alpha virtual orbital,
                        #  only beta virtual orbitals, and we initialize
                        #  the array with one element. See the regression
                        #  QVGXLLKOCUKJST-UHFFFAOYAJmult3Fixed.out
                        #  donated by Gregory Magoon (gmagoon).
                        if (hasattr(self, "homos")):
                            # Extend the array to two elements
                            # 'HOMO' indexes the HOMO in the arrays
                            self.homos.resize([2])
                            self.homos[1] = i-1
                        else:
                            # 'HOMO' indexes the HOMO in the arrays
                            self.homos = numpy.array([i-1], "i")
                    parts = line[17:].split()
                    for x in parts:
                        self.mosyms[1].append(self.normalisesym(x.strip('()')))
                        i += 1
                    line = next(inputfile)

        # Alpha/Beta electron eigenvalues.
        if line[1:6] == "Alpha" and line.find("eigenvalues") >= 0:

            # For counterpoise fragments, skip these lines.
            if self.counterpoise != 0: return

            # For ONIOM calcs, ignore this section in order to bypass assertion failure.
            if self.oniom: return

            self.updateprogress(inputfile, "Eigenvalues", self.fupdate)
            self.moenergies = [[]]
            HOMO = -2

            while line.find('Alpha') == 1:
                if line.split()[1] == "virt." and HOMO == -2:

                    # If there aren't any symmetries, this is a good way to find the HOMO.
                    # Also, check for consistency if homos was already parsed.
                    HOMO = len(self.moenergies[0])-1
                    if hasattr(self, "homos"):
                        assert HOMO == self.homos[0]
                    else:
                        self.homos = numpy.array([HOMO], "i")

                # Convert to floats and append to moenergies, but sometimes Gaussian
                #  doesn't print correctly so test for ValueError (bug 1756789).
                part = line[28:]
                i = 0
                while i*10+4 < len(part):
                    s = part[i*10:(i+1)*10]
                    try:
                        x = self.float(s)
                    except ValueError:
                        x = numpy.nan
                    self.moenergies[0].append(utils.convertor(x, "hartree", "eV"))
                    i += 1
                line = next(inputfile)

            # If, at this point, self.homos is unset, then there were not
            # any alpha virtual orbitals
            if not hasattr(self, "homos"):
                HOMO = len(self.moenergies[0])-1
                self.homos = numpy.array([HOMO], "i")
            

            if line.find('Beta') == 2:
                self.moenergies.append([])

            HOMO = -2
            while line.find('Beta') == 2:
                if line.split()[1] == "virt." and HOMO == -2:

                    # If there aren't any symmetries, this is a good way to find the HOMO.
                    # Also, check for consistency if homos was already parsed.
                    HOMO = len(self.moenergies[1])-1
                    if len(self.homos) == 2:
                        assert HOMO == self.homos[1]
                    else:
                        self.homos.resize([2])
                        self.homos[1] = HOMO

                part = line[28:]
                i = 0
                while i*10+4 < len(part):
                    x = part[i*10:(i+1)*10]
                    self.moenergies[1].append(utils.convertor(self.float(x), "hartree", "eV"))
                    i += 1
                line = next(inputfile)

            self.moenergies = [numpy.array(x, "d") for x in self.moenergies]
            
        # Gaussian Rev <= B.0.3 (?)
        # AO basis set in the form of general basis input:
        #  1 0
        # S   3 1.00       0.000000000000
        #      0.7161683735D+02  0.1543289673D+00
        #      0.1304509632D+02  0.5353281423D+00
        #      0.3530512160D+01  0.4446345422D+00
        # SP   3 1.00       0.000000000000
        #      0.2941249355D+01 -0.9996722919D-01  0.1559162750D+00
        #      0.6834830964D+00  0.3995128261D+00  0.6076837186D+00
        #      0.2222899159D+00  0.7001154689D+00  0.3919573931D+00
        if line[1:16] == "AO basis set in":
        
            # For counterpoise fragment calcualtions, skip these lines.
            if self.counterpoise != 0: return
        
            self.gbasis = []
            line = next(inputfile)
            while line.strip():
                gbasis = []
                line = next(inputfile)
                while line.find("*")<0:
                    temp = line.split()
                    symtype = temp[0]
                    numgau = int(temp[1])
                    gau = []
                    for i in range(numgau):
                        temp = list(map(self.float, next(inputfile).split()))
                        gau.append(temp)
                        
                    for i, x in enumerate(symtype):
                        newgau = [(z[0], z[i+1]) for z in gau]
                        gbasis.append((x, newgau))
                    line = next(inputfile) # i.e. "****" or "SP ...."
                self.gbasis.append(gbasis)
                line = next(inputfile) # i.e. "20 0" or blank line

        # Start of the IR/Raman frequency section.
        # Caution is advised here, as additional frequency blocks
        #   can be printed by Gaussian (with slightly different formats),
        #   often doubling the information printed.
        # See, for a non-standard exmaple, regression Gaussian98/test_H2.log
        if line[1:14] == "Harmonic freq":

            self.updateprogress(inputfile, "Frequency Information", self.fupdate)
            removeold = False

            # The whole block should not have any blank lines.
            while line.strip() != "":

                # The line with indices
                if line[1:15].strip() == "" and line[15:22].strip().isdigit():
                    freqbase = int(line[15:22])
                    if freqbase == 1 and hasattr(self, 'vibfreqs'):
                        # This is a reparse of this information
                        removeold = True

                # Lines with symmetries and symm. indices begin with whitespace.
                if line[1:15].strip() == "" and not line[15:22].strip().isdigit():

                    if not hasattr(self, 'vibsyms'):
                        self.vibsyms = []
                    syms = line.split()
                    self.vibsyms.extend(syms)
            
                if line[1:15] == "Frequencies --":
                
                    if not hasattr(self, 'vibfreqs'):
                        self.vibfreqs = []
                        
                    if removeold: # This is a reparse, so throw away the old info
                        if hasattr(self, "vibsyms"):
                            # We have already parsed the vibsyms so don't throw away!
                            self.vibsyms = self.vibsyms[-len(line[15:].split()):]
                        if hasattr(self, "vibirs"):
                            self.vibirs = []
                        if hasattr(self, 'vibfreqs'):
                            self.vibfreqs = []
                        if hasattr(self, 'vibramans'):
                            self.vibramans = []
                        if hasattr(self, 'vibdisps'):
                            self.vibdisps = []
                        removeold = False
                        
                    freqs = [self.float(f) for f in line[15:].split()]
                    self.vibfreqs.extend(freqs)
            
                if line[1:15] == "IR Inten    --":
                
                    if not hasattr(self, 'vibirs'):
                        self.vibirs = []
                    irs = [self.float(f) for f in line[15:].split()]
                    self.vibirs.extend(irs)

                if line[1:15] == "Raman Activ --":
                
                    if not hasattr(self, 'vibramans'):
                        self.vibramans = []
                    ramans = [self.float(f) for f in line[15:].split()]
                    self.vibramans.extend(ramans)
                
                # Block with displacement should start with this.
                if line.strip().split()[0:3] == ["Atom", "AN", "X"]:
                    if not hasattr(self, 'vibdisps'):
                        self.vibdisps = []
                    disps = []
                    for n in range(self.natom):
                        line = next(inputfile)
                        numbers = [float(s) for s in line[10:].split()]
                        N = len(numbers) // 3
                        if not disps:
                            for n in range(N):
                                disps.append([])
                        for n in range(N):
                            disps[n].append(numbers[3*n:3*n+3])
                    self.vibdisps.extend(disps)
                
                line = next(inputfile)

# Below is the old code for the IR/Raman frequency block, can probably be removed.
#            while len(line[:15].split()) == 0:
#                self.logger.debug(line)
#                self.vibsyms.extend(line.split()) # Adding new symmetry
#                line = inputfile.next()
#                # Read in frequencies.
#                freqs = [self.float(f) for f in line.split()[2:]]
#                self.vibfreqs.extend(freqs)
#                line = inputfile.next()
#                line = inputfile.next()
#                line = inputfile.next()
#                irs = [self.float(f) for f in line.split()[3:]]
#                self.vibirs.extend(irs)
#                line = inputfile.next() # Either the header or a Raman line
#                if line.find("Raman") >= 0:
#                    if not hasattr(self, "vibramans"):
#                        self.vibramans = []
#                    ramans = [self.float(f) for f in line.split()[3:]]
#                    self.vibramans.extend(ramans)
#                    line = inputfile.next() # Depolar (P)
#                    line = inputfile.next() # Depolar (U)
#                    line = inputfile.next() # Header
#                line = inputfile.next() # First line of cartesian displacement vectors
#                p = [[], [], []]
#                while len(line[:15].split()) > 0:
#                    # Store the cartesian displacement vectors
#                    broken = map(float, line.strip().split()[2:])
#                    for i in range(0, len(broken), 3):
#                        p[i/3].append(broken[i:i+3])
#                    line = inputfile.next()
#                self.vibdisps.extend(p[0:len(broken)/3])
#                line = inputfile.next() # Should be the line with symmetries
#            self.vibfreqs = numpy.array(self.vibfreqs, "d")
#            self.vibirs = numpy.array(self.vibirs, "d")
#            self.vibdisps = numpy.array(self.vibdisps, "d")
#            if hasattr(self, "vibramans"):
#                self.vibramans = numpy.array(self.vibramans, "d")
                
        # Electronic transitions.
        if line[1:14] == "Excited State":
        
            if not hasattr(self, "etenergies"):
                self.etenergies = []
                self.etoscs = []
                self.etsyms = []
                self.etsecs = []
            # Need to deal with lines like:
            # (restricted calc)
            # Excited State   1:   Singlet-BU     5.3351 eV  232.39 nm  f=0.1695
            # (unrestricted calc) (first excited state is 2!)
            # Excited State   2:   ?Spin  -A      0.1222 eV 10148.75 nm  f=0.0000
            # (Gaussian 09 ZINDO)
            # Excited State   1:      Singlet-?Sym    2.5938 eV  478.01 nm  f=0.0000  <S**2>=0.000
            p = re.compile(":(?P<sym>.*?)(?P<energy>-?\d*\.\d*) eV")
            groups = p.search(line).groups()
            self.etenergies.append(utils.convertor(self.float(groups[1]), "eV", "cm-1"))
            self.etoscs.append(self.float(line.split("f=")[-1].split()[0]))
            self.etsyms.append(groups[0].strip())
            
            line = next(inputfile)

            p = re.compile("(\d+)")
            CIScontrib = []
            while line.find(" ->") >= 0: # This is a contribution to the transition
                parts = line.split("->")
                self.logger.debug(parts)
                # Has to deal with lines like:
                #       32 -> 38         0.04990
                #      35A -> 45A        0.01921
                frommoindex = 0 # For restricted or alpha unrestricted
                fromMO = parts[0].strip()
                if fromMO[-1] == "B":
                    frommoindex = 1 # For beta unrestricted
                fromMO = int(p.match(fromMO).group())-1 # subtract 1 so that it is an index into moenergies
                
                t = parts[1].split()
                tomoindex = 0
                toMO = t[0]
                if toMO[-1] == "B":
                    tomoindex = 1
                toMO = int(p.match(toMO).group())-1 # subtract 1 so that it is an index into moenergies

                percent = self.float(t[1])
                # For restricted calculations, the percentage will be corrected
                # after parsing (see after_parsing() above).
                CIScontrib.append([(fromMO, frommoindex), (toMO, tomoindex), percent])
                line = next(inputfile)
            self.etsecs.append(CIScontrib)

# Circular dichroism data (different for G03 vs G09)

# G03

## <0|r|b> * <b|rxdel|0>  (Au), Rotatory Strengths (R) in
## cgs (10**-40 erg-esu-cm/Gauss)
##       state          X           Y           Z     R(length)
##         1         0.0006      0.0096     -0.0082     -0.4568
##         2         0.0251     -0.0025      0.0002     -5.3846
##         3         0.0168      0.4204     -0.3707    -15.6580
##         4         0.0721      0.9196     -0.9775     -3.3553

# G09

## 1/2[<0|r|b>*<b|rxdel|0> + (<0|rxdel|b>*<b|r|0>)*]
## Rotatory Strengths (R) in cgs (10**-40 erg-esu-cm/Gauss)
##       state          XX          YY          ZZ     R(length)     R(au)
##         1        -0.3893     -6.7546      5.7736     -0.4568     -0.0010
##         2       -17.7437      1.7335     -0.1435     -5.3845     -0.0114
##         3       -11.8655   -297.2604    262.1519    -15.6580     -0.0332

        if (line[1:52] == "<0|r|b> * <b|rxdel|0>  (Au), Rotatory Strengths (R)" or
            line[1:50] == "1/2[<0|r|b>*<b|rxdel|0> + (<0|rxdel|b>*<b|r|0>)*]"):

            self.etrotats = []
            next(inputfile) # Units
            headers = next(inputfile) # Headers
            Ncolms = len(headers.split())
            line = next(inputfile)
            parts = line.strip().split()
            while len(parts) == Ncolms:
                try:
                    R = self.float(parts[4])
                except ValueError:
                    # nan or -nan if there is no first excited state
                    # (for unrestricted calculations)
                    pass
                else:
                    self.etrotats.append(R)
                line = next(inputfile)
                temp = line.strip().split()
                parts = line.strip().split()                
            self.etrotats = numpy.array(self.etrotats, "d")

        # Number of basis sets functions.
        # Has to deal with lines like:
        #  NBasis =   434 NAE=    97 NBE=    97 NFC=    34 NFV=     0
        # and...
        #  NBasis = 148  MinDer = 0  MaxDer = 0
        # Although the former is in every file, it doesn't occur before
        #   the overlap matrix is printed.
        if line[1:7] == "NBasis" or line[4:10] == "NBasis":

            # For counterpoise fragment, skip these lines.
            if self.counterpoise != 0: return

            # For ONIOM calcs, ignore this section in order to bypass assertion failure.
            if self.oniom: return

            # If nbasis was already parsed, check if it changed.
            nbasis = int(line.split('=')[1].split()[0])
            if hasattr(self, "nbasis"):
                assert nbasis == self.nbasis
            else:
                self.nbasis = nbasis
                
        # Number of linearly-independent basis sets.
        if line[1:7] == "NBsUse":

            # For counterpoise fragment, skip these lines.
            if self.counterpoise != 0: return

            # For ONIOM calcs, ignore this section in order to bypass assertion failure.
            if self.oniom: return

            nmo = int(line.split('=')[1].split()[0])
            if hasattr(self, "nmo"):
                assert nmo == self.nmo
            else:
                self.nmo = nmo

        # For AM1 calculations, set nbasis by a second method,
        #   as nmo may not always be explicitly stated.
        if line[7:22] == "basis functions, ":
        
            nbasis = int(line.split()[0])
            if hasattr(self, "nbasis"):
                assert nbasis == self.nbasis
            else:
                self.nbasis = nbasis

        # Molecular orbital overlap matrix.
        # Has to deal with lines such as:
        #   *** Overlap ***
        #   ****** Overlap ******
        # Note that Gaussian sometimes drops basis functions,
        #  causing the overlap matrix as parsed below to not be
        #  symmetric (which is a problem for population analyses, etc.)
        if line[1:4] == "***" and (line[5:12] == "Overlap"
                                 or line[8:15] == "Overlap"):

            # Ensure that this is the main calc and not a fragment
            if self.counterpoise != 0: return

            self.aooverlaps = numpy.zeros( (self.nbasis, self.nbasis), "d")
            # Overlap integrals for basis fn#1 are in aooverlaps[0]
            base = 0
            colmNames = next(inputfile)
            while base < self.nbasis:
                 
                self.updateprogress(inputfile, "Overlap", self.fupdate)
                        
                for i in range(self.nbasis-base): # Fewer lines this time
                    line = next(inputfile)
                    parts = line.split()
                    for j in range(len(parts)-1): # Some lines are longer than others
                        k = float(parts[j+1].replace("D", "E"))
                        self.aooverlaps[base+j, i+base] = k
                        self.aooverlaps[i+base, base+j] = k
                base += 5
                colmNames = next(inputfile)
            self.aooverlaps = numpy.array(self.aooverlaps, "d")                    

        # Molecular orbital coefficients (mocoeffs).
        # Essentially only produced for SCF calculations.
        # This is also the place where aonames and atombasis are parsed.
        if line[5:35] == "Molecular Orbital Coefficients" or line[5:41] == "Alpha Molecular Orbital Coefficients" or line[5:40] == "Beta Molecular Orbital Coefficients":

            # If counterpoise fragment, return without parsing orbital info
            if self.counterpoise != 0: return
            # Skip this for ONIOM calcs
            if self.oniom: return

            if line[5:40] == "Beta Molecular Orbital Coefficients":
                beta = True
                if self.popregular:
                    return
                    # This was continue before refactoring the parsers.
                    #continue # Not going to extract mocoeffs
                # Need to add an extra array to self.mocoeffs
                self.mocoeffs.append(numpy.zeros((self.nmo, self.nbasis), "d"))
            else:
                beta = False
                self.aonames = []
                self.atombasis = []
                mocoeffs = [numpy.zeros((self.nmo, self.nbasis), "d")]

            base = 0
            self.popregular = False
            for base in range(0, self.nmo, 5):
                
                self.updateprogress(inputfile, "Coefficients", self.fupdate)
                         
                colmNames = next(inputfile)   

                if not colmNames.split():
                    self.logger.warning("Molecular coefficients header found but no coefficients.")
                    break;

                if base == 0 and int(colmNames.split()[0]) != 1:
                    # Implies that this is a POP=REGULAR calculation
                    # and so, only aonames (not mocoeffs) will be extracted
                    self.popregular = True
                symmetries = next(inputfile)
                eigenvalues = next(inputfile)
                for i in range(self.nbasis):
                                   
                    line = next(inputfile)
                    if i == 0:
                        # Find location of the start of the basis function name
                        start_of_basis_fn_name = line.find(line.split()[3]) - 1
                    if base == 0 and not beta: # Just do this the first time 'round
                        parts = line[:start_of_basis_fn_name].split()
                        if len(parts) > 1: # New atom
                            if i > 0:
                                self.atombasis.append(atombasis)
                            atombasis = []
                            atomname = "%s%s" % (parts[2], parts[1])
                        orbital = line[start_of_basis_fn_name:20].strip()
                        self.aonames.append("%s_%s" % (atomname, orbital))
                        atombasis.append(i)

                    part = line[21:].replace("D", "E").rstrip()
                    temp = [] 
                    for j in range(0, len(part), 10):
                        temp.append(float(part[j:j+10]))
                    if beta:
                        self.mocoeffs[1][base:base + len(part) / 10, i] = temp
                    else:
                        mocoeffs[0][base:base + len(part) / 10, i] = temp

                if base == 0 and not beta: # Do the last update of atombasis
                    self.atombasis.append(atombasis)
                if self.popregular:
                    # We now have aonames, so no need to continue
                    break
            if not self.popregular and not beta:
                self.mocoeffs = mocoeffs

        # Natural Orbital Coefficients (nocoeffs) - alternative for mocoeffs.
        # Most extensively formed after CI calculations, but not only.
        # Like for mocoeffs, this is also where aonames and atombasis are parsed.
        if line[5:33] == "Natural Orbital Coefficients":

            self.aonames = []
            self.atombasis = []
            nocoeffs = numpy.zeros((self.nmo, self.nbasis), "d")

            base = 0
            self.popregular = False
            for base in range(0, self.nmo, 5):
                
                self.updateprogress(inputfile, "Coefficients", self.fupdate)
                         
                colmNames = next(inputfile)   
                if base == 0 and int(colmNames.split()[0]) != 1:
                    # Implies that this is a POP=REGULAR calculation
                    # and so, only aonames (not mocoeffs) will be extracted
                    self.popregular = True

                # No symmetry line for natural orbitals.
                # symmetries = inputfile.next()
                eigenvalues = next(inputfile)

                for i in range(self.nbasis):
                                   
                    line = next(inputfile)

                    # Just do this the first time 'round.
                    if base == 0:

                        # Changed below from :12 to :11 to deal with Elmar Neumann's example.
                        parts = line[:11].split()
                        # New atom.
                        if len(parts) > 1:
                            if i > 0:
                                self.atombasis.append(atombasis)
                            atombasis = []
                            atomname = "%s%s" % (parts[2], parts[1])
                        orbital = line[11:20].strip()
                        self.aonames.append("%s_%s" % (atomname, orbital))
                        atombasis.append(i)

                    part = line[21:].replace("D", "E").rstrip()
                    temp = [] 

                    for j in range(0, len(part), 10):
                        temp.append(float(part[j:j+10]))

                    nocoeffs[base:base + len(part) / 10, i] = temp

                # Do the last update of atombasis.
                if base == 0:
                    self.atombasis.append(atombasis)

                # We now have aonames, so no need to continue.
                if self.popregular:
                    break

            if not self.popregular:
                self.nocoeffs = nocoeffs

        # For FREQ=Anharm, extract anharmonicity constants
        if line[1:40] == "X matrix of Anharmonic Constants (cm-1)":
            Nvibs = len(self.vibfreqs)
            self.vibanharms = numpy.zeros( (Nvibs, Nvibs), "d")

            base = 0
            colmNames = next(inputfile)
            while base < Nvibs:
                 
                for i in range(Nvibs-base): # Fewer lines this time
                    line = next(inputfile)
                    parts = line.split()
                    for j in range(len(parts)-1): # Some lines are longer than others
                        k = float(parts[j+1].replace("D", "E"))
                        self.vibanharms[base+j, i+base] = k
                        self.vibanharms[i+base, base+j] = k
                base += 5
                colmNames = next(inputfile)

        # Pseudopotential charges.
        if line.find("Pseudopotential Parameters") > -1:

            dashes = next(inputfile)
            label1 = next(inputfile)
            label2 = next(inputfile)
            dashes = next(inputfile)

            line = next(inputfile)
            if line.find("Centers:") < 0:
                return
                # This was continue before parser refactoring.
                # continue

# Needs to handle code like the following:
#
#  Center     Atomic      Valence      Angular      Power                                                       Coordinates
#  Number     Number     Electrons     Momentum     of R      Exponent        Coefficient                X           Y           Z
# ===================================================================================================================================
# Centers:   1
# Centers:  16
# Centers:  21 24
# Centers:  99100101102
#    1         44           16                                                                      -4.012684 -0.696698  0.006750
#                                      F and up 
#                                                     0      554.3796303       -0.05152700                

            centers = []
            while line.find("Centers:") >= 0:
                temp = line[10:]
                for i in range(0, len(temp)-3, 3):
                    centers.append(int(temp[i:i+3]))
                line = next(inputfile)
            centers.sort() # Not always in increasing order
            
            self.coreelectrons = numpy.zeros(self.natom, "i")

            for center in centers:
                front = line[:10].strip()
                while not (front and int(front) == center):
                    line = next(inputfile)
                    front = line[:10].strip()
                info = line.split()
                self.coreelectrons[center-1] = int(info[1]) - int(info[2])
                line = next(inputfile)

        # This will be printed for counterpoise calcualtions only.
        # To prevent crashing, we need to know which fragment is being considered.
        # Other information is also printed in lines that start like this.
        if line[1:14] == 'Counterpoise:':
        
            if line[42:50] == "fragment":
                self.counterpoise = int(line[51:54])

        # This will be printed only during ONIOM calcs; use it to set a flag
        # that will allow assertion failures to be bypassed in the code.
        if line[1:7] == "ONIOM:":
            self.oniom = True

        if (line[1:24] == "Mulliken atomic charges" or
            line[1:22] == "Lowdin Atomic Charges"):
            if not hasattr(self, "atomcharges"):
                self.atomcharges = {}
            ones = next(inputfile)
            charges = []
            nline = next(inputfile)
            while not "Sum of" in nline:
                charges.append(float(nline.split()[2]))
                nline = next(inputfile)
            if "Mulliken" in line:
                self.atomcharges["mulliken"] = charges
            else:
                self.atomcharges["lowdin"] = charges



if __name__ == "__main__":
    import doctest, gaussianparser, sys

    if len(sys.argv) == 1:
        doctest.testmod(gaussianparser, verbose=False)

    if len(sys.argv) >= 2:
        parser = gaussianparser.Gaussian(sys.argv[1])
        data = parser.parse()

    if len(sys.argv) > 2:
        for i in range(len(sys.argv[2:])):
            if hasattr(data, sys.argv[2 + i]):
                print(getattr(data, sys.argv[2 + i]))

