"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"


import re

import numpy

import utils
import logfileparser


class Gaussian(logfileparser.Logfile):
    """A Gaussian 98/03 log file."""

    def __init__(self, *args):

        # Call the __init__ method of the superclass
        super(Gaussian, self).__init__(logname="Gaussian", *args)
        
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
        for k,v in greek:
            if label.startswith(k):
                tmp = label[len(k):]
                label = v
                if tmp:
                    label = v + "." + tmp
        
        ans = label.replace("U", "u").replace("G", "g") 
        return ans

    def before_parsing(self):

        # Used to index self.scftargets[].
        SCFRMS, SCFMAX, SCFENERGY = range(3)
        # Flag that indicates whether it has reached the end of a geoopt.
        self.optfinished = False
        # Flag for identifying Coupled Cluster runs.
        self.coupledcluster = False
        # Fragment number for counterpoise calculations (normally zero).
        self.counterpoise = 0

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""
        
        # Number of atoms.
        if line[1:8] == "NAtoms=":

            self.updateprogress(inputfile, "Attributes", self.fupdate)
                    
            natom = int(line.split()[1])
            if hasattr(self, "natom"):
                assert self.natom == natom
            else:
                # I wonder whether this code will ever be executed.
                self.natom = natom

        # Catch message about completed optimization.
        if line[1:23] == "Optimization completed":
            self.optfinished = True
        
        # Extract the atomic numbers and coordinates from the input orientation,
        #   in the event the standard orientation isn't available.
        if line.find("Input orientation") > -1 or line.find("Z-Matrix orientation") > -1:

            # If this is a counterpoise calculation, this output means that
            #   the supermolecule is now being considered, so we can set:
            self.counterpoise = 0

            self.updateprogress(inputfile, "Attributes", self.cupdate)
                    
            self.inputcoords = []
            self.inputatoms = []
            
            hyphens = inputfile.next()
            colmNames = inputfile.next()
            colmNames = inputfile.next()
            hyphens = inputfile.next()
            
            atomcoords = []
            line = inputfile.next()
            while line != hyphens:
                broken = line.split()
                self.inputatoms.append(int(broken[1]))
                atomcoords.append(map(float, broken[3:6]))
                line = inputfile.next()

            self.inputcoords.append(atomcoords)

            if not hasattr(self, "natom"):
                self.atomnos = numpy.array(self.inputatoms, 'i')
                self.natom = len(self.atomnos)

        # Extract the atomic numbers and coordinates of the atoms.
        if not self.optfinished and line[25:45] == "Standard orientation":

            self.updateprogress(inputfile, "Attributes", self.cupdate)

            # If this is a counterpoise calculation, this output means that
            #   the supermolecule is now being considered, so we can set:
            self.counterpoise = 0

            if not hasattr(self, "atomcoords"):
                self.atomcoords = []
            
            hyphens = inputfile.next()
            colmNames = inputfile.next()
            colmNames = inputfile.next()
            hyphens = inputfile.next()
            
            atomnos = []
            atomcoords = []
            line = inputfile.next()
            while line != hyphens:
                broken = line.split()
                atomnos.append(int(broken[1]))
                atomcoords.append(map(float, broken[3:6]))
                line = inputfile.next()
            self.atomcoords.append(atomcoords)
            if not hasattr(self, "natom"):
                self.atomnos = numpy.array(atomnos, 'i')
                self.natom = len(self.atomnos)

        # Find the targets for SCF convergence (QM calcs).
        if line[1:44] == 'Requested convergence on RMS density matrix':

            if not hasattr(self, "scftargets"):
                self.scftargets = []

            scftargets = []
            # The RMS density matrix.
            scftargets.append(self.float(line.split('=')[1].split()[0]))
            line = inputfile.next()
            # The MAX density matrix.
            scftargets.append(self.float(line.strip().split('=')[1][:-1]))
            line = inputfile.next()
            # For G03, there's also the energy (not for G98).
            if line[1:10] == "Requested":
                scftargets.append(self.float(line.strip().split('=')[1][:-1]))

            self.scftargets.append(scftargets)

        # Extract SCF convergence information (QM calcs).
        if line[1:10] == 'Cycle   1':
                    
            if not hasattr(self, "scfvalues"):
                self.scfvalues = []

            scfvalues = []
            line = inputfile.next()
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
                    line = inputfile.next()
                # May be interupted by EOF.
                except StopIteration:
                    break

            self.scfvalues.append(scfvalues)

        # Extract SCF convergence information (AM1 calcs).
        if line[1:4] == 'It=':
                    
            self.scftargets = numpy.array([1E-7], "d") # This is the target value for the rms
            self.scfvalues = [[]]

            line = inputfile.next()
            while line.find(" Energy") == -1:
            
                if self.progress:
                    step = inputfile.tell()
                    if step != oldstep:
                        self.progress.update(step, "AM1 Convergence")
                        oldstep = step
                        
                parts = line.strip().split()
                self.scfvalues[0].append(self.float(parts[-1][:-1]))
                line = inputfile.next()

        # Note: this needs to follow the section where 'SCF Done' is used
        #   to terminate a loop when extracting SCF convergence information.
        if line[1:9] == 'SCF Done':

            if not hasattr(self, "scfenergies"):
                self.scfenergies = []

            self.scfenergies.append(utils.convertor(self.float(line.split()[4]), "hartree", "eV"))

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
            line = inputfile.next()
            if line[34:43] == "UMP4(SDQ)":
              mp4energy = self.float(line.split("=")[2])
              line = inputfile.next()
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
        #  but the output is the same as when parsing for mpenergies.
        # First turn on flag for Coupled Cluster runs.
        if line[1:23] == "Coupled Cluster theory" or line[1:8] == "CCSD(T)":

            self.coupledcluster = True
            if not hasattr(self, "ccenergies"):
                self.ccenergies = []

        # Now read the consecutive correlated energies when ,
        #  but append only the last one to ccenergies.
        # Only the highest level energy is appended - ex. CCSD(T), not CCSD.
        if self.coupledcluster and line[27:35] == "E(CORR)=":
            self.ccenergy = self.float(line.split()[3])
        if self.coupledcluster and line[1:9] == "CCSD(T)=":
            self.ccenergy = self.float(line.split()[1])
        # Append when leaving link 913
        if self.coupledcluster and line[1:16] == "Leave Link  913":
            self.ccenergies.append(utils.convertor(self.ccenergy, "hartree", "eV"))

        # Geometry convergence information.
        if line[49:59] == 'Converged?':

            if not hasattr(self, "geotargets"):
                self.geovalues = []
                self.geotargets = numpy.array([0.0, 0.0, 0.0, 0.0], "d")

            newlist = [0]*4
            for i in range(4):
                line = inputfile.next()
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

        # Charge and multiplicity.
        # If counterpoise correction is used, multiple lines match.
        # The first one contains charge/multiplicity of the whole molecule.:
        #   Charge =  0 Multiplicity = 1 in supermolecule
        #   Charge =  0 Multiplicity = 1 in fragment  1.
        #   Charge =  0 Multiplicity = 1 in fragment  2.
        if line[1:7] == 'Charge' and line.find("Multiplicity")>=0:

            parts = line.split()
            if len(parts) == 6 or parts[-1] == 'supermolecule':
                self.charge = int(parts[2])
                self.mult = int(parts[5])

        # Orbital symmetries.
        if line[1:20] == 'Orbital symmetries:' and not hasattr(self, "mosyms"):

            # For counterpoise fragments, skip these lines.
            if self.counterpoise != 0: return

            self.updateprogress(inputfile, "MO Symmetries", self.fupdate)
                    
            self.mosyms = [[]]
            line = inputfile.next()
            unres = False
            if line.find("Alpha Orbitals") == 1:
                unres = True
                line = inputfile.next()
            i = 0
            while len(line) > 18 and line[17] == '(':
                if line.find('Virtual') >= 0:
                    self.homos = numpy.array([i-1], "i") # 'HOMO' indexes the HOMO in the arrays
                parts = line[17:].split()
                for x in parts:
                    self.mosyms[0].append(self.normalisesym(x.strip('()')))
                    i += 1 
                line = inputfile.next()
            if unres:
                line = inputfile.next()
                # Repeat with beta orbital information
                i = 0
                self.mosyms.append([])
                while len(line) > 18 and line[17] == '(':
                    if line.find('Virtual')>=0:
                        self.homos.resize([2]) # Extend the array to two elements
                        self.homos[1] = i-1 # 'HOMO' indexes the HOMO in the arrays
                    parts = line[17:].split()
                    for x in parts:
                        self.mosyms[1].append(self.normalisesym(x.strip('()')))
                        i += 1
                    line = inputfile.next()

        # Alpha/Beta electron eigenvalues.
        if line[1:6] == "Alpha" and line.find("eigenvalues") >= 0:

            # For counterpoise fragments, skip these lines.
            if self.counterpoise != 0: return

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

                part = line[28:]
                i = 0
                while i*10+4 < len(part):
                    x = part[i*10:(i+1)*10]
                    self.moenergies[0].append(utils.convertor(self.float(x), "hartree", "eV"))
                    i += 1
                line = inputfile.next()            

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
                line = inputfile.next()

            self.moenergies = [numpy.array(x, "d") for x in self.moenergies]
            
        if line[1:14] == "AO basis set ":
            ## Gaussian Rev <= B.0.3
            self.gbasis = []
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
            line = inputfile.next()
            while line.strip():
                gbasis = []
                line = inputfile.next()
                while line.find("*")<0:
                    temp = line.split()
                    symtype = temp[0]
                    numgau = int(temp[1])
                    gau = []
                    for i in range(numgau):
                        temp = map(self.float, inputfile.next().split())
                        gau.append(temp)
                        
                    for i,x in enumerate(symtype):
                        newgau = [(z[0],z[i+1]) for z in gau]
                        gbasis.append( (x,newgau) )
                    line = inputfile.next() # i.e. "****" or "SP ...."
                self.gbasis.append(gbasis)
                line = inputfile.next() # i.e. "20 0" or blank line

        # Start of the IR/Raman frequency section.
        # Caution is advised here, as additional frequency blocks
        #   can be printed by Gaussian (with slightly different formats),
        #   often doubling the information printed.
        # See, for a non-standard exmaple, regression Gaussian98/test_H2.log
        if line[1:14] == "Harmonic freq":

            self.updateprogress(inputfile, "Frequency Information", self.fupdate)

            # The whole block should not have any blank lines.
            while line.strip() != "":

                # Lines with symmetries and symm. indices begin with whitespace.
                if line[1:15].strip() == "" and not line[21].isdigit():

                    if not hasattr(self, 'vibsyms'):
                        self.vibsyms = []
                    syms = line.split()
                    self.vibsyms.extend(syms)
            
                if line[1:15] == "Frequencies --":
                
                    if not hasattr(self, 'vibfreqs'):
                        self.vibfreqs = []
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
                # Remember, it is possible to have less than three columns!
                # There should be as many lines as there are atoms.
                if line[1:29] == "Atom AN      X      Y      Z":
                
                    if not hasattr(self, 'vibdisps'):
                        self.vibdisps = []
                    disps = []
                    for n in range(self.natom):
                        line = inputfile.next()
                        numbers = [float(s) for s in line[10:].split()]
                        N = len(numbers) / 3
                        if not disps:
                            for n in range(N):
                                disps.append([])
                        for n in range(N):
                            disps[n].append(numbers[n:n+3])
                    self.vibdisps.extend(disps)
                
                line = inputfile.next()

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
            parts = line[36:].split()
            self.etenergies.append(utils.convertor(self.float(parts[0]), "eV", "cm-1"))
            self.etoscs.append(self.float(parts[4].split("=")[1]))
            self.etsyms.append(line[21:36].strip())
            
            line = inputfile.next()

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
                # Changed: parse the actual coefficient
                #sqr = percent**2*2 # The fractional contribution of this CI
                #if percent < 0:
                #    sqr = -sqr
                #CIScontrib.append([(fromMO, frommoindex), (toMO, tomoindex), sqr])
                CIScontrib.append([(fromMO, frommoindex), (toMO, tomoindex), percent])
                line = inputfile.next()
            self.etsecs.append(CIScontrib)

        # Circular dichroism data.
        if line[1:52] == "<0|r|b> * <b|rxdel|0>  (Au), Rotatory Strengths (R)":

            self.etrotats = []
            inputfile.next()
            inputfile.next()
            line = inputfile.next()
            parts = line.strip().split()
            while len(parts) == 5:
                try:
                    R = self.float(parts[-1])
                except ValueError:
                    # nan or -nan if there is no first excited state
                    # (for unrestricted calculations)
                    pass
                else:
                    self.etrotats.append(R)
                line = inputfile.next()
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

            # If nmo was already parsed, check if it changed.
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
        if line[1:4] == "***" and (line[5:12] == "Overlap"
                                 or line[8:15] == "Overlap"):

            self.aooverlaps = numpy.zeros( (self.nbasis, self.nbasis), "d")
            # Overlap integrals for basis fn#1 are in aooverlaps[0]
            base = 0
            colmNames = inputfile.next()
            while base < self.nbasis:
                 
                self.updateprogress(inputfile, "Overlap", self.fupdate)
                        
                for i in range(self.nbasis-base): # Fewer lines this time
                    line = inputfile.next()
                    parts = line.split()
                    for j in range(len(parts)-1): # Some lines are longer than others
                        k = float(parts[j+1].replace("D", "E"))
                        self.aooverlaps[base+j, i+base] = k
                        self.aooverlaps[i+base, base+j] = k
                base += 5
                colmNames = inputfile.next()
            self.aooverlaps = numpy.array(self.aooverlaps, "d")                    

        # Molecular orbital coefficients (mocoeffs).
        # Essentially only produced for SCF calculations.
        # This is also the place where aonames and atombasis are parsed.
        if line[5:35] == "Molecular Orbital Coefficients" or line[5:41] == "Alpha Molecular Orbital Coefficients" or line[5:40] == "Beta Molecular Orbital Coefficients":

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
                         
                colmNames = inputfile.next()   
                if base==0 and int(colmNames.split()[0])!=1:
                    # Implies that this is a POP=REGULAR calculation
                    # and so, only aonames (not mocoeffs) will be extracted
                    self.popregular = True
                symmetries = inputfile.next()
                eigenvalues = inputfile.next()
                for i in range(self.nbasis):
                                   
                    line = inputfile.next()
                    if base == 0 and not beta: # Just do this the first time 'round
                        # Changed below from :12 to :11 to deal with Elmar Neumann's example
                        parts = line[:11].split()
                        if len(parts) > 1: # New atom
                            if i>0:
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
                         
                colmNames = inputfile.next()   
                if base==0 and int(colmNames.split()[0])!=1:
                    # Implies that this is a POP=REGULAR calculation
                    # and so, only aonames (not mocoeffs) will be extracted
                    self.popregular = True

                # No symmetry line for natural orbitals.
                # symmetries = inputfile.next()
                eigenvalues = inputfile.next()

                for i in range(self.nbasis):
                                   
                    line = inputfile.next()

                    # Just do this the first time 'round.
                    if base == 0:

                        # Changed below from :12 to :11 to deal with Elmar Neumann's example.
                        parts = line[:11].split()
                        # New atom.
                        if len(parts) > 1:
                            if i>0:
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

        # Pseudopotential charges.
        if line.find("Pseudopotential Parameters") > -1:

            dashes = inputfile.next()
            label1 = inputfile.next()
            label2 = inputfile.next()
            dashes = inputfile.next()

            line = inputfile.next()
            if line.find("Centers:") < 0:
                return
                # This was continue before parser refactoring.
                # continue

            centers = line.split()[1:]
            
            self.coreelectrons = numpy.zeros(self.natom, "i")

            for center in centers:
                while line[:10].find(center) < 0:
                    line = inputfile.next()
                
                info = line.split()
                self.coreelectrons[int(center)-1] = int(info[1]) - int(info[2])

        # This will be printed for counterpoise calcualtions only.
        # To prevent crashing, we need to know which fragment is being considered.
        # Other information is also printed in lines that start like this.
        if line[1:14] == 'Counterpoise:':
        
            if line[42:50] == "fragment":
                self.counterpoise = int(line[51:54])


if __name__ == "__main__":
    import doctest, gaussianparser
    doctest.testmod(gaussianparser, verbose=False)
