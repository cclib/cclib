"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

import re
import Numeric
import random # For sometimes running the progress updater
import utils
import logfileparser

class Gaussian(logfileparser.Logfile):
    """A Gaussian 98/03 log file"""
    SCFRMS, SCFMAX, SCFENERGY = range(3) # Used to index self.scftargets[]
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
        >>> labels = ['A1', 'AG', 'A1G', "SG", "PI", "PHI", "DLTA"]
        >>> map(sym, labels)
        ['A1', 'Ag', 'A1g', 'sigma', 'pi', 'phi', 'delta']
        """
        _greek = {'SG':'sigma', 'PI':'pi', 'PHI':'phi', 'DLTA':'delta'}
        if label in _greek:
            return _greek[label]

        ans = label.replace("U", "u").replace("G", "g") 
        return ans

    def parse(self, fupdate=0.05, cupdate=0.002):
        """Extract information from the logfile."""
        inputfile = utils.openlogfile(self.filename)
        
        if self.progress:
            
            inputfile.seek(0, 2) #go to end of file
            nstep = inputfile.tell()
            inputfile.seek(0)
            self.progress.initialize(nstep)
            oldstep = 0
            
        optfinished = False # Flag that indicates whether it has reached the end of a geoopt
        
        for line in inputfile:
            
            if self.progress and random.random() < cupdate:
                
                step = inputfile.tell()
                if step != oldstep:
                    self.progress.update(step, "Unsupported Information")
                    oldstep = step
                
            if line[1:8] == "NAtoms=":
# Find the number of atoms
                if self.progress and random.random() < fupdate:
                    step = inputfile.tell()
                    if step != oldstep:
                        self.progress.update(step, "Attributes")
                        oldstep = step
                        
                natom = int(line.split()[1])
                if hasattr(self, "natom"):
                    assert self.natom == natom
                else:
                    # I wonder whether this code will ever be executed
                    self.natom = natom
                    self.logger.info("Creating attribute natom: %d" % self.natom)                    
            
            if line[1:23] == "Optimization completed":
                optfinished = True
            
            if line.find("Input orientation") > -1 or line.find("Z-Matrix orientation") > -1:
# Extract the atomic numbers and coordinates in the event standard orientation isn't available

                if self.progress and random.random() < cupdate:
                    step = inputfile.tell()
                    if step != oldstep:
                        self.progress.update(step, "Attributes")
                        oldstep = step
                        
                inputcoords = []
                inputatoms = []
                
                hyphens = inputfile.next()
                colmNames = inputfile.next()
                colmNames = inputfile.next()
                hyphens = inputfile.next()
                
                atomnos = []
                atomcoords = []
                line = inputfile.next()
                while line != hyphens:
                    broken = line.split()
                    inputatoms.append(int(broken[1]))
                    atomcoords.append(map(float, broken[3:6]))
                    line = inputfile.next()
                inputcoords.append(atomcoords)
                if not hasattr(self, "natom"):
                    self.atomnos = Numeric.array(inputatoms, 'i')
                    self.logger.info("Creating attribute atomnos[]")
                    self.natom = len(self.atomnos)
                    self.logger.info("Creating attribute natom: %d" % self.natom)                

            if not optfinished and line[25:45] == "Standard orientation":
# Extract the atomic numbers and coordinates of the atoms
                if self.progress and random.random() < cupdate:
                    step = inputfile.tell()
                    if step != oldstep:
                        self.progress.update(step, "Attributes")
                        oldstep = step
                        
                if not hasattr(self, "atomcoords"):
                    self.logger.info("Creating attribute atomcoords[]")
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
                    self.atomnos = Numeric.array(atomnos, 'i')
                    self.logger.info("Creating attribute atomnos[]")
                    self.natom = len(self.atomnos)
                    self.logger.info("Creating attribute natom: %d" % self.natom)

            if line[1:44] == 'Requested convergence on RMS density matrix':
# Find the targets for SCF convergence (QM calcs)                
                if not hasattr(self, "scftargets"):
                    self.logger.info("Creating attribute scftargets[]")
                    self.scftargets = []
                scftargets = []
                # The RMS density matrix
                scftargets.append(self.float(line.split('=')[1].split()[0]))
                line = inputfile.next()
                # The MAX density matrix
                scftargets.append(self.float(line.strip().split('=')[1][:-1]))
                line = inputfile.next()
                # For G03, there's also the energy (not for G98)
                if line[1:10] == "Requested":
                    scftargets.append(self.float(line.strip().split('=')[1][:-1]))
                self.scftargets.append(scftargets)

            if line[1:10] == 'Cycle   1':
# Extract SCF convergence information (QM calcs)
                        
                if not hasattr(self, "scfvalues"):
                    self.logger.info("Creating attribute scfvalues")
                    self.scfvalues = []
                scfvalues = []
                line = inputfile.next()
                while line.find("SCF Done") == -1:
                
                    if self.progress and random.random() < fupdate:
                        step = inputfile.tell()
                        if step != oldstep:
                            self.progress.update(step, "QM Convergence")
                            oldstep = step                
                          
                    if line.find(' E=') == 0:
                        self.logger.debug(line)
                    if line.find(" RMSDP") == 0:
#  RMSDP=3.74D-06 MaxDP=7.27D-05 DE=-1.73D-07 OVMax= 3.67D-05
# or
#  RMSDP=1.13D-05 MaxDP=1.08D-04              OVMax= 1.66D-04
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
                    except StopIteration: # May be interupted by EOF
                        break
                self.scfvalues.append(scfvalues)

            if line[1:4] == 'It=':
# Extract SCF convergence information (AM1 calcs)
                        
                self.logger.info("Creating attributes scftargets, scfvalues")
                self.scftargets = Numeric.array([1E-7], "f") # This is the target value for the rms
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

            if line[1:9] == 'SCF Done':
# Note: this needs to follow the section where 'SCF Done' is used to terminate
# a loop when extracting SCF convergence information
                if not hasattr(self, "scfenergies"):
                    self.logger.info("Creating attribute scfenergies[]")
                    self.scfenergies = []
                self.scfenergies.append(utils.convertor(self.float(line.split()[4]), "hartree", "eV"))

            if line[49:59] == 'Converged?':
# Extract Geometry convergence information
                if not hasattr(self, "geotargets"):
                    self.logger.info("Creating attributes geotargets[], geovalues[[]]")
                    self.geovalues = []
                    self.geotargets = Numeric.array([0.0, 0.0, 0.0, 0.0], "f")
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

            if line[1:19] == 'Orbital symmetries' and not hasattr(self, "mosyms"):
# Extracting orbital symmetries
                if self.progress and random.random() < fupdate:
                    step = inputfile.tell()
                    if step != oldstep:
                        self.progress.update(step, "MO Symmetries")
                        oldstep = step
                        
                self.logger.info("Creating attribute mosyms[[]]")
                self.mosyms = [[]]
                line = inputfile.next()
                unres = False
                if line.find("Alpha Orbitals") == 1:
                    unres = True
                    line = inputfile.next()
                i = 0
                while len(line) > 18 and line[17] == '(':
                    if line.find('Virtual') >= 0:
                        self.homos = Numeric.array([i-1], "i") # 'HOMO' indexes the HOMO in the arrays
                        self.logger.info("Creating attribute homos[]")
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

            if line[1:6] == "Alpha" and line.find("eigenvalues") >= 0:
# Extract the alpha electron eigenvalues
                if self.progress and random.random() < fupdate:
                    step = inputfile.tell()
                    if step != oldstep:
                        self.progress.update(step, "Eigenvalues")
                        oldstep = step
                        
                self.logger.info("Creating attribute moenergies[[]]")
                self.moenergies = [[]]
                HOMO = -2
                while line.find('Alpha') == 1:
                    if line.split()[1] == "virt." and HOMO == -2:
                        # If there aren't any symmetries,
                        # this is a good way to find the HOMO
                        HOMO = len(self.moenergies[0])-1
                        if hasattr(self, "homos"):
                            assert HOMO == self.homos[0]
                        else:
                            self.logger.info("Creating attribute homos[]")
                            self.homos = Numeric.array([HOMO], "i")
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
                        # If there aren't any symmetries,
                        # this is a good way to find the HOMO
                        HOMO = len(self.moenergies[1])-1
                        if len(self.homos) == 2:
                            # It already has a self.homos (with the Alpha value)
                            # but does it already have a Beta value?
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
                self.moenergies = Numeric.array(self.moenergies, "f")                    
            if line[1:14] == "AO basis set ":
                ## Gaussian Rev <= B.0.3
                self.logger.info("Creating attribute gbasis")
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

            if line[1:14] == "Harmonic freq":
# Start of the IR/Raman frequency section
                if self.progress and random.random() < fupdate:
                    step = inputfile.tell()
                    if step != oldstep:
                        self.progress.update(step, "Frequency Information")
                        oldstep = step
                        
                self.vibsyms = []
                self.vibirs = []
                self.vibfreqs = []
                self.vibcarts = []
                self.logger.info("Creating attribute vibsyms[], vibfreqs[]")
                self.logger.info("Creating attribute vibirs[], vibcarts[]")                
                line = inputfile.next()
                while len(line[:15].split()) > 0:
                    # Get past the three/four line title of the columns
                    line = inputfile.next()
                line = inputfile.next() # The line with symmetries
                while len(line[:15].split()) == 0:
                    self.logger.debug(line)
                    self.vibsyms.extend(line.split()) # Adding new symmetry
                    line = inputfile.next()
                    self.vibfreqs.extend(map(self.float, line[15:].split())) # Adding new frequencies
                    line = inputfile.next()
                    line = inputfile.next()
                    line = inputfile.next()
                    self.vibirs.extend(map(self.float, line[15:].split())) # Adding IR intensities
                    line = inputfile.next() # Either the header or a Raman line
                    if line.find("Raman") >= 0:
                        if not hasattr(self, "vibramans"):
                            self.vibramans = []
                            self.logger.info("Creating attribute vibramans[]")
                        self.vibramans.extend(map(self.float, line[15:].split())) # Adding Raman intensities
                        line = inputfile.next() # Depolar (P)
                        line = inputfile.next() # Depolar (U)
                        line = inputfile.next() # Header
                    line = inputfile.next() # First line of cartesian displacement vectors
                    p = [[], [], []]
                    while len(line[:15].split()) > 0:
                        # Store the cartesian displacement vectors
                        broken = map(float, line.strip().split()[2:])
                        for i in range(0, len(broken), 3):
                            p[i/3].append(broken[i:i+3])
                        line = inputfile.next()
                    self.vibcarts.extend(p[0:len(broken)/3])
                    line = inputfile.next() # Should be the line with symmetries
                self.vibfreqs = Numeric.array(self.vibfreqs, "f")
                self.vibirs = Numeric.array(self.vibirs, "f")
                self.vibcarts = Numeric.array(self.vibcarts, "f")
                if hasattr(self, "vibramans"):
                    self.vibramans = Numeric.array(self.vibramans, "f")
                    
            if line[1:14] == "Excited State":
# Extract the electronic transitions
                if not hasattr(self, "etenergies"):
                    self.etenergies = []
                    self.etoscs = []
                    self.etsyms = []
                    self.etsecs = []
                    self.logger.info("Creating attributes etenergies[], etoscs[], etsyms[], etsecs[]")
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
                    sqr = percent**2*2 # The fractional contribution of this CI
                    if percent < 0:
                        sqr = -sqr
                    CIScontrib.append([(fromMO, frommoindex), (toMO, tomoindex), sqr])
                    line = inputfile.next()
                self.etsecs.append(CIScontrib)

            if line[1:52] == "<0|r|b> * <b|rxdel|0>  (Au), Rotatory Strengths (R)":
# Extract circular dichroism data
                self.etrotats = []
                self.logger.info("Creating attribute etrotats[]")
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
                self.etrotats = Numeric.array(self.etrotats, "f")

            if line[1:7] == "NBasis" or line[4:10] == "NBasis":
# Extract the number of basis sets
                nbasis = int(line.split('=')[1].split()[0])
                # Has to deal with lines like:
                #  NBasis =   434 NAE=    97 NBE=    97 NFC=    34 NFV=     0
                #     NBasis = 148  MinDer = 0  MaxDer = 0
                # Although the former is in every file, it doesn't occur before
                # the overlap matrix is printed
                if hasattr(self, "nbasis"):
                    assert nbasis == self.nbasis
                else:
                    self.nbasis = nbasis
                    self.logger.info("Creating attribute nbasis: %d" % self.nbasis)
                    
            if line[1:7] == "NBsUse":
# Extract the number of linearly-independent basis sets
                nmo = int(line.split('=')[1].split()[0])
                if hasattr(self, "nmo"):
                    assert nmo == self.nmo
                else:
                    self.nmo = nmo
                    self.logger.info("Creating attribute nmo: %d" % self.nmo)

            if line[7:22] == "basis functions, ":
# For AM1 calculations, set nbasis by a second method
# (nmo may not always be explicitly stated)
                nbasis = int(line.split()[0])
                if hasattr(self, "nbasis"):
                    assert nbasis == self.nbasis
                else:
                    self.nbasis = nbasis
                    self.logger.info("Creating attribute nbasis: %d" % self.nbasis)

            if line[1:4] == "***" and (line[5:12] == "Overlap"
                                     or line[8:15] == "Overlap"):
# Extract the molecular orbital overlap matrix
                # Has to deal with lines such as:
                #  *** Overlap ***
                #  ****** Overlap ******
                self.logger.info("Creating attribute aooverlaps[x, y]")
                self.aooverlaps = Numeric.zeros( (self.nbasis, self.nbasis), "float")
                # Overlap integrals for basis fn#1 are in aooverlaps[0]
                base = 0
                colmNames = inputfile.next()
                while base < self.nbasis:
                     
                    if self.progress and random.random() < fupdate:
                        step = inputfile.tell()
                        if step != oldstep:
                            self.progress.update(step, "Overlap")
                            oldstep = step
                            
                    for i in range(self.nbasis-base): # Fewer lines this time
                        line = inputfile.next()
                        parts = line.split()
                        for j in range(len(parts)-1): # Some lines are longer than others
                            k = float(parts[j+1].replace("D", "E"))
                            self.aooverlaps[base+j, i+base] = k
                            self.aooverlaps[i+base, base+j] = k
                    base += 5
                    colmNames = inputfile.next()
                self.aooverlaps = Numeric.array(self.aooverlaps, "f")                    


            if line[5:35] == "Molecular Orbital Coefficients" or line[5:41] == "Alpha Molecular Orbital Coefficients" or line[5:40] == "Beta Molecular Orbital Coefficients":
                if line[5:40] == "Beta Molecular Orbital Coefficients":
                    beta = True
                    # Need to add an extra dimension to self.mocoeffs
                    self.mocoeffs = Numeric.resize(self.mocoeffs, (2, nmo, nbasis))
                else:
                    beta = False
                    self.logger.info("Creating attributes aonames[], mocoeffs[][]")
                    self.aonames = []
                    self.mocoeffs = Numeric.zeros((1, nmo, nbasis), "float")

                base = 0
                for base in range(0, nmo, 5):
                    
                    if self.progress:
                        step = inputfile.tell()
                        if step != oldstep and random.random() < fupdate:
                            self.progress.update(step, "Coefficients")
                            oldstep = step
                             
                    colmNames = inputfile.next()
                    symmetries = inputfile.next()
                    eigenvalues = inputfile.next()
                    for i in range(nbasis):
                    
                    
                        line = inputfile.next()
                        if base == 0 and not beta: # Just do this the first time 'round
                            # Changed below from :12 to :11 to deal with Elmar Neumann's example
                            parts = line[:11].split()
                            if len(parts) > 1: # New atom
                                atomname = "%s%s" % (parts[2], parts[1])
                            orbital = line[11:20].strip()
                            self.aonames.append("%s_%s" % (atomname, orbital))

                        part = line[21:].replace("D", "E").rstrip()
                        temp = [] 
                        for j in range(0, len(part), 10):
                            temp.append(float(part[j:j+10]))
                        if beta:
                            self.mocoeffs[1, base:base + len(part) / 10, i] = temp
                        else:
                            self.mocoeffs[0, base:base + len(part) / 10, i] = temp

            if line.find("Pseudopotential Parameters") > -1:
#parse pseudopotential charges

                dashes = inputfile.next()
                label1 = inputfile.next()
                label2 = inputfile.next()
                dashes = inputfile.next()

                line = inputfile.next()
                if line.find("Centers:") < 0:
                    continue

                centers = line.split()[1:]
                
                self.logger.info("Creating attribute coreelectrons[]")
                self.coreelectrons = Numeric.zeros(self.natom)

                for center in centers:
                    while line[:10].find(center) < 0:
                        line = inputfile.next()
                    
                    info = line.split()
                    self.coreelectrons[int(center)-1] = int(info[1]) - int(info[2])

        inputfile.close()

        if self.progress:
            self.progress.update(nstep, "Done")
            
        _toarray = ['geovalues', 'scfenergies', 'scftargets', 'atomcoords',
                     'etenergies', 'etoscs']
        for attr in _toarray:
            if hasattr(self, attr):
                setattr(self, attr, Numeric.array(getattr(self, attr), 'f'))

        if hasattr(self, "scfvalues"):
            self.scfvalues = [Numeric.array(x, "f") for x in self.scfvalues]

        if not hasattr(self,"atomcoords"):
            self.logger.info("Creating attribute atomcoords[]")
            self.atomcoords = Numeric.array(inputcoords, 'f')

        if not hasattr(self,"coreelectrons"):
            self.logger.info("Creating attribute coreelectrons[]")
            self.coreelectrons = Numeric.zeros(self.natom)

        self.parsed = True

if __name__ == "__main__":
    import doctest, gaussianparser
    doctest.testmod(gaussianparser, verbose=False)
