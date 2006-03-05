"""
cclib is a parser for computational chemistry log files.

See http://cclib.sf.net for more information.

Copyright (C) 2006 Noel O'Boyle and Adam Tenderholt

 This program is free software; you can redistribute and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY, without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

Contributions (monetary as well as code :-) are encouraged.
"""
import math,sys,logging,copy,re,os,time # How many of these are necessary?
import Numeric

def convertor(value,fromunits,tounits):
    """Convert from one set of units to another.

    >>> convertor(8,"eV","cm-1")
    64000
    """
    _convertor = {"eV_to_cm-1": lambda x: x*8065.6,
                  "nm_to_cm-1": lambda x: 1e7/x,
                  "cm-1_to_nm": lambda x: 1e7/x}

    return _convertor["%s_to_%s" % (fromunits,tounits)] (value)

class PeriodicTable(object):
    """Allows conversion between element name and atomic no.

    >>> t = PeriodicTable()
    >>> t.element[6]
    'C'
    >>> t.number['C']
    6
    """
    def __init__(self):
        self.element = [None,"H","He","Li","Be","B","C","N","O","F","Ne"]
        self.number = {}
        for i in range(1,len(self.element)):
            self.number[self.element[i]] = i

class Logfile(object):
    """Abstract class that contains the methods that act on data
    parsed from various types of logfile."""
    def __init__(self,filename):
        self.filename = filename

    def float(number):
        """Convert a string to a float avoiding the problem with Ds."""
        number = number.replace("D","E")
        return float(number)

class G03(Logfile):
    """A Gaussian 03 log file
    
    Attributes:
     filename -- the name of the log file
     logger -- a logging object
     NAtoms -- the number of atoms in the molecule
     atomicNo[] -- the atomic numbers of the atoms

     scfenergy[] -- the SCF energies

    Class Methods:
     float(a) -- convert a string to a float

    Methods:
     parse() -- extract general info from the logfile
    """
    SCFRMS,SCFMAX,SCFENERGY = range(3) # Used to index self.scftarget[]
    def __init__(self,filename):

        # Set up the logger...
        # Note: all loggers with the same name share the logger.
        # Here, loggers with different filenames are different
        self.logger = logging.getLogger('G03.%s' % self.filename)
        self.logger.setLevel(logging.INFO)
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("[%(name)s %(levelname)s] %(message)s"))
        self.logger.addHandler(handler)

    def __str__(self):
        """Return a string representation of the object."""
        return "Gaussian 03 log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'G03("%s")' % (self.filename)

    def parse(self):
        """Extract general info from the logfile.
        
        Creates the following instance attributes:
         NAtoms, atomicNo[], 
         scfProgress[], SCF_target_rms, SCF_target_energy, SCF_target_max,
         geoProgress[], scfenergy[]
         evalue [[]]
        """
        inputfile = open(self.filename,"r")
        for line in inputfile:
            
            if line[1:8]=="NAtoms=":
# Find the number of atoms
                NAtoms = int(line.split()[1])
                if hasattr(self,"NAtoms"):
                    assert self.NAtoms==NAtoms
                else:
                    # I wonder whether this code will ever be executed
                    self.NAtoms = NAtoms
                    self.logger.info("Creating attribute NAtoms: %d" % self.NAtoms)                    
            
            if not hasattr(self,"atomicNo") and (line.find("Z-Matrix orientation")>=0
                                                 or line[25:45]=="Standard orientation"
                                                 or line[26:43]=="Input orientation"):
# Extract the atomic numbers of the atoms
                self.logger.info("Creating attribute atomicNo[]")
                self.atomicNo = []
                hyphens = inputfile.next()
                colmNames = inputfile.next(); colmNames = inputfile.next()
                hyphens = inputfile.next()
                line = inputfile.next()
                while line!=hyphens:
                    self.atomicNo.append(int(line.split()[1]))
                    line = inputfile.next()
                NAtoms = len(self.atomicNo)
                if hasattr(self,"NAtoms"):
                    assert self.NAtoms==NAtoms
                else:
                    self.NAtoms = NAtoms
                    self.logger.info("Creating attribute NAtoms: %d" % self.NAtoms)


# Find the targets for the SCF convergence (QM calcs)
# We assume that the targets don't change, although it's
# easy enough to store all of the targets
            if line[1:44]=='Requested convergence on RMS density matrix':
                if not hasattr(self,"scftarget"):
                    self.logger.info("Creating attribute scftarget[]")
                self.scftarget = [None]*3
                self.scftarget[G03.SCFRMS] = self.float(line.split('=')[1].split()[0])
            if line[1:44]=='Requested convergence on MAX density matrix':
                self.scftarget[G03.SCFMAX] = self.float(line.strip().split('=')[1][:-1])
            if line[1:44]=='Requested convergence on             energy':
                self.scftarget[G03.SCFENERGY] = self.float(line.strip().split('=')[1][:-1])

            if line[1:10]=='Cycle   1':
# Extract SCF convergence information (QM calcs)
                if not hasattr(self,"scfvalue"):
                    self.logger.info("Creating attribute scfvalue[[]]")
                    self.scfvalue = []
                newlist = [ [] for x in self.scftarget ]
                line = inputfile.next()
                while line.find("SCF Done")==-1:
                    if line.find(' E=')==0:
                        self.logger.debug(line)
                    if line.find(" RMSDP")==0:
                        parts = line.split()
                        newlist[G03.SCFRMS].append(self.float(parts[0].split('=')[1]))
                        newlist[G03.SCFMAX].append(self.float(parts[1].split('=')[1]))
                        energy = 1.0
                        if len(parts)>4:
                            energy = parts[2].split('=')[1]
                            if energy=="":
                                energy = self.float(parts[3])
                            else:
                                energy = self.float(energy)
                        # I moved the following line back a TAB to see the effect
                        # (it was originally part of the above "if len(parts)")
                        newlist[G03.SCFENERGY].append(energy)
                        self.logger.debug(line)
                    try:
                        line = inputfile.next()
                    except StopIteration: # May be interupted by EOF
                        break
                self.scfvalue.append(newlist)

            if line[1:4]=='It=':
# Extract SCF convergence information (AM1 calcs)
                self.logger.info("Creating attributes scftarget[],scfvalue[[]]")
                self.scftarget = [1E-7] # This is the target value for the rms
                self.scfvalue = [[]]
                line = inputfile.next()
                while line.find(" Energy")==-1:
                    self.logger.debug(line)
                    parts = line.strip().split()
                    self.scfvalue[0].append(self.float(parts[-1][:-1]))
                    line = inputfile.next()

            if line[1:9]=='SCF Done':
# Note: this needs to follow the section where 'SCF Done' is used to terminate
# a loop when extract SCF convergence information
                self.logger.debug(line)
                self.logger.debug("SCF Done")
                if hasattr(self,"scfenergy"):
                    self.scfenergy.append(line.split()[4])
                else:
                    self.scfenergy = [line.split()[4]]
                    self.logger.info("Creating attribute scfenergy[]")

            if line[49:59]=='Converged?':
# Extract Geometry convergence information
                if not hasattr(self,"geotarget"):
                    self.logger.info("Creating attributes geotarget[],geovalue[[]]")
                    self.geovalue = []
                    self.geotarget = [None]*4
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
                    self.geotarget[i] = self.float(parts[3])
                self.geovalue.append(newlist)

            if line[1:19]=='Orbital symmetries' and not hasattr(self,"orbsym"):
# Extracting orbital symmetries
                self.logger.info("Creating attribute orbsym[[]]")
                self.orbsym = [[]]
                line = inputfile.next()
                unres = False
                if line.find("Alpha Orbitals")==1:
                    unres = True
                    line = inputfile.next()
                i = 0
                while len(line)>18 and line[17]=='(':
                    if line.find('Virtual')>=0:
                        self.HOMO = [i-1] # 'HOMO' indexes the HOMO in the arrays
                        self.logger.info("Creating attribute HOMO[]")
                    parts = line[17:].split()
                    for x in parts:
                        self.orbsym[0].append(x.strip('()'))
                        i+= 1 
                    line = inputfile.next()
                if unres:
                    line = inputfile.next()
                    # Repeat with beta orbital information
                    i = 0
                    self.orbsym.append([])
                    while len(line)>18 and line[17]=='(':
                        if line.find('Virtual')>=0:
                            self.HOMO.append(i-1) # 'HOMO' indexes the HOMO in the arrays
                        parts = line[17:].split()
                        for x in parts:
                            self.orbsym[1].append(x.strip('()'))
                            i+= 1
                        line = inputfile.next()

            if line[1:6]=="Alpha" and line.find("eigenvalues")>=0:
# Extract the alpha electron eigenvalues
                self.logger.info("Creating attribute evalue[[]]")
                self.evalue = [[]]
                HOMO = -2
                while line.find('Alpha')==1:
                    if line.split()[1]=="virt." and HOMO==-2:
                        # If there aren't any symmetries,
                        # this is a good way to find the HOMO
                        HOMO = len(self.evalue[0])-1
                        if hasattr(self,"HOMO"):
                            assert HOMO==self.HOMO[0]
                        else:
                            self.logger.info("Creating attribute HOMO[]")
                            self.HOMO = [HOMO]
                    part = line[28:]
                    i = 0
                    while i*10+4<len(part):
                        x = part[i*10:(i+1)*10]
                        self.evalue[0].append(self.float(x)*27.2114) # from a.u. (hartrees) to eV
                        i += 1
                    line = inputfile.next()            
                if line.find('Beta')==2:
                    self.evalue.append([])
                HOMO = -2
                while line.find('Beta')==2:
                    if line.split()[1]=="virt." and HOMO==-2:
                        # If there aren't any symmetries,
                        # this is a good way to find the HOMO
                        HOMO = len(self.evalue[1])-1
                        if len(self.HOMO)==2:
                            # It already has a self.HOMO (with the Alpha value)
                            # but does it already have a Beta value?
                            assert HOMO==self.HOMO[1]
                        else:
                             self.HOMO.append(HOMO)
                    part = line[28:]
                    i = 0
                    while i*10+4<len(part):
                        x = part[i*10:(i+1)*10]
                        self.evalue[1].append(self.float(x)*27.2114) # from a.u. (hartrees) to eV
                        i += 1
                    line = inputfile.next()   

            if line[1:14]=="Harmonic freq":
# Start of the IR/Raman frequency section
                self.vibsym = []
                self.ir = []
                self.vibfreq = []
                self.logger.info("Creating attribute vibsym[]")
                self.logger.info("Creating attribute vibfreq[]")
                self.logger.info("Creating attribute ir[]")                
                line = inputfile.next()
                while len(line[:15].split())>0:
                    # Get past the three/four line title of the columns
                    line = inputfile.next()
                line = inputfile.next() # The line with symmetries
                while len(line[:15].split())==0:
                    self.logger.debug(line)
                    self.vibsym.extend(line.split()) # Adding new symmetry
                    line = inputfile.next()
                    self.vibfreq.extend(map(self.float,line[15:].split())) # Adding new frequencies
                    [inputfile.next() for i in [0,1]] # Skip two lines
                    line = inputfile.next()
                    self.ir.extend(map(self.float,line[15:].split())) # Adding IR intensities
                    line = inputfile.next()
                    if line.find("Raman")>=0:
                        if not hasattr(self,"raman"):
                            self.raman = []
                            self.logger.info("Creating attribute raman[]")
                        line = inputfile.next()
                        self.raman.extend(map(self.float,line[15:].split())) # Adding Raman intensities
                    line = inputfile.next()
                    while len(line[:15].split())>0:
                        line = inputfile.next()
                    line = inputfile.next() # Should be the line with symmetries
                    
            if line[1:14]=="Excited State":
# Extract the electronic transitions
                if not hasattr(self,"etenergy"):
                    self.etenergy = []
                    self.etwavelen = []
                    self.etosc = []
                    self.etsym = []
                    self.etcis = []
                    self.logger.info("Creating attributes etenergy[], etwavelen[], etosc[], etsym[], etcis[]")
                # Need to deal with lines like:
                # (restricted calc)
                # Excited State   1:   Singlet-BU     5.3351 eV  232.39 nm  f=0.1695
                # (unrestricted calc) (first excited state is 2!)
                # Excited State   2:   ?Spin  -A      0.1222 eV 10148.75 nm  f=0.0000
                parts = line[36:].split()
                self.logger.debug(parts)
                self.etenergy.append(convertor(self.float(parts[0]),"eV","cm-1"))
                self.etwavelen.append(self.float(parts[2]))
                self.etosc.append(self.float(parts[4].split("=")[1]))
                self.etsym.append(line[21:36].split())
                
                line = inputfile.next()

                p = re.compile("(\d+)")
                CIScontrib = []
                while line.find(" ->")>=0: # This is a contribution to the transition
                    parts = line.split("->")
                    self.logger.debug(parts)
                    # Has to deal with lines like:
                    #       32 -> 38         0.04990
                    #      35A -> 45A        0.01921
                    frommoindex = 0 # For restricted or alpha unrestricted
                    fromMO = parts[0].strip()
                    if fromMO[-1]=="B":
                        frommoindex = 1 # For beta unrestricted
                    fromMO = int(p.match(fromMO).group()) # extract the number
                    
                    t = parts[1].split()
                    tomoindex = 0
                    toMO = t[0]
                    if toMO[-1]=="B":
                        tomoindex = 1
                    toMO = int(p.match(toMO).group())

                    percent = self.float(t[1])
                    sqr = percent**2*2 # The fractional contribution of this CI
                    if percent<0:
                        sqr = -sqr
                    CIScontrib.append([(fromMO,frommoindex),(toMO,tomoindex),sqr])
                    line = inputfile.next()
                self.etcis.append(CIScontrib)


            if line[1:52]=="<0|r|b> * <b|rxdel|0>  (Au), Rotatory Strengths (R)":
# Extract circular dichroism data
                self.rotatory = []
                self.logger.info("Creating attribute rotatory[]")
                inputfile.next()
                inputfile.next()
                line = inputfile.next()
                parts = line.strip().split()
                while len(parts)==5:
                    try:
                        R = self.float(parts[-1])
                    except ValueError:
                        # nan or -nan if there is no first excited state
                        # (for unrestricted calculations)
                        pass
                    else:
                        self.rotatory.append(R)
                    line = inputfile.next()
                    temp = line.strip().split()
                    parts = line.strip().split()                

            if line[1:7]=="NBasis" or line[4:10]=="NBasis":
# Extract the number of basis sets
                NBasis = int(line.split('=')[1].split()[0])
                # Has to deal with lines like:
                #  NBasis=   434 NAE=    97 NBE=    97 NFC=    34 NFV=     0
                #     NBasis = 148  MinDer = 0  MaxDer = 0
                # Although the former is in every file, it doesn't occur before
                # the overlap matrix is printed
                if hasattr(self,"NBasis"):
                    assert NBasis==self.NBasis
                else:
                    self.NBasis = NBasis
                    self.logger.info("Creating attribute NBasis: %d" % self.NBasis)
                    
            if line[1:7]=="NBsUse":
# Extract the number of linearly-independent basis sets
                NBsUse = int(line.split('=')[1].split()[0])
                if hasattr(self,"NBsUse"):
                    assert NBsUse==self.NBsUse
                else:
                    self.NBsUse = NBsUse
                    self.logger.info("Creating attribute NBsUse: %d" % self.NBsUse)

            if line[7:22]=="basis functions,":
# For AM1 calculations, set NBasis by a second method
# (NBsUse may not always be explicitly stated)
                    NBasis = int(line.split()[0])
                    if hasattr(self,"NBasis"):
                        assert NBasis==self.NBasis
                    else:
                        self.NBasis = NBasis
                        self.logger.info("Creating attribute NBasis: %d" % self.NBasis)

            if line[1:4]=="***" and (line[5:12]=="Overlap"
                                             or line[8:15]=="Overlap"):
# Extract the molecular orbital overlap matrix
                # Has to deal with lines such as:
                #  *** Overlap ***
                #  ****** Overlap ******
                self.logger.info("Creating attribute overlap[x,y]")
                import time; oldtime = time.time()
                self.overlap = Numeric.zeros( (self.NBasis,self.NBasis), "float")
                # Overlap integrals for basis fn#1 are in overlap[0]
                base = 0
                colmNames = inputfile.next()
                while base<self.NBasis:
                    for i in range(self.NBasis-base): # Fewer lines this time
                        line = inputfile.next()
                        parts = line.split()
                        for j in range(len(parts)-1): # Some lines are longer than others
                            k = float(parts[j].replace("D","E"))
                            self.overlap[base+j,i+base] = k
                            self.overlap[i+base,base+j] = k
                    base += 5
                    colmNames = inputfile.next()
                self.logger.info("Took %f seconds" % (time.time()-oldtime))

            if line[5:35]=="Molecular Orbital Coefficients" or line[5:41]=="Alpha Molecular Orbital Coefficients" or line[5:40]=="Beta Molecular Orbital Coefficients":
                import time; oldtime = time.time()
                if line[5:40]=="Beta Molecular Orbital Coefficients":
                    beta = True
                    # Need to add an extra dimension to self.mocoeff
                    self.mocoeff = Numeric.resize(self.mocoeff,(2,NBsUse,NBasis))
                else:
                    beta = False
                    self.logger.info("Creating attributes orbitals[], mocoeff[][]")
                    self.orbitals = []
                    self.mocoeff = Numeric.zeros((NBsUse,NBasis),"float")

                base = 0
                for base in range(0,NBsUse,5):
                    colmNames = inputfile.next()
                    symmetries = inputfile.next()
                    eigenvalues = inputfile.next()
                    for i in range(NBasis):
                        line = inputfile.next()
                        if base==0 and not beta: # Just do this the first time 'round
                            # Changed below from :12 to :11 to deal with Elmar Neumann's example
                            parts = line[:11].split()
                            if len(parts)>1: # New atom
                                atomname = "%s%s" % (parts[2],parts[1])
                            orbital = line[11:20].strip()
                            self.orbitals.append("%s_%s" % (atomname,orbital))

                        part = line[21:].replace("D","E").rstrip()
                        temp = [] 
                        for j in range(0,len(part),10):
                            temp.append(float(part[j:j+10]))
                        if beta:
                            self.mocoeff[1,base:base+len(part)/10,i] = temp
                        else:
                            self.mocoeff[base:base+len(part)/10,i] = temp
                self.logger.info("Took %f seconds" % (time.time()-oldtime))

                 
        inputfile.close()

# Note to self: Needs to be added to the main parser
    def extractTrajectory(self):
        """Extract trajectory information from a Gaussian logfile."""
        inputfile = open(self.filename,"r")
        self.traj = []
        self.trajSummary = []
        for line in inputfile:
            if line.find(" Cartesian coordinates:")==0:
                coords = []
                for i in range(self.NAtoms):
                    line = inputfile.next()
                    parts = line.strip().split()
                    # Conversion from a.u. to Angstrom
                    coords.append([ self.float(x)*0.5292 for x in [parts[3],parts[5],parts[7]] ])
                self.traj.append(coords)
            if line==" Trajectory summary\n":
                # self.trajSummaryHeader = inputfile.next().strip().split()
                header = inputfile.next()
                line = inputfile.next()
                while line.find("Max Error")==-1:
                    parts = line.strip().split("  ")
                    self.trajSummary.append(parts)
                    line = inputfile.next()
        inputfile.close()
        assert len(self.traj)==len(self.trajSummary)
        
if __name__=="__main__":
    import doctest,parser
    doctest.testmod(parser,verbose=False)
