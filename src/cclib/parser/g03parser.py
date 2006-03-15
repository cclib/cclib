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
import re,time
import Numeric
import random # For sometimes running the progress updater
from logfileparser import Logfile # import the superclass

class G03(Logfile):
    """A Gaussian 98/03 log file"""
    SCFRMS,SCFMAX,SCFENERGY = range(3) # Used to index self.scftargets[]
    def __init__(self,*args):

        # Call the __init__ method of the superclass
        super(G03, self).__init__(*args,logname="G03")
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Gaussian 03 log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'G03("%s")' % (self.filename)

    def parse(self):
        """Extract information from the logfile."""
        inputfile = open(self.filename,"r")
        
        if self.progress:
            
            inputfile.seek(0,2) #go to end of file
            nstep=inputfile.tell()
            inputfile.seek(0)
            self.progress.initialize(nstep)
            oldstep=0
            
        for line in inputfile:
            
            if self.progress and random.random()<0.05:
                
                step = inputfile.tell()
                if step!=oldstep:
                    self.progress.update(step)
                    oldstep = step
                
            if line[1:8]=="NAtoms=":
# Find the number of atoms
                natom = int(line.split()[1])
                if hasattr(self,"natom"):
                    assert self.natom==natom
                else:
                    # I wonder whether this code will ever be executed
                    self.natom = natom
                    self.logger.info("Creating attribute natom: %d" % self.natom)                    
            
            if not hasattr(self,"atomnos") and (line.find("Z-Matrix orientation")>=0
                                                 or line[25:45]=="Standard orientation"
                                                 or line[26:43]=="Input orientation"):
# Extract the atomic numbers of the atoms
                self.logger.info("Creating attribute atomnos[]")
                self.atomnos = []
                hyphens = inputfile.next()
                colmNames = inputfile.next(); colmNames = inputfile.next()
                hyphens = inputfile.next()
                line = inputfile.next()
                while line!=hyphens:
                    self.atomnos.append(int(line.split()[1]))
                    line = inputfile.next()
                natom = len(self.atomnos)
                if hasattr(self,"natom"):
                    assert self.natom==natom
                else:
                    self.natom = natom
                    self.logger.info("Creating attribute natom: %d" % self.natom)


# Find the targets for the SCF convergence (QM calcs)
# We assume that the targets don't change, although it's
# easy enough to store all of the targets
            if line[1:44]=='Requested convergence on RMS density matrix':
                if not hasattr(self,"scftargets"):
                    self.logger.info("Creating attribute scftargets[]")
                self.scftargets = [None]*3
                self.scftargets[G03.SCFRMS] = self.float(line.split('=')[1].split()[0])
            if line[1:44]=='Requested convergence on MAX density matrix':
                self.scftargets[G03.SCFMAX] = self.float(line.strip().split('=')[1][:-1])
            if line[1:44]=='Requested convergence on             energy':
                self.scftargets[G03.SCFENERGY] = self.float(line.strip().split('=')[1][:-1])

            if line[1:10]=='Cycle   1':
# Extract SCF convergence information (QM calcs)
                if not hasattr(self,"scfvalues"):
                    self.logger.info("Creating attribute scfvalues[[]]")
                    self.scfvalues = []
                newlist = [ [] for x in self.scftargets ]
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
                self.scfvalues.append(newlist)

            if line[1:4]=='It=':
# Extract SCF convergence information (AM1 calcs)
                self.logger.info("Creating attributes scftargets[],scfvalues[[]]")
                self.scftargets = [1E-7] # This is the target value for the rms
                self.scfvalues = [[]]
                line = inputfile.next()
                while line.find(" Energy")==-1:
                    self.logger.debug(line)
                    parts = line.strip().split()
                    self.scfvalues[0].append(self.float(parts[-1][:-1]))
                    line = inputfile.next()

            if line[1:9]=='SCF Done':
# Note: this needs to follow the section where 'SCF Done' is used to terminate
# a loop when extract SCF convergence information
                self.logger.debug(line)
                self.logger.debug("SCF Done")
                if hasattr(self,"scfenergies"):
                    self.scfenergies.append(line.split()[4])
                else:
                    self.scfenergies = [line.split()[4]]
                    self.logger.info("Creating attribute scfenergies[]")

            if line[49:59]=='Converged?':
# Extract Geometry convergence information
                if not hasattr(self,"geotargets"):
                    self.logger.info("Creating attributes geotargets[],geovalues[[]]")
                    self.geovalues = []
                    self.geotargets = [None]*4
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

            if line[1:19]=='Orbital symmetries' and not hasattr(self,"mosyms"):
# Extracting orbital symmetries
                self.logger.info("Creating attribute mosyms[[]]")
                self.mosyms = [[]]
                line = inputfile.next()
                unres = False
                if line.find("Alpha Orbitals")==1:
                    unres = True
                    line = inputfile.next()
                i = 0
                while len(line)>18 and line[17]=='(':
                    if line.find('Virtual')>=0:
                        self.homos = [i-1] # 'HOMO' indexes the HOMO in the arrays
                        self.logger.info("Creating attribute homos[]")
                    parts = line[17:].split()
                    for x in parts:
                        self.mosyms[0].append(x.strip('()'))
                        i+= 1 
                    line = inputfile.next()
                if unres:
                    line = inputfile.next()
                    # Repeat with beta orbital information
                    i = 0
                    self.mosyms.append([])
                    while len(line)>18 and line[17]=='(':
                        if line.find('Virtual')>=0:
                            self.homos.append(i-1) # 'HOMO' indexes the HOMO in the arrays
                        parts = line[17:].split()
                        for x in parts:
                            self.mosyms[1].append(x.strip('()'))
                            i+= 1
                        line = inputfile.next()

            if line[1:6]=="Alpha" and line.find("eigenvalues")>=0:
# Extract the alpha electron eigenvalues
                self.logger.info("Creating attribute moenergies[[]]")
                self.moenergies = [[]]
                HOMO = -2
                while line.find('Alpha')==1:
                    if line.split()[1]=="virt." and HOMO==-2:
                        # If there aren't any symmetries,
                        # this is a good way to find the HOMO
                        HOMO = len(self.moenergies[0])-1
                        if hasattr(self,"homos"):
                            assert HOMO==self.homos[0]
                        else:
                            self.logger.info("Creating attribute homos[]")
                            self.homos = [HOMO]
                    part = line[28:]
                    i = 0
                    while i*10+4<len(part):
                        x = part[i*10:(i+1)*10]
                        self.moenergies[0].append(self.float(x)*27.2114) # from a.u. (hartrees) to eV
                        i += 1
                    line = inputfile.next()            
                if line.find('Beta')==2:
                    self.moenergies.append([])
                HOMO = -2
                while line.find('Beta')==2:
                    if line.split()[1]=="virt." and HOMO==-2:
                        # If there aren't any symmetries,
                        # this is a good way to find the HOMO
                        HOMO = len(self.moenergies[1])-1
                        if len(self.homos)==2:
                            # It already has a self.homos (with the Alpha value)
                            # but does it already have a Beta value?
                            assert HOMO==self.homos[1]
                        else:
                             self.homos.append(HOMO)
                    part = line[28:]
                    i = 0
                    while i*10+4<len(part):
                        x = part[i*10:(i+1)*10]
                        self.moenergies[1].append(self.float(x)*27.2114) # from a.u. (hartrees) to eV
                        i += 1
                    line = inputfile.next()   

            if line[1:14]=="Harmonic freq":
# Start of the IR/Raman frequency section
                self.vibsyms = []
                self.vibirs = []
                self.vibfreqs = []
                self.logger.info("Creating attribute vibsyms[]")
                self.logger.info("Creating attribute vibfreqs[]")
                self.logger.info("Creating attribute vibirs[]")                
                line = inputfile.next()
                while len(line[:15].split())>0:
                    # Get past the three/four line title of the columns
                    line = inputfile.next()
                line = inputfile.next() # The line with symmetries
                while len(line[:15].split())==0:
                    self.logger.debug(line)
                    self.vibsyms.extend(line.split()) # Adding new symmetry
                    line = inputfile.next()
                    self.vibfreqs.extend(map(self.float,line[15:].split())) # Adding new frequencies
                    [inputfile.next() for i in [0,1]] # Skip two lines
                    line = inputfile.next()
                    self.vibirs.extend(map(self.float,line[15:].split())) # Adding IR intensities
                    line = inputfile.next()
                    if line.find("Raman")>=0:
                        if not hasattr(self,"vibramans"):
                            self.vibramans = []
                            self.logger.info("Creating attribute vibramans[]")
                        line = inputfile.next()
                        self.vibramans.extend(map(self.float,line[15:].split())) # Adding Raman intensities
                    line = inputfile.next()
                    while len(line[:15].split())>0:
                        line = inputfile.next()
                    line = inputfile.next() # Should be the line with symmetries
                    
            if line[1:14]=="Excited State":
# Extract the electronic transitions
                if not hasattr(self,"etenergy"):
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
                self.etenergies.append(convertor(self.float(parts[0]),"eV","cm-1"))
                self.etoscs.append(self.float(parts[4].split("=")[1]))
                self.etsyms.append(line[21:36].split())
                
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
                self.etsecs.append(CIScontrib)


            if line[1:52]=="<0|r|b> * <b|rxdel|0>  (Au), Rotatory Strengths (R)":
# Extract circular dichroism data
                self.etrotats = []
                self.logger.info("Creating attribute etrotats[]")
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
                        self.etrotats.append(R)
                    line = inputfile.next()
                    temp = line.strip().split()
                    parts = line.strip().split()                

            if line[1:7]=="NBasis" or line[4:10]=="NBasis":
# Extract the number of basis sets
                nbasis = int(line.split('=')[1].split()[0])
                # Has to deal with lines like:
                #  NBasis =   434 NAE=    97 NBE=    97 NFC=    34 NFV=     0
                #     NBasis = 148  MinDer = 0  MaxDer = 0
                # Although the former is in every file, it doesn't occur before
                # the overlap matrix is printed
                if hasattr(self,"nbasis"):
                    assert nbasis==self.nbasis
                else:
                    self.nbasis= nbasis
                    self.logger.info("Creating attribute nbasis: %d" % self.nbasis)
                    
            if line[1:7]=="NBsUse":
# Extract the number of linearly-independent basis sets
                nindep = int(line.split('=')[1].split()[0])
                if hasattr(self,"nindep"):
                    assert nindep==self.nindep
                else:
                    self.nindep = nindep
                    self.logger.info("Creating attribute nindep: %d" % self.nindep)

            if line[7:22]=="basis functions,":
# For AM1 calculations, set nbasis by a second method
# (nindep may not always be explicitly stated)
                    nbasis = int(line.split()[0])
                    if hasattr(self,"nbasis"):
                        assert nbasis==self.nbasis
                    else:
                        self.nbasis = nbasis
                        self.logger.info("Creating attribute nbasis: %d" % self.nbasis)

            if line[1:4]=="***" and (line[5:12]=="Overlap"
                                             or line[8:15]=="Overlap"):
# Extract the molecular orbital overlap matrix
                # Has to deal with lines such as:
                #  *** Overlap ***
                #  ****** Overlap ******
                self.logger.info("Creating attribute aooverlaps[x,y]")
                # oldtime = time.time()
                self.aooverlaps = Numeric.zeros( (self.nbasis,self.nbasis), "float")
                # Overlap integrals for basis fn#1 are in aooverlaps[0]
                base = 0
                colmNames = inputfile.next()
                while base<self.nbasis:
                    for i in range(self.nbasis-base): # Fewer lines this time
                        line = inputfile.next()
                        parts = line.split()
                        for j in range(len(parts)-1): # Some lines are longer than others
                            k = float(parts[j].replace("D","E"))
                            self.aooverlaps[base+j,i+base] = k
                            self.aooverlaps[i+base,base+j] = k
                    base += 5
                    colmNames = inputfile.next()
                # self.logger.info("Took %f seconds" % (time.time()-oldtime))

            if line[5:35]=="Molecular Orbital Coefficients" or line[5:41]=="Alpha Molecular Orbital Coefficients" or line[5:40]=="Beta Molecular Orbital Coefficients":
                # oldtime = time.time()
                if line[5:40]=="Beta Molecular Orbital Coefficients":
                    beta = True
                    # Need to add an extra dimension to self.mocoeffs
                    self.mocoeffs = Numeric.resize(self.mocoeffs,(2,nindep,nbasis))
                else:
                    beta = False
                    self.logger.info("Creating attributes aonames[], mocoeffs[][]")
                    self.aonames = []
                    self.mocoeffs = Numeric.zeros((nindep,nbasis),"float")

                base = 0
                for base in range(0,nindep,5):
                    colmNames = inputfile.next()
                    symmetries = inputfile.next()
                    eigenvalues = inputfile.next()
                    for i in range(nbasis):
                        line = inputfile.next()
                        if base==0 and not beta: # Just do this the first time 'round
                            # Changed below from :12 to :11 to deal with Elmar Neumann's example
                            parts = line[:11].split()
                            if len(parts)>1: # New atom
                                atomname = "%s%s" % (parts[2],parts[1])
                            orbital = line[11:20].strip()
                            self.aonames.append("%s_%s" % (atomname,orbital))

                        part = line[21:].replace("D","E").rstrip()
                        temp = [] 
                        for j in range(0,len(part),10):
                            temp.append(float(part[j:j+10]))
                        if beta:
                            self.mocoeffs[1,base:base+len(part)/10,i] = temp
                        else:
                            self.mocoeffs[base:base+len(part)/10,i] = temp
                # self.logger.info("Took %f seconds" % (time.time()-oldtime))
                 
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
                for i in range(self.natom):
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
