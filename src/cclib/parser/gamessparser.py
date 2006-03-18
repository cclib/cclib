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
from logfileparser import Logfile,convertor

class GAMESS(Logfile):
    """A GAMESS log file."""
    SCFRMS,SCFMAX,SCFENERGY = range(3) # Used to index self.scftargets[]
    def __init__(self,*args):

        # Call the __init__ method of the superclass
        super(GAMESS, self).__init__(logname="GAMESS",*args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "GAMESS log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'GAMESS("%s")' % (self.filename)

    def parse(self):
        """Extract information from the logfile."""
        inputfile = open(self.filename,"r")
        
        if self.progress:
            
            inputfile.seek(0,2) #go to end of file
            nstep=inputfile.tell()
            inputfile.seek(0)
            self.progress.initialize(nstep)
            oldstep=0


        endofopt = False
            
        for line in inputfile:
            
            if self.progress and random.random()<0.05:
                
                step = inputfile.tell()
                if step!=oldstep:
                    self.progress.update(step)
                    oldstep = step

            if line.find("OPTTOL")>=0:
                # Two possibilities:
                #           OPTTOL = 1.000E-04          RMIN   = 1.500E-03
                # INPUT CARD> $STATPT OPTTOL=0.0001 NSTEP=100 $END
                if not hasattr(self,"geotargets"):
                    self.logger.info("Creating attribute geotargets[]")
                    temp = line.split()
                    for i,x in enumerate(temp):
                        if x.find("OPTTOL")>=0:
                            if x=="OPTTOL":
                                opttol = float(temp[i+2])
                            else:
                                opttol = float(x.split('=')[1])
                            self.geotargets = Numeric.array([opttol,3./opttol])
                            
            if line.find("EQUILIBRIUM GEOMETRY LOCATED")>=0:
# This is necessary if a frequency calculation follows a geometry optimisation
                endofopt = True
                
            if not endofopt and line.find("FINAL")==1:
                if not hasattr(self,"scfenergies"):
                    self.logger.info("Creating attribute scfenergies[]")
                    self.scfenergies = []
# Here is an example from Neil Berry:
#  FINAL R-B3LYP ENERGY IS     -382.0507446475 AFTER  10 ITERATIONS
# but in some cases the energy can be in position [3] not [4] so let's
# take the number after the "IS"
                temp = line.split()
                self.scfenergies.append(temp[temp.index("IS")+1])

            if not endofopt and line.find("MAXIMUM GRADIENT")>0:
                if not hasattr(self,"geovalues"):
                    self.logger.info("Creating attribute geovalues[]")
                    self.geovalues = []
                temp = line.strip().split()
                self.geovalues.append([float(temp[3]),float(temp[7])])



            if line.find("DENSITY CONV=")==5:
                self.scftargets[0] = float(line.strip().split()[-1])
                
            if line.find("ITER EX DEM")==1:
# This is the section with the SCF information                
                if not hasattr(self,"scfvalues"):
                    self.logger.info("Creating attribute scfvalues")
                    self.scfvalues = []
                line = inputfile.next()
                den = []
                while line!='\n':
# The SCF information is terminated by a blank line                    
                    try:
                        temp = int(line[0:4])
                    except ValueError:
# Occurs for:
#  * * *   INITIATING DIIS PROCEDURE   * * *
#  CONVERGED TO SWOFF, SO DFT CALCULATION IS NOW SWITCHED ON
#  DFT CODE IS SWITCHING BACK TO THE FINER GRID
                        pass
                    else:
                        den.append(float(line.split()[5]))
                    line = inputfile.next()
                self.scfvalues.append(den)

            if line.find("NORMAL COORDINATE ANALYSIS IN THE HARMONIC APPROXIMATION")>=0:
                # Start of the frequency section
                self.logger.info("Creating attributes vibfreqs, vibirs")
                self.vibfreqs = []
                self.vibirs = []

                # Need to get past the list of atomic weights
                hyphens = inputfile.next()
                blank = inputfile.next()
                line = inputfile.next()
                blank = inputfile.next()
                line = inputfile.next()
                numAtom = 0
                while line!="\n":
                    numAtom += 1
                    line = inputfile.next()

                # Print out the following lines which may contain some useful info:
                # e.g. WARNING, MODE 7 HAS BEEN CHOSEN AS A VIBRATION
                line = inputfile.next()
                while line.find("FREQUENCIES IN CM**-1")==-1:
                    line = inputfile.next()
                line = inputfile.next()

                blank = inputfile.next()
                freqNo = inputfile.next()
                while freqNo.find("SAYVETZ")==-1:
                    freq = inputfile.next().strip().split()
                    self.vibfreqs.extend(map(float,freq[1:]))
                    reducedMass = inputfile.next()
                    irIntensity = inputfile.next().strip().split()
                    self.vibirs.extend(map(float,irIntensity[2:]))
                    blank = inputfile.next()
                    # Skip XYZ data for each atom plus
                    # the Sayvetz stuff at the end
                    for j in range(numAtom*3+10):
                        line = inputfile.next()
                    blank = inputfile.next()
                    freqNo = inputfile.next()
                self.vibfreqs = Numeric.array(self.vibfreqs,"f")
                self.vibirs = Numeric.array(self.vibirs,"f")


            if line.find("EIGENVECTORS")==10 or line.find("MOLECULAR OBRITALS")==10:
                # The details returned come from the *final* report of evalues and
                # the last list of symmetries in the log file
                # This is fine for GeoOpt and SP, but may be weird for TD and Freq(?)
                
                # Take the last one of either in the file
                if not hasattr(self,"moenergies"):
                    self.logger.info("Creating attributes moenergies, mosyms")
                self.moenergies = []
                self.mosyms = []
                if not hasattr(self,"nindep"):
                    self.logger.info("Creating attribute nindep with default value")
                    self.nindep = self.nbasis
                self.mocoeffs = Numeric.zeros((1,self.nindep,self.nbasis),"f")
                line = inputfile.next()
                blank = inputfile.next() # blank line
                base = 0
                while line.find("END OF RHF")==-1:
                    line = inputfile.next() # Eigenvector no
                    line = inputfile.next()
                    self.moenergies.extend(map(float,line.split()))
                    line = inputfile.next()
                    self.mosyms.extend(line.split())
                    line = inputfile.next()
                    i=0
                    while line!=blank and line.find("END OF RHF")==-1:
                        if base==0: # Just do this the first time 'round
                            atomno=int(line.split()[2])-1
                            # atomorb[atomno].append(int(line.split()[0])-1)
                            # What's the story with the previous line?
                        temp = line[15:] # Strip off the crud at the start
                        j = 0
                        while j*11+4<len(temp):
                            self.mocoeffs[0,base+j,i] = float(temp[j*11:(j+1)*11])
                            j+=1
                        line = inputfile.next()
                        i+=1
                    base+=5
                self.moenergies = Numeric.array(self.moenergies,"f")

            if line.find("NUMBER OF OCCUPIED ORBITALS")>=0:
                if not hasattr(self,"homos"):
                    self.logger.info("Creating attribute homos")
                    temp = line.strip().split('=')
                    self.homos = Numeric.array(int(temp[-1])-1,"d")

            if line.find("TOTAL NUMBER OF ATOMS")==1:
                self.logger.info("Creating attribute natom")
                self.natom = int(line.split()[-1])
                
            if line.find("NUMBER OF CARTESIAN GAUSSIAN BASIS")==1 or line.find("TOTAL NUMBER OF BASIS FUNCTIONS")==1:
                # The first is from Julien's Example and the second is from Alexander's
                # I think it happens if you use a polar basis function instead of a cartesian one
                self.logger.info("Creating attribute nbasis")
                self.nbasis = int(line.split()[-1])
                    
            elif line.find("TOTAL NUMBER OF MOS IN VARIATION SPACE")==1:
                # Note that this line is not always present, so by default
                # NBsUse is set equal to NBasis (see below).
                self.logger.info("Creating attribute nindep")
                self.indep = int(line.split()[-1])
                
            elif line.find("OVERLAP MATRIX")==0 or line.find("OVERLAP MATRIX")==1:
                # The first is for PC-GAMESS, the second for GAMESS
                # Read 1-electron overlap matrix
                self.logger.info("Creating attribute aooverlaps")
                self.aooverlaps = Numeric.zeros((self.nbasis,self.nbasis), "f")
                base = 0
                while base<self.nbasis:
                    blank = inputfile.next()
                    line = inputfile.next() # Basis fn number
                    blank = inputfile.next()
                    for i in range(self.nbasis-base): # Fewer lines each time
                        line = inputfile.next()
                        temp = line.split()
                        for j in range(4,len(temp)):
                            self.aooverlaps[base+j-4,i+base] = float(temp[j])
                            self.aooverlaps[i+base,base+j-4] = float(temp[j])
                    base+=5

        inputfile.close()

        if not hasattr(self,"geotargets"):
            self.logger.info("Creating attribute geotargets[] with default values")
            opttol = 1e-4
            self.geotargets = Numeric.array([opttol,3./opttol])
        if not hasattr(self,"scftargets"):
            self.logger.info("Creating attribute scftargets[] with default values")
            self.scftargets = Numeric.array([1e-5])
        if hasattr(self,"geovalues"): self.geovalues = Numeric.array(self.geovalues,"f")
        if not hasattr(self,"nindep"):
            self.logger.info("Creating attribute nindep with default value")
            self.nindep = self.nbasis


        
if __name__=="__main__":
    import doctest,parser
    doctest.testmod(parser,verbose=False)
