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

class Jaguar(Logfile):
    """A Jaguar output file"""

    def __init__(self,*args):

        # Call the __init__ method of the superclass
        super(Jaguar, self).__init__(logname="Jaguar",*args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Jaguar output file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Jaguar("%s")' % (self.filename)

    def parse(self,fupdate=0.05,cupdate=0.002):
        """Extract information from the logfile."""
        inputfile = open(self.filename,"r")
        
        if self.progress:
            
            inputfile.seek(0,2) #go to end of file
            nstep=inputfile.tell()
            inputfile.seek(0)
            self.progress.initialize(nstep)
            oldstep=0
            
        for line in inputfile:
            
            if self.progress and random.random()<cupdate:
                
                step = inputfile.tell()
                if step!=oldstep:
                    self.progress.update(step,"Unsupported Information")
                    oldstep = step

            if line[0:4]=="etot":
# Get SCF convergence information
                if not hasattr(self,"scfvalues"):
                    self.scfvalues = []
                    self.logger.info("Creating attribute: scfvalues")
                values = []
                while line[0:4]=="etot":
                    if line[39:47].strip():
                        denergy = float(line[39:47])
                    else:
                        denergy = 0 # Should really be greater than target value
                                    # or should we just ignore the values in this line
                    ddensity = float(line[48:56])
                    maxdiiserr = float(line[57:65])
                    values.append([denergy,ddensity,maxdiiserr])
                    line = inputfile.next()
                self.scfvalues.append(values)

            if line[1:5]=="SCFE":
# Get the energy of the molecule
                if not hasattr(self,"scfenergies"):
                    self.logger.info("Creating attribute scfenergies")
                    self.scfenergies = []
                temp = line.strip().split()
                self.scfenergies.append(float(temp[temp.index("hartrees")-1]))

            if line[2:28]=="geometry optimization step":
# Get Geometry Opt convergence information
                if not hasattr(self,"geovalues"):
                    self.geovalues = []
                    self.geotargets = Numeric.zeros(4,"float")
                    self.logger.info("Creating attributes: geovalues,geotargets")
                blank = inputfile.next()
                blank = inputfile.next()
                line = inputfile.next()
                i = 0
                values = []
                while line!=blank:
                    if line[41]=="(":
                        # A new geo convergence value
                        values.append(float(line[26:37]))
                        self.geotargets[i] = float(line[43:54])
                        i+=1
                    line = inputfile.next()
                self.geovalues.append(values)

            if line[2:33]=="Orbital energies/symmetry label":
# Get MO Energies and symmetrys
                if not hasattr(self,"moenergies"):
                    self.logger.info("Creating attributes: moenergies, mosyms")
                self.mosyms = [[]]
                self.moenergies = [[]]
                line = inputfile.next()
                while line.strip():
                    temp = line.strip().split()
                    for i in range(0,len(temp),2):
                        self.moenergies[0].append(convertor(float(temp[i]),"hartree","eV"))
                        self.mosyms[0].append(temp[i+1])
                    line = inputfile.next()
                self.moenergies = Numeric.array(self.moenergies,"f")

            if line[1:28]=="number of occupied orbitals":
                if not hasattr(self,"homos"):
                    self.logger.info("Creating attribute: homos")
                self.homos = Numeric.array([float(line.strip().split()[-1])-1],"i")

            if line[2:27]=="number of basis functions":
                if not hasattr(self,"nbasis"):
                    self.logger.info("Creating attribute: nbasis")
                self.nbasis = float(line.strip().split()[-1])

            if line[2:23]=="start of program freq":
# IR stuff
                self.logger.info("Creating attribute: vibfreqs")
                self.vibfreqs = []
                blank = inputfile.next()
                line = inputfile.next(); line = inputfile.next()
                blank = inputfile.next()
                
                freqs = inputfile.next()
                while freqs!=blank:
                    temp = freqs.strip().split()
                    self.vibfreqs.extend(map(float,temp[1:]))
                    temp = inputfile.next().strip().split()
                    if temp[0]=="symmetries": # May go straight from frequencies to reduced mass
                        if not hasattr(self,"vibsyms"):
                            self.logger.info("Creating attributes: vibsyms, vibirs")
                            self.vibsyms = []
                            self.vibirs = []
                        self.vibsyms.extend(map(self.normalisesym,temp[1:]))
                        temp = inputfile.next().strip().split()                                
                        self.vibirs.extend(map(float,temp[1:]))
                        reducedmass = inputfile.next()
                    line = inputfile.next()
                    while line!=blank: # Read the cartesian displacements
                        line = inputfile.next()
                    freqs = inputfile.next()
                self.vibfreqs = Numeric.array(self.vibfreqs)
                if hasattr(self,"vibirs"):
                    self.vibirs = Numeric.array(self.vibirs)

        inputfile.close()

##        if hasattr(self,"scfvalues"):
##            self.scfvalues = Numeric.array(self.scfvalues,"f")
        if hasattr(self,"scfenergies"):
            self.scfenergies = Numeric.array(self.scfenergies,"f")
        self.parsed = True
        
if __name__=="__main__":
    import doctest,jaguarparser
    doctest.testmod(jaguarparser,verbose=False)
