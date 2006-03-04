"""
GaussSum is a parser for computational chemistry log files.

See http://gausssum.sf.net for more information.

Copyright (C) 2005 Noel O'Boyle <baoilleach@gmail.com>

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
import math,sys,logging,copy,re,os,time
import Numeric
from utils import *
from gnupy import Gnuplot

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

    def tidyevalue(self):
        """Create a namedlist of orbital number, orbital name and evalue."""

        if not hasattr(self,"orbitalname"):
            self.calcorbitalname()

        numbers = range(len(self.evalue[0]),0,-1)
        orbitalname = copy.deepcopy(self.orbitalname)

        for x in orbitalname:
            x.reverse()
        evalue = copy.deepcopy(self.evalue)
        for y in evalue:
            y.reverse()
        
        if len(self.HOMO)==1:
            # Restricted
            colA = Namedlist("Orbital no.",numbers)
            colB = Namedlist("Name",orbitalname[0])
            colC = Namedlist("Eigenvalue (eV)",evalue[0])
            t = colA + colB + colC
        else:
            # Unrestricted
            colA = Namedlist("Orbital no.",numbers)
            colB = Namedlist("Alpha Name",orbitalname[0])
            colC = Namedlist("Alpha Eigenvalue (eV?)",evalue[0])
            colD = Namedlist("Beta Name",orbitalname[1])
            colE = Namedlist("Beta Eigenvalue (eV?)",evalue[1])
            t = colA + colB + colC + colD + colE

        return t

    def calcorbitalname(self):
        """Calculate the name of each orbital in terms of HOMO/LUMO."""

        self.logger.info("Calculating orbitalname [[]]")
        self.orbitalname = []
        for homo in self.HOMO:
            names = []
            for i in range(0,len(self.evalue[0])): # i runs over orbital indices
                if i==homo:
                    name = "HOMO"
                elif i==homo+1:
                    name = "LUMO"
                elif i>homo:
                    name = "L + %d" % (i-homo-1)
                else:
                    name = "H - %d" % (homo-i)
                names.append(name)
            self.orbitalname.append(names)
                
    
    def calcscfprogress(self):
        """Calculate the progress of the self-consistent field convergence.

        Requires:
         scfvalue
         scftarget
        """
        self.logger.info("Calculating scfprogress [[]]")
        answer = []
        for geom in self.scfvalue:
            thisgeom = []
            for scfcycle in range(len(geom[0])):
                progress = 0
                for i in range(len(self.scftarget)):
                    target = self.scftarget[i]
                    value = geom[i][scfcycle]
                    if value>target:
                        progress += math.log(value/target)
                thisgeom.append(progress)
            answer.append(thisgeom)
        self.scfprogress = answer

    def tidyscfprogress(self):
        """Print out the progress of the self-consistent field convergence."""

        df = Namedlist("Step",range(1,len(self.scfprogress[-1])+1))
        df += Namedlist("Progress",self.scfprogress[-1])

        return df

    def plotscfprogress(self,filename):
        """Plot the progress of the SCF optimisation."""
        self.logger.info("Plotting the progress of the SCF convergence")

        g = Gnuplot(filename)        

        g.data( zip(range(1,len(self.scfprogress[-1])+1),self.scfprogress[-1]),"notitle with lines")
        g.data( zip(range(1,len(self.scfprogress[-1])+1),self.scfprogress[-1]),"notitle with points")

        g.commands( "set xrange [%d:%d]" % (1,len(self.scfprogress[-1])),
                    "set yrange [0:]",
                    "set title 'Progress of the SCF convergence'",
                    "set xlabel 'Step number",
                    "set ylabel 'Progress'" )

        g.plot()
        

    def tidygeoprogress(self):
        """Print out the progress of the geometry optimisation."""

        df = Namedlist("Step",range(1,len(self.geoprogress)+1))
        df += Namedlist("Progress",self.geoprogress)

        return df
                   
    def calcgeoprogress(self):
        """Calculate the progress of the geometry optimisation.

        Requires:
         geovalue
         geotarget
        """
        self.logger.info("Calculating geoprogress []")        
        answer = []
        for geom in self.geovalue:
            progress = 0
            for i in range(len(self.geotarget)):
                target = self.geotarget[i]
                value = geom[i]
                if value>target:
                    progress += math.log(value/target)
            answer.append(progress)
        self.geoprogress = answer

    def plotgeoprogress(self,filename):
        """Plot the progress of the geometry optimisation."""
        self.logger.info("Plotting the progress of the geometry optimisation")

        g = Gnuplot(filename)

        g.data( zip(range(1,len(self.geoprogress)+1),self.geoprogress),"notitle with lines")
        g.data( zip(range(1,len(self.geoprogress)+1),self.geoprogress),"notitle with points")

        g.commands(  "set xrange [%d:%d]" % (1,len(self.geoprogress)),
                     "set yrange [0:]",
                     'set title "Progress of the geometry optimisation"',
                     'set xlabel "Step number"',
                     'set ylabel "Progress"' )
                    
        g.plot()


    def tidyirspectrum(self):
        """Print out the IR spectrum."""

        df = Namedlist("Energy (cm-1)",self.irspectrum.xvalues)
        df += Namedlist("Intensity",self.irspectrum.spectrum[:,0])

        return df

    def plotirspectrum(self,filename):
        """Plot the IR spectrum."""

        self.logger.info("Plotting the IR spectrum")
        
        start = self.irspectrum.start
        end = self.irspectrum.end

        g = Gnuplot(filename)

        g.data( zip(self.irspectrum.xvalues,self.irspectrum.spectrum[:,0]),"notitle with lines")

        g.commands( "set xlabel 'Frequency (cm-1)'",
                    "set ylabel 'IR activity'",
                    "set xrange [%f:%f] reverse" % (start,end),
                    "set yrange [*:*] reverse",
                    "set key bottom",
                    "set title 'IR spectrum'" )
        
        g.plot()

    def tidyramanspectrum(self):
        """Print out the Raman spectrum."""

        df = Namedlist("Energy (cm-1)",self.ramanspectrum.xvalues)
        df += Namedlist("Intensity",self.ramanspectrum.spectrum[:,0])

        return df

    def plotramanspectrum(self,filename):
        """Plot the Raman spectrum."""

        self.logger.info("Plotting the Raman spectrum")
        
        start = self.ramanspectrum.start
        end = self.ramanspectrum.end

        g = Gnuplot(filename)

        g.data( zip(self.ramanspectrum.xvalues,self.ramanspectrum.spectrum[:,0]),"notitle with lines")

        g.commands( "set xlabel 'Frequency (cm-1)'",
                    "set ylabel 'Raman activity'",
                    "set title 'Raman spectrum'" )
        
        g.plot()

    def tidyuvspectrum(self):
        """Print out the UV spectrum."""

        df = Namedlist("Energy (cm-1)",self.uvspectrum.xvalues)
        df += Namedlist("Wavelength (nm)",
                        [convertor(x,"cm-1","nm") for x in self.uvspectrum.xvalues])
        df += Namedlist("Oscillator strength", self.uvspectrum.spectrum[:,0])

        return df    


    def calcuvspectrum(self):
        """Convolute the UV spectrum.

        Requires:
         etenergy
         etosc
        """
        self.logger.info("Calculating uvspectrum 'GaussianSpectrum'")        
        # UV spectrum is Gaussian in *energy* not *wavelength*
        start = convertor(self.pref['uvvis.start'],"nm","cm-1")
        end = convertor(self.pref['uvvis.end'],"nm","cm-1")
        
        t = GaussianSpectrum(start,end,self.pref['uvvis.numpoints'],
                             ( self.etenergy,[[x*2.174e8/self.pref['uvvis.fwhm'] for x in self.etosc]] ),
                             self.pref['uvvis.fwhm'])
        self.uvspectrum = t

    def plotuvspectrum(self,filename):
        """Plot the UV spectrum."""

        self.logger.info("Plotting the UV spectrum")
        wavelen = [convertor(x,"cm-1","nm") for x in self.uvspectrum.xvalues]
        start = convertor(self.uvspectrum.start,"cm-1","nm")
        end = convertor(self.uvspectrum.end,"cm-1","nm")

        g = Gnuplot(filename)

        g.data( zip(wavelen,self.uvspectrum.spectrum[:,0]),"notitle with lines")
        g.data( zip(self.etwavelen,self.etosc),"notitle axes x1y2 with impulses")

        g.commands(  "set xrange [%f:%f]" % (start,end),
                     "set ytics nomirror",
                     "set y2tics",
                     'set title "UV spectrum"',
                     'set ylabel "Epsilon"',
                     'set y2label "Oscillator strength"' )

        g.plot()
    

    def calcirspectrum(self):
        """Convolute the IR spectrum.

        Requires:
         vibfreq
         ir
        """
        self.logger.info("Calculating irspectrum 'Spectrum'")                
        # IR spectrum is Lorentzian in energy
        start = self.pref['ir.start']
        end = self.pref['ir.end']
        
        t = Spectrum(start,end,self.pref['ir.numpoints'],
                     [zip(self.vibfreq,self.ir)],
                     self.pref['ir.fwhm'], lorentzian)
        self.irspectrum = t

    def calcramanspectrum(self):
        """Convolute the Raman spectrum.

        Requires:
         vibfreq
         raman
        """
        self.logger.info("Calculating ramanspectrum 'Spectrum'")                
        # Raman spectrum is Lorentzian in energy
        start = self.pref['raman.start']
        end = self.pref['raman.end']
        
        t = Spectrum(start,end,self.pref['raman.numpoints'],
                     [zip(self.vibfreq,self.raman)],
                     self.pref['raman.fwhm'], lorentzian)
        self.ramanspectrum = t


    def calcpdos(self):
        """Calculate the contribution of various groups to each orbital.

        Requires:
         mocoeff
         overlap

        Optional:
         groups (defaults to "allorbitals")
        """
        self.logger.info("Calculating pdos")
        # contrib is a matrix of the contributions of each basis fn to each MO
        # dot() is matrix multiplication
        #  * is tensor multiplication(?): corresponding entries are multiplied
        contrib = self.mocoeff * Numeric.dot(self.mocoeff,self.overlap)

        if not hasattr(self,"groups"):
            self.groups = Groups(self.orbitals,type="allorbitals")

        groups = self.groups.groups
        groupnames = groups.keys()
        data = Numeric.zeros( (self.NBsUse,len(groups)), "float")

        for i,groupname in enumerate(groupnames):
            for basisfn in groups[groupname]:
                data[:,i] += contrib[:,basisfn]

        self.pdos = data

    def tidypdos(self):
        """Create a nice dataframe for the pdos."""

        if not hasattr(self,"orbitalname"):
            self.calcorbitalname()
        orbitalname = copy.deepcopy(self.orbitalname)
        for x in orbitalname:
            x.reverse()
        evalue = copy.deepcopy(self.evalue)
        for y in evalue:
            y.reverse()
            
        groupnames = self.groups.groups.keys()
        numbers = range(len(self.evalue[0]),0,-1)            
        colA = Namedlist("Orbital number",numbers)

        if len(self.HOMO)==1:
            # Restricted
            colB = Namedlist("Name",orbitalname[0])
            colC = Namedlist("Eigenvalue (eV)",evalue[0])
            t = colA + colB + colC
            for i in range(len(groupnames)):
                groupdata = list(self.pdos[:,i])
                groupdata.reverse()
                t += Namedlist(groupnames[i],groupdata)
        else:
            # Unrestricted
            pass

        return t

    def tidydosspectrum(self):
        """Print the dos spectrum nicely."""

        df = Namedlist("Energy (eV)",self.dosspectrum[0].xvalues)
        for i,spectrum in enumerate(self.dosspectrum):
            name = "DOS"
            if len(self.HOMO)==2:
                name = ["Alpha DOS","Beta DOS"] [i]
            df += Namedlist(name, spectrum.spectrum[:,i])

        return df

    def plotdosspectrum(self,filename):
        """Plot the DOS spectrum."""

        self.logger.info("Plotting the DOS spectrum")
        
        start = self.dosspectrum[0].start
        end = self.dosspectrum[0].end

        g = Gnuplot(filename)

        g.data( zip(self.dosspectrum[0].xvalues,self.dosspectrum[0].spectrum[:,0]),'title "DOS spectrum"')
        # HOMO holds the index of the HOMO in evalue
        realorbs = [ (x,-1) for x in self.evalue[0][:self.HOMO[0]+1] ]
        virtorbs = [ (x,-1) for x in self.evalue[0][self.HOMO[0]+1:] ]
        g.data(realorbs,'title "Occupied orbitals" with impulses')
        g.data(virtorbs,'title "Virtual orbitals" with impulses')
        
        g.commands(  "set data style lines",
                     "set xlabel 'Energy (eV)'",
                     "set xrange [%f:%f]" % (start,end),
                     "set yrange [-1:*]",
                     'set title "Density of states spectrum"' )
        
        g.plot()
   

    def calcdosspectrum(self):
        """Convolute the density of states spectrum."""
        self.logger.info("Calculating dosspectrum")
        self.dosspectrum = []

        for i in range(len(self.HOMO)):
             t = GaussianSpectrum(self.pref['dos.start'],self.pref['dos.end'],self.pref['dos.numpoints'],
                                 ( self.evalue[i], [[1]*len(self.evalue[i])]),
                                 self.pref['dos.fwhm'])
             self.dosspectrum.append(t)

    def tidypdosspectrum(self):
        """Make a nice print out of the pdosspectrum."""

        groups = self.groups.groups
        groupnames = groups.keys()

        df = Namedlist("Energy (eV)",self.pdosspectrum[0].xvalues)
        for i,spectrum in enumerate(self.pdosspectrum):
            prefix = ""
            if len(self.HOMO)==2:
                prefix = ["Alpha","Beta"] [i]
            df += Namedlist([prefix+x for x in groupnames], spectrum.spectrum)

        return df

    def plotpdosspectrum(self,filename,stacked=False):
        """Plot the PDOS spectrum."""

        self.logger.info("Plotting the PDOS spectrum")
        
        start = self.pdosspectrum[0].start
        end = self.pdosspectrum[0].end

        g = Gnuplot(filename)

        tot = Numeric.zeros( len(self.pdosspectrum[0].xvalues),"double" )
        for i,v in enumerate(self.groups.groups.keys()):
            g.data( zip(self.pdosspectrum[0].xvalues,tot+self.pdosspectrum[0].spectrum[:,i]),                    'title "%s"' % v)
            if stacked:
                tot += self.pdosspectrum[0].spectrum[:,i]
        
        # HOMO holds the index of the HOMO in evalue
        realorbs = [ (x,-1) for x in self.evalue[0][:self.HOMO[0]+1] ]
        virtorbs = [ (x,-1) for x in self.evalue[0][self.HOMO[0]+1:] ]
        g.data(realorbs,'title "Occupied orbitals" with impulses')
        g.data(virtorbs,'title "Virtual orbitals" with impulses')
        
        g.commands( "set data style lines",
                    "set xlabel 'Energy (eV)'",
                    "set xrange [%f:%f]" % (start,end),
                    "set yrange [-1:*]",
                    'set title "Partial density of states spectrum"' )
        g.plot()

        
    def calcpdosspectrum(self):
        """Convolute the partial density of states spectrum."""
        self.logger.info("Calculating pdosspectrum")
        
        self.pdosspectrum = []
        groups = self.groups.groups
        groupnames = groups.keys()

        if len(self.HOMO)==1:
            # Restricted
            peaks = ( self.evalue[0],[self.pdos[:,x] for x in range(len(groups))] )
            t = GaussianSpectrum(self.pref['dos.start'],self.pref['dos.end'],self.pref['dos.numpoints'],
                         peaks,self.pref['dos.fwhm'])
            self.pdosspectrum.append(t)
        else:
            # Unrestricted
            peaks = ( self.evalue[0],[self.pdos[0][:,x] for x in range(len(groups))] )
            t = GaussianSpectrum(self.pref['dos.start'],self.pref['dos.end'],self.pref['dos.numpoints'],
                         peaks,self.pref['dos.fwhm'])
            self.pdosspectrum.append(t)
            peaks = ( self.evalue[1],[self.pdos[1][:,x] for x in range(len(groups))] )
            t = GaussianSpectrum(self.pref['dos.start'],self.pref['dos.end'],self.pref['dos.numpoints'],
                         peaks,self.pref['dos.fwhm'])
            self.pdosspectrum.append(t)
            

        
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

        self.filename = filename

        # Set up the preferences (with the default values)
        self.pref = Preferences()
        
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
                self.scftarget[G03.SCFRMS] = G03.float(line.split('=')[1].split()[0])
            if line[1:44]=='Requested convergence on MAX density matrix':
                self.scftarget[G03.SCFMAX] = G03.float(line.strip().split('=')[1][:-1])
            if line[1:44]=='Requested convergence on             energy':
                self.scftarget[G03.SCFENERGY] = G03.float(line.strip().split('=')[1][:-1])

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
                        newlist[G03.SCFRMS].append(G03.float(parts[0].split('=')[1]))
                        newlist[G03.SCFMAX].append(G03.float(parts[1].split('=')[1]))
                        energy = 1.0
                        if len(parts)>4:
                            energy = parts[2].split('=')[1]
                            if energy=="":
                                energy = G03.float(parts[3])
                            else:
                                energy = G03.float(energy)
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
                    self.scfvalue[0].append(G03.float(parts[-1][:-1]))
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
                        value = G03.float(parts[2])
                    except ValueError:
                        self.logger.error("Problem parsing the value for geometry optimisation: %s is not a number." % parts[2])
                    else:
                        newlist[i] = value
                    self.geotarget[i] = G03.float(parts[3])
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
                        self.evalue[0].append(G03.float(x)*27.2114) # from a.u. (hartrees) to eV
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
                        self.evalue[1].append(G03.float(x)*27.2114) # from a.u. (hartrees) to eV
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
                    self.vibfreq.extend(map(G03.float,line[15:].split())) # Adding new frequencies
                    [inputfile.next() for i in [0,1]] # Skip two lines
                    line = inputfile.next()
                    self.ir.extend(map(G03.float,line[15:].split())) # Adding IR intensities
                    line = inputfile.next()
                    if line.find("Raman")>=0:
                        if not hasattr(self,"raman"):
                            self.raman = []
                            self.logger.info("Creating attribute raman[]")
                        line = inputfile.next()
                        self.raman.extend(map(G03.float,line[15:].split())) # Adding Raman intensities
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
                self.etenergy.append(convertor(G03.float(parts[0]),"eV","cm-1"))
                self.etwavelen.append(G03.float(parts[2]))
                self.etosc.append(G03.float(parts[4].split("=")[1]))
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

                    percent = G03.float(t[1])
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
                        R = G03.float(parts[-1])
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

        
#    @staticmethod
    def float(number):
        """Convert a string to a float."""
        try:
            ans = float(number)
        except ValueError:
            number = "E".join(number.split("D"))
            ans = float(number)
        return ans
    float = staticmethod(float)
    
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
                    coords.append([ G03.float(x)*0.5292 for x in [parts[3],parts[5],parts[7]] ])
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
        
    def writeTrajectory(self,outputfilename):
        """Write trajectory information to a quasi-PDB file."""
        outputfile = open(outputfilename,"w")
        pt = PeriodicTable()
        for i in range(len(self.traj)):
            for j in range(self.NAtoms):
                line = "HETATM %4d %2s   UNK            % .3f  % .3f  % .3f  0.00  0.00          %2s\n" % (j+1,pt.element[self.atomicNo[j]],self.traj[i][j][0],self.traj[i][j][1],self.traj[i][j][2],pt.element[self.atomicNo[j]])
                outputfile.write(line)
            if hasattr(self,"connect"):
                for k,v in self.connect[i].iteritems():
                    line = "CONECT %4d" % (k+1)
                    for value in v:
                        line += " %4d" % (value+1)
                    outputfile.write(line+"\n")
                    
            outputfile.write("END\n")
        outputfile.close()

    def writeTrajectorySummary(self,outputfilename):
        """Write trajectory summary information to TSF."""
        outputfile = open(outputfilename,"w")
        header = "\t".join(["Time (fs)","Kinetic (au)","Potent (au)","Delta E (au)","Delta A (h-bar)"])
        outputfile.write(header+"\n")
        for x in self.trajSummary:
            outputfile.write("\t".join(x)+"\n")
        outputfile.close()

if __name__=="__main__":
    import doctest,parser
    doctest.testmod(parser,verbose=False)
