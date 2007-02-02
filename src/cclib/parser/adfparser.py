"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

import Numeric
import random # For sometimes running the progress updater
import utils
import logfileparser

class ADF(logfileparser.Logfile):
    """An ADF log file"""
    SCFCNV, SCFCNV2 = range(2) #used to index self.scftargets[]
    maxelem, norm = range(2) # used to index scf.values
    def __init__(self, *args):

        # Call the __init__ method of the superclass
        super(ADF, self).__init__(logname="ADF", *args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "ADF log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'ADF("%s")' % (self.filename)

    def normalisesym(self, label):
        """Use standard symmetry labels instead of ADF labels.

        To normalise:
        (1) any periods are removed (except in the case of greek letters)
        (2) XXX is replaced by X, and a " added.
        (3) XX is replaced by X, and a ' added.
        (4) The greek letters Sigma, Pi, Delta and Phi are replaced by
            their lowercase equivalent.

        >>> sym = ADF("dummyfile").normalisesym
        >>> labels = ['A','s','A1','A1.g','Sigma','Pi','Delta','Phi','Sigma.g','A.g','AA','AAA','EE1','EEE1']
        >>> map(sym,labels)
        ['A', 's', 'A1', 'A1g', 'sigma', 'pi', 'delta', 'phi', 'sigma.g', 'Ag', "A'", 'A"', "E1'", 'E1"']
        """
        greeks = ['Sigma', 'Pi', 'Delta', 'Phi']
        for greek in greeks:
            if label.startswith(greek):
                return label.lower()
            
        ans = label.replace(".", "")
        l = len(ans)
        if l > 1 and ans[0] == ans[1]: # Python only tests the second condition if the first is true
            if l > 2 and ans[1] == ans[2]:
                ans = ans.replace(ans[0]*3, ans[0]) + '"'
            else:
                ans = ans.replace(ans[0]*2, ans[0]) + "'"
        return ans
        
    def normalisedegenerates(self, label, num):
        """Generate a string used for matching degenerate orbital labels

        To normalise:
        (1) if label is E or T, return label:num
        (2) if label is P or D, look up in dict, and return answer
        """

        ndict = { 'P': {0:"P:x", 1:"P:y", 2:"P:z"},\
                  'D': {0:"D:z2", 1:"D:x2-y2", 2:"D:xy", 3:"D:xz", 4:"D:yz"}}
                  
        if label == 'P' or label == 'D':
            return ndict[label][num]

        else:
            return "%s:%i"%(label,num+1)

    def extract(self, inputfile, fupdate=0.05, cupdate=0.002):
        """Extract information from the file object inputfile."""
        
        oldstep = 0
            
        # Used to avoid extracting the final geometry twice in a GeoOpt
        NOTFOUND, GETLAST, NOMORE = range(3)
        finalgeometry = NOTFOUND 
        
        # Used for calculating the scftarget (variables names taken from the ADF manual)
        accint = SCFconv = sconv2 = None
        
        # keep track of nosym and unrestricted case to parse Energies since it doens't have an all Irreps section
        nosymflag = False
        unrestrictedflag = False
        
        for line in inputfile:
            
            if self.progress and random.random() < cupdate:
                step = inputfile.tell()
                if step != oldstep:
                    self.progress.update(step, "Unsupported Information")
                    oldstep = step
                
            if line.find("INPUT FILE") >= 0:
#check to make sure we aren't parsing Create jobs
                while line:
                
                    if self.progress and random.random() < fupdate:
                        step = inputfile.tell()
                        #if step!=oldstep:
                        self.progress.update(step, "Unsupported Information")
                        oldstep = step
                  
                    if line.find("INPUT FILE") >= 0:
                        line2 = inputfile.next()
                    else:
                        line2 = None

                    if line2 and line2.find("Create") < 0:
                        break
                            
                    line = inputfile.next()
            
            if line[1:10] == "Symmetry:":
                info = line.split()
                if info[1] == "NOSYM":
                    nosymflag = True

            if line[4:13] == 'Molecule:':
                info = line.split()
                if info[1] == 'UNrestricted':
                    unrestrictedflag = True

            if line[1:6] == "ATOMS":
# Find the number of atoms and their atomic numbers
# Also extract the starting coordinates (for a GeoOpt anyway)
                if self.progress and random.random() < cupdate:
                    step = inputfile.tell()
                    if step != oldstep:
                        self.progress.update(step, "Attributes")
                        oldstep = step
                
                self.logger.info("Creating attribute atomnos[], atomcoords[], coreelectrons[]")
                self.atomnos = []
                self.atomcoords = []
                self.coreelectrons = []
                
                underline = inputfile.next()  #clear pointless lines
                label1 = inputfile.next()     # 
                label2 = inputfile.next()     #
                line = inputfile.next()
                atomcoords = []
                while len(line)>2: #ensure that we are reading no blank lines
                    info = line.split()
                    element = info[1].split('.')[0]
                    self.atomnos.append(self.table.number[element])
                    atomcoords.append(map(float, info[2:5]))
                    self.coreelectrons.append(int(float(info[5]) - float(info[6])))
                    line = inputfile.next()
                self.atomcoords.append(atomcoords)
                
                self.natom = len(self.atomnos)
                self.logger.info("Creating attribute natom: %d" % self.natom)
                self.atomnos = Numeric.array(self.atomnos, "i")
                
            if line[1:10] == "FRAGMENTS":
                header = inputfile.next()

                self.frags = []
                self.fragnames = []

                line = inputfile.next()
                while len(line) > 2: #ensure that we are reading no blank lines
                    info = line.split()
                    
                    if len(info) == 7: #fragment name is listed here
                        self.fragnames.append("%s_%s"%(info[1],info[0]))
                        self.frags.append([])
                        self.frags[-1].append(int(info[2]) - 1)

                    elif len(info) == 5: #add atoms into last fragment
                        self.frags[-1].append(int(info[0]) - 1)

                    line = inputfile.next()

            if line[1:22] == "S C F   U P D A T E S":
# find targets for SCF convergence

                if not hasattr(self,"scftargets"):
                    self.logger.info("Creating attribute scftargets[]")
                    self.scftargets = []

                #underline, blank, nr
                for i in range(3):
                    inputfile.next()

                line = inputfile.next()
                SCFconv = float(line.split()[-1])
                line = inputfile.next()
                sconv2 = float(line.split()[-1])
              
            if line[1:11] == "CYCLE    1":
              
                if self.progress and random.random() < fupdate:
                    step = inputfile.tell()
                    if step != oldstep:
                        self.progress.update(step, "QM Convergence")
                        oldstep = step
              
                newlist = []
                line = inputfile.next()

                if not hasattr(self,"geovalues"):
                    # This is the first SCF cycle
                    self.scftargets.append([sconv2*10, sconv2])
                elif finalgeometry in [GETLAST, NOMORE]:
                    # This is the final SCF cycle
                    self.scftargets.append([SCFconv*10, SCFconv])
                else:
                    # This is an intermediate SCF cycle
                    oldscftst = self.scftargets[-1][1]
                    grdmax = self.geovalues[-1][1]
                    scftst = max(SCFconv, min(oldscftst, grdmax/30, 10**(-accint)))
                    self.scftargets.append([scftst*10, scftst])
                        
                while line.find("SCF CONVERGED") == -1 and line.find("SCF not fully converged, result acceptable") == -1 and line.find("SCF NOT CONVERGED") == -1:
                    if line[4:12] == "SCF test":
                        if not hasattr(self, "scfvalues"):
                            self.logger.info("Creating attribute scfvalues")
                            self.scfvalues = []
                                                
                        info = line.split()
                        newlist.append([float(info[4]), abs(float(info[6]))])
                    try:
                        line = inputfile.next()
                    except StopIteration: #EOF reached?
                        self.logger.warning("SCF did not converge, so attributes may be missing")
                        break            

                if line.find("SCF not fully converged, result acceptable") > 0:
                    self.logger.warning("SCF not fully converged, results acceptable")

                if line.find("SCF NOT CONVERGED") > 0:
                    self.logger.warning("SCF did not converge! moenergies and mocoeffs are unreliable")

                if hasattr(self, "scfvalues"):
                    self.scfvalues.append(newlist)
              
            if line[51:65] == "Final Geometry":
                finalgeometry = GETLAST
            
            if line[1:24] == "Coordinates (Cartesian)" and finalgeometry in [NOTFOUND, GETLAST]:
                # Get the coordinates from each step of the GeoOpt
                if not hasattr(self, "atomcoords"):
                    self.logger.info("Creating attribute atomcoords")
                    self.atomcoords = []
                equals = inputfile.next()
                blank = inputfile.next()
                title = inputfile.next()
                title = inputfile.next()
                hyphens = inputfile.next()

                atomcoords = []
                line = inputfile.next()
                while line != hyphens:
                    atomcoords.append(map(float, line.split()[5:8]))
                    line = inputfile.next()
                self.atomcoords.append(atomcoords)
                if finalgeometry == GETLAST: # Don't get any more coordinates
                    finalgeometry = NOMORE

            if line[1:27] == 'Geometry Convergence Tests':
# Extract Geometry convergence information
                if not hasattr(self, "geotargets"):
                    self.logger.info("Creating attributes geotargets[], geovalues[[]]")
                    self.geovalues = []
                    self.geotargets = Numeric.array([0.0, 0.0, 0.0, 0.0, 0.0], "f")
                if not hasattr(self, "scfenergies"):
                    self.logger.info("Creating attribute scfenergies[]")
                    self.scfenergies = []
                equals = inputfile.next()
                blank = inputfile.next()
                line = inputfile.next()
                temp = inputfile.next().strip().split()
                self.scfenergies.append(utils.convertor(float(temp[-1]), "hartree", "eV"))
                for i in range(6):
                    line = inputfile.next()
                values = []
                for i in range(5):
                    temp = inputfile.next().split()
                    self.geotargets[i] = float(temp[-3])
                    values.append(float(temp[-4]))
                self.geovalues.append(values)
 
            if line[1:27] == 'General Accuracy Parameter':
                # Need to know the accuracy of the integration grid to
                # calculate the scftarget...note that it changes with time
                accint = float(line.split()[-1])
            
            if line.find('Orbital Energies, per Irrep and Spin') > 0 and not hasattr(self, "mosyms") and nosymflag and not unrestrictedflag:
#Extracting orbital symmetries and energies, homos for nosym case
#Should only be for restricted case because there is a better text block for unrestricted and nosym
                
                self.logger.info("Creating attribute mosyms[[]]")
                self.mosyms = [[]]

                self.logger.info("Creating attribute moenergies[[]]")
                self.moenergies = [[]]
                 
                underline = inputfile.next()
                header = inputfile.next()
                underline = inputfile.next()
                label = inputfile.next()
                line = inputfile.next()

                info = line.split()

                if not info[0] == '1':
                    self.logger.warning("MO info up to #%s is missing" % info[0])
 
                #handle case where MO information up to a certain orbital are missing
                while int(info[0]) != len(self.moenergies[0]):
                    self.moenergies[0].append(99999)
                    self.mosyms[0].append('A')
          
                homoA = None

                while len(line) > 10:
                    info = line.split()
                    self.mosyms[0].append('A')
                    self.moenergies[0].append(utils.convertor(float(info[2]), 'hartree', 'eV'))
                    if info[1] == '0.000' and not hasattr(self, 'homos'):
                        self.logger.info("Creating attribute homos[]")
                        self.homos = [len(self.moenergies[0]) - 2]
                    line = inputfile.next()

                self.moenergies = [Numeric.array(self.moenergies[0], "f")]
                self.homos = Numeric.array(self.homos, "i")

            if line[1:29] == 'Orbital Energies, both Spins' and not hasattr(self, "mosyms") and nosymflag and unrestrictedflag:
#Extracting orbital symmetries and energies, homos for nosym case
#should only be here if unrestricted and nosym

                self.logger.info("Creating attribute mosymms[[]]")
                self.mosyms = [[], []]

                self.logger.info("Creating attribute moenergies[[]]")
                moenergies = [[], []]

                underline = inputfile.next()
                blank = inputfile.next()
                header = inputfile.next()
                underline = inputfile.next()
                line = inputfile.next()

                homoa = 0
                homob = None

                while len(line) > 5:
                    info = line.split()
                    if info[2] == 'A': 
                        self.mosyms[0].append('A')
                        moenergies[0].append(utils.convertor(float(info[4]), 'hartree', 'eV'))
                        if info[3] != '0.00':
                            homoa = len(moenergies[0]) - 1
                    elif info[2] == 'B':
                        self.mosyms[1].append('A')
                        moenergies[1].append(utils.convertor(float(info[4]), 'hartree', 'eV'))
                        if info[3] != '0.00':
                            homob = len(moenergies[1]) - 1
                    else:
                        print "Error reading line: %s" % line

                    line = inputfile.next()

                self.moenergies = [Numeric.array(x, "f") for x in moenergies]
                self.logger.info("Creating attribute homos[]")
                self.homos = Numeric.array([homoa, homob], "i")


            if line[1:29] == 'Orbital Energies, all Irreps' and not hasattr(self, "mosyms"):
#Extracting orbital symmetries and energies, homos
                self.logger.info("Creating attribute mosyms[[]]")
                self.mosyms = [[]]
                symlist = {}

                self.logger.info("Creating attribute moenergies[[]]")
                self.moenergies = [[]]
                
                underline = inputfile.next()
                blank = inputfile.next()
                header = inputfile.next()
                underline2 = inputfile.next()
                line = inputfile.next()
                
                homoa = None
                homob = None

                multiple = {'E':2, 'T':3, 'P':3, 'D':5}
  
                while line.strip():
                    info = line.split()
                    if len(info) == 5: #this is restricted
                        count = multiple.get(info[0][0],1)
                        for repeat in range(count): # i.e. add E's twice, T's thrice
                            self.mosyms[0].append(self.normalisesym(info[0]))
                            self.moenergies[0].append(utils.convertor(float(info[3]), 'hartree', 'eV'))

                            sym = info[0]
                            if count > 1: # add additional sym label
                                sym = self.normalisedegenerates(info[0],repeat)

                            try:
                                symlist[sym][0].append(len(self.moenergies[0])-1)
                            except KeyError:
                                symlist[sym]=[[]]
                                symlist[sym][0].append(len(self.moenergies[0])-1)

                        if info[2] == '0.00' and not hasattr(self, 'homos'):
                            self.logger.info("Creating attribute homos[]")
                            self.homos = [len(self.moenergies[0]) - (count + 1)] #count, because need to handle degenerate cases
                        line = inputfile.next()
                    elif len(info) == 6: #this is unrestricted
                        if len(self.moenergies) < 2: #if we don't have space, create it
                            self.moenergies.append([])
                            self.mosyms.append([])
                        count = multiple.get(info[0][0], 1)
                        if info[2] == 'A':
                            for repeat in range(count): # i.e. add E's twice, T's thrice
                                self.mosyms[0].append(self.normalisesym(info[0]))
                                self.moenergies[0].append(utils.convertor(float(info[4]), 'hartree', 'eV'))
                                
                                sym = info[0]
                                if count > 1: #add additional sym label
                                    sym = self.normalisedegenerates(info[0],repeat)
                                
                                try:
                                    symlist[sym][0].append(len(self.moenergies[0])-1)
                                except KeyError:
                                    symlist[sym]=[[],[]]
                                    symlist[sym][0].append(len(self.moenergies[0])-1)

                            if info[3] == '0.00' and homoa == None:
                                homoa = len(self.moenergies[0]) - (count + 1) #count because degenerate cases need to be handled
                                
                        if info[2] == 'B':
                            for repeat in range(count): # i.e. add E's twice, T's thrice
                                self.mosyms[1].append(self.normalisesym(info[0]))
                                self.moenergies[1].append(utils.convertor(float(info[4]), 'hartree', 'eV'))

                                sym = info[0]
                                if count > 1: #add additional sym label
                                    sym = self.normalisedegenerates(info[0],repeat)
                                
                                try:
                                    symlist[sym][1].append(len(self.moenergies[1])-1)
                                except KeyError:
                                    symlist[sym]=[[],[]]
                                    symlist[sym][1].append(len(self.moenergies[1])-1)

                            if info[3] == '0.00' and homob == None:
                                homob = len(self.moenergies[1]) - (count + 1)
                                
                        line = inputfile.next()
                        
                    else: #different number of lines
                        print "Error", info
      
                if len(info) == 6: #still unrestricted, despite being out of loop
                    self.logger.info("Creating attribute homos[]")
                    self.homos = [homoa, homob]
    
                self.moenergies = [Numeric.array(x, "f") for x in self.moenergies]
                self.homos = Numeric.array(self.homos, "i")

            if line[1:28] == "Vibrations and Normal Modes":
                # Section on extracting vibdisps
                # Also contains vibfreqs, but these are extracted in the
                # following section (see below)
                self.logger.info("Creating attribute vibdisps")
                self.vibdisps = []
                equals = inputfile.next()
                blank = inputfile.next()
                header = inputfile.next()
                header = inputfile.next()
                blank = inputfile.next()
                blank = inputfile.next()

                freqs = inputfile.next()
                while freqs.strip()!="":
                    minus = inputfile.next()
                    p = [ [], [], [] ]
                    for i in range(len(self.atomnos)):
                        broken = map(float, inputfile.next().split()[1:])
                        for j in range(0, len(broken), 3):
                            p[j/3].append(broken[j:j+3])
                    self.vibdisps.extend(p[:(len(broken)/3)])
                    blank = inputfile.next()
                    blank = inputfile.next()
                    freqs = inputfile.next()
                self.vibdisps = Numeric.array(self.vibdisps, "f")
  
            if line[1:24] == "List of All Frequencies":
# Start of the IR/Raman frequency section
                if self.progress and random.random() < fupdate:
                    step = inputfile.tell()
                    if step != oldstep:
                        self.progress.update(step, "Frequency Information")
                        oldstep = step
                         
#                 self.vibsyms = [] # Need to look into this a bit more
                self.vibirs = []
                self.vibfreqs = []
#                 self.logger.info("Creating attribute vibsyms[]")
                self.logger.info("Creating attribute vibfreqs[], vibirs[]")
                for i in range(8):
                    line = inputfile.next()
                line = inputfile.next().strip()
                while line:
                    temp = line.split()
                    self.vibfreqs.append(float(temp[0]))                    
                    self.vibirs.append(float(temp[2])) # or is it temp[1]?
                    line = inputfile.next().strip()
                self.vibfreqs = Numeric.array(self.vibfreqs, "f")
                self.vibirs = Numeric.array(self.vibirs, "f")
                if hasattr(self, "vibramans"):
                    self.vibramans = Numeric.array(self.vibramans, "f")


#******************************************************************************************************************8
#delete this after new implementation using smat, eigvec print,eprint?
            if line[1:49] == "Total nr. of (C)SFOs (summation over all irreps)":
# Extract the number of basis sets
                self.nbasis = int(line.split(":")[1].split()[0])
                self.logger.info("Creating attribute nbasis: %i" % self.nbasis)
                   
   # now that we're here, let's extract aonames
   
                self.logger.info("Creating attribute fonames[]")
                self.fonames = []
                   
                blank = inputfile.next()
                note = inputfile.next()
                symoffset = 0
                blank = inputfile.next(); blank = inputfile.next(); blank = inputfile.next()
                
                nosymreps = []
                while len(self.fonames) < self.nbasis:
                            
                    sym = inputfile.next()
                    line = inputfile.next()
                    num = int(line.split(':')[1].split()[0])
                    nosymreps.append(num)
                       
                    #read until line "--------..." is found
                    while line.find('-----') < 0:
                        line = inputfile.next()
                       
                    line = inputfile.next() # the start of the first SFO

                    while len(self.fonames) < symoffset + num:
                        info = line.split()
                          
                        #index0 index1 occ2 energy3/4 fragname5 coeff6 orbnum7 orbname8 fragname9
                        orbname = info[8]
                        orbital = info[7] + orbname.replace(":", "")
                          
                        fragname = info[5]
                        frag = fragname + info[9]
                          
                        coeff = float(info[6])

                        line = inputfile.next()
                        while line.strip() and not line[:7].strip(): # while it's the same SFO
                            # i.e. while not completely blank, but blank at the start
                            info = line[43:].split()
                            if len(info)>0: # len(info)==0 for the second line of dvb_ir.adfout
                                frag += "+" + fragname + info[-1]
                                coeff = float(info[-4])
                                if coeff < 0:
                                    orbital += '-' + info[-3] + info[-2].replace(":", "")
                                else:
                                    orbital += '+' + info[-3] + info[-2].replace(":", "")
                            line = inputfile.next()
                        # At this point, we are either at the start of the next SFO or at
                        # a blank line...the end

                        self.fonames.append("%s_%s" % (frag, orbital))
                    symoffset += num
                    
                    # blankline blankline
                    inputfile.next(); inputfile.next()
                    
                    
            if line[1:32] == "S F O   P O P U L A T I O N S ,":
#Extract overlap matrix

                self.logger.info("Creating attribute fooverlaps[x, y]")
                self.fooverlaps = Numeric.zeros((self.nbasis, self.nbasis), "float")
                
                symoffset = 0
                
                for nosymrep in nosymreps:
                            
                    line = inputfile.next()
                    while line.find('===') < 10: #look for the symmetry labels
                        line = inputfile.next()
                    #blank blank text blank col row
                    for i in range(6):
                        inputfile.next()
                    
                    base = 0
                    
                    while base < nosymrep: #have we read all the columns?
                          
                        for i in range(nosymrep - base):
                        
                            if self.progress:
                                step = inputfile.tell()
                                if step != oldstep and random.random() < fupdate:
                                    self.progress.update(step, "Overlap")
                                    oldstep = step
                                
                            line = inputfile.next()
                            parts = line.split()[1:]
                            
                            for j in range(len(parts)):
                                k = float(parts[j])
                                self.fooverlaps[base + symoffset + j, base + symoffset +i] = k
                                self.fooverlaps[base + symoffset + i, base + symoffset + j] = k
                              
                        #blank, blank, column
                        for i in range(3):
                            inputfile.next()
                        
                        base += 4
                                        
                    symoffset += nosymrep
                    base = 0
                        
            if line[48:67] == "SFO MO coefficients":
#extract MO coefficients              
                #read stars and three blank lines
                inputfile.next()
                inputfile.next()
                inputfile.next()
                inputfile.next()

                self.logger.info("Creating attribute mocoeffs: array[]")
                
                line = inputfile.next()
                
                if line.find("***** SPIN 1 *****") > 0:
                    beta = 1
                    self.mocoeffs = [Numeric.zeros((self.nbasis, self.nbasis), "float")
                                     for x in range(2)]
                    
                    #get rid of two blank lines and symmetry label
                    inputfile.next()
                    inputfile.next()
                    line = inputfile.next()
                    sym = line.split()[1]
                    
                else:
                    beta = 0
                    self.mocoeffs = [Numeric.zeros((self.nbasis, self.nbasis), "float")]
                    sym = line.split()[1]
                    
                #get rid of 12 lines of text
                for i in range(10):
                    inputfile.next()
                  
                for spin in range(beta + 1):
                    symoffset = 0
                    base = 0
                    
                    if spin == 1:
                        #read spin, blank, blank, blank, symlabel, text, underline, blank
                        for i in range(4):
                            inputfile.next()
                        line = inputfile.next()
                        sym = line.split()[1]
                        for i in range(3):
                            inputfile.next()
                        
                    while symoffset + base < self.nbasis:
                
                        line = inputfile.next()
                        if len(line) < 3:
                            symoffset += base
                            base = 0
                            #print symoffset
                            
                        monumbers = line.split()
                        #print monumbers
                        #get rid of unneeded lines
                        line = inputfile.next()
                        info = line.split()
                        if info[0] == "occup:":
                            inputfile.next()
                        elif info[0] == "===":
                            sym = line.split()[1]
                            inputfile.next()

                        if nosymflag:
                            #aolist = range(len(self.moenergies[spin]))
                            aolist = range(self.nbasis)
                        else:
                            aolist = symlist[sym][spin]
                          
                        row = 0
                        oldindex = 0
                        line = inputfile.next()
                        while len(line) > 5:
                           
                            if self.progress:
                                step = inputfile.tell()
                                if step != oldstep and random.random() < fupdate:
                                    self.progress.update(step, "Coefficients")
                                    oldstep = step
            
                            cols = line.split()
                            for i in range(len(cols[1:])):
                                #self.mocoeffs[spin,row+symoffset,i+symoffset+base]=float(cols[i+1])
                                self.mocoeffs[spin][aolist[i+base], row + symoffset] = float(cols[i + 1])
                                
                            line = inputfile.next()
                            row += 1
                        base += len(cols[1:])
                        
        
if __name__ == "__main__":
    import doctest, adfparser
    doctest.testmod(adfparser, verbose=False)
