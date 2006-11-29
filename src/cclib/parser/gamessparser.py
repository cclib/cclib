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

class GAMESS(logfileparser.Logfile):
    """A GAMESS log file."""
    SCFRMS, SCFMAX, SCFENERGY = range(3) # Used to index self.scftargets[]
    def __init__(self, *args):

        # Call the __init__ method of the superclass
        super(GAMESS, self).__init__(logname="GAMESS", *args)

    def __str__(self):
        """Return a string representation of the object."""
        return "GAMESS log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'GAMESS("%s")' % (self.filename)

    def normalisesym(self, label):
        """Normalise the symmetries used by GAMESS.

        To normalise, two rules need to be applied:
        (1) Occurences of U/G in the 2/3 position of the label
            must be lower-cased
        (2) Two single quotation marks must be replaced by a double

        >>> t = GAMESS("dummyfile").normalisesym
        >>> labels = ['A', 'A1', 'A1G', "A'", "A''", "AG"]
        >>> answers = map(t, labels)
        >>> print answers
        ['A', 'A1', 'A1g', "A'", 'A"', 'Ag']
        """
        if label[1:] == "''":
            end = '"'
        else:
            end = label[1:].replace("U", "u").replace("G", "g")
        return label[0] + end

    def normalise_aonames(self, listoflines):
        """Normalise the aonames attribute to agree with the other parsers.

        We want this to work even if there are 1000+ atoms. Our only assumption
        is that all of the relevant information is in the first 17 characters
        of the line.
        
        >>> t = GAMESS("dummyfile")
        >>> data = ['    5  C  1  S   ', '    6  C  1  S   ',\
                    '    7  C  1  S   ', '   56  C  1XXXX  ',\
                    '  100  C  2  S   ', '    1  SI  1  S  ' ]
        >>> print t.normalise_aonames(data)
        ['C1_S', 'C1_S', 'C1_S', 'C1_XXXX', 'C2_S', 'Si1_S']
        """
        p = re.compile("(\d+)\s*([A-Z][A-Z]?)\s*(\d+)\s*([A-Z]+)")
        ans = []
        oldatom = "0"
        for line in listoflines:
            m = p.search(line.strip())
            assert m, "Cannot pick out the aoname from this information: %s" % line
            g = m.groups()
            aoname = "%s%s_%s" % (g[1].capitalize(), g[2], g[3])
            oldatom = g[2]
            ans.append(aoname)
            
        return ans
    
    def parse(self):
        """Extract information from the logfile."""
        inputfile = utils.openlogfile(self.filename)
        
        if self.progress:
            
            inputfile.seek(0, 2) #go to end of file
            nstep = inputfile.tell()
            inputfile.seek(0)
            self.progress.initialize(nstep)
            oldstep = 0


        firststdorient = True # Used to decide whether to wipe the atomcoords clean
        geooptfinished = False # Used to avoid extracting the final geometry twice
            
        for line in inputfile:
            
            if self.progress and random.random() < 0.05:
                
                step = inputfile.tell()
                if step != oldstep:
                    self.progress.update(step)
                    oldstep = step

            if line.find("OPTTOL") >= 0:
                # Two possibilities:
                #           OPTTOL = 1.000E-04          RMIN   = 1.500E-03
                # INPUT CARD> $STATPT OPTTOL=0.0001 NSTEP=100 $END
                if not hasattr(self, "geotargets"):
                    self.logger.info("Creating attribute geotargets[]")
                    temp = line.split()
                    for i, x in enumerate(temp):
                        if x.find("OPTTOL") >= 0:
                            if x == "OPTTOL":
                                opttol = float(temp[i + 2])
                            else:
                                opttol = float(x.split('=')[1])
                            self.geotargets = Numeric.array([opttol, 3. / opttol])
                            
            if line.find("FINAL") == 1:
                if not hasattr(self, "scfenergies"):
                    self.logger.info("Creating attribute scfenergies[]")
                    self.scfenergies = []
# Has to deal with such lines as:
#  FINAL R-B3LYP ENERGY IS     -382.0507446475 AFTER  10 ITERATIONS
#  FINAL ENERGY IS     -379.7594673378 AFTER   9 ITERATIONS
# ...so take the number after the "IS"
                temp = line.split()
                self.scfenergies.append(utils.convertor(float(temp[temp.index("IS") + 1]), "hartree", "eV"))

            if line.find("MAXIMUM GRADIENT") > 0:
                if not hasattr(self, "geovalues"):
                    self.logger.info("Creating attribute geovalues[]")
                    self.geovalues = []
                temp = line.strip().split()
                self.geovalues.append([float(temp[3]), float(temp[7])])

            if line[11:50] == "ATOMIC                      COORDINATES":
                # This is the input orientation, which is the only data available for
                # SP calcs, but which should be overwritten by the standard orientation
                # values, which is the only information available for all geoopt cycles.
                if not hasattr(self, "atomcoords"):
                    self.logger.info("Creating attribute atomcoords, atomnos")
                    self.atomcoords = []
                    self.atomnos = []
                line = inputfile.next()
                atomcoords = []
                atomnos = []
                line = inputfile.next()
                while line.strip():
                    temp = line.strip().split()
                    atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom") for x in temp[2:5]])
                    atomnos.append(int(round(float(temp[1])))) # Don't use the atom name as this is arbitary
                    line = inputfile.next()
                self.atomnos = Numeric.array(atomnos, "i")
                self.atomcoords.append(atomcoords)

            if line[12:40] == "EQUILIBRIUM GEOMETRY LOCATED":
                # Prevent extraction of the final geometry twice
                geooptfinished = True
            
            if line[1:29] == "COORDINATES OF ALL ATOMS ARE" and not geooptfinished:
                # This is the standard orientation, which is the only coordinate
                # information available for all geometry optimisation cycles.
                # The input orientation will be overwritten if this is a geometry optimisation
                # We assume that a previous Input Orientation has been found and
                # used to extract the atomnos
                if firststdorient:
                    firststdorient = False
                    # Wipes out the single input coordinate at the start of the file
                    self.atomcoords = []
                    
                line = inputfile.next()
                hyphens = inputfile.next()

                atomcoords = []
                line = inputfile.next()                
                while line.strip():
                    temp = line.strip().split()
                    atomcoords.append(map(float, temp[2:5]))
                    line = inputfile.next()
                self.atomcoords.append(atomcoords)
            
            if line.rstrip()[-15:] == "SCF CALCULATION":
                # This is the section with the SCF information
                line = inputfile.next()
                while line.find("ITER EX") < 0:
                    if line.find("DENSITY CONV=") >= 0 or line.find("DENSITY MATRIX CONV=") >= 0:
# Needs to deal with lines like:
# (GAMESS VERSION = 12 DEC 2003)
#     DENSITY MATRIX CONV=  2.00E-05  DFT GRID SWITCH THRESHOLD=  3.00E-04
# (GAMESS VERSION = 22 FEB 2006)
#           DENSITY MATRIX CONV=  1.00E-05
# (PC GAMESS version 6.2, Not DFT?)
#     DENSITY CONV=  1.00E-05
                        index = line.find("DENSITY CONV=")
                        if index < 0:
                            index = line.find("DENSITY MATRIX CONV=")
                            index += len("DENSITY MATRIX CONV=")
                        else:
                            index += len("DENSITY CONV=")
                        scftarget = float(line[index:].split()[0])
                    line = inputfile.next()

                if not hasattr(self, "scftargets"):
                    self.logger.info("Creating attribute scftargets")
                    self.scftargets = []
                self.scftargets.append([scftarget])

                if not hasattr(self,"scfvalues"):
                    self.logger.info("Creating attribute scfvalues")
                    self.scfvalues = []
                line = inputfile.next()
                values = []
                while line.strip():
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
                        values.append([float(line.split()[5])])
                    line = inputfile.next()
                self.scfvalues.append(values)

            if line.find("NORMAL COORDINATE ANALYSIS IN THE HARMONIC APPROXIMATION") >= 0:
# GAMESS has...
# MODES 1 TO 6 ARE TAKEN AS ROTATIONS AND TRANSLATIONS.
#
#     FREQUENCIES IN CM**-1, IR INTENSITIES IN DEBYE**2/AMU-ANGSTROM**2,
#     REDUCED MASSES IN AMU.
#
#                          1           2           3           4           5
#       FREQUENCY:        52.49       41.45       17.61        9.23       10.61  
#    REDUCED MASS:      3.92418     3.77048     5.43419     6.44636     5.50693
#    IR INTENSITY:      0.00013     0.00001     0.00004     0.00000     0.00003

# whereas PC-GAMESS has...
# MODES 1 TO 6 ARE TAKEN AS ROTATIONS AND TRANSLATIONS.
#
#     FREQUENCIES IN CM**-1, IR INTENSITIES IN DEBYE**2/AMU-ANGSTROM**2
#
#                          1           2           3           4           5
#       FREQUENCY:         5.89        1.46        0.01        0.01        0.01  
#    IR INTENSITY:      0.00000     0.00000     0.00000     0.00000     0.00000

# If Raman is present we have (for PC-GAMESS)...
# MODES 1 TO 6 ARE TAKEN AS ROTATIONS AND TRANSLATIONS.
#
#     FREQUENCIES IN CM**-1, IR INTENSITIES IN DEBYE**2/AMU-ANGSTROM**2
#     RAMAN INTENSITIES IN ANGSTROM**4/AMU, DEPOLARIZATIONS ARE DIMENSIONLESS
#
#                          1           2           3           4           5
#       FREQUENCY:         5.89        1.46        0.04        0.03        0.01  
#    IR INTENSITY:      0.00000     0.00000     0.00000     0.00000     0.00000
# RAMAN INTENSITY:       12.675       1.828       0.000       0.000       0.000
#  DEPOLARIZATION:        0.750       0.750       0.124       0.009       0.750

                self.logger.info("Creating attributes vibfreqs, vibirs")
                self.vibfreqs = []
                self.vibirs = []
                self.logger.info("Creating attributes vibdisps")                
                self.vibdisps = []

                # Need to get past the list of atomic weights
                hyphens = inputfile.next()
                blank = inputfile.next()
                line = inputfile.next()
                blank = inputfile.next()
                line = inputfile.next()
                numAtom = 0
                while line.strip():
                    numAtom += 1
                    line = inputfile.next()

                line = inputfile.next()
                while line.find("FREQUENCIES IN CM**-1") == -1:
                    line = inputfile.next()
                while line != blank:
                    line = inputfile.next()
                
                freqNo = inputfile.next()
                while freqNo.find("SAYVETZ") == -1:
                    freq = inputfile.next().strip().split()
                    self.vibfreqs.extend(map(float, freq[1:]))
                    line = inputfile.next()
                    if line.find("REDUCED") >= 0: # skip the reduced mass (not always present)
                        line = inputfile.next()
                    irIntensity = line.strip().split()
                    self.vibirs.extend(map(float, irIntensity[2:]))
                    line = inputfile.next()
                    if line.find("RAMAN") >= 0:
                        if not hasattr(self,"vibramans"):
                            self.logger.info("Creating attribute vibramans")
                            self.vibramans = []
                        ramanIntensity = line.strip().split()
                        self.vibramans.extend(map(float, ramanIntensity[2:]))
                        depolar = inputfile.next()
                        line = inputfile.next()
                    assert line == blank

                    # Extract the Cartesian displacement vectors
                    p = [ [], [], [], [], [] ]
                    for j in range(numAtom):
                        q = [ [], [], [], [], [] ]
                        for k in range(3): # x, y, z
                            line = inputfile.next()[21:]
                            broken = map(float, line.split())
                            for l in range(len(broken)):
                                q[l].append(broken[l])
                        for k in range(len(broken)):
                            p[k].append(q[k])
                    self.vibdisps.extend(p[:len(broken)])

                    # Skip the Sayvetz stuff at the end
                    for j in range(10):
                        line = inputfile.next()
                    blank = inputfile.next()
                    freqNo = inputfile.next()
                self.vibfreqs = Numeric.array(self.vibfreqs, "f")
                self.vibirs = Numeric.array(self.vibirs, "f")
                self.vibdisps = Numeric.array(self.vibdisps, "f")

            if line[5:21] == "ATOMIC BASIS SET":
                if not hasattr(self, "gbasis"):
                    self.logger.info("Creating attribute gbasis")
                    self.gbasis = []
                line = inputfile.next()
                while line.find("SHELL")<0:
                    line = inputfile.next()
                blank = inputfile.next()
                atomname = inputfile.next()
                while line.find("TOTAL NUMBER")<0:
                    gbasis = [] # Stores basis sets on one atom
                    blank = inputfile.next()
                    line = inputfile.next()
                    while len(line.split())!=1 and line.find("TOTAL NUMBER")<0:
                        coeff = {}
                        # coefficients and symmetries for a block of rows
                        while line.strip():
                            temp = line.strip().split()
                            sym = temp[1]
                            assert sym in ['S', 'P', 'D', 'F', 'G', 'L']
                            if sym == "L": # L refers to SP
                                if len(temp)==6: # GAMESS US
                                    coeff.setdefault("S", []).append( (float(temp[3]), float(temp[4])) )
                                    coeff.setdefault("P", []).append( (float(temp[3]), float(temp[5])) )
                                else: # PC GAMESS
                                    assert temp[6][-1] == temp[9][-1] == ')'
                                    coeff.setdefault("S", []).append( (float(temp[3]), float(temp[6][:-1])) )
                                    coeff.setdefault("P", []).append( (float(temp[3]), float(temp[9][:-1])) )
                            else:
                                if len(temp)==5: # GAMESS US
                                    coeff.setdefault(sym, []).append( (float(temp[3]), float(temp[4])) )
                                else: # PC GAMESS
                                    assert temp[6][-1] == ')'
                                    coeff.setdefault(sym, []).append( (float(temp[3]), float(temp[6][:-1])) )
                            line = inputfile.next()
# either a blank or a continuation of the block
                        if sym == "L":
                            gbasis.append( ('S', coeff['S']))
                            gbasis.append( ('P', coeff['P']))
                        else:
                            gbasis.append( (sym, coeff[sym]))
                        line = inputfile.next()
# either the start of the next block or the start of a new atom or
# the end of the basis function section
                    self.gbasis.append(gbasis)

            if line.find("EIGENVECTORS") == 10 or line.find("MOLECULAR OBRITALS") == 10:
                # The details returned come from the *final* report of evalues and
                # the last list of symmetries in the log file
                # This is fine for GeoOpt and SP, but may be weird for TD and Freq(?)
                
                # Take the last one of either in the file
                if not hasattr(self, "moenergies"):
                    self.logger.info("Creating attributes moenergies, mosyms")
                self.moenergies = [[]]
                self.mosyms = [[]]
                if not hasattr(self, "nmo"):
                    self.logger.info("Creating attribute nmo with default value")
                    self.nmo = self.nbasis
                if not hasattr(self, "mocoeffs"):
                    self.logger.info("Creating attribute mocoeffs")
                self.mocoeffs = Numeric.zeros((1, self.nmo, self.nbasis), "f")
                line = inputfile.next()
                for base in range(0, self.nmo, 5):
                    blank = inputfile.next()
                    line = inputfile.next() # Eigenvector no
                    line = inputfile.next()
                    self.moenergies[0].extend([utils.convertor(float(x), "hartree", "eV") for x in line.split()])
                    line = inputfile.next()
                    self.mosyms[0].extend(map(self.normalisesym, line.split()))
                    for i in range(self.nbasis):
                        line = inputfile.next()
                        # if base==0: # Just do this the first time 'round
                            # atomno=int(line.split()[2])-1
                            # atomorb[atomno].append(int(line.split()[0])-1)
                            # What's the story with the previous line?
                        temp = line[15:] # Strip off the crud at the start
                        j = 0
                        while j*11+4 < len(temp):
                            self.mocoeffs[0, base+j, i] = float(temp[j * 11:(j + 1) * 11])
                            j += 1
                line = inputfile.next()
                if line.find("END OF RHF") == -1: # implies unrestricted
# If it's restricted we have
#  ...... END OF RHF CALCULATION ......
# If it's unrestricted we have...
#
#  ----- BETA SET ----- 
#
#          ------------
#          EIGENVECTORS
#          ------------
#
#                      1          2          3          4          5

                    self.mocoeffs.resize((2, self.nmo, self.nbasis))
                    self.moenergies.append([])
                    self.mosyms.append([])
                    for i in range(5):
                        line = inputfile.next()
                    for base in range(0, self.nmo, 5):
                        blank = inputfile.next()
                        line = inputfile.next() # Eigenvector no
                        line = inputfile.next()
                        self.moenergies[1].extend([utils.convertor(float(x), "hartree", "eV") for x in line.split()])
                        line = inputfile.next()
                        self.mosyms[1].extend(map(self.normalisesym, line.split()))
                        for i in range(self.nbasis):
                            line = inputfile.next()
                            temp = line[15:] # Strip off the crud at the start
                            j = 0
                            while j * 11 + 4 < len(temp):
                                self.mocoeffs[1, base+j, i] = float(temp[j * 11:(j + 1) * 11])
                                j += 1
                    line = inputfile.next()
                assert line.find("END OF") >= 0
                self.moenergies = Numeric.array(self.moenergies, "f")

            if line.find("NUMBER OF OCCUPIED ORBITALS") >= 0:
                if not hasattr(self," homos"):
                    self.logger.info("Creating attribute homos")
                homos = [int(line.split()[-1])-1]
                line = inputfile.next()
                homos.append(int(line.split()[-1])-1)
                # Note that we cannot trust this self.homos until we come to
                # a line that contains the phrase:
                # "SYMMETRIES FOR INITAL GUESS ORBITALS FOLLOW"
                # which either is followed by "ALPHA" or "BOTH"
                # at which point we can say for certain that it is an
                # un/restricted calculations
                self.homos = Numeric.array(homos, "i")

            if line.find("SYMMETRIES FOR INITIAL GUESS ORBITALS FOLLOW") >= 0:
                # Not unrestricted, so lop off the second index
                if line.find("BOTH SET(S)") >= 0:
                    self.homos = Numeric.resize(self.homos, [1])

            if line.find("TOTAL NUMBER OF ATOMS") == 1:
                self.logger.info("Creating attribute natom")
                self.natom = int(line.split()[-1])
                
            if line.find("NUMBER OF CARTESIAN GAUSSIAN BASIS") == 1 or line.find("TOTAL NUMBER OF BASIS FUNCTIONS") == 1:
                # The first is from Julien's Example and the second is from Alexander's
                # I think it happens if you use a polar basis function instead of a cartesian one
                self.logger.info("Creating attribute nbasis")
                self.nbasis = int(line.strip().split()[-1])
                    
            elif line.find("SPHERICAL HARMONICS KEPT IN THE VARIATION SPACE") >= 0:
                # Note that this line is present if ISPHER=1, e.g. for C_bigbasis
                if not hasattr(self, "nmo"):
                    self.logger.info("Creating attribute nmo")
                self.nmo = int(line.strip().split()[-1])
                
            elif line.find("TOTAL NUMBER OF MOS IN VARIATION SPACE") == 1:
                # Note that this line is not always present, so by default
                # NBsUse is set equal to NBasis (see below).
                if not hasattr(self, "nmo"):
                    self.logger.info("Creating attribute nmo")
                self.nmo = int(line.split()[-1])

            elif line.find("OVERLAP MATRIX") == 0 or line.find("OVERLAP MATRIX") == 1:
                # The first is for PC-GAMESS, the second for GAMESS
                # Read 1-electron overlap matrix
                if not hasattr(self, "aooverlaps"):
                    self.logger.info("Creating attribute aooverlaps, aonames")
                    self.aooverlaps = Numeric.zeros((self.nbasis, self.nbasis), "f")
                    self.aonames = []
                else:
                    self.logger.info("Reading additional aooverlaps...")
                base = 0
                aonames = []
                while base < self.nbasis:
                    blank = inputfile.next()
                    line = inputfile.next() # Basis fn number
                    blank = inputfile.next()
                    for i in range(self.nbasis - base): # Fewer lines each time
                        line = inputfile.next()
                        temp = line.split()
                        if base == 0: # Only do this for the first block
                            aonames.append(line[:17])
                        for j in range(4, len(temp)):
                            self.aooverlaps[base+j-4, i+base] = float(temp[j])
                            self.aooverlaps[i+base, base+j-4] = float(temp[j])
                    base += 5
                self.aonames = self.normalise_aonames(aonames)

        inputfile.close()

        if not hasattr(self, "geotargets"):
            self.logger.info("Creating attribute geotargets[] with default values")
            opttol = 1e-4
            self.geotargets = Numeric.array([opttol, 3. / opttol])
        if hasattr(self, "scftargets"):
            self.scftargets = Numeric.array(self.scftargets, "f")
        if hasattr(self, "scfvalues"):
            self.scfvalues = [Numeric.array(x, "f") for x in self.scfvalues]
        if hasattr(self,"geovalues"): self.geovalues = Numeric.array(self.geovalues, "f")
        if hasattr(self, "atomcoords"):
            self.atomcoords = Numeric.array(self.atomcoords, "f")
        if not hasattr(self, "nmo"):
            self.logger.info("Creating attribute nmo with default value")
            self.nmo = self.nbasis
        if not hasattr(self,"coreelectrons"):
            self.coreelectrons = Numeric.zeros(self.natom, "i")

        self.parsed = True


        
if __name__ == "__main__":
    import doctest, gamessparser
    doctest.testmod(gamessparser, verbose=False)
