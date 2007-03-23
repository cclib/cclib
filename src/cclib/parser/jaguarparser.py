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

class Jaguar(logfileparser.Logfile):
    """A Jaguar output file"""

    def __init__(self, *args):

        # Call the __init__ method of the superclass
        super(Jaguar, self).__init__(logname="Jaguar", *args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Jaguar output file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Jaguar("%s")' % (self.filename)

    def normalisesym(self, label):
        """Normalise the symmetries used by Jaguar.

        To normalise, three rules need to be applied:
        (1) To handle orbitals of E symmetry, retain everything before the /
        (2) Replace two p's by "
        (2) Replace any remaining single p's by '

        >>> t = Jaguar("dummyfile").normalisesym
        >>> labels = ['A', 'A1', 'Ag', 'Ap', 'App', "A1p", "A1pp", "E1pp/Ap"]
        >>> answers = map(t, labels)
        >>> print answers
        ['A', 'A1', 'Ag', "A'", 'A"', "A1'", 'A1"', 'E1"']
        """
        ans = label.split("/")[0].replace("pp", '"').replace("p", "'")
        return ans

    def extract(self, inputfile, fupdate=0.05, cupdate=0.002):
        """Extract information from the file object inputfile."""

        oldstep = 0

        geoopt = False # Is this a GeoOpt? Needed for SCF targets/values.
            
        for line in inputfile:
            
            if self.progress and random.random() < cupdate:
                
                step = inputfile.tell()
                if step != oldstep:
                    self.progress.update(step, "Unsupported Information")
                    oldstep = step

            if line[0:4] == "etot":
# Get SCF convergence information
                if not hasattr(self, "scfvalues"):
                    self.scfvalues = []
                    self.logger.info("Creating attribute: scfvalues,scftargets")
                    self.scftargets = [[5E-5, 5E-6]]
                values = []
                while line[0:4] == "etot":
# Jaguar 4.2
# etot   1  N  N  0  N  -382.08751886450           2.3E-03  1.4E-01
# etot   2  Y  Y  0  N  -382.27486023153  1.9E-01  1.4E-03  5.7E-02
# Jaguar 6.5
# etot   1  N  N  0  N    -382.08751881733           2.3E-03  1.4E-01
# etot   2  Y  Y  0  N    -382.27486018708  1.9E-01  1.4E-03  5.7E-02
                    temp = line.split()[7:]
                    if len(temp)==3:
                        denergy = float(temp[0])
                    else:
                        denergy = 0 # Should really be greater than target value
                                    # or should we just ignore the values in this line
                    ddensity = float(temp[-2])
                    maxdiiserr = float(temp[-1])
                    if not geoopt:
                        values.append([denergy, ddensity])
                    else:
                        values.append([ddensity])
                    line = inputfile.next()
                self.scfvalues.append(values)

            # Hartree-Fock energy after SCF
            if line[1:18] == "SCFE: SCF energy:":
                if not hasattr(self, "scfenergies"):
                    self.logger.info("Creating attribute scfenergies")
                    self.scfenergies = []
                temp = line.strip().split()
                scfenergy = float(temp[temp.index("hartrees") - 1])
                scfenergy = utils.convertor(scfenergy, "hartree", "eV")
                self.scfenergies.append(scfenergy)

            # Energy after LMP2 correction
            if line[1:18] == "Total LMP2 Energy":
                if not hasattr(self, "mpenergies"):
                    self.logger.info("Creating attribute mpenergies")
                    self.mpenergies = [[]]
                lmp2energy = float(line.split()[-1])
                lmp2energy = utils.convertor(lmp2energy, "hartree", "eV")
                self.mpenergies[-1].append(lmp2energy)

            if line[2:14] == "new geometry" or line[1:21] == "Symmetrized geometry" or line.find("Input geometry") > 0:
# Get the atom coordinates
                if not hasattr(self, "atomcoords"):
                    self.logger.info("Creating attributes: atomcoords, atomnos, natom")
                if not hasattr(self, "atomcoords") or line[1:21] == "Symmetrized geometry":
                    # Wipe the "Input geometry" if "Symmetrized geometry" present
                    self.atomcoords = []
                p = re.compile("(\D+)\d+") # One/more letters followed by a number
                atomcoords = []
                atomnos = []
                angstrom = inputfile.next()
                title = inputfile.next()
                line = inputfile.next()
                while line.strip():
                    temp = line.split()
                    element = p.findall(temp[0])[0]
                    atomnos.append(self.table.number[element])
                    atomcoords.append(map(float, temp[1:]))
                    line = inputfile.next()
                self.atomcoords.append(atomcoords)
                self.atomnos = Numeric.array(atomnos, "i")
                self.natom = len(atomcoords)

            if line[2:24] == "start of program geopt":
                if not geoopt:
                    # Need to keep only the RMS density change info
                    # if this is a geoopt
                    self.scftargets = [[self.scftargets[0][0]]]
                    if hasattr(self, "scfvalues"):
                        self.scfvalues[0] = [[x[0]] for x in self.scfvalues[0]]
                    geoopt = True
                else:
                    self.scftargets.append([5E-5])

            if line[2:28] == "geometry optimization step":
# Get Geometry Opt convergence information
                if not hasattr(self, "geovalues"):
                    self.geovalues = []
                    self.geotargets = Numeric.zeros(5, "d")
                    self.logger.info("Creating attributes: geovalues,geotargets")
                gopt_step = int(line.split()[-1])
                blank = inputfile.next()
                blank = inputfile.next()
                line = inputfile.next()
                values = []
                target_index = 0                
                if gopt_step == 1:
                    # The first optimization step does not produce an energy change
                    values.append(0.0)
                    target_index = 1
                while line != blank:
                    if line[41] == "(":
                        # A new geo convergence value
                        values.append(float(line[26:37]))
                        self.geotargets[target_index] = float(line[43:54])
                        target_index += 1
                    line = inputfile.next()
                self.geovalues.append(values)

            if line.find("number of occupied orbitals") > 0:
# Get number of MOs
                occs = int(line.split()[-1])
                line = inputfile.next()
                virts = int(line.split()[-1])
                self.nmo = occs + virts
                self.homos = Numeric.array([occs-1], "i")

                unrestrictedflag = False

            if line.find("number of alpha occupied orb") > 0:
# Get number of MOs for an unrestricted calc

                aoccs = int(line.split()[-1])
                line = inputfile.next()
                avirts = int(line.split()[-1])
                line = inputfile.next()
                boccs = int(line.split()[-1])
                line = inputfile.next()
                bvirt = int(line.split()[-1])

                self.nmo = aoccs + avirts
                self.homos = Numeric.array([aoccs-1,boccs-1], "i")
                unrestrictedflag = True
                
            if line[2:33] == "Orbital energies/symmetry label":
# Get MO Energies and symmetrys
                if not hasattr(self,"moenergies"):
                    self.logger.info("Creating attributes: moenergies, mosyms")
                self.mosyms = [[]]
                self.moenergies = [[]]
                line = inputfile.next()
                while line.strip():
                    temp = line.strip().split()
                    for i in range(0, len(temp), 2):
                        self.moenergies[0].append(utils.convertor(float(temp[i]), "hartree", "eV"))
                        self.mosyms[0].append(self.normalisesym(temp[i+1]))
                    line = inputfile.next()
                self.moenergies = [Numeric.array(self.moenergies[0], "d")]

            if line.find("Orbital energies:") == 2:
# Get MO Energies
# Jaguar 6.0
                if not hasattr(self,"moenergies"):
                    self.logger.info("Creating attribute moenergies")
                    self.moenergies = [[]]
                    line = inputfile.next()
                    while line.strip():
                        temp = line.strip().split()
                        for i in range(len(temp)):
                            self.moenergies[0].append(utils.convertor(float(temp[i]), "hartree", "eV"))
                        line = inputfile.next()
                    self.moenergies = [Numeric.array(self.moenergies[0], "d")]

            if line.find("Alpha Orbital energies:") == 2:
# Get alpha MO Energies
# Jaguar 6.0
                if not hasattr(self,"moenergies"):
                    self.logger.info("Creating attribute moenergies")
                    self.moenergies = [[],[]]
                    line = inputfile.next()
                    while line.strip():
                        temp = line.strip().split()
                        for i in range(len(temp)):
                            self.moenergies[0].append(utils.convertor(float(temp[i]), "hartree", "eV"))
                        line = inputfile.next()

                    blank = inputfile.next()
                    homo = inputfile.next()
                    lumo = inputfile.next()
                    blank = inputfile.next()
                    title = inputfile.next()

                    line = inputfile.next()
                    while line.strip():
                        temp = line.strip().split()
                        for i in range(len(temp)):
                            self.moenergies[1].append(utils.convertor(float(temp[i]), "hartree", "eV"))
                        line = inputfile.next()

                    self.moenergies = [Numeric.array(x, "d") for x in self.moenergies]

            if line.find("Occupied + virtual Orbitals- final wvfn") > 0:
                
                blank = inputfile.next()
                stars = inputfile.next()
                blank = inputfile.next()
                blank = inputfile.next()
                
                if not hasattr(self,"mocoeffs"):
                    self.logger.info("Creating mocoeffs")
                    if unrestrictedflag:
                        spin = 2
                    else:
                        spin = 1

                    self.mocoeffs = []
                    
                
                aonames = []
                lastatom = "X"

                offset = 0

                for s in range(spin):
                    mocoeffs = Numeric.zeros((len(self.moenergies[s]), self.nbasis), "d")

                    if s == 1: #beta case
                        stars = inputfile.next()
                        blank = inputfile.next()
                        title = inputfile.next()
                        blank = inputfile.next()
                        stars = inputfile.next()
                        blank = inputfile.next()
                        blank = inputfile.next()

                    for k in range(0,len(self.moenergies[s]),5):

                        numbers = inputfile.next()
                        eigens = inputfile.next()
                        line = inputfile.next()

                        for i in range(self.nbasis):

                            info = line.split()

                            if not hasattr(self,"aonames"):
                                if lastatom != info[1]:
                                    scount = 1
                                    pcount = 3
                                    dcount = 6 #six d orbitals in Jaguar

                                if info[2] == 'S':
                                    aonames.append("%s_%i%s"%(info[1], scount, info[2]))
                                    scount += 1
                            
                                if info[2] == 'X' or info[2] == 'Y' or info[2] == 'Z':
                                    aonames.append("%s_%iP%s"%(info[1], pcount / 3, info[2]))
                                    pcount += 1
                            
                                if info[2] == 'XX' or info[2] == 'YY' or info[2] == 'ZZ' or \
                                   info[2] == 'XY' or info[2] == 'XZ' or info[2] == 'YZ':

                                    aonames.append("%s_%iD%s"%(info[1], dcount / 6, info[2]))
                                    dcount += 1

                                lastatom = info[1]

                            for j in range(len(info[3:])):
                                mocoeffs[j+k,i] = float(info[3+j])

                            line = inputfile.next()

                        if not hasattr(self,"aonames"):
                            self.aonames = aonames

                        offset += 5
                    self.mocoeffs.append(mocoeffs)
                            
                            
            if line[2:6] == "olap":
                if line[6]=="-":
                    continue # avoid "olap-dev"
                self.logger.info("Creating attribute aooverlaps")
                self.aooverlaps = Numeric.zeros((self.nbasis, self.nbasis), "d")

                for i in range(0, self.nbasis, 5):
                    blank = inputfile.next()
                    header = inputfile.next()
                    for j in range(i, self.nbasis):
                        temp = map(float, inputfile.next().split()[1:])
                        self.aooverlaps[j, i:(i+len(temp))] = temp
                        self.aooverlaps[i:(i+len(temp)), j] = temp
                
            if line[1:28] == "number of occupied orbitals":
                if not hasattr(self, "homos"):
                    self.logger.info("Creating attribute: homos")
                self.homos = Numeric.array([float(line.strip().split()[-1])-1], "i")

            if line[2:27] == "number of basis functions":
                if not hasattr(self, "nbasis"):
                    self.logger.info("Creating attribute: nbasis")
                self.nbasis = int(line.strip().split()[-1])

            if line[2:23] == "start of program freq":
# IR stuff
                self.logger.info("Creating attribute: vibfreqs, vibdisps")
                self.vibfreqs = []
                self.vibdisps = []
                blank = inputfile.next()
                line = inputfile.next()
                line = inputfile.next()
                forceconstants = False
                if line.find("force constants")>=0:
                    forceconstants = True
                    # Could handle this differently if a problem in future
                line = inputfile.next()
                while line.strip():
                    line = inputfile.next()
                
                freqs = inputfile.next()
                while freqs.strip():
                    temp = freqs.strip().split()
                    self.vibfreqs.extend(map(float, temp[1:]))
                    temp = inputfile.next().strip().split()
                    if temp[0] == "symmetries": # May go straight from frequencies to reduced mass
                        if not hasattr(self, "vibsyms"):
                            self.logger.info("Creating attributes: vibsyms, vibirs")
                            self.vibsyms = []
                            self.vibirs = []
                        self.vibsyms.extend(map(self.normalisesym, temp[1:]))
                        temp = inputfile.next().strip().split()                                
                        self.vibirs.extend(map(float, temp[1:]))
                        reducedmass = inputfile.next()
                        if forceconstants:
                            forceconst = inputfile.next()
                    line = inputfile.next()

                    q = [ [] for i in range(7) ] # Hold up to 7 lists of triplets
                    ## while line.strip(): # Read the cartesian displacements
                    for x in range(len(self.atomnos)):
                        p = [ [] for i in range(7) ] # Hold up to 7 triplets
                        for i in range(3):
                            broken = map(float, line.strip().split()[2:])
                            for j in range(len(broken)):
                                p[j].append(broken[j])
                            line = inputfile.next()
                        for i in range(len(broken)):
                            q[i].append(p[i])

                    self.vibdisps.extend(q[:len(broken)])
                    freqs = inputfile.next()
                self.vibfreqs = Numeric.array(self.vibfreqs, "d")
                self.vibdisps = Numeric.array(self.vibdisps, "d")
                if hasattr(self, "vibirs"):
                    self.vibirs = Numeric.array(self.vibirs, "d")
                    
            # Parse excited state output (for CIS calculations).
            if line[2:15] == "Excited State":
                if not hasattr(self, "etenergies"):
                    self.etenergies = []
                    self.logger.info("Creating attribute: etenergies")
                if not hasattr(self, "etoscs"):
                    self.etoscs = []
                    self.logger.info("Creating attribute: etoscs")
                etenergy = float(line.split()[3])
                self.etenergies.append(etenergy)
                while line[2:21] != "Oscillator strength":
                    line = inputfile.next()
                strength = float(line.split()[-1])
                self.etoscs.append(strength)

        
if __name__ == "__main__":
    import doctest, jaguarparser
    doctest.testmod(jaguarparser, verbose=False)
