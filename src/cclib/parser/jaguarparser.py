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

    def parse(self, fupdate=0.05, cupdate=0.002):
        """Extract information from the logfile."""
        inputfile = open(self.filename, "r")
        
        if self.progress:
            
            inputfile.seek(0, 2) #go to end of file
            nstep = inputfile.tell()
            inputfile.seek(0)
            self.progress.initialize(nstep)
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

            if line[1:5] == "SCFE":
# Get the energy of the molecule
                if not hasattr(self, "scfenergies"):
                    self.logger.info("Creating attribute scfenergies")
                    self.scfenergies = []
                temp = line.strip().split()
                self.scfenergies.append(utils.convertor(float(temp[temp.index("hartrees") - 1]), "hartree", "eV"))

            if line[2:14] == "new geometry" or line[1:21] == "Symmetrized geometry":
# Get the atom coordinates
                if not hasattr(self, "atomcoords"):
                    self.logger.info("Creating attributes: atomcoords, atomnos, natom")
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
                    geotargets = Numeric.zeros(4, "f")
                    i = 0
                    self.logger.info("Creating attributs: geovalues,geotargets")
                blank = inputfile.next()
                blank = inputfile.next()
                line = inputfile.next()
                values = []
                while line != blank:
                    if line[41] == "(":
                        # A new geo convergence value
                        values.append(float(line[26:37]))
                        if not hasattr(self,"geotargets"):
                            geotargets[i] = float(line[43:54])
                            i += 1
                    line = inputfile.next()
                self.geovalues.append(values)
                if not hasattr(self, "geotargets"):
                    self.geotargets = geotargets

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
                        self.mosyms[0].append(temp[i+1])
                    line = inputfile.next()
                self.moenergies = Numeric.array(self.moenergies, "f")
            
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
                self.logger.info("Creating attribute: vibfreqs")
                self.vibfreqs = []
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
                    while line.strip(): # Read the cartesian displacements
                        line = inputfile.next()
                    freqs = inputfile.next()
                self.vibfreqs = Numeric.array(self.vibfreqs)
                if hasattr(self, "vibirs"):
                    self.vibirs = Numeric.array(self.vibirs)

        inputfile.close()

        if hasattr(self, "scfvalues"):
            self.scfvalues = [Numeric.array(x, "f") for x in self.scfvalues]
        if hasattr(self, "scftargets"):
            self.scftargets = Numeric.array(self.scftargets, "f")
        if hasattr(self, "scfenergies"):
            self.scfenergies = Numeric.array(self.scfenergies, "f")
        if hasattr(self, "atomcoords"):
            self.atomcoords = Numeric.array(self.atomcoords, "f")
        self.parsed = True
        
if __name__ == "__main__":
    import doctest, jaguarparser
    doctest.testmod(jaguarparser, verbose=False)
