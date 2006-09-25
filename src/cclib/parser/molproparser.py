import re
import Numeric
import random # For sometimes running the progress updater
import utils
import logfileparser

class Molpro(logfileparser.Logfile):
    """Molpro file parser"""
    SCFMAX,SCFENERGY = range(2) # Used to index self.scftargets[]
    def __init__(self,*args):
        # Call the __init__ method of the superclass
        super(Molpro, self).__init__(logname="Molpro",*args)

    def __str__(self):
        """Return a string representation of the object."""
        return "Molpro log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Molpro("%s")' % (self.filename)

    def normalisesym(self,label):
        """Normalise the symmetries used by Molpro."""
        ans = label.replace("`", "'").replace("``", "''")
        return ans

    def normalise_aonames(self,listoflines):
        return

    def parse(self, fupdate=0.05, cupdate=0.002):
        """Extract information from the logfile."""
        inputfile = open(self.filename,"r")
        if self.progress:
            inputfile.seek(0,2) #go to end of file
            nstep=inputfile.tell()
            inputfile.seek(0)
            self.progress.initialize(nstep)
            oldstep=0

        optfinished = False
            
        for line in inputfile:
            if self.progress and random.random() < cupdate:
                step = inputfile.tell()
                if step!=oldstep:
                    self.progress.update(step)
                    oldstep = step
            if line[1:23] == "CONVERGENCE THRESHOLDS":
                if not hasattr(self, "scftargets"):
                    self.logger.info("Creating attribute scftargets[]")
                self.scftargets = []
                # The MAX density matrix
                # The MAX Energy
                scftargets = map(float, line.split()[2::2])
                self.scftargets.append(scftargets)
            if line[1:23] == "NUMBER OF CONTRACTIONS":
                nbasis = int(line.split()[3])
                if hasattr(self, "nbasis"):
                    assert nbasis == self.nbasis
                else:
                    self.nbasis = nbasis
                    self.logger.info("Creating attribute nbasis: %d" % self.nbasis)

            if line[1:19] == "ATOMIC COORDINATES":
                if not hasattr(self,"atomcoords"):
                    self.logger.info("Creating attribute atomcoords, atomnos")
                    self.atomcoords = []
                    self.atomnos = []
                line = inputfile.next()
                line = inputfile.next()
                line = inputfile.next()
                atomcoords = []
                atomnos = []
                line = inputfile.next()
                while line.strip():
                    temp = line.strip().split()
                    atomcoords.append([utils.convertor(float(x),"bohr","Angstrom") for x in temp[3:6]]) #bohrs to angs
                    atomnos.append(int(round(float(temp[2]))))
                    line = inputfile.next()
                self.atomnos = Numeric.array(atomnos,"i")
                self.atomcoords.append(atomcoords)
                self.logger.info("Creating attribute natom")
                self.natom = len(self.atomnos)
            if line[1:13] == "Normal Modes" and not hasattr(self,"vibfreqs"):
                self.logger.info("Creating attributes vibsyms, vibfreqs, vibirs")
                self.vibfreqs = []
                self.vibirs = []
                self.vibsyms = []
                line = inputfile.next()
                line = inputfile.next()
                for i in line.split():
                    try: map(int, i)
                    except: self.vibsyms.append(i)
                while line.strip():
                    if line[1:12] == "Wavenumbers":
                        freqs = line.strip().split()[2:]
                        self.vibfreqs.extend(map(float,freqs))
                    if line[1:21] == "Intensities [km/mol]":
                        irIntensity = line.strip().split()[2:]
                        self.vibirs.extend(map(float,irIntensity))
                    line = inputfile.next()
                    if not line.strip():
                        line = inputfile.next()
                        if not line.strip(): break
                        else: self.vibsyms += line.split()[1::2]
            if line[1:20] == "Normal Modes of low" and hasattr(self,"vibfreqs"):
                line = inputfile.next()
                line = inputfile.next()
                self.lowfreqs = []
                self.lowirs = []
                while line.strip():
                    if line[1:12] == "Wavenumbers":
                        freqs = line.strip().split()[2:]
                        self.lowfreqs.extend(map(float,freqs))
                    if line[1:21] == "Intensities [km/mol]":
                        irIntensity = line.strip().split()[2:]
                        self.lowirs.extend(map(float,irIntensity))
                    line = inputfile.next()
                    if not line.strip():
                        line = inputfile.next()
                        if not line.strip(): break
                        else: continue
                self.vibfreqs = self.lowfreqs + self.vibfreqs
                self.vibirs = self.lowirs + self.vibirs
            
            if line[1:16] == "Force Constants":
                self.logger.info("Creating attribute hessian")
                self.hessian = []
                line = inputfile.next()
                hess = []
                tmp = []
                while line.strip():
                    try: map(float, line.strip().split()[2:])
                    except: 
                        line = inputfile.next()
                    line.strip().split()[1:]
                    hess.extend([map(float,line.strip().split()[1:])])
                    line = inputfile.next()
                line = 0
                while (line==0) or (len(hess[0]) > 1):
                    tmp.append(hess.pop(0))
                    line += 1
                k = 5
                while len(hess) != 0:
                    tmp[k] += hess.pop(0)
                    k += 1
                    if (len(tmp[k-1]) == line): break
                    if k >= line: k = len(tmp[-1])
                for l in tmp: self.hessian += l

        inputfile.close()

        if self.progress:
            self.progress.update(nstep, "Done")
        _toarray = ['atomcoords','scftargets', 'vibfreqs', 'vibirs', 'hessian']
        _normalise = ['vibsyms','mosyms']
        for attr in _toarray:
            if hasattr(self, attr):
                setattr(self, attr, Numeric.array(getattr(self, attr), 'f'))
        for attr in _normalise:
            if hasattr(self, attr):
                for i in xrange(len(self.vibsyms)):
                    self.vibsyms[i] = self.normalisesym(self.vibsyms[i])

        self.parsed = True

if __name__=="__main__":
    import doctest,molproparser
    doctest.testmod(molproparser,verbose=False)
