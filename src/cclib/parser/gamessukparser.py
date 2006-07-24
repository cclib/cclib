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

class GAMESSUK(logfileparser.Logfile):
    """A GAMESS UK log file"""
    SCFRMS, SCFMAX, SCFENERGY = range(3) # Used to index self.scftargets[]
    def __init__(self, *args):

        # Call the __init__ method of the superclass
        super(GAMESSUK, self).__init__(logname="GAMESSUK", *args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "GAMESS UK log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'GAMESSUK("%s")' % (self.filename)
    
    def normalisesym(self, label):
        """Use standard symmetry labels instead of GAMESS UK labels.

        >>> t = GAMESSUK("dummyfile.txt")
        >>> labels = ['au', 'bg', "a'", "a''"]
        >>> answer = [t.normalisesym(x) for x in labels]
        >>> answer
        ['Au', 'Bg', "A'", "A''"]
        """
        ans = label.replace("a", "A").replace("b", "B") 
        
        return ans

    def parse(self, fupdate=0.05, cupdate=0.002):
        """Extract information from the logfile."""
        inputfile = open(self.filename, "r")
        
        if self.progress:
            
            inputfile.seek(0, 2) #go to end of file
            nstep = inputfile.tell()
            inputfile.seek(0)
            self.progress.initialize(nstep)
            oldstep = 0

        for line in inputfile:
            
            if self.progress and random.random() < cupdate:
                
                step = inputfile.tell()
                if step != oldstep:
                    self.progress.update(step, "Unsupported Information")
                    oldstep = step

            if line[1:22] == "total number of atoms":
                if not hasattr(self, "natom"):
                    self.natom = int(line.split()[-1])
                    self.logger.info("Creating attribute natom: %d" % self.natom)

            if line[3:44] == "convergence threshold in optimization run":
                # Assuming that this is only found in the case of OPTXYZ
                # (i.e. an optimization in Cartesian coordinates)
                if not hasattr(self, "geotargets"):
                    self.logger.info("Creating attributes geotargets")
                self.geotargets = [float(line.split()[-2])]

            if line[32:61] == "largest component of gradient":
                # This is the geotarget in the case of OPTXYZ
                if not hasattr(self, "geovalues"):
                    self.logger.info("Creating attribute geovalues")
                    self.geovalues = []
                self.geovalues.append(float(line.split()[4]))

            if line[37:49] == "convergence?":
                # Get the geovalues and geotargets for OPTIMIZE
                if not hasattr(self, "geovalues"):
                    self.logger.info("Creating attribute geotargets, geovalues")
                    self.geovalues = []
                    self.geotargets = []
                geotargets = []
                geovalues = []
                for i in range(4):
                    temp = line.split()
                    geovalues.append(float(temp[2]))
                    if not self.geotargets:
                        geotargets.append(float(temp[-2]))
                    line = inputfile.next()
                self.geovalues.append(geovalues)
                if not self.geotargets:
                    self.geotargets = geotargets
            
            if line[40:59] == "nuclear coordinates":
                if not hasattr(self, "atomcoords"):
                    self.logger.info("Creating attribute atomcoords[], atomnos[]")
                    self.atomcoords = []
                    self.atomnos = []
                    
                asterisk = inputfile.next()
                blank = inputfile.next()
                colmname = inputfile.next()
                equals = inputfile.next()

                atomcoords = []
                atomnos = []
                line = inputfile.next()
                while line != equals:
                    temp = line.strip().split()
                    atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom") for x in temp[0:2]])
                    if not self.atomnos:
                        atomnos.append(int(float(temp[3])))
                        
                    line = inputfile.next()

                self.atomcoords.append(atomcoords)
                if not self.atomnos:
                    self.atomnos = atomnos

            if line[1:32] == "total number of basis functions":
                self.logger.info("Creating attribute nbasis")
                self.nbasis = int(line.split()[-1])
            
            if line[1:36] == "number of occupied orbitals (alpha)":
                if not hasattr(self, "homos"):
                    self.logger.info("Creating attribute homos")
                alpha = int(line.split()[-1])-1
                beta = int(inputfile.next().split()[-1])-1
                self.homos = Numeric.array([alpha, beta], "i")

            if line[37:69] == "s-matrix over gaussian basis set":
                self.logger.info("Creating attribute aooverlaps")
                self.aooverlaps = Numeric.zeros((self.nbasis, self.nbasis), "d")

                minus = inputfile.next()
                blank = inputfile.next()
                for i in range(0,self.nbasis,12):
                    blank = inputfile.next()
                    blank = inputfile.next()
                    header = inputfile.next()
                    blank = inputfile.next()
                    blank = inputfile.next()

                    for j in range(self.nbasis):
                        temp = map(float, inputfile.next().split()[1:])
                        self.aooverlaps[j,(0+i):(len(temp)+i)] = temp
                

            if line[3:27] == "Wavefunction convergence":
                self.logger.info("Creating attribute scftargets")
                scftarget = float(line.split()[-2])
                self.scftargets = []

            if line[3:11] == "SCF TYPE":
                type = line.split()[-2]
                assert type in ['rhf', 'uhf'] # what about rohf?
                if type != 'uhf':
                    self.homos = Numeric.array(self.homos[0])

            if line[15:31] == "convergence data":
                if not hasattr(self, "scfvalues"):
                    self.logger.info("Creating attribute scfvalues")
                    self.scfvalues = []
                self.scftargets.append([scftarget]) # Assuming it does not change over time
                while line[1:10] != "="*9:
                    line = inputfile.next()
                line = inputfile.next()
                while line[1:10] != "="*9: # May be two or three lines (unres)
                    line = inputfile.next()
                
                scfvalues = []
                line = inputfile.next()
                while line.strip():
                    temp = line.split()
                    scfvalues.append([float(temp[5])])
                    line = inputfile.next()
                self.scfvalues.append(scfvalues)   

            if line[10:22] == "total energy":
                if not hasattr(self, "scfenergies"):
                    self.logger.info("Creating attribute scfenergies")
                    self.scfenergies = []
                scfenergy = utils.convertor(float(line.split()[-1]), "hartree", "eV")
                self.scfenergies.append(scfenergy)

            if line[31:50] == "SYMMETRY ASSIGNMENT":
                # (Need IPRINT SCF to get the mosyms for every step of a geoopt)
                if not hasattr(self, "mosyms"):
                    self.logger.info("Creating attribute mosyms")
                    self.mosyms = []
                equals = inputfile.next()
                title = inputfile.next()
                title = inputfile.next()
                equals = inputfile.next()

                mosyms = []
                line = inputfile.next()
                while line != equals:
                    temp = line[91:].strip().split()
                    for i in range(1,len(temp),2):
                        mosyms.append(self.normalisesym(temp[i]))
                    line = inputfile.next()
                assert len(mosyms) == self.nmo
                self.mosyms.append(mosyms)

            if line[2:12] == "m.o. irrep":
                ########## eigenvalues ###########
                # This section appears once at the start of a geo-opt and once at the end
                # unless IPRINT SCF is used (when it appears at every step in addition)
                # This means that the energies of the final geometry appear twice in the
                # case of IPRINT SCF (there's no way around this [unless we test
                # for IPRINT SCF, which I don't like to do], since the alternative would
                # be for the final geometry not to appear at all in the case of no
                # IPRINT SCF, which I think would be a bad idea)
                if not hasattr(self, "moenergies"):
                    self.logger.info("Creating attribute moenergies, nmo")
                    self.moenergies = []

                title = inputfile.next()
                equals = inputfile.next()

                moenergies = []
                line = inputfile.next()
                while line != equals:
                    temp = line.strip().split()
                    moenergies.append(float(temp[3]))
                    line = inputfile.next()
                self.nmo = len(moenergies)
                self.moenergies.append(moenergies)
                
        if self.progress:
            self.progress.update(nstep, "Done")

        _toarray = ['atomcoords', 'geotargets', 'geovalues', 'scftargets',
                    'scfenergies']
        for attr in _toarray:
            if hasattr(self, attr):
                setattr(self, attr, Numeric.array(getattr(self, attr), 'f'))

        if hasattr(self, "scfvalues"):
            self.scfvalues = [Numeric.array(x, "f") for x in self.scfvalues]

             
if __name__ == "__main__":
    import doctest
    doctest.testmod()
