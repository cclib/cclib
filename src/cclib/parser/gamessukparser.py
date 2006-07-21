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
        """
        pass

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
                
        if self.progress:
            self.progress.update(nstep, "Done")

        _toarray = ['atomcoords']
        for attr in _toarray:
            if hasattr(self, attr):
                setattr(self, attr, Numeric.array(getattr(self, attr), 'f'))

            
if __name__ == "__main__":
    import doctest, gamessukparser
    doctest.testmod()
