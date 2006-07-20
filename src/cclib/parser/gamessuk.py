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
        """Use standard symmetry labels instead of Gaussian labels.

        To normalise:
        (1) If label is one of [SG, PI, PHI, DLTA], replace by [sigma, pi, phi, delta]
        (2) replace any G or U by their lowercase equivalent

        >>> sym = Gaussian("dummyfile").normalisesym
        >>> labels = ['A1', 'AG', 'A1G', "SG", "PI", "PHI", "DLTA"]
        >>> map(sym, labels)
        ['A1', 'Ag', 'A1g', 'sigma', 'pi', 'phi', 'delta']
        """
        _greek = {'SG':'sigma', 'PI':'pi', 'PHI':'phi', 'DLTA':'delta'}
        if label in _greek:
            return _greek[label]

        ans = label.replace("U", "u").replace("G", "g") 
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
            
        optfinished = False # Flag that indicates whether it has reached the end of a geoopt
        
        for line in inputfile:
            
            if self.progress and random.random() < cupdate:
                
                step = inputfile.tell()
                if step != oldstep:
                    self.progress.update(step, "Unsupported Information")
                    oldstep = step
                
                self.scfenergies.append(utils.convertor(self.float(line.split()[4]), "hartree", "eV"))

            if line[49:59] == 'Converged?':
# Extract Geometry convergence information
                if not hasattr(self, "geotargets"):
                    self.logger.info("Creating attributes geotargets[], geovalues[[]]")
                    self.geovalues = []
                    self.geotargets = Numeric.array([0.0, 0.0, 0.0, 0.0], "f")
                newlist = [0]*4

        if self.progress:
            self.progress.update(nstep, "Done")
            
if __name__ == "__main__":
    import doctest, gamessukparser
    doctest.testmod()
