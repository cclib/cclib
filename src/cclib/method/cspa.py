"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

import Numeric
import random # For sometimes running the progress updater
from population import Population

class CSPA(Population):
    """The C-squared population analysis"""
    def __init__(self, *args):

        # Call the __init__ method of the superclass
        super(CSPA, self).__init__(logname="CSPA", *args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "CSPA of" % (self.parser)

    def __repr__(self):
        """Return a representation of the object."""
        return 'CSPA("%s")' % (self.parser)
    
    def calculate(self, indices=None, fupdate=0.05):
        """Perform a C^2 population analysis given the results of a parser"""
    
        if not self.parser.parsed:
            self.parser.parse()

#do we have the needed info in the parser?
        if not hasattr(self.parser, "mocoeffs"):
            self.logger.error("Missing mocoeffs")
            return False
        if not hasattr(self.parser, "nbasis"):
            self.logger.error("Missing nbasis")
            return False
        if not hasattr(self.parser, "homos"):
            self.logger.error("Missing homos")
            return False
#end attributes check

        self.logger.info("Creating attribute aoresults: array[3]")
        unrestricted = (len(self.parser.mocoeffs)==2)
        nmocoeffs = len(self.parser.mocoeffs[0])
        nbasis = self.parser.nbasis

        #determine number of steps, and whether process involves beta orbitals
        nstep = nmocoeffs
        if unrestricted:
            self.aoresults = Numeric.zeros([2, nmocoeffs, nbasis], "f")
            nstep += nmocoeffs
        else:
            self.aoresults = Numeric.zeros([1, nmocoeffs, nbasis], "f")

        #intialize progress if available
        if self.progress:
            self.progress.initialize(nstep)

        step = 0
        for spin in range(len(self.parser.mocoeffs)):

            for i in range(nmocoeffs):

                if self.progress and random.random() < fupdate:
                    self.progress.update(step, "C^2 Population Analysis")

                submocoeffs = self.parser.mocoeffs[spin][i]
                scale = Numeric.innerproduct(submocoeffs, submocoeffs)
                tempcoeffs = Numeric.multiply(submocoeffs, submocoeffs)
                self.aoresults[spin][i] = Numeric.divide(tempcoeffs, scale)

                step += 1

        if self.progress:
            self.progress.update(nstep, "Done")

        retval = super(CSPA, self).partition(indices)

        if not retval:
            self.logger.error("Error in partitioning results")
            return False

#create array for mulliken charges
        self.logger.info("Creating fragcharges: array[1]")
        size = len(self.fragresults[0][0])
        self.fragcharges = Numeric.zeros([size], "f")
        
        for spin in range(len(self.fragresults)):

            for i in range(self.parser.homos[spin] + 1):

                temp = Numeric.reshape(self.fragresults[spin][i], (size,))
                self.fragcharges = Numeric.add(self.fragcharges, temp)
        
        if not unrestricted:
            self.fragcharges = Numeric.multiply(self.fragcharges, 2)

        return True

if __name__ == "__main__":
    import doctest, cspa
    doctest.testmod(cspa, verbose=False)
