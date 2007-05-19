"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

import random # For sometimes running the progress updater

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy
    numpy.inner = numpy.innerproduct

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
        nbasis = self.parser.nbasis

        #determine number of steps, and whether process involves beta orbitals
        self.aoresults = []
        alpha = len(self.parser.mocoeffs[0])
        self.aoresults.append(numpy.zeros([alpha, nbasis], "d"))
        nstep = alpha
        if unrestricted:
            beta = len(self.parser.mocoeffs[1])
            self.aoresults.append(numpy.zeros([beta, nbasis], "d"))
            nstep += beta

        #intialize progress if available
        if self.progress:
            self.progress.initialize(nstep)

        step = 0
        for spin in range(len(self.parser.mocoeffs)):

            for i in range(len(self.parser.mocoeffs[spin])):

                if self.progress and random.random() < fupdate:
                    self.progress.update(step, "C^2 Population Analysis")

                submocoeffs = self.parser.mocoeffs[spin][i]
                scale = numpy.inner(submocoeffs, submocoeffs)
                tempcoeffs = numpy.multiply(submocoeffs, submocoeffs)
                tempvec = tempcoeffs/scale
                self.aoresults[spin][i] = numpy.divide(tempcoeffs, scale).astype("d")

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
        self.fragcharges = numpy.zeros([size], "d")
        
        for spin in range(len(self.fragresults)):

            for i in range(self.parser.homos[spin] + 1):

                temp = numpy.reshape(self.fragresults[spin][i], (size,))
                self.fragcharges = numpy.add(self.fragcharges, temp)
        
        if not unrestricted:
            self.fragcharges = numpy.multiply(self.fragcharges, 2)

        return True

if __name__ == "__main__":
    import doctest, cspa
    doctest.testmod(cspa, verbose=False)
