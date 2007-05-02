"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

import random # For sometimes running the progress updater
import logging

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

from calculationmethod import Method


class Density(Method):
    """Calculate the density matrix"""
    def __init__(self, parser, progress=None, loglevel=logging.INFO,
                 logname="Density"):

        # Call the __init__ method of the superclass
        super(Density, self).__init__(parser, progress, loglevel, logname)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Density matrix of" % (self.parser)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Density matrix("%s")' % (self.parser)
    
    def calculate(self, fupdate=0.05):
        """Calculate the density matrix given the results of a parser"""
    
        if not self.parser.parsed:
            self.parser.parse()

#do we have the needed info in the parser?
        if not hasattr(self.parser, "mocoeffs"): 
            self.logger.error("Parser missing mocoeffs")
            return False
        if not hasattr(self.parser,"nbasis"):
            self.logger.error("Parser missing nbasis")
            return False
        if not hasattr(self.parser,"homos"):
            self.logger.error("Parser missing homos")
            return False
#end attribute check

        self.logger.info("Creating attribute density: array[3]")
        size = self.parser.nbasis
        unrestricted = (len(self.parser.mocoeffs) == 2)

        #determine number of steps, and whether process involves beta orbitals
        nstep = self.parser.homos[0] + 1
        if unrestricted:
            self.density = numpy.zeros([2, size, size], "d")
            nstep += self.parser.homos[1] + 1
        else:
            self.density = numpy.zeros([1, size, size], "d")

        #intialize progress if available
        if self.progress:
            self.progress.initialize(nstep)

        step = 0
        for spin in range(len(self.parser.mocoeffs)):

            for i in range(self.parser.homos[spin] + 1):

                if self.progress and random.random() < fupdate:
                    self.progress.update(step, "Density Matrix")

                col = numpy.reshape(self.parser.mocoeffs[spin][i], (size, 1))
                colt = numpy.reshape(col, (1, size))

                tempdensity = numpy.dot(col, colt)
                self.density[spin] = numpy.add(self.density[spin],
                                                 tempdensity)

                step += 1

        if not unrestricted: #multiply by two to account for second electron
            self.density[0] = numpy.add(self.density[0], self.density[0])

        if self.progress:
            self.progress.update(nstep, "Done")

        return True #let caller know we finished density
