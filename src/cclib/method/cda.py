"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 453 $"

import random # For sometimes running the progress updater

import numpy

from fragments import FragmentAnalysis


class CDA(FragmentAnalysis):
    """The Charge decomposition analysis"""
    def __init__(self, *args):

        # Call the __init__ method of the superclass
        super(FragmentAnalysis, self).__init__(logname="CDA", *args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "CDA of" % (self.parser)

    def __repr__(self):
        """Return a representation of the object."""
        return 'CDA("%s")' % (self.parser)
    
    def calculate(self, fragments, cupdate=0.05):
        """Perform a charge decomposition analysis."""
    
        retval = super(CDA, self).calculate(fragments, cupdate)
        if not retval:
            return False

        #at this point, there should be a mocoeffs and fooverlaps in analogy to a parser

        donations = []
        bdonations = []
        repulsions = []

        if len(self.mocoeffs) == 2:
            occs = 1
        else:
            occs = 2

#intialize progress if available
        nstep = self.parser.homos[0]
        if len(self.parser.homos) == 2:
            nstep += self.parser.homos[1]

        if self.progress:
            self.progress.initialize(nstep)

#begin method
        step = 0
        for spin in range(len(self.mocoeffs)):

            size = len(self.mocoeffs[spin])
            homo = self.parser.homos[spin]

            if len(fragments[0].homos) == 2:
                homoa = fragments[0].homos[spin]
            else:
                homoa = fragments[0].homos[0]

            if len(fragments[1].homos) == 2:
                homob = fragments[1].homos[spin]
            else:
                homob = fragments[1].homos[0]

            offset = fragments[0].nbasis

            self.logger.info("Creating donations, bdonations, and repulsions: array[]")
            donations.append(numpy.zeros(size, "d"))
            bdonations.append(numpy.zeros(size, "d"))
            repulsions.append(numpy.zeros(size, "d"))

            for i in range(self.parser.homos[spin] + 1):

#calculate donation for each MO
                for k in range(0, homoa + 1):
                    for n in range(offset + homob + 1, self.parser.nbasis):
                        donations[spin][i] += 2 * occs * self.mocoeffs[spin][i,k] \
                                                * self.mocoeffs[spin][i,n] * self.fooverlaps[k][n]

                for l in range(offset, offset + homob + 1):
                    for m in range(homoa + 1, offset):
                        bdonations[spin][i] += 2 * occs * self.mocoeffs[spin][i,l] \
                                                * self.mocoeffs[spin][i,m] * self.fooverlaps[l][m]


                for k in range(0, homoa + 1):
                    for m in range(offset, offset+homob + 1):
                        repulsions[spin][i] += 2 * occs * self.mocoeffs[spin][i,k] \
                                                * self.mocoeffs[spin][i, m] * self.fooverlaps[k][m]

                step += 1
                if self.progress and random.random() < cupdate:
                    self.progress.update(step, "Charge Decomposition Analysis...")

        if self.progress:
            self.progress.update(nstep, "Done.")

        self.donations = donations
        self.bdonations = bdonations
        self.repulsions = repulsions

        return True
