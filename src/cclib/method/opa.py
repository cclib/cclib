"""
cclib is a parser for computational chemistry log files.

See http://cclib.sf.net for more information.

Copyright (C) 2006 Noel O'Boyle and Adam Tenderholt

 This program is free software; you can redistribute and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY, without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

Contributions (monetary as well as code :-) are encouraged.
"""

__revision__ = "$Revision$"

import Numeric
import random # For sometimes running the progress updater
from calculationmethod import Method

class OPA(Method):
    """The overlap population analysis"""
    def __init__(self, *args):

        # Call the __init__ method of the superclass
        super(OPA, self).__init__(logname="OPA", *args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "OPA of" % (self.parser)

    def __repr__(self):
        """Return a representation of the object."""
        return 'OPA("%s")' % (self.parser)
    
    def calculate(self, indices=None, fupdate=0.05):
        """Perform an overlap population analysis given the results of a parser"""
    
        if not self.parser.parsed:
            self.parser.parse()

        #do we have the needed info in the parser?
        if not hasattr(self.parser, "mocoeffs") \
          and not ( hasattr(self.parser, "aooverlaps") \
                    or hasattr(self.parser, "fooverlaps") ) \
          and not hasattr(self.parser, "nbasis"):
            self.logger.error("Missing mocoeffs, aooverlaps/fooverlaps or nbasis")
            return False #let the caller of function know we didn't finish

        unrestricted = (len(self.parser.mocoeffs) == 2)
        nmocoeffs = len(self.parser.mocoeffs[0])
        nbasis = self.parser.nbasis

        if not indices:
#build list of groups of orbitals in each atom for atomresults
            if hasattr(self.parser, "aonames"):
                names = self.parser.aonames
            elif hasattr(self.parser, "foonames"):
                names = self.parser.fonames

            atoms = []
            indices = []

            name = names[0].split('_')[0]
            atoms.append(name)
            indices.append([0])

            for i in range(1, len(names)):
                name = names[i].split('_')[0]
                try:
                    index = atoms.index(name)
                except ValueError: #not found in atom list
                    atoms.append(name)
                    indices.append([i])
                else:
                    indices[index].append(i)

        #determine number of steps, and whether process involves beta orbitals
        nstep = nmocoeffs
        nfrag = len(indices) #nfrag
        self.logger.info("Creating attribute results: array[4]")
        if unrestricted:
            self.results = Numeric.ones([2, nmocoeffs, nfrag, nfrag], "f")
            nstep += nmocoeffs
        else:
            self.results=Numeric.ones([1, nmocoeffs, nfrag, nfrag], "f")

        #intialize progress if available
        if self.progress:
            self.progress.initialize(nstep)

        step = 0
        for spin in range(len(self.parser.mocoeffs)):

            for i in range(nmocoeffs):

                if self.progress and random.random() < fupdate:
                    self.progress.update(step, "Overlap Population Analysis")

                # OP_{AB,i} = \sum_{a in A} \sum_{b in B} 2 c_{ai} c_{bi} S_{ab}
                #       = \sum_{a in A} c_{ai} \sum_{b in B} c_{bi} S_{ab}
                #       = \sum_{a in A} c_{ai} \sub_{b in B} C_{Bi}
                #          where C_{Bi} = C(i) times S(a)

                ci = self.parser.mocoeffs[spin][i]

                if hasattr(self.parser, "aooverlaps"):
                    temp = Numeric.multiply(ci, self.parser.aooverlaps)
                elif hasattr(self.parser,"fooverlaps"):
                    temp = Numeric.multiply(ci, self.parser.fooverlaps)

                for A in range(len(indices)-1):

                    for B in range(A+1, len(indices)):

                        tempb = 0
                        for b in indices[B]:
                            tempb += temp[:, b]

                        tempresult = 0
                        for a in indices[A]:
                            tempresult += ci[a] * tempb[a]

                        self.results[spin][i][A][B] = 2 * tempresult
                        self.results[spin][i][B][A] = 2 * tempresult

                step += 1

        if self.progress:
            self.progress.update(nstep, "Done")

        return True

if __name__ == "__main__":
    import doctest, opa
    doctest.testmod(opa, verbose=False)
