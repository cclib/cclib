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
import re,time
import Numeric
import random # For sometimes running the progress updater
from calculationmethod import Method

class CSPA(Method):
    """The C-squared population analysis"""
    def __init__(self,*args):

        # Call the __init__ method of the superclass
        super(CSPA, self).__init__(logname="Density",*args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Density matrix of" % (self.parser)

    def __repr__(self):
        """Return a representation of the object."""
        return 'CSPA("%s")' % (self.parser)
    
    def calculate(self,fupdate=0.05):
        """Perform a C-squared population analysis given the results of a parser"""
    
        if not self.parser.parsed:
            self.parser.parse()

        #do we have the needed info in the parser?
        if not hasattr(self.parser,"mocoeffs") \
          and not hasattr(self.parser,"nbasis"):
            self.logger.error("Missing mocoeffs or nbasis")
            return False #let the caller of function know we didn't finish

        self.logger.info("Creating attribute aoresults: array[3]")
        unrestricted=(len(self.parser.mocoeffs)==2)
        nmocoeffs=len(self.parser.mocoeffs[0])
        nbasis=self.parser.nbasis

        #determine number of steps, and whether process involves beta orbitals
        nstep=nmocoeffs
        if unrestricted:
            self.aoresults=Numeric.zeros([2,nmocoeffs,nbasis],"f")
            nstep+=nmocoeffs
        else:
            self.aoresults=Numeric.zeros([1,nmocoeffs,nbasis],"f")

        #intialize progress if available
        if self.progress:
            self.progress.initialize(nstep)

        step=0
        for spin in range(len(self.parser.mocoeffs)):

            for i in range(nmocoeffs):

                if self.progress and random.random()<fupdate:
                    self.progress.update(step,"C^2 Population Analysis")

                submocoeffs=self.parser.mocoeffs[spin][i]
                scale=Numeric.innerproduct(submocoeffs,submocoeffs)
                tempcoeffs=Numeric.multiply(submocoeffs,submocoeffs)
                self.aoresults[spin][i]=Numeric.divide(tempcoeffs,scale)

                step+=1

        if self.progress:
            self.progress.update(nstep,"Done")

if __name__=="__main__":
    import doctest,g03parser
    doctest.testmod(g03parser,verbose=False)
