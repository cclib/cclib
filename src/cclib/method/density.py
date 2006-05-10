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

class Density(Method):
    """Calculate the density matrix"""
    def __init__(self,*args):

        # Call the __init__ method of the superclass
        super(Density, self).__init__(logname="Density",*args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Density matrix of" % (self.parser)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Density matrix("%s")' % (self.parser)
    
    def calculate(self,fupdate=0.05,cupdate=0.002):
        """Calculate the density matrix given the results of a parser"""
    
        if not self.parser.parsed:
            self.parser.parse()

        #do we have the needed info in the parser?
        if not hasattr(self.parser,"mocoeffs") \
          and not hasattr(self.parser,"aooverlaps") \
          and not hasattr(self.parser,"nbasis") \
          and not hasattr(self.parser,"homos"):
            self.logger.error("Missing mocoeffs or aooverlaps")
            return False #let the caller of function know we didn't finish

        self.logger.info("Creating attribute density: array[3]")
        size=self.parser.nbasis
        unrestricted=(len(self.parser.mocoeffs)==2)

        if unrestricted:
            self.density=Numeric.zeros([2,size,size],"f")
        else:
            self.density=Numeric.zeros([1,size,size],"f")

        for spin in range(len(self.parser.mocoeffs)):

            for i in range(self.parser.homos[spin]+1):
                col=Numeric.reshape(self.parser.mocoeffs[spin][i],(size,1))
                colt=Numeric.reshape(col,(1,size))

                tempdensity=Numeric.matrixmultiply(col,colt)
                self.density[spin]=Numeric.add(self.density[spin],tempdensity)

        if not unrestricted: #multiply by two to account for second electron
            self.density[0]=Numeric.add(self.density[0],self.density[0])

if __name__=="__main__":
    import doctest,g03parser
    doctest.testmod(g03parser,verbose=False)
