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
from density import Density

class MBO(Density):
    """Calculate the density matrix"""
    def __init__(self,*args):

        # Call the __init__ method of the superclass
        super(MBO, self).__init__(logname="MBO",*args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Mayer's bond order of" % (self.parser)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Mayer\'s bond order("%s")' % (self.parser)
    
    def calculate(self,indices=None,fupdate=0.05):
        """Calculate Mayer's bond orders given the results of a parser"""
    
        if not self.parser.parsed:
            self.parser.parse()

        super(MBO,self).calculate(fupdate)

        #do we have the needed info in the parser?
        if not ( hasattr(self.parser,"aooverlaps") or hasattr(self.parser,"fooverlaps")):
            self.logger.error("Missing overlap matrix")
            return False #let the caller of function know we didn't finish

        if not indices:
#build list of groups of orbitals in each atom for atomresults
            if hasattr(self.parser,"aonames"):
                names=self.parser.aonames
            elif hasattr(self.parser,"foonames"):
                names=self.parser.fonames

            atoms=[]
            indices=[]

            name=names[0].split('_')[0]
            atoms.append(name)
            indices.append([0])

            for i in range(1,len(names)):
                name=names[i].split('_')[0]
                try:
                    index=atoms.index(name)
                except ValueError: #not found in atom list
                    atoms.append(name)
                    indices.append([i])
                else:
                    indices[index].append(i)
#done building list

        self.logger.info("Creating attribute fragresults: array[3]")
        size=len(indices)
        unrestricted=(len(self.parser.mocoeffs)==2)

        #determine number of steps, and whether process involves beta orbitals
        PS=[]
        PS.append(Numeric.matrixmultiply(self.density[0],self.parser.aooverlaps))
        nstep=size**2 #approximately quadratic in size
        if unrestricted:
            self.fragresults=Numeric.zeros([2,size,size],"f")
            PS.append(Numeric.matrixmultiply(self.density[1],self.parser.aooverlaps))
        else:
            self.fragresults=Numeric.zeros([1,size,size],"f")

        #intialize progress if available
        if self.progress:
            self.progress.initialize(nstep)

        step=0
        for i in range(len(indices)):

            if self.progress and random.random() < fupdate:
                self.progress.update(step,"Mayer's Bond Order")

            for j in range(i+1,len(indices)):

                tempsumA=0
                tempsumB=0
                
                for a in indices[i]:

                    for b in indices[j]:

                        tempsumA+=2*PS[0][a][b]*PS[0][b][a]
                        if unrestricted:
                            tempsumB+=2*PS[1][a][b]*PS[1][b][a]

                self.fragresults[0][i,j]=tempsumA
                self.fragresults[0][j,i]=tempsumA

                if unrestricted:
                    self.fragresults[1][i,j]=tempsumB
                    self.fragresults[1][j,i]=tempsumB

        if self.progress:
            self.progress.update(nstep,"Done")

        return True

if __name__=="__main__":
    import doctest,g03parser
    doctest.testmod(g03parser,verbose=False)
