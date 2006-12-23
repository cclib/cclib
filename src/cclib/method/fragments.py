"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 238 $"

import Numeric
import LinearAlgebra
import random # For sometimes running the progress updater
from calculationmethod import *

class FragmentAnalysis(Method):
    """Convert a molecule's basis functions from atomic-based to fragment MO-based"""
    def __init__(self, parser, progress=None, loglevel=logging.INFO,
                 logname="FragmentAnalysis of"):

        # Call the __init__ method of the superclass
        super(FragmentAnalysis, self).__init__(parser, progress, loglevel, logname)
        self.parsed = False
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Fragment molecule basis of" % (self.parser)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Fragment molecular basis("%s")' % (self.parser)

    def calculate(self, fragments, cupdate=0.05):

        nFragBasis = 0
        nFragAlpha = 0
        nFragBeta = 0
        self.fonames = []

        unrestricted = ( len(self.parser.mocoeffs) == 2 )

        self.logger.info("Creating attribute fonames[]")

#collect basis info on the fragments
        for j in range(len(fragments)):
            nFragBasis += fragments[j].nbasis
            nFragAlpha += fragments[j].homos[0] + 1
            if unrestricted:
                nFragBeta += fragments[j].homos[0] + 1 #assume restricted fragments
            
            #assign fonames based on fragment name and MO number
            for i in range(fragments[j].nbasis):
                if hasattr(fragments[j],"name"):
                    self.fonames.append("%s_%i"%(fragments[j].name,i))
                else:
                    self.fonames.append("noname%i_%i"%(j,i))

        nBasis = self.parser.nbasis
        nAlpha = self.parser.homos[0] + 1
        if unrestricted:
            nBeta = self.parser.homos[1] + 1

        if nBasis != nFragBasis:
            self.logger.error("Basis functions don't match")
            return

        if nAlpha != nFragAlpha:
            self.logger.error("Alpha electrons don't match")
            return

        if unrestricted and nBeta != nFragBeta:
            self.logger.error("Beta electrons don't match")
            return

        blockMatrix = Numeric.zeros((nBasis,nBasis),"float")
        pos = 0

        #build up block-diagonal matrix from fragment mocoeffs
        #need to switch ordering from [mo,ao] to [ao,mo]
        for i in range(len(fragments)):
            size = fragments[i].nbasis
            blockMatrix[pos:pos+size,pos:pos+size] = Numeric.transpose(fragments[i].mocoeffs[0])
            pos += size
        
        #invert and mutliply to result in fragment MOs as basis
        iBlockMatrix = LinearAlgebra.inverse(blockMatrix) 
        self.mocoeffs = [Numeric.transpose(Numeric.matrixmultiply(iBlockMatrix, Numeric.transpose(self.parser.mocoeffs[0])))]
        if unrestricted:
            self.mocoeffs.append(Numeric.transpose(Numeric.matrixmultiply(iBlockMatrix, \
                                    Numeric.transpose(self.parser.mocoeffs[1]))))
        
        self.logger.info("Creating mocoeffs in new fragment MO basis: mocoeffs[]")

        self.parsed = True
        self.nbasis = nBasis
        self.homos = self.parser.homos
