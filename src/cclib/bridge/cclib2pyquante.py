"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

from PyQuante.Molecule import Molecule

def makepyquante(atomcoords, atomnos):
    """Create a PyQuante Molecule.

    >>> import numpy
    >>> from PyQuante.hartree_fock import hf
    >>> atomnos = numpy.array([1,8,1],"i")
    >>> a = numpy.array([[-1,1,0],[0,0,0],[1,1,0]],"f")
    >>> pyqmol = pyquante(a,atomnos)
    >>> en,orbe,orbs = hf(pyqmol)
    >>> print en
    -73.8001234204
    """
# The only thing missing is charge, but this is also missing
# from cclib...things to do

    return Molecule("notitle", zip(atomnos, atomcoords), units = "Angstrom")

if __name__ == "__main__":
    import doctest, cclib2pyquante
    doctest.testmod(cclib2pyquante)
