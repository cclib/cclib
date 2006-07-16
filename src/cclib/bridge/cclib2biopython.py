"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

from Bio.PDB.Atom import Atom
from cclib.parser.utils import PeriodicTable

def makebiopython(atomcoords, atomnos):
    """Create a list of BioPython Atoms.

    This creates a list of BioPython Atoms suitable for use
    by Bio.PDB.Superimposer, for example.

    >>> import Numeric
    >>> from Bio.PDB.Superimposer import Superimposer
    >>> atomnos = Numeric.array([1,8,1],"i")
    >>> a = Numeric.array([[-1,1,0],[0,0,0],[1,1,0]],"f")
    >>> b = Numeric.array([[1.1,2,0],[1,1,0],[2,1,0]],"f")
    >>> si = Superimposer()
    >>> si.set_atoms(biopython(a,atomnos),biopython(b,atomnos))
    >>> print si.rms
    0.29337859596
    """
    pt = PeriodicTable()
    bioatoms = []
    for coords, atomno in zip(atomcoords, atomnos):
        bioatoms.append(Atom(pt.element[atomno], coords, 0, 0, 0, 0, 0))
    return bioatoms

if __name__ == "__main__":
    import doctest, biopython
    doctest.testmod(biopython)
