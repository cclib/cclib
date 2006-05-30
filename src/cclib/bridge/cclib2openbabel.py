import pyopenbabel as pob
import openbabel as ob
from cclib.parser.utils import PeriodicTable

def makeopenbabel(atomcoords,atomnos):
    """Create a list of pyopenbabel molecules.

    >>> import Numeric
    >>> atomnos = Numeric.array([1,8,1],"i")
    >>> a = Numeric.array([[-1,1,0],[0,0,0],[1,1,0]],"f")
    >>> pyOBmol = pyopenbabel(a,atomnos)
    >>> print pyOBmol.write("inchi").strip()
    InChI=1/H2O/h1H2
    """
# The only thing missing is charge, but this is also missing
# from cclib...things to do
    obmol = ob.OBMol()
    for coords,atomno in zip(atomcoords,atomnos):
        obatom = ob.OBAtom()
        obatom.SetAtomicNum(atomno)
        obatom.SetVector(*coords)
        obmol.AddAtom(obatom)
    obmol.ConnectTheDots()
    obmol.PerceiveBondOrders()
    return pob.Molecule(obmol)
    
if __name__=="__main__":
    import doctest,cclib2openbabel
    doctest.testmod(cclib2openbabel)

