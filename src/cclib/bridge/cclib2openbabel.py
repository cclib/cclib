"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision$"

import openbabel as ob

def makeopenbabel(atomcoords, atomnos):
    """Create an Open Babel molecule.

    >>> import Numeric, openbabel
    >>> atomnos = Numeric.array([1,8,1],"i")
    >>> coords = Numeric.array([[-1.,1.,0.],[0.,0.,0.],[1.,1.,0.]])
    >>> obmol = makeopenbabel(coords, atomnos)
    >>> obconversion = openbabel.OBConversion()
    >>> formatok = obconversion.SetOutFormat("inchi")
    >>> print obconversion.WriteString(obmol).strip()
    InChI=1/H2O/h1H2
    """
# The only thing missing is charge, but this is also missing
# from cclib...things to do
    obmol = ob.OBMol()
    for i in range(len(atomnos)):
        # Note that list(atomcoords[i]) is not equivalent!!!
        coords = atomcoords[i].tolist()
        atomno = atomnos[i]
        obatom = ob.OBAtom()
        obatom.SetAtomicNum(atomno)
        obatom.SetVector(*coords)
        obmol.AddAtom(obatom)
    obmol.ConnectTheDots()
    obmol.PerceiveBondOrders()
    return obmol

if __name__ == "__main__":
    import doctest
    doctest.testmod()
