import Numeric
import copy     # Needed for copying objects

from PyQuante.CGBF import CGBF

from cclib.bridge import makepyquante
from cclib.parser.utils import convertor

from pyvtk import *
from pyvtk.DataSetAttr import *

class Volume(object):
    """Represent a volume in space.

    Attributes:
       origin -- the bottom left hand corner of the volume :-)
       data -- a Numeric array of values for each point in the volume
       spacing -- a tuple of (dx, dy, dz)

    """
    def __init__(self, origin, data, spacing):
        self.origin = origin
        self.data = data
        self.spacing = spacing

    def corner(self):
        topcorner = []
        for i in range(3):
            topcorner.append(self.origin[i] + 
                             self.spacing[i] * (self.data.shape[i]-1))
        return tuple(topcorner)
    
    def __str__(self):
        """Return a string representation."""
        return "Volume %s to %s (density: %s)" % (self.origin, self.corner(),
                                                  self.spacing)

    def write(self, filename, format="VTK"):
        """Write the volume to file."""

        if format.upper() not in ["VTK"]:
            raise "Format must be VTK"

        ranges = (Numeric.arange(self.data.shape[2]),
                  Numeric.arange(self.data.shape[1]),
                  Numeric.arange(self.data.shape[0]))
        v = VtkData(RectilinearGrid(*ranges), "Test",
                    PointData(Scalars(self.data.flat)))
        outputfile = open(filename, "w")
        outputfile.write(v.to_string())
        outputfile.close()

def wavefunction(coords, mocoeffs, bfs, volume):
    """Calculate the magnitude of the wavefunction at every point in a volume.
    
    Attributes:
        coords -- the coordinates of the atoms
        mocoeffs -- mocoeffs for one eigenvalue
        bfs -- typically gbasis (although eventually STOs too)
        volume -- a template Volume object (will not be altered)
    """
    mymol = makepyquante(coords, [0 for x in coords])

    sym2powerlist = {
        'S' : [(0,0,0)],
        'P' : [(1,0,0),(0,1,0),(0,0,1)],
        'D' : [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(0,1,1),(1,0,1)],
        'F' : [(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),
               (0,3,0),(0,2,1),(0,1,2), (0,0,3)]
        }

    bfs = []
    for i,atom in enumerate(mymol):
        bs = a.gbasis[i]
        for sym,prims in bs:
            for power in sym2powerlist[sym]:
                bf = CGBF(atom.pos(),power)
                for expnt,coef in prims:
                    bf.add_primitive(expnt,coef)
                bf.normalize()
                bfs.append(bf)

    wavefn = copy.copy(volume)
    wavefn.data = Numeric.zeros( wavefn.data.shape, "d")

    conversion = convertor(1,"bohr","Angstrom")
    x = Numeric.arange(wavefn.origin[0], wavefn.corner()[0]+wavefn.spacing[0], wavefn.spacing[0]) / conversion
    y = Numeric.arange(wavefn.origin[1], wavefn.corner()[1]+wavefn.spacing[1], wavefn.spacing[1]) / conversion
    z = Numeric.arange(wavefn.origin[2], wavefn.corner()[2]+wavefn.spacing[2], wavefn.spacing[2]) / conversion

    for bs in range(len(bfs)):
        data = Numeric.zeros( wavefn.data.shape, "d")
        for i in range(wavefn.data.shape[0]):
            for j in range(wavefn.data.shape[1]):
                for k in range(wavefn.data.shape[2]):
                    data[i, j, k] = bfs[bs].amp(x[i], y[j], z[k])
        data = data * mocoeffs[bs]
        wavefn.data += data
    
    return wavefn

if __name__=="__main__":
    import psyco
    psyco.full() # Down from 3.5 to 1.5

    from cclib.parser import guesstype
    import logging
    a = guesstype("../../../data/Gaussian/basicGaussian03/dvb_sp_basis.log")
    a.logger.setLevel(logging.ERROR)
    a.parse()
    
    b = guesstype("../../../data/Gaussian/basicGaussian03/dvb_sp.out")
    b.logger.setLevel(logging.ERROR)
    b.parse()

    vol = Volume( (-2.5,-5,-1.5), Numeric.zeros((20,40,12)), spacing=(0.25,0.25,0.25) )
    homowavefn = wavefunction(b.atomcoords[0], b.mocoeffs[0,b.homos[0]], a.gbasis, vol)
    homowavefn.write("cubefile.vtk")
    
