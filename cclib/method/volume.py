# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculation methods related to volume based on cclib data."""

from __future__ import print_function
import copy

import numpy

from cclib.bridge import cclib2pyquante
from cclib.parser.utils import convertor
from cclib.parser.utils import find_package

""" In the dictionary sym2powerlist below, each element is a list that contain the combinations of
    powers that are applied to x, y, and z in the equation for the gaussian primitives --
    \psi (x, y, z) = x^a * y^b * z^c * exp(-\lambda * r^2)
"""
sym2powerlist = {
    "S": [(0, 0, 0)],
    "P": [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    "D": [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (0, 1, 1), (1, 0, 1)],
    "F": [
        (3, 0, 0),
        (2, 1, 0),
        (2, 0, 1),
        (1, 2, 0),
        (1, 1, 1),
        (1, 0, 2),
        (0, 3, 0),
        (0, 2, 1),
        (0, 1, 2),
        (0, 0, 3),
    ],
}


_found_pyquante = find_package("PyQuante")
if _found_pyquante:
    from PyQuante.CGBF import CGBF

    def getbfs(ccdata):
        pymol = cclib2pyquante.makepyquante(ccdata)

        bfs = []
        for i, atom in enumerate(pymol):
            bs = ccdata.gbasis[i]
            for sym, prims in bs:
                for power in sym2powerlist[sym]:
                    bf = CGBF(atom.pos(), power)
                    for expnt, coef in prims:
                        bf.add_primitive(expnt, coef)
                    bf.normalize()
                    bfs.append(bf)
        return bfs

    # Small wrapper PyQuante & pyquante2 function that evaluates basis function on a given point
    # Used in both `wavefunction` and `electrondensity`
    def pyamp(bfs, bs, x, y, z):
        return bfs[bs].amp(x, y, z)


_found_pyquante2 = find_package("pyquante2")
if _found_pyquante2:
    from pyquante2 import cgbf

    def getbfs(ccdata):
        pymol = cclib2pyquante.makepyquante(ccdata)

        bfs = []
        for i, atom in enumerate(pymol):  # `atom` is instance of pyquante2.geo.atom.atom class.
            basis = ccdata.gbasis[i]  # `basis` is basis coefficients stored in ccData.
            for sym, primitives in basis:
                for power in sym2powerlist[sym]:  # `sym` is S, P, D, F and is used as key here.
                    exponentlist = []
                    coefficientlist = []

                    for exponents, coefficients in primitives:
                        exponentlist.append(exponents)
                        coefficientlist.append(coefficients)

                    basisfunction = cgbf(
                        atom.atuple()[1:4], powers=power, exps=exponentlist, coefs=coefficientlist,
                    )
                    basisfunction.normalize()
                    bfs.append(basisfunction)

        return bfs

    # Small wrapper PyQuante & pyquante2 function that evaluates basis function on a given point
    # Used in both `wavefunction` and `electrondensity`
    def pyamp(bfs, bs, x, y, z):
        return bfs[bs](x, y, z)  # 1D numpy array with size 1 is returned from __call__ here.


_found_pyvtk = find_package("pyvtk")
if _found_pyvtk:
    from pyvtk import *
    from pyvtk.DataSetAttr import *


def _check_pyquante():
    if (not _found_pyquante) and (not _found_pyquante2):
        raise ImportError("You must install `pyquante2` or `PyQuante` to use this function.")


def _check_pyvtk(found_pyvtk):
    if not found_pyvtk:
        raise ImportError("You must install `pyvtk` to use this function.")


_check_pyquante()


class Volume(object):
    """Represent a volume in space.

    Required parameters:
       origin -- the bottom left hand corner of the volume
       topcorner -- the top right hand corner
       spacing -- the distance between the points in the cube

    Attributes:
       data -- a NumPy array of values for each point in the volume
               (set to zero at initialisation)
       numpts -- the numbers of points in the (x,y,z) directions
    """

    def __init__(self, origin, topcorner, spacing):

        self.origin = numpy.asarray(origin, dtype=float)
        self.topcorner = numpy.asarray(topcorner, dtype=float)
        self.spacing = numpy.asarray(spacing, dtype=float)
        self.numpts = []
        for i in range(3):
            self.numpts.append(int((self.topcorner[i] - self.origin[i]) / self.spacing[i] + 1))
        self.data = numpy.zeros(tuple(self.numpts), "d")

    def __str__(self):
        """Return a string representation."""
        return "Volume %s to %s (density: %s)" % (self.origin, self.topcorner, self.spacing,)

    def write(self, filename, fformat="Cube"):
        """Write the volume to a file."""

        fformat = fformat.lower()

        writers = {
            "vtk": self.writeasvtk,
            "cube": self.writeascube,
        }

        if fformat not in writers:
            raise RuntimeError("File format must be either VTK or Cube")

        writers[fformat](filename)

    def writeasvtk(self, filename):
        _check_pyvtk(_found_pyvtk)
        ranges = (
            numpy.arange(self.data.shape[2]),
            numpy.arange(self.data.shape[1]),
            numpy.arange(self.data.shape[0]),
        )
        v = VtkData(
            RectilinearGrid(*ranges),
            "Test",
            PointData(Scalars(self.data.ravel(), "from cclib", "default")),
        )
        v.tofile(filename)

    def integrate(self):
        boxvol = (
            self.spacing[0]
            * self.spacing[1]
            * self.spacing[2]
            * convertor(1, "Angstrom", "bohr") ** 3
        )
        return sum(self.data.ravel()) * boxvol

    def integrate_square(self):
        boxvol = (
            self.spacing[0]
            * self.spacing[1]
            * self.spacing[2]
            * convertor(1, "Angstrom", "bohr") ** 3
        )
        return sum(self.data.ravel() ** 2) * boxvol

    def writeascube(self, filename):
        # Remember that the units are bohr, not Angstroms
        def convert(x):
            return convertor(x, "Angstrom", "bohr")

        ans = []
        ans.append("Cube file generated by cclib")
        ans.append("")
        format = "%4d%12.6f%12.6f%12.6f"
        origin = [convert(x) for x in self.origin]
        ans.append(format % (0, origin[0], origin[1], origin[2]))
        ans.append(format % (self.data.shape[0], convert(self.spacing[0]), 0.0, 0.0))
        ans.append(format % (self.data.shape[1], 0.0, convert(self.spacing[1]), 0.0))
        ans.append(format % (self.data.shape[2], 0.0, 0.0, convert(self.spacing[2])))
        line = []
        for i in range(self.data.shape[0]):
            for j in range(self.data.shape[1]):
                for k in range(self.data.shape[2]):
                    line.append(scinotation(self.data[i, j, k]))
                    if len(line) == 6:
                        ans.append(" ".join(line))
                        line = []
                if line:
                    ans.append(" ".join(line))
                    line = []
        with open(filename, "w") as outputfile:
            outputfile.write("\n".join(ans))


def scinotation(num):
    """Write in scientific notation."""
    ans = "%10.5E" % num
    broken = ans.split("E")
    exponent = int(broken[1])
    if exponent < -99:
        return "  0.000E+00"
    if exponent < 0:
        sign = "-"
    else:
        sign = "+"
    return ("%sE%s%s" % (broken[0], sign, broken[1][-2:])).rjust(12)


def getGrid(vol):
    """Helper function that returns (x, y, z), each of which are numpy array of the values that
       correspond to grid points.
       
    Input:
       vol -- Volume object (will not be altered)
       """
    conversion = convertor(1, "bohr", "Angstrom")
    gridendpt = vol.topcorner + 0.5 * vol.spacing
    x = numpy.arange(vol.origin[0], gridendpt[0], vol.spacing[0]) / conversion
    y = numpy.arange(vol.origin[1], gridendpt[1], vol.spacing[1]) / conversion
    z = numpy.arange(vol.origin[2], gridendpt[2], vol.spacing[2]) / conversion
    
    return (x, y, z)


def wavefunction(ccdata, volume, mocoeffs):
    """Calculate the magnitude of the wavefunction at every point in a volume.
    
    Inputs:
        ccdata -- ccData object
        volume -- Volume object (will not be altered)
        mocoeffs -- molecular orbital to use for calculation; i.e. ccdata.mocoeffs[0][3]
    
    Output:
        Volume object with wavefunction at each grid point stored in data attribute
    """
    bfs = getbfs(ccdata)

    wavefn = copy.copy(volume)
    wavefn.data = numpy.zeros(wavefn.data.shape, "d")

    x, y, z = getGrid(wavefn)

    # PyQuante & pyquante2
    for bs in range(len(bfs)):
        data = numpy.zeros(wavefn.data.shape, "d")
        for i, xval in enumerate(x):
            for j, yval in enumerate(y):
                for k, zval in enumerate(z):
                    data[i, j, k] = pyamp(bfs, bs, xval, yval, zval)
        data *= mocoeffs[bs]
        wavefn.data += data

    return wavefn


def electrondensity_spin(ccdata, volume, mocoeffslist):
    """Calculate the magnitude of the electron density at every point in a volume for either up or down spin

    Inputs:
        ccdata -- ccData object
        volume -- Volume object (will not be altered)
        mocoeffslist -- list of molecular orbital to calculate electron density from;
                        i.e. [ccdata.mocoeffs[0][1:2]]
    
    Output:
        Volume object with wavefunction at each grid point stored in data attribute

    Attributes:
        coords -- the coordinates of the atoms
        mocoeffs -- mocoeffs for all of the occupied eigenvalues
        gbasis -- gbasis from a parser object
        volume -- a template Volume object (will not be altered)

    Note: mocoeffs is a list of NumPy arrays. The list will be of length 1.
    """
    assert len(mocoeffslist) == 1, "mocoeffslist input to the function should have length of 1."

    bfs = getbfs(ccdata)

    density = copy.copy(volume)
    density.data = numpy.zeros(density.data.shape, "d")

    x, y, z = getGrid(density)

    # For occupied orbitals
    # `mocoeff` and `gbasis` in ccdata object is ordered in a way `homos` can specify which orbital
    # is the highest lying occupied orbital in mocoeff and gbasis.
    for mocoeffs in mocoeffslist:
        for mocoeff in mocoeffs:
            wavefn = numpy.zeros(density.data.shape, "d")
            for bs in range(len(bfs)):
                data = numpy.zeros(density.data.shape, "d")
                for i, xval in enumerate(x):
                    for j, yval in enumerate(y):
                        tmp = []
                        for zval in z:
                            tmp.append(pyamp(bfs, bs, xval, yval, zval))
                        data[i, j, :] = tmp
                data *= mocoeff[bs]
                wavefn += data
            density.data += wavefn ** 2

    return density


def electrondensity(ccdata, volume, mocoeffslist):
    """Calculate the magnitude of the total electron density at every point in a volume.

    Inputs:
        ccdata -- ccData object
        volume -- Volume object (will not be altered)
        mocoeffslist -- list of molecular orbital to calculate electron density from;
                        i.e. [ccdata.mocoeffs[0][1:2]]
    
    Output:
        Volume object with wavefunction at each grid point stored in data attribute

    Attributes:
        coords -- the coordinates of the atoms
        mocoeffs -- mocoeffs for all of the occupied eigenvalues
        gbasis -- gbasis from a parser object
        volume -- a template Volume object (will not be altered)

    Note: mocoeffs is a list of NumPy arrays. The list will be of length 1
          for restricted calculations, and length 2 for unrestricted.
    """

    if len(mocoeffslist) == 2:
        return electrondensity_spin(ccdata, volume, [mocoeffslist[0]]) + electrondensity_spin(
            ccdata, volume, [mocoeffslist[1]]
        )
    else:
        edens = electrondensity_spin(ccdata, volume, [mocoeffslist[0]])
        edens.data *= 2
        return edens
