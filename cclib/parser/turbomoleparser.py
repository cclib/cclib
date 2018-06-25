# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Turbomole output files"""

from __future__ import print_function

import re

import numpy

from cclib.parser import logfileparser
from cclib.parser import utils


class Turbomole(logfileparser.Logfile):
    """A Turbomole log file"""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(Turbomole, self).__init__(logname="Turbomole", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "Turbomole output file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Turbomole("%s")' % (self.filename)

    # This still needs to be implemented.
    def normalisesym(self, label):
        """Normalise the symmetries used by Turbomole."""

    def before_parsing(self):
        self.geoopt = False # Is this a GeoOpt? Needed for SCF targets/values.

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        ## This information is in the control file.
        #   $rundimensions
        #   dim(fock,dens)=1860
        #   natoms=20
        #   nshell=40
        #   nbf(CAO)=60
        #   nbf(AO)=60
        #   dim(trafo[SAO<-->AO/CAO])=60
        #   rhfshells=1
        if line[3:10]=="natoms=":
            self.natom=int(line[10:])

        if line[3:11] == "nbf(AO)=":
            nmo = int(line.split('=')[1])
            self.set_attribute('nbasis', nmo)
            self.set_attribute('nmo', nmo)

        ## Atomic coordinated in job.last:
        #              +--------------------------------------------------+
        #              | Atomic coordinate, charge and isotop information |
        #              +--------------------------------------------------+
        #
        #
        #              atomic coordinates              atom shells charge pseudo isotop
        #    -2.69176330   -0.00007129   -0.44712612    c      3    6.000    0     0
        #    -1.69851645   -0.00007332    2.06488947    c      3    6.000    0     0
        #     0.92683848   -0.00007460    2.49592179    c      3    6.000    0     0
        #     2.69176331   -0.00007127    0.44712612    c      3    6.000    0     0
        #     1.69851645   -0.00007331   -2.06488947    c      3    6.000    0     0
        #...
        #    -7.04373606    0.00092244    2.74543891    h      1    1.000    0     0
        #    -9.36352819    0.00017229    0.07445322    h      1    1.000    0     0
        #    -0.92683849   -0.00007461   -2.49592179    c      3    6.000    0     0
        #    -1.65164853   -0.00009927   -4.45456858    h      1    1.000    0     0
        if 'Atomic coordinate, charge and isotop information' in line:
            self.skip_lines(inputfile,['border','b','b'])
            line = next(inputfile)
            if 'atomic coordinates' in line:
                if not hasattr(self, 'atomnos'):
                    self.atomnos = []

                if not hasattr(self, 'atomcoords'):
                    self.atomcoords = []

                atomcoords = []
                atomnos = []
                line = next(inputfile)
               
                while len(line) > 2:
                    atomnos.append(line.split()[3].upper())
                    atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom") for x in line.split()[0:3]])
                    line = next(inputfile)

                self.atomcoords.append(atomcoords)

                if self.atomnos == []:
                    self.atomnos = [utils.PeriodicTable().number[atomelement] for atomelement in atomnos]

        ## If we are unable to find coordinates in the job file, we will look for them
        ## in the coord file.
        #$coord
        #   -2.69176330280845     -0.00007129445712     -0.44712612093731      c
        #   -1.69851644615005     -0.00007332282157      2.06488947265450      c
        #    0.92683848474542     -0.00007460039817      2.49592179180606      c
        #    2.69176330586455     -0.00007127328394      0.44712611937145      c
        #    1.69851645303760     -0.00007330575699     -2.06488946767951      c
        #...
        #   -7.04373605585274      0.00092243787879      2.74543890990978      h
        #   -9.36352819434217      0.00017228707094      0.07445321997922      h
        #   -0.92683848797451     -0.00007460625018     -2.49592178685491      c
        #   -1.65164852640697     -0.00009927259342     -4.45456857763148      h
        #$redundant
        if line[0:6] == "$coord":
            if '$user' not in next(inputfile) and '$coordinate' not in line:
                if not hasattr(self, 'atomcoords'): 
                    self.atomcoords = []
                    if not hasattr(self, 'atomnos'):
                        self.atomnos = []

                    atomcoords = []
                    atomnos = []
                    line = next(inputfile)

                    while line[0] != "$":
                        atomnos.append(line.split()[3].upper())
                        atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom") for x in line.split()[0:3]])
                        line = next(inputfile)

                    self.atomcoords.append(atomcoords)
                    if self.atomnos == []:
                        self.atomnos = [utils.PeriodicTable().number[atomelement] for atomelement in atomnos]


    def after_parsing(self):
        pass