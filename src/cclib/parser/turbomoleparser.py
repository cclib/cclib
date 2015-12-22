# -*- coding: utf-8 -*-
#
# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Parser for Turbomole output files"""


from __future__ import print_function
import re

import numpy

from . import logfileparser
from . import utils


class Turbomole(logfileparser.Logfile):
    """A Turbomole output file"""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(Turbomole, self).__init__(logname="Turbomole", *args, **kwargs)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Turbomole log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Turbomole("%s")' % (self.filename)

    def normalisesym(self, label):
        """Normalise the symmetries used by Turbomole."""
        return ans

    def before_parsing(self):
        self.geoopt = False # Is this a GeoOpt? Needed for SCF targets/values.

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        if line[14:32] == "atomic coordinates":
            atomcoords = []
            atomnos = []

            line = inputfile.next()
           
            while len(line) > 2:
                split_line = line.split()
                atsym = split_line[3].capitalize()
                atomnos.append(self.table.number[atsym])
                atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom")
                                    for x in split_line[0:3]])
                line = inputfile.next()

            if not hasattr(self,"atomcoords"):
                self.atomcoords = []

            self.atomcoords.append(atomcoords)
            self.atomnos = numpy.array(atomnos, dtype=int)
            self.natom = len(self.atomnos)

        if line[26:49] == "a o f o r c e - program":
            self.vibirs = []
            self.vibfreqs = []
            self.vibsyms = []
            self.vibdisps = []

        if line[12:26] == "ATOMIC WEIGHTS":
        # begin parsing atomic weights
            self.vibmasses=[]
            line=inputfile.next() # lines =======
            line=inputfile.next() # notes
            line=inputfile.next() # start reading
            split_line=line.split()
            while(len(split_line) > 0):
                self.vibmasses.append(float(split_line[2]))
                line=inputfile.next()
                split_line=line.split()

        # Block with displacements should start with this.
        #        mode               1        2        3        4        5        6        
        #                                                                                 
        #      frequency           0.00     0.00     0.00     0.00     0.00     0.00      
        #                                                                                 
        #      symmetry                                                                   
        #                                                                                 
        #         IR                -        -        -        -        -        -        
        # |dDIP/dQ|   (a.u.)     0.0000   0.0000   0.0000   0.0000   0.0000   0.0000      
        # intensity (km/mol)       0.00     0.00     0.00     0.00     0.00     0.00      
        # intensity (  %   )       0.00     0.00     0.00     0.00     0.00     0.00      
        #                                                                                 
        #        RAMAN              -        -        -        -        -        -        
        #                                                                                 
        #   1   c           x   0.03470 -0.00288  0.16407  0.04551 -0.16878 -0.02291      
        #                   y  -0.04449 -0.00240  0.12840  0.04538  0.15081 -0.04529      
        #                   z   0.17586 -0.10961 -0.01476  0.05942  0.05105 -0.00434      
        #   2   c           x   0.02365 -0.00250  0.14854  0.06901 -0.15407 -0.00293      
        #                   y  -0.01436 -0.00342  0.17074 -0.01869  0.11069 -0.09976      
        #                   z   0.24320  0.11146 -0.03720  0.11432  0.06397 -0.05581      
        # ...
        #  20   h           x   0.02122 -0.00242  0.14513  0.07418 -0.15083  0.00146      
        #                   y  -0.09180 -0.00081  0.06191  0.14599  0.21381  0.04026      
        #                   z   0.00000 -0.44974 -0.01306  0.00000  0.00000  0.00374      
        #                                                                                 
        # reduced mass(g/mol)     4.932    3.431    5.953    5.397    5.971    5.078      
        #                                                                                 
        #                                                                                 
        #        mode               7        8        9       10       11       12        
        if line[5:14] == "frequency":
            split_line=line.replace("i","-").split()
            freqs = [self.float(f) for f in split_line[1:]]
            self.vibfreqs.extend(freqs)
            line=inputfile.next() # blank
            line=inputfile.next() # symmetry
            syms=line.split()
            self.vibsyms.extend(syms[1:])
            line=inputfile.next() # blank
            line=inputfile.next() # IR (YES/NO)
            line=inputfile.next() # |dDIP/dQ|   (a.u.)
            line=inputfile.next() # intensity (km/mol)
            split_line=line.split()
            irs = [self.float(f) for f in split_line[2:]]
            self.vibirs.extend(irs)
            line=inputfile.next() # intensity (  %   )
            line=inputfile.next() # blank
            line=inputfile.next() # RAMAN (YES/NO)
            line=inputfile.next() # blank
            # now read in displacements
            disps = [] # temporary 
            for n in range(3*self.natom):
                line = next(inputfile)
                numbers = [float(s) for s in line[20:].split()]
                atomindex = int(n/3)
                numbermodes = len(numbers)
                if not disps: # i.e. if n == 0 initialize disps
                    for mode in range(numbermodes):
                        # For each mode, make list of list [atom][coord_index]
                        disps.append([[] for x in range(0,self.natom)]) 
                for mode in range(numbermodes): 
                    disps[mode][atomindex].append(numbers[mode])
            self.vibdisps.extend(disps)

    def after_parsing(self):

        # delete all frequencies that correspond to translations or rotations
        if hasattr(self,"vibfreqs"):
            n_modes = 3*self.natom-6
            self.vibfreqs = self.vibfreqs[-n_modes:]
            self.vibdisps = self.vibdisps[-n_modes:]
            self.vibirs = self.vibirs[-n_modes:]
            self.vibsyms = self.vibsyms[-n_modes:]

#EOF
