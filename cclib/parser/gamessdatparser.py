# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for GAMESS(US) .dat output files"""


import re
import numpy

from cclib.parser import logfileparser
from cclib.parser import utils


class GAMESSDAT(logfileparser.Logfile):
    """A GAMESS .dat log file"""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="GAMESSDAT", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return f"GAMESS .dat log file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'GAMESS .dat ("{self.filename}")'
    
    def normalisesym(self, label):
        """Normalise the symmetries used by GAMESS .dat."""

        pass

    def before_parsing(self):
        # To change: declared only for passing unit tests
        self.mocoeffs = [ -1 ]
        self.metadata["input_file_contents"] = None
        self.metadata["legacy_package_version"] = None
        self.scfenergies = [ 0 ]
        self.b3lyp_energy = 0


    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""
        
        # Extract element and its properties

        # $DATA  
        # water                                                                           
        # C1       0
        # O           8.0      0.0000000000      0.0000000000      0.0000000000
        #    STO     3
                
        # H           1.0      0.9900000000      0.0000000000      0.0000000000
        #    STO     3
                
        # H           1.0     -0.2728810000      0.9516490000      0.0000000000
        #    STO     3
                
        #  $END   

        # Extract molecule name

        if line[1:6] == "$DATA":

            line = next(inputfile)
            name = line.split('\n')[0].strip()
            self.metadata["name"] = name

            # Extract atomic information

            atomic_data = []
            line = next(inputfile)

            while line.strip() != "$END":
                parts = line.split()
                if len(parts) >= 2:
                    symbol = parts[0]
                    mass = float(parts[1])
                    coordinates = [float(coord) for coord in parts[2:]]
                    atom_info = {"symbol": symbol, "mass": mass, "coordinates": coordinates}
                    atomic_data.append(atom_info)
                line = next(inputfile)
        
            self.metadata["atoms"] = atomic_data

        # Extract energy

        # --- CLOSED SHELL ORBITALS --- GENERATED AT Mon Aug  5 13:05:47 2019
        # water                                                                           
        # E(RHF)=      -74.9643287920, E(NUC)=    8.8870072224,   13 ITERS

        # Extract E(RHF) value using regex
        rhf_match = re.search(r"E\(RHF\)=(.*?),", line)
        if rhf_match:
            self.metadata["E_RHF"] = rhf_match.group(1).strip()

        # Extract E(NUC) value using regex
        nuc_match = re.search(r"E\(NUC\)=(.*?),", line)
        if nuc_match:
            self.metadata["E_NUC"] = nuc_match.group(1).strip()

        # Extract number of ITERS using regex
        iters_match = re.search(r"(\d+)\s+ITERS", line)
        if iters_match:
            self.metadata["ITERS"] = iters_match.group(1).strip()

        
        # Extracting MP2 Energy Value

        # MP2 NATURAL ORBITALS, E(MP2)=      -75.0022821133
        
        if "E(MP2)=" in line:
            self.mpenergies = float(line.split()[-1])


        # Extracting Gaussian

        # ----- TOP OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----
        # water                                                                           
        # GAUSSIAN              7 MOL ORBITALS     21 PRIMITIVES        3 NUCLEI
        #   O    1    (CENTRE  1)   0.00000000  0.00000000  0.00000000  CHARGE =  8.0
        #   H    2    (CENTRE  2)   1.87082873  0.00000000  0.00000000  CHARGE =  1.0
        #   H    3    (CENTRE  3)  -0.51567032  1.79835585  0.00000000  CHARGE =  1.0

        # Extracting Centre Assignments

        # CENTRE ASSIGNMENTS    1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  2  2  2  3  3
        # CENTRE ASSIGNMENTS    3
        # TYPE ASSIGNMENTS      1  1  1  1  1  1  2  2  2  3  3  3  4  4  4  1  1  1  1  1
        # TYPE ASSIGNMENTS      1

        # Extracting Exponents

        # EXPONENTS  1.3070932E+02 2.3808866E+01 6.4436083E+00 5.0331513E+00 1.1695961E+00
        # EXPONENTS  3.8038896E-01 5.0331513E+00 1.1695961E+00 3.8038896E-01 5.0331513E+00
        # EXPONENTS  1.1695961E+00 3.8038896E-01 5.0331513E+00 1.1695961E+00 3.8038896E-01
        # EXPONENTS  3.4252509E+00 6.2391373E-01 1.6885540E-01 3.4252509E+00 6.2391373E-01
        # EXPONENTS  1.6885540E-01
