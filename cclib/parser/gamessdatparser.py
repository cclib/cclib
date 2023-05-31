# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
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
        super().__init__(logname="GAMESS .dat", *args, **kwargs)

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
        pass

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

        # TODO


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

        # Extract vectors

        #  $VEC   
        #  1  1 9.94202990E-01 2.59157151E-02 2.40311554E-03 3.18904296E-03 0.00000000E+00
        #  1  2-5.62726478E-03-5.62726567E-03
        #  2  1-2.34217935E-01 8.45881798E-01 7.04411127E-02 9.34785483E-02-0.00000000E+00
        #  2  2 1.56449309E-01 1.56449336E-01
        #  3  1-1.17687509E-08 6.36837896E-08 4.81820843E-01-3.63078078E-01 0.00000000E+00
        #  3  2 4.46376801E-01-4.46376807E-01
        #  4  1 1.00458159E-01-5.21395067E-01 4.65965488E-01 6.18357071E-01 0.00000000E+00
        #  4  2 2.89063958E-01 2.89063907E-01
        #  5  1-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00 1.00000000E+00
        #  5  2-0.00000000E+00-0.00000000E+00
        #  6  1-1.28350522E-01 8.32525679E-01 4.40905546E-01 5.85101037E-01-0.00000000E+00
        #  6  2-7.75800880E-01-7.75800504E-01
        #  7  1 3.90260255E-08-2.79856684E-07 7.79855187E-01-5.87663268E-01 0.00000000E+00
        #  7  2-8.08915389E-01 8.08915850E-01
        #  $END   

        # TODO

        # Extracting Moments at Point

        # MOMENTS AT POINT    1 X,Y,Z=  0.075831  0.100631  0.000000
        # MP2 NATURAL ORBITALS, E(MP2)=      -75.0022821133

        if line.startswith("MOMENTS AT POINT"):
            coordinates_pattern = r"X,Y,Z=\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)"
            coordinates_match = re.match(coordinates_pattern, line)
            
            if coordinates_match:
                self.metadata["Coordinates"] = {
                    "X": float(coordinates_match.group(1)),
                    "Y": float(coordinates_match.group(2)),
                    "Z": float(coordinates_match.group(3))
                }

        # Extracting Dipole

        # DIPOLE       1.007144  1.336525  0.000000

        if line.startswith("DIPOLE"):
            dipole_pattern = r"DIPOLE\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)"
            dipole_match = re.match(dipole_pattern, line)
            
            if dipole_match:
                self.metadata["Dipole"] = {
                    "X": float(dipole_match.group(1)),
                    "Y": float(dipole_match.group(2)),
                    "Z": float(dipole_match.group(3))
                }
        
        # Extracting MP2 Energy Value

        # MP2 NATURAL ORBITALS, E(MP2)=      -75.0022821133
        
        if "E(MP2)=" in line:
            energy_pattern = r"E\(MP2\)=\s+([\d.-]+)"
            energy_match = re.search(energy_pattern, line)
            
            if energy_match:
                self.metadata["Energy(MP2)"] = float(energy_match.group(1))
