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

        while "E (RHF)" not in line:
            line = next(inputfile)

        # Extract E(RHF) value

        if line.startswith("E(RHF)="):
            rhf_value = float(line.split("=")[1].strip())
            self.scfenergies = [ rhf_value ]

        # Extract E(NUC) value 

        if "E(NUC)=" in line:
            nuc_value = float(line.split("=")[1].strip())
            self.metadata["E_NUC"] = nuc_value

        # Extract number of ITERS 

        if "ITERS" in line:
            iters_value = int(line.split()[0])
            self.metadata["ITERS"] = iters_value

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

        # Extract vector information 

        if line[0:5] == "$VEC":

            self.metadata["vectors"] = []

            while "$END" not in line:
                line = next(inputfile).strip()
                fields = line.split()
                vector_info = {
                    "atom_index": int(fields[0]),
                    "component": int(fields[1]),
                    "values": [float(value) for value in fields[2:]]
                }
                self.metadata["vectors"].append(vector_info)

        # Extracting Population Analysis

        #  POPULATION ANALYSIS
        # O            8.31989  -0.31989   8.22116  -0.22116
        # H            0.84006   0.15994   0.88942   0.11058
        # H            0.84006   0.15994   0.88942   0.11058

        if "POPULATION ANALYSIS" in line:

            self.metadata["population"] = []

            line = next(inputfile)
            
            while line[0] != ' ':
                fields = line.split()
                atom_info = {
                    "atom": fields[0],
                    "charge": float(fields[1]),
                    "spin": float(fields[2]),
                    "net_charge": float(fields[3]),
                    "net_spin": float(fields[4])
                }
                self.metadata["population"].append(atom_info)
                line = next(inputfile).strip()

        # Extracting Moments at Point

        # MOMENTS AT POINT    1 X,Y,Z=  0.075831  0.100631  0.000000
        # MP2 NATURAL ORBITALS, E(MP2)=      -75.0022821133

        # if "MOMENTS AT POINT" in line:
        #     coordinates_pattern = r"X,Y,Z=\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)"
        #     coordinates_match = re.match(coordinates_pattern, line)
            
        #     if coordinates_match:
        #         self.moments = [
        #             float(coordinates_match.group(1)),
        #             float(coordinates_match.group(2)),
        #             float(coordinates_match.group(3))
        #         ]

        # Extracting Dipole

        # DIPOLE       1.007144  1.336525  0.000000

        # if "DIPOLE" in line:
        #     dipole_pattern = r"DIPOLE\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)"
        #     dipole_match = re.match(dipole_pattern, line)
            
        #     if dipole_match:
        #         self.metadata["Dipole"] = {
        #             "X": float(dipole_match.group(1)),
        #             "Y": float(dipole_match.group(2)),
        #             "Z": float(dipole_match.group(3))
        #         }
        
        # Extracting MP2 Energy Value

        # MP2 NATURAL ORBITALS, E(MP2)=      -75.0022821133
        
        # if "E(MP2)=" in line:
        #     energy_pattern = r"E\(MP2\)=\s+([\d.-]+)"
        #     energy_match = re.search(energy_pattern, line)
            
        #     if energy_match:
        #         self.mpenergies = float(energy_match.group(1))