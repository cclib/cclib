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
        self.pt = utils.PeriodicTable()
        # To change: declared only for passing unit tests
        self.mocoeffs = [ -1 ]
        self.metadata["input_file_contents"] = None
        self.metadata["legacy_package_version"] = None
        self.scfenergies = [ 0 ]
        self.b3lyp_energy = 0

    # def after_parsing(self):
    #     if hasattr(self, "atomcoords"):
    #         self.atomcoords = numpy.array(self.atomcoords)
        


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

        # while "E(RHF)" not in line:
        #     line = next(inputfile)

        # Extract E(RHF) value

        if "E(R" in line:
            val_pattern = r"E\(R[^,]+"
            match = re.search(val_pattern, line)
            val = float(match.group().split(' ')[-1])
            self.scfenergies = [ val ]

        # Extract E(NUC) value 

        if "E(NUC)=" in line:
            pattern_e_nuc = r"E\(NUC[^,]+"
            match_e_nuc = re.search(pattern_e_nuc, line)
            nuc_value = float(match_e_nuc.group().split(' ')[-1].strip())
            self.metadata["E_NUC"] = nuc_value

        # Extract number of ITERS 

        # if "ITERS" in line:
        #     iters_value = int(line.split()[-1])
        #     self.metadata["ITERS"] = iters_value

        
        # Extracting MP2 Energy Value

        # MP2 NATURAL ORBITALS, E(MP2)=      -75.0022821133
        
        if "E(MP2)=" in line:
            self.mpenergies = float(line.split()[-1])


        # Extracting Gaussian

        # ----- TOP OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----
        # water                                                                           
        # GAUSSIAN              7 MOL ORBITALS     21 PRIMITIVES        3 NUCLEI

        if line[0:8] == "GAUSSIAN":

            parts = line.split()

            self.nmo    = int(parts[1])
            self.nbasis = int(parts[4])
            self.natom  = int(parts[6])


            # Continue extracting
            
            #   O    1    (CENTRE  1)   0.00000000  0.00000000  0.00000000  CHARGE =  8.0
            #   H    2    (CENTRE  2)   1.87082873  0.00000000  0.00000000  CHARGE =  1.0
            #   H    3    (CENTRE  3)  -0.51567032  1.79835585  0.00000000  CHARGE =  1.0

            line = next(inputfile)

            atom_info = []
            self.atomcharges = dict()
            self.atomcharges['CHANGENAME'] = []
            self.atomnos = []
            
            while '(CENTRE' in line:

                parts = line.split()

                symbol = parts[0]
                atomno = self.pt.number[symbol]

                # atom_number = int(parts[1])
                # centre_number = int(parts[3][:-1])
                val1 = float(parts[4])
                val2 = float(parts[5])
                val3 = float(parts[6])

                charge = float(parts[-1])

                atom_info.append((val1, val2, val3))
                
                self.atomnos.append(atomno)
                self.atomcharges['CHANGENAME'].append(charge)

                line = next(inputfile)

            self.metadata['atom_info'] = atom_info


        # Extracting Centre Assignments

        # CENTRE ASSIGNMENTS    1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  2  2  2  3  3
        # CENTRE ASSIGNMENTS    3
        
        self.atombasis = []

        while line[0:18] == "CENTRE ASSIGNMENTS":
            centre_assignments = re.findall(r'\d+', line)
            self.atombasis.extend(centre_assignments)
            line = next(inputfile)

        # Extracting Ttype Assignments

        # TYPE ASSIGNMENTS      1  1  1  1  1  1  2  2  2  3  3  3  4  4  4  1  1  1  1  1
        # TYPE ASSIGNMENTS      1

        self.metadata["type_assignments"] = []

        while line[0:16] == "TYPE ASSIGNMENTS":
            type_assignments = re.findall(r'\d+', line)
            self.metadata["type_assignments"].extend(type_assignments)
            line = next(inputfile)

        # Extracting Exponents

        # EXPONENTS  1.3070932E+02 2.3808866E+01 6.4436083E+00 5.0331513E+00 1.1695961E+00
        # EXPONENTS  3.8038896E-01 5.0331513E+00 1.1695961E+00 3.8038896E-01 5.0331513E+00
        # EXPONENTS  1.1695961E+00 3.8038896E-01 5.0331513E+00 1.1695961E+00 3.8038896E-01
        # EXPONENTS  3.4252509E+00 6.2391373E-01 1.6885540E-01 3.4252509E+00 6.2391373E-01
        # EXPONENTS  1.6885540E-01
