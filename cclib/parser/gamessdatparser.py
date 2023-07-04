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

from cclib.parser.utils import PeriodicTable

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
        self.pt = PeriodicTable()
        # To change: declared only for passing unit tests

    def after_parsing(self):
        if hasattr(self, "atomcoords"):
            len_coords = len(self.atomcoords)
            self.atomcoords = numpy.reshape(self.atomcoords, (1, len_coords, 3))


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

            self.atommasses = []
            self.atomcoords = []

            line = next(inputfile)

            while line.strip() != "$END":
                parts = line.split()
                if len(parts) > 2:

                    symbol      = parts[0]
                    mass        = float(parts[1])
                    coordinates = [float(coord) for coord in parts[2:]]

                    assert len(coordinates) == 3

                    self.atommasses.append(mass)
                    self.atomcoords.append(coordinates)
                    
                line = next(inputfile)


        # Extract energy

        # --- CLOSED SHELL ORBITALS --- GENERATED AT Mon Aug  5 13:05:47 2019
        # water                                                                           
        # E(RHF)=      -74.9643287920, E(NUC)=    8.8870072224,   13 ITERS

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

        # Extract vector information and populate mocoeffs.
        # Also extract nbasis from here.

        if line[1:5] == "$VEC":

            self.nbasis = None

            self.mocoeffs = []

            while "$END" not in line:
                
                line = next(inputfile)

                if '$END' in line: return

                moc_line = line.replace('-', ' -').replace('E -', 'E-').strip()
                mocoeff = [float(m) for m in moc_line.split()[1:]]
                atom_number = line.split()[0]
                line_number = moc_line.split()[1]

                if not self.mocoeffs:
                    self.mocoeffs.append(mocoeff)

                elif atom_number == str(len(self.mocoeffs)):
                    self.mocoeffs[-1].extend(mocoeff)

                elif len(mocoeff) > 0:
                    self.mocoeffs.append(mocoeff)
            
                if atom_number.isdigit(): 
                    self.nbasis = int(atom_number)
                    self.nmo = int(line_number) * len(mocoeff)
                


            self.mocoeffs = [ self.mocoeffs ]

        
        # Extracting MP2 Energy Value

        # MP2 NATURAL ORBITALS, E(MP2)=      -75.0022821133
        
        if "E(MP2)=" in line:
            self.mpenergies = float(line.split()[-1])

        #  POPULATION ANALYSIS
        # C            6.00435  -0.00435   5.99623   0.00377
        # C            6.00435  -0.00435   5.99623   0.00377
        # C            6.07658  -0.07658   6.04383  -0.04383
        # C            6.07658  -0.07658   6.04383  -0.04383
        # C            6.07696  -0.07696   6.04409  -0.04409
        # C            6.07696  -0.07696   6.04409  -0.04409
        # H            0.92217   0.07783   0.95790   0.04210
        # H            0.92217   0.07783   0.95790   0.04210
        # H            0.92084   0.07916   0.95636   0.04364
        # H            0.92084   0.07916   0.95636   0.04364
        # C            6.07622  -0.07622   6.03889  -0.03889
        # C            6.07622  -0.07622   6.03889  -0.03889
        # H            0.92346   0.07654   0.95907   0.04093
        # H            0.92346   0.07654   0.95907   0.04093
        # C            6.15454  -0.15454   6.09150  -0.09150
        # C            6.15454  -0.15454   6.09150  -0.09150
        # H            0.92403   0.07597   0.95707   0.04293
        # H            0.92403   0.07597   0.95707   0.04293
        # H            0.92086   0.07914   0.95507   0.04493
        # H            0.92086   0.07914   0.95507   0.04493
        #  MOMENTS AT POINT    1 X,Y,Z= -0.000000  0.000000  0.000000
        #  DIPOLE       0.000000  0.000000  0.000000
        # ----- TOP OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----

        # Extract Moments and Dipole

        if line[1:17] == 'MOMENTS AT POINT':
            self.moments = [[ float(moment) for moment in line.split()[-3:] ]]

            line = next(inputfile)

            self.moments.append([ float(dipole) for dipole in line.split()[-3:] ])


        # Extracting Gaussian

        # ----- TOP OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----
        # water                                                                           
        # GAUSSIAN              7 MOL ORBITALS     21 PRIMITIVES        3 NUCLEI

        if line[0:8] == "GAUSSIAN":

            parts = line.split()

            # self.nmo    = int(parts[1])
            self.natom  = int(parts[6])

            # Continue extracting
            
            #   O    1    (CENTRE  1)   0.00000000  0.00000000  0.00000000  CHARGE =  8.0
            #   H    2    (CENTRE  2)   1.87082873  0.00000000  0.00000000  CHARGE =  1.0
            #   H    3    (CENTRE  3)  -0.51567032  1.79835585  0.00000000  CHARGE =  1.0

            line = next(inputfile)

            self.atomnos = []
            self.atomcoords = []

            
            while '(CENTRE' in line:

                parts = line.split()

                symbol = parts[0]
                atomno = self.pt.number[symbol]

                # atom_number = int(parts[1])
                # centre_number = int(parts[3][:-1])
                coords = [ float(n) for n in parts[4:7] ]

                charge = float(parts[-1])

                self.atomnos.append(atomno)
                self.atomcoords.append(coords)

                line = next(inputfile)


        # Extracting Centre Assignments

        # CENTRE ASSIGNMENTS    1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  2  2  2  3  3
        # CENTRE ASSIGNMENTS    3
        
        if line[0:18] == "CENTRE ASSIGNMENTS":
            self.atombasis = []
            current_number = 1
            start_num, end_num = 0, 0
            num = 1

            while line[0:18] == 'CENTRE ASSIGNMENTS':
                numbers = [int(num) for num in line.split()[2:]]  # Adjusted the index to skip the first number
                for num in numbers:
                    if num != current_number:
                        self.atombasis.append(list(range(start_num, end_num)))
                        start_num = end_num
                        current_number = num
                    
                    end_num += 1

                line = next(inputfile)

            if start_num > 0:
                self.atombasis.append(list(range(start_num, end_num)))

        # Extracting Type Assignments

        # TYPE ASSIGNMENTS      1  1  1  1  1  1  2  2  2  3  3  3  4  4  4  1  1  1  1  1
        # TYPE ASSIGNMENTS      1

        self.metadata["type_assignments"] = []

        while line[0:16] == "TYPE ASSIGNMENTS":
            type_assignments = [ int(n) for n in line.split()[2:] ]
            self.metadata["type_assignments"].extend(type_assignments)
            line = next(inputfile)
