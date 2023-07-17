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
        super().__init__(logname="NBO", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return f"NBO log file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'NBO ("{self.filename}")'
    
    def normalisesym(self, label):
        """Normalise the symmetries used by NBO."""

        pass

    def before_parsing(self):

        pass

    def after_parsing(self):
        
        pass

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        ''' NATURAL POPULATIONS:  Natural atomic orbital occupancies

        NAO Atom No lang   Type(AO)    Occupancy      Energy
        -------------------------------------------------------
        1    O  1  s      Cor( 1s)     1.99998     -20.54367
        2    O  1  s      Val( 2s)     1.74925      -0.99433
        3    O  1  s      Ryd( 3s)     0.00135       1.41731
        4    O  1  px     Val( 2p)     1.43668      -0.30356
        5    O  1  px     Ryd( 3p)     0.00246       1.32660
        6    O  1  py     Val( 2p)     1.99548      -0.49006
        7    O  1  py     Ryd( 3p)     0.00095       1.20360
        8    O  1  pz     Val( 2p)     1.73514      -0.41093
        9    O  1  pz     Ryd( 3p)     0.00034       1.24945
        10    O  1  dxy    Ryd( 3d)     0.00000       3.11873
        11    O  1  dxz    Ryd( 3d)     0.00359       3.69068
        12    O  1  dyz    Ryd( 3d)     0.00121       3.07541
        13    O  1  dx2y2  Ryd( 3d)     0.00101       3.43979
        14    O  1  dz2    Ryd( 3d)     0.00133       3.14194

        15    H  2  s      Val( 1s)     0.52980       0.25861
        16    H  2  s      Ryd( 2s)     0.00136       0.45874
        17    H  2  px     Ryd( 2p)     0.00209       2.60593
        18    H  2  py     Ryd( 2p)     0.00118       1.94640
        19    H  2  pz     Ryd( 2p)     0.00117       2.28805

        20    H  3  s      Val( 1s)     0.52980       0.25861
        21    H  3  s      Ryd( 2s)     0.00136       0.45874
        22    H  3  px     Ryd( 2p)     0.00209       2.60593
        23    H  3  py     Ryd( 2p)     0.00118       1.94640
        24    H  3  pz     Ryd( 2p)     0.00117       2.28805'''

        if 'NAO Atom No lang   Type(AO)    Occupancy      Energy' in line:

            line = next(inputfile)
            line = next(inputfile)

            while 'Summary of Natural Population Analysis:' not in line:
                if len(line.strip()) <= 0:
                    line = next(inputfile)
                    continue
                    
                natural_population = line.split()

                nao       = int(natural_population[0])
                atom      = natural_population[1]
                no        = int(natural_population[2])
                lang      = natural_population[3]
                type_ao   = natural_population[4] + natural_population[5]
                occupancy = float(natural_population[6])
                energy    = float(natural_population[7])

                # TODO append to attibutes
                # self.append_attribute('coreelectrons', core)
                # 
                # self.set_attribute('charge', sum(natural_charge))

                line = next(inputfile)
                    

        ''' Summary of Natural Population Analysis:

                                            Natural Population
                    Natural    ---------------------------------------------
        Atom No    Charge        Core      Valence    Rydberg      Total
        --------------------------------------------------------------------
            O  1   -0.92878      1.99998     6.91655    0.01225     8.92878
            H  2    0.46439      0.00000     0.52980    0.00581     0.53561
            H  3    0.46439      0.00000     0.52980    0.00581     0.53561
        ====================================================================
        * Total *  0.00000      1.99998     7.97616    0.02386    10.00000'''

        if '  Atom No    Charge' in line:
    
            line = next(inputfile)
            line = next(inputfile)
            
            while '==============' not in line:
                population_analysis = line.split()
                
                print(population_analysis)

                atom = population_analysis[0]
                no = int(population_analysis[1])
                natural_charge = float(population_analysis[2])
                core = float(population_analysis[3])
                valence = float(population_analysis[4])
                rydberg = float(population_analysis[5])
                total = float(population_analysis[6])
                                
                # TODO append to attibutes
                #             
                # self.append_attribute('naos', nao)

                line = next(inputfile)

        

