# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for CFOUR output files"""

import datetime
import re

from cclib.parser import data, logfileparser, utils
from cclib.parser.logfileparser import StopParsing

import numpy as np


class CFOUR(logfileparser.Logfile):
    """A CFOUR v2.1 log file."""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="CFOUR", *args, **kwargs)

    def __str__(self):
        #return a string representation of the object
        return f"CFOUR log file {self.filename}"

    def __repr__(self):
        #return a representation of the object
        return f'CFOUR("{self.filename}")'

    def normalisesym(self, label):
        pass

    def before_parsing(self):
        #set to true so that atomic numbers are parsed on the first block of coordinates
        self.set_attribute('first_coord_block',True)
        #set the atomic coordinates to []
        self.set_attribute('atomcoords',[])
        #set the list of scf energies to []
        self.set_attribute('scfenergies',[])
        #set mo energies list to []
        self.set_attribute('moenergies',[])
        #set mo symmetries list to []
        self.set_attribute('mosyms',[])

    def after_parsing(self):
        #change atomic coordinates to a numpy array
        self.atomcoords=np.array(self.atomcoords)
        #change scf energies to a numpy array
        self.scfenergies=np.array(self.scfenergies)
        #get the number of atoms
        if len(self.atomcoords)>=1:
            self.set_attribute('natom',len(self.atomcoords[-1]))
        #get the number of molecular orbitals
        if len(self.moenergies)>=1:
            self.set_attribute('nmo',len(self.moenergies[0]))

    def extract(self, inputfile, line):
        #set package metadata to CFOUR
        self.metadata['package']='CFOUR'
        #get the version of CFOUR
        if 'Version' in line:
            self.metadata['package_version']=line.split()[1]
        #get the name of the basis set used
        if 'BASIS                IBASIS' in line:
            self.metadata['basis_set']=line.split()[2]
        #get whether the reference is unrestricted or not
        if 'REFERENCE            IREFNC' in line:
            self.metadata['unrestricted']=(True if line.split()[2][0]=='U' else False)
        #get the number of atomic orbitals in the basis
        if ('There are' in line) and ('functions in the AO basis.' in line):
            self.set_attribute('nbasis',int(line.split()[2]))
        #get the net charge of the system
        if 'CHARGE               ICHRGE' in line:
            self.set_attribute('charge',utils.float(line.split()[2]))
        #get the spin multiplicity of the system
        if 'MULTIPLICTY          IMULTP' in line:
            self.set_attribute('mult',int(line.split()[2]))
        #get the coordinates at each step in a geometry optimization
        #if this is the first time parsing a block of coordinates also get the atomic numbers
        '''
         ----------------------------------------------------------------
                     Coordinates used in calculation (QCOMP)
         ----------------------------------------------------------------
          Z-matrix   Atomic            Coordinates (in bohr)
           Symbol    Number           X              Y              Z
         ----------------------------------------------------------------
             C         6        -0.00000000     0.00000000     1.23064327
             O         8         0.00000000     0.00000000    -0.92327590
         ----------------------------------------------------------------
        '''
        if 'Coordinates used in calculation (QCOMP)' in line:
            line=next(inputfile)
            line=next(inputfile)
            line=next(inputfile)
            line=next(inputfile)
            line=next(inputfile)
            atomnos=[]
            temp_atomcoords=[]
            while not '----------------------------------------------------------------' in line:
                if self.first_coord_block:
                    atomnos.append(int(line.strip().split()[1]))
                temp_atomcoords.append([utils.convertor(utils.float(line.strip().split()[2]),'bohr','Angstrom'),utils.convertor(utils.float(line.strip().split()[3]),'bohr','Angstrom'),utils.convertor(utils.float(line.strip().split()[4]),'bohr','Angstrom')])
                line=next(inputfile)
            if self.first_coord_block:
                self.set_attribute('atomnos',atomnos)
            self.atomcoords.append(temp_atomcoords)
            self.first_coord_block=False
        #get the scf energy at each step in a geometry optimization
        if 'E(SCF)=' in line:
            self.scfenergies.append(utils.float(line.split()[1]))
        #get alpha mo energies of the last ran scf method
        if 'ORBITAL EIGENVALUES (ALPHA)  (1H = 27.2113834 eV)' in line:
            line=next(inputfile)
            line=next(inputfile)
            line=next(inputfile)
            line=next(inputfile)
            self.moenergies=[]
            self.mosyms=[]
            alpha_moenergies=[]
            alpha_mosyms=[]
            while not (('VSCF finished.' in line)or('ORBITAL EIGENVALUES ( BETA)  (1H = 27.2113834 eV)' in line)):
                if ('+++++' in line) or (line.strip()==''):
                    line=next(inputfile)
                    continue
                alpha_moenergies.append(utils.float(line.split()[2]))
                alpha_mosyms.append(line.split()[5])
                line=next(inputfile)
            self.moenergies.append(np.array(alpha_moenergies))
            self.mosyms.append(alpha_mosyms)
        #get beta mo energies of the last ran scf method if an unrestricted reference is used
        if 'ORBITAL EIGENVALUES ( BETA)  (1H = 27.2113834 eV)' in line:
            line=next(inputfile)
            line=next(inputfile)
            line=next(inputfile)
            line=next(inputfile)
            beta_moenergies=[]
            beta_mosyms=[]
            while not ('VSCF finished.' in line):
                if ('+++++' in line) or (line.strip()==''):
                    line=next(inputfile)
                    continue
                beta_moenergies.append(utils.float(line.split()[2]))
                beta_mosyms.append(line.split()[5])
                line=next(inputfile)
            self.moenergies.append(np.array(beta_moenergies))
            self.mosyms.append(beta_mosyms)




