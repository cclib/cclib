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
<<<<<<< HEAD
        #CFOUR uses 1g, 1u, 2g, 2u,... for E1g, E1u, E2g, E2u,...
        try:
            label_int=int(label[0])
            return 'E'+label
        except:
            if len(label)>=2:
                #CFOUR uses '' for "
                if "''" in label:
                    label.replace("''",'"')
                #CFOUR uses SG for sigma, PI for pi, DE for delta, and PH for phi
                if 'SG'==label[:2]:
                    if (len(label)==2)or((len(label)==3)and((label[2]=='+')or(label[2]=='-'))):
                        label='sigma'
=======
        # CFOUR uses A'' instead of A"
        label.replace("''",'"')
        # CFOUR uses 1g, 1u, 2g, 2u,... for E1g, E1u, E2g, E2u,...
        try:
            label_int = int(label[0])
            return 'E'+ label
        except:
            if len(label)>=2:
                # CFOUR uses SG for sigma, PI for pi, DE for delta, and PH for phi
                if "SG" == label[:2]:
                    if len(label) == 2:
                        label = "sigma"
>>>>>>> ec752d40 (Improved normalizing symmetry labels in CFOUR and its corresponding test)
                    else:
                        label='sigma.'+label[2]
                if 'PI'==label[:2]:
                    if len(label)==2:
                        label='pi'
                    else:
                        label='pi.'+label[2:]
                if 'DE'==label[:2]:
                    if len(label)==2:
                        label='delta'
                    else:
                        label='delta.'+label[2:]
                if 'PH'==label[:2]:
                    if len(label)==2:
                        label='phi'
                    else:
                        label='phi.'+label[2:]
            return label

    def before_parsing(self):
        #set package metadata to CFOUR
        self.metadata['package']='CFOUR'
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
        #set to True so that alpha MOs are parsed first
        self.set_attribute('alpha_mos_to_parse',True)
        #set temp alpha mo coeffs to []
        self.set_attribute('temp_alpha_mocoeffs',[])
        #set temp beta mo coeffs to []
        self.set_attribute('temp_beta_mocoeffs',[])
        #set mo coeffs to []
        self.set_attribute('mocoeffs',[])
        #should mo coeffs be reset the next time an mo is parsed
        self.set_attribute('mocoeffs_should_be_reset',True)
        #set to True so that the first time parsing mo coeffs it also parses ao names
        self.set_attribute('parse_aonames',True)
        #set cc energies to []
        self.set_attribute('ccenergies',[])

    def after_parsing(self):
        #change atomic coordinates to a numpy array
        self.atomcoords=np.array(self.atomcoords)
        #change scf energies to a numpy array
        self.scfenergies=np.array(self.scfenergies)
        #get the number of atoms
        if len(self.atomcoords)>=1:
            self.set_attribute('natom',len(self.atomnos))
        #get the number of atomic orbitals in the basis
        if hasattr(self,'aonames'):
            self.set_attribute('nbasis',len(self.aonames))
        #get the number of molecular orbitals
        if len(self.moenergies)>=1:
            self.set_attribute('nmo',len(self.moenergies[0]))

    def extract(self, inputfile, line):
        #get the version of CFOUR
        if 'Version' in line:
            self.metadata['package_version']=line.split()[1]
        #get the name of the basis set used
        if 'BASIS                IBASIS' in line:
            self.metadata['basis_set']=line.split()[2]
        #get whether the reference is unrestricted or not
        if 'REFERENCE            IREFNC' in line:
            self.metadata['unrestricted']=(True if line.split()[2][0]=='U' else False)
        #get full point group
        if 'The full molecular point group is' in line:
            self.metadata['symmetry_detected']=line.split()[6]
        #get used point group
        if 'The computational point group is' in line:
            self.metadata['symmetry_used']=line.split()[5]
        #get the net charge of the system
        if 'CHARGE               ICHRGE' in line:
            self.set_attribute('charge',utils.float(line.split()[2]))
        #get the spin multiplicity of the system
        if 'MULTIPLICTY          IMULTP' in line:
            self.set_attribute('mult',int(line.split()[2]))
        if 'A miracle has come to pass. The CC iterations have converged.' in line:
            cc_lines=[]
            while not '@CHECKOUT-I,' in line:
                cc_lines.append(line)
                line=next(inputfile)
            if cc_lines[-1].split()[-1]=='a.u.':
                self.ccenergies.append(utils.float(cc_lines[-1].split()[-2]))
            else:
                self.ccenergies.append(utils.float(cc_lines[-1].split()[-1]))
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
            atomic_symbols=[]
            temp_atomcoords=[]
            while not '----------------------------------------------------------------' in line:
                if self.first_coord_block:
                    atomic_number=int(line.strip().split()[1])
                    atomic_symbol=line.strip().split()[0]
                    if atomic_number>0:
                        atomnos.append(atomic_number)
                        atomic_symbols.append(atomic_symbol)
                temp_atomcoords.append([utils.convertor(utils.float(line.strip().split()[2]),'bohr','Angstrom'),utils.convertor(utils.float(line.strip().split()[3]),'bohr','Angstrom'),utils.convertor(utils.float(line.strip().split()[4]),'bohr','Angstrom')])
                line=next(inputfile)
            if self.first_coord_block:
                self.set_attribute('atomnos',atomnos)
                self.set_attribute('atomic_symbols',atomic_symbols)
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
            while not (('VSCF finished.' in line)or('ORBITAL EIGENVALUES ( BETA)  (1H = 27.2113834 eV)' in line)or('SCF failed to converge in' in line)):
                if ('+++++' in line) or (line.strip()==''):
                    line=next(inputfile)
                    continue
                alpha_moenergies.append(utils.float(line.split()[2]))
                alpha_mosyms.append(self.normalisesym(line.split()[4]))
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
            while not (('VSCF finished.' in line)or('SCF failed to converge in' in line)):
                if ('+++++' in line) or (line.strip()==''):
                    line=next(inputfile)
                    continue
                beta_moenergies.append(utils.float(line.split()[2]))
                beta_mosyms.append(self.normalisesym(line.split()[4]))
                line=next(inputfile)
            self.moenergies.append(np.array(beta_moenergies))
            self.mosyms.append(beta_mosyms)
        #add on to the molecular orbital coefficeints
        if len(line.split())>=2:
            if 'Symmetry'==line.split()[0]:
                line=next(inputfile)
                if '--------------------------------------------------------------------------------' in line:
                    if self.mocoeffs_should_be_reset:
                        self.mocoeffs=[]
                        self.mocoeffs_should_be_reset=False
                    line=next(inputfile)
                    len_split_line=len(line.split())
                    if self.alpha_mos_to_parse:
                        if self.parse_aonames:
                            aonames=[]
                            atombasis=[]
                            start_atom=0
                            sub_amount=0
                            basis_index=0
                        for i in range(len_split_line-2):
                            self.temp_alpha_mocoeffs.append([])
                        while not '----------' in line:
                            split_line=line.split()
                            if self.parse_aonames:
                                if int(split_line[0])-start_atom>1:
                                    sub_amount+=((int(split_line[0])-start_atom)-1)
                                    start_atom=int(split_line[0])
                                    atombasis.append([])
                                    subshell_number=0
                                elif int(split_line[0])-start_atom==1:
                                    start_atom=int(split_line[0])
                                    atombasis.append([])
                                    subshell_number=0
                                if split_line[1]=='S':
                                    subshell_number+=1
                                if len(split_line[1])==2:
                                    if split_line[1][1]=='X':
                                        subshell_number+=1
                                if len(split_line[1])>=3:
                                    if split_line[1][1]=='X':
                                        try:
                                            x_power=int(split_line[1][2:])
                                            subshell_number+=1
                                        except:
                                            pass
                                aonames.append('{}{}_{}{}'.format(self.atomic_symbols[(int(split_line[0])-sub_amount)-1],int(split_line[0])-sub_amount,subshell_number,split_line[1]))
                                atombasis[-1].append(basis_index)
                                basis_index+=1
                            for i in range(len_split_line-2):
                                self.temp_alpha_mocoeffs[i-(len_split_line-2)].append(utils.float(split_line[i+2]))
                            line=next(inputfile)
                        if self.parse_aonames:
                            self.set_attribute('aonames',aonames)
                            self.set_attribute('atombasis',atombasis)
                            self.parse_aonames=False
                        line=next(inputfile)
                        if ('-----------' in line)or(''==line.strip()):
                            self.alpha_mos_to_parse=False
                            self.mocoeffs.append(self.temp_alpha_mocoeffs)
                            self.temp_alpha_mocoeffs=[]
                            if (''==line.strip()):
                                self.alpha_mos_to_parse=True
                                self.mocoeffs_should_be_reset=True
                    else:
                        for i in range(len_split_line-2):
                            self.temp_beta_mocoeffs.append([])
                        while not '----------------' in line:
                            split_line=line.split()
                            for i in range(len_split_line-2):
                                self.temp_beta_mocoeffs[i-(len_split_line-2)].append(utils.float(split_line[i+2]))
                            line=next(inputfile)
                        line=next(inputfile)
                        if (''==line.strip()):
                            self.alpha_mos_to_parse=True
                            self.mocoeffs.append(self.temp_beta_mocoeffs)
                            self.temp_beta_mocoeffs=[]
                            self.mocoeffs_should_be_reset=True







