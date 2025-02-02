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
        # CFOUR uses A'' instead of A"
        label = label.replace("''", '"')
        label = label.replace("+","")
        # CFOUR uses 1g, 1u, 2g, 2u,... for E1g, E1u, E2g, E2u,...
        try:
            label_int = int(label[0])
            return 'E'+ label
        except:
            if len(label)>=2:
                # CFOUR uses SG for sigma, PI for pi, DE for delta, and PH for phi
                if "u" == label:
                    label = "Eu"
                if "g" == label:
                    label = "Eg"
                if "SG" == label[:2]:
                    if len(label) == 2:
                        label = "sigma"
                    else:
                        label = "sigma." + label[2]
                if "PI" == label[:2]:
                    if len(label) == 2:
                        label = "pi"
                    else:
                        label = "pi." + label[2]
                if "DE" == label[:2]:
                    if len(label) == 2:
                        label = "delta"
                    else:
                        label = "delta." + label[2]
                if "PH" == label[:2]:
                    if len(label) == 2:
                        label = "phi"
                    else:
                        label = "phi." + label[2]
            return label

    def before_parsing(self):
        # set package metadata to CFOUR
        self.metadata["package"] = "CFOUR"
        # geting atomic number and symbol is different for 1 atom
        self.set_attribute("only_one_atom", False)
        # set to true so that atomic numbers are parsed on the first block of coordinates
        self.set_attribute("first_coord_block", True)
        # set the atomic coordinates to []
        self.set_attribute("atomcoords", [])
        # set the list of scf energies to []
        self.set_attribute("scfenergies", [])
        # set mo energies list to []
        self.set_attribute("moenergies", [])
        # set mo symmetries list to []
        self.set_attribute("mosyms", [])
        # set to True so that alpha MOs are parsed first
        self.set_attribute("alpha_mos_to_parse", True)
        # set temp alpha mo coeffs to []
        self.set_attribute("temp_alpha_mocoeffs", [])
        # set temp beta mo coeffs to []
        self.set_attribute("temp_beta_mocoeffs", [])
        # set mo coeffs to []
        self.set_attribute("mocoeffs", [])
        # should mo coeffs be reset the next time an mo is parsed
        self.set_attribute("mocoeffs_should_be_reset", True)
        # set to True so that the first time parsing mo coeffs it also parses ao names
        self.set_attribute("parse_aonames", True)
        # set cc energies to []
        self.set_attribute("ccenergies", [])
        # dict of ECP labels
        self.set_attribute("ecp_labels", {})
        # core electron dict
        self.set_attribute("core_electron_dict", {})
        # set excitation energies to []
        self.set_attribute("etenergies", [])
        # set etoscs to []
        self.set_attribute("etoscs", [])
        # set etsecs to []
        self.set_attribute("etsecs", [])
        # set etsyms to []
        self.set_attribute("etsyms", [])
        # set sym numbering to {}
        self.set_attribute("sym_numbering", {})
        # set current symmetry to "0"
        self.set_attribute("curr_sym", "0")

    def after_parsing(self):
        # change atomic coordinates to a numpy array
        self.atomcoords = np.array(self.atomcoords)
        # change scf energies to a numpy array
        self.scfenergies = np.array(self.scfenergies)
        # get the number of atoms
        if len(self.atomcoords) >= 1:
            self.set_attribute("natom", len(self.atomnos))
        # get the number of atomic orbitals in the basis
        if hasattr(self, "aonames"):
            self.set_attribute("nbasis", len(self.aonames))
        # get the number of molecular orbitals
        if len(self.moenergies) >= 1:
            self.set_attribute("nmo", len(self.moenergies[0]))
        # get core electrons
        for i in self.ecp_labels.keys():
            if i in self.core_electron_dict.keys():
                for j in self.ecp_labels[i]:
                    self.coreelectrons[j] = self.core_electron_dict[i]
        # sort etenergies, etoscs, etsecs, and etsyms
        sort_inds=np.argsort(self.etenergies)
        temp_etenergies=[]
        temp_etoscs=[]
        temp_etsecs=[]
        temp_etsyms=[]
        for i in sort_inds:
            temp_etenergies.append(self.etenergies[i])
            if self.estate_prop_on:
                temp_etoscs.append(self.etoscs[i])
            temp_etsecs.append(self.etsecs[i])
            temp_etsyms.append(self.etsyms[i])
        self.etenergies=temp_etenergies
        self.etoscs=temp_etoscs
        self.etsecs=temp_etsecs
        self.etsyms=temp_etsyms

    def extract(self, inputfile, line):
        # get the version of CFOUR
        if "Version" in line:
            self.metadata["package_version"] = line.split()[1]
        # get the name of the basis set used
        if "BASIS                IBASIS" in line:
            self.metadata["basis_set"] = line.split()[2]
        # get calc_level
        if "CALCLEVEL            ICLLVL" in line:
            self.set_attribute("calc_level",line.strip().split()[2])
        # get excited_states_method
        if "EXCITE               IEXCIT" in line:
            if not line.strip().split()[2]=="NONE":
                self.metadata["excited_states_method"]=line.strip().split()[2]+"-"+self.calc_level
        # get whether the reference is unrestricted or not
        if "REFERENCE            IREFNC" in line:
            if line.split()[2][0] == "U":
                self.metadata["unrestricted"] = True
                self.set_attribute('homos',np.array([0,0]))
            else:
                if line.split()[2] == "ROHF":
                    self.set_attribute("homos", np.array([0, 0]))
                else:
                    self.set_attribute("homos", np.array([0]))
                self.metadata["unrestricted"] = False
        # get full point group
        if "The full molecular point group is" in line:
            self.metadata["symmetry_detected"] = line.split()[6]
        # get used point group
        if "The computational point group is" in line:
            self.metadata["symmetry_used"] = line.split()[5]
        # get the net charge of the system
        if "CHARGE               ICHRGE" in line:
            self.set_attribute("charge", utils.float(line.split()[2]))
        # get estate_prop state
        if "ESTATE_PROP          IEXPRP" in line:
            if line.strip().split()[2]=="OFF":
                self.set_attribute("estate_prop_on",False)
            else:
                self.set_attribute("estate_prop_on",True)
        # get the spin multiplicity of the system
        if "MULTIPLICTY          IMULTP" in line:
            self.set_attribute("mult", int(line.split()[2]))
        # get coupled cluster energy
        if "A miracle has come to pass. The CC iterations have converged." in line:
            cc_lines = []
            while "@CHECKOUT-I," not in line:
                if not line.strip()=='':
                    cc_lines.append(line)
                line = next(inputfile)
            if cc_lines[-1].split()[-1] == "a.u.":
                self.ccenergies.append(utils.float(cc_lines[-1].split()[-2]))
            else:
                self.ccenergies.append(utils.float(cc_lines[-1].split()[-1]))
        # get coefficients and exponents of the gaussian basis set
        if "ATOM                 EXPONENT      COEFFICIENTS" in line:
            atom_index = {}
            line = next(inputfile)
            should_parse = True
            last_line = ""
            temp_basis_info = []
            gbasis = []
            self.ecp_labels={}
            for i in range(len(self.atomic_symbols)):
                gbasis.append([])
            hashtag_in_last_line = False
            first_iter = True
            while True:
                if line.strip()=='':
                    line=next(inputfile)
                    continue
                line=line.replace('+','')
                split_line=line.strip().split()
                line_length=len(split_line)
                if not ((hashtag_in_last_line) or ('#' in line) or (len(last_line.strip().split())==line_length)):
                    break
                if '#' in line:
                    if (not first_iter) and should_parse:
                        for i in temp_basis_info:
                            gbasis[atom_index[curr_atom]].append((curr_ang_mom, i))
                    first_iter = False
                    hashtag_in_last_line = True
                    if len(split_line[0])==1:
                        curr_atom = split_line[0] + split_line[1] + split_line[2]
                        ecp_label = split_line[0] + split_line[1]
                    else:
                        curr_atom = split_line[0] + split_line[1]
                        ecp_label = split_line[0]
                    symbol_len = 0
                    for i in curr_atom:
                        if not i=="#":
                            symbol_len+=1
                        else:
                            break
                    if curr_atom not in atom_index.keys():
                        for i in range(len(self.atomic_symbols)):
                            if self.atomic_symbols[i] == curr_atom[:symbol_len] and (
                                i not in atom_index.values()
                            ):
                                if ecp_label in self.ecp_labels.keys():
                                    self.ecp_labels[ecp_label] += [i]
                                else:
                                    self.ecp_labels[ecp_label] = [i]
                                atom_index[curr_atom] = i
                                break
                    if len(split_line[-1]) == 1:
                        if split_line[-1] == "S":
                            curr_ang_mom = "S"
                            should_parse = True
                        elif split_line[-1] == "X":
                            curr_ang_mom = "P"
                            should_parse = True
                        else:
                            should_parse = False
                    elif len(split_line[-1]) == 2:
                        if split_line[-1] == "XX":
                            curr_ang_mom = "D"
                            should_parse = True
                        else:
                            should_parse=False
                    else:
                        if (
                            (int(split_line[-1][-3]) > 2)
                            and (int(split_line[-1][-2]) == 0)
                            and (int(split_line[-1][-1]) == 0)
                        ):
                            curr_ang_mom = split_line[-1][:-3]
                            should_parse = True
                        else:
                            should_parse=False
                if should_parse and (not ('#' in line)):
                    if hashtag_in_last_line:
                        temp_basis_info=[]
                        for i in range(line_length-2):
                            temp_basis_info.append([])
                    for i in range(line_length-2):
                        if not float(split_line[2+i])==0.:
                            temp_basis_info[i].append((float(split_line[1]),float(split_line[2+i])))
                if hashtag_in_last_line and (not ('#' in line)):
                    hashtag_in_last_line=False
                last_line=line
                line=next(inputfile)
            if should_parse:
                for i in temp_basis_info:
                    gbasis[atom_index[curr_atom]].append((curr_ang_mom, i))
            self.set_attribute("gbasis", gbasis)
        # exception for only one atom
        if "1 entries found in Z-matrix" in line:
            self.only_one_atom=True
        if ("NUCLEAR CHARGE:" in line) and self.only_one_atom:
            self.set_attribute("atomnos", [int(line.strip().split()[2])])
        if ("NUCLEAR COORDINATES (IN A.U.) ARE :" in line) and self.only_one_atom:
            line = next(inputfile)
            line = next(inputfile)
            self.set_attribute("atomic_symbols", [line.strip().split()[0]])
            self.atomcoords.append([[float(line.strip().split()[2]),float(line.strip().split()[3]),float(line.strip().split()[4])]])
        # get the coordinates at each step in a geometry optimization
        # if this is the first time parsing a block of coordinates also get the atomic numbers
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
                self.set_attribute("atomnos", atomnos)
                self.set_attribute("coreelectrons", np.zeros(len(atomnos)))
                self.set_attribute("atomic_symbols", atomic_symbols)
            self.atomcoords.append(temp_atomcoords)
            self.first_coord_block = False
        # get core electrons in each atoms ECP
        if "ECP PARAMETERS FOR ATOM" in line:
            num_ce_index = line.strip().split()[-1]
            line = next(inputfile)
            let_ce_index = line.strip().split(":")[0]
            ce_index = let_ce_index + "#" + num_ce_index
            while "NCORE =" not in line:
                line = next(inputfile)
            self.core_electron_dict[ce_index] = int(line.strip().split()[2])
        # get the scf energy at each step in a geometry optimization
        if "E(SCF)=" in line:
            self.scfenergies.append(utils.float(line.split()[1]))
        # get alpha mo energies of the last ran scf method
        if "ORBITAL EIGENVALUES (ALPHA)  (1H = 27.2113834 eV)" in line:
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            self.moenergies = []
            self.mosyms = []
            alpha_moenergies = []
            alpha_mosyms = []
            while not (
                ("VSCF finished." in line)
                or ("ORBITAL EIGENVALUES ( BETA)  (1H = 27.2113834 eV)" in line)
                or ("SCF failed to converge in" in line)
            ):
                if ("+++++" in line) or (line.strip() == ""):
                    if "+++++" in line:
                        self.homos[0] = int(last_line.strip().split()[0])-1
                    line = next(inputfile)
                    continue
                alpha_moenergies.append(utils.float(line.split()[2]))
                alpha_mosyms.append(self.normalisesym(line.split()[4]))
                if self.normalisesym(line.split()[5]) not in self.sym_numbering.values():
                    self.sym_numbering[line.split()[6][1:-1]]=self.normalisesym(line.split()[5])
                last_line=line
                line = next(inputfile)
            self.moenergies.append(np.array(alpha_moenergies))
            self.mosyms.append(alpha_mosyms)
        # get beta mo energies of the last ran scf method if an unrestricted reference is used
        if "ORBITAL EIGENVALUES ( BETA)  (1H = 27.2113834 eV)" in line:
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            beta_moenergies = []
            beta_mosyms = []
            while not (("VSCF finished." in line) or ("SCF failed to converge in" in line)):
                if ("+++++" in line) or (line.strip() == ""):
                    if "+++++" in line:
                        self.homos[1] = int(last_line.strip().split()[0])-1
                    line = next(inputfile)
                    continue
                beta_moenergies.append(utils.float(line.split()[2]))
                beta_mosyms.append(self.normalisesym(line.split()[4]))
                last_line=line
                line = next(inputfile)
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
                            self.temp_beta_mocoeffs = []
                            self.mocoeffs_should_be_reset = True
        # change curr_sym
        if "Beginning symmetry block" in line:
            self.curr_sym=line.strip().split()[3][:-1]
        # get excitation energies
        if 'Converged eigenvalue:' in line:
            temp_etenergy=float(line.strip().split()[2])
            if self.estate_prop_on:
                num_dash_lines=0
                parse_etsecs=False
                keep_parse=True
                temp_etsecs=[]
                mult="Singlet"
                line = next(inputfile)
                while "Right Transition Moment" not in line:
                    if "Converged eigenvalue:" in line:
                        temp_etenergy = float(line.strip().split()[2])
                        parse_etsecs = True
                    if parse_etsecs and (
                        "--------------------------------------------------------------------------------"
                        in line
                    ):
                        keep_parse = False
                    if parse_etsecs and keep_parse and num_dash_lines == 2:
                        i=int(line.strip().split()[0]) if int(line.strip().split()[0])==0 else int(line.strip().split()[0])-1
                        j=int(line.strip().split()[1]) if int(line.strip().split()[1])==0 else int(line.strip().split()[1])-1
                        a=int(line.strip().split()[2]) if int(line.strip().split()[2])==0 else int(line.strip().split()[2])-1
                        b=int(line.strip().split()[3]) if int(line.strip().split()[3])==0 else int(line.strip().split()[3])-1
                        temp_etsecs.append(((i, j), (a, b), np.float64(line.strip().split()[4])))
                    if (
                        parse_etsecs
                        and keep_parse
                        and num_dash_lines == 3
                        and ("The state is a" in line)
                    ):
                        mult = line.strip().split()[4]
                        if mult == "triplet":
                            mult = "Triplet"
                        else:
                            mult='Singlet'
                    if parse_etsecs and ('--------------------------------------------------------------------------------' in line):
                        num_dash_lines+=1
                        keep_parse=True
                    line = next(inputfile)
                self.etenergies.append(temp_etenergy)
                self.etsecs.append(temp_etsecs)
                self.etsyms.append('{}-{}'.format(mult,self.sym_numbering[self.curr_sym]))
            else:
                self.etenergies.append(temp_etenergy)
                line = next(inputfile)
                line = next(inputfile)
                line = next(inputfile)
                line = next(inputfile)
                line = next(inputfile)
                line = next(inputfile)
                line = next(inputfile)
                temp_etsecs = []
                while (
                    "--------------------------------------------------------------------------------"
                    not in line
                ):
                    i=int(line.strip().split()[0]) if int(line.strip().split()[0])==0 else int(line.strip().split()[0])-1
                    j=int(line.strip().split()[1]) if int(line.strip().split()[1])==0 else int(line.strip().split()[1])-1
                    a=int(line.strip().split()[2]) if int(line.strip().split()[2])==0 else int(line.strip().split()[2])-1
                    b=int(line.strip().split()[3]) if int(line.strip().split()[3])==0 else int(line.strip().split()[3])-1
                    temp_etsecs.append(((i, j), (a, b), np.float64(line.strip().split()[4])))
                    line = next(inputfile)
                self.etsecs.append(temp_etsecs)
                self.etsyms.append('Singlet-{}'.format(self.sym_numbering[self.curr_sym]))
        # get etoscs
        if 'Norm of oscillator strength :' in line:
            self.etoscs.append(float(line.strip().split()[-1]))
