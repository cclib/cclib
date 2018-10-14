# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Turbomole output files."""

from __future__ import print_function

import re

import numpy

from cclib.parser import logfileparser
from cclib.parser import utils

class AtomBasis:
    def __init__(self, atname, basis_name, inputfile):
        self.symmetries=[]
        self.coefficients=[]
        self.atname=atname
        self.basis_name=basis_name

        self.parse_basis(inputfile)

    def parse_basis(self, inputfile):
        i=0
        line=inputfile.next()

        while(line[0]!="*"):
            (nbasis_text, symm)=line.split()
            self.symmetries.append(symm)

            nbasis=int(nbasis_text)
            coeff_arr=numpy.zeros((nbasis, 2), float)

            for j in range(0, nbasis, 1):
                line=inputfile.next()
                (e1_text, e2_text)=line.split()
                coeff_arr[j][0]=float(e1_text)
                coeff_arr[j][1]=float(e2_text)

            self.coefficients.append(coeff_arr)
            line=inputfile.next()

class Turbomole(logfileparser.Logfile):
    """A Turbomole log file."""

    def __init__(self, *args, **kwargs):
        super(Turbomole, self).__init__(logname="Turbomole", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "Turbomole output file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Turbomole("%s")' % (self.filename)

    def normalisesym(self, label):
        """Normalise the symmetries used by Turbomole."""
        raise NotImplementedError('Now yet implemented for Turbomole.')

    def before_parsing(self):
        self.geoopt = False # Is this a GeoOpt? Needed for SCF targets/values.
        self.periodic_table = utils.PeriodicTable()

    @staticmethod
    def split_molines(inline):
        """Splits the lines containing mocoeffs (each of length 20)
        and converts them to float correctly.
        """
        line = inline.replace("D", "E")
        f1 = line[0:20]
        f2 = line[20:40]
        f3 = line[40:60]
        f4 = line[60:80]

        if(len(f4) > 1):
            return [float(f1), float(f2), float(f3), float(f4)]
        if(len(f3) > 1):
            return [float(f1), float(f2), float(f3)]
        if(len(f2) > 1):
            return [float(f1), float(f2)]
        if(len(f1) > 1):
            return [float(f1)]

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        ## This information is in the control file.
        #   $rundimensions
        #   dim(fock,dens)=1860
        #   natoms=20
        #   nshell=40
        #   nbf(CAO)=60
        #   nbf(AO)=60
        #   dim(trafo[SAO<-->AO/CAO])=60
        #   rhfshells=1
        if line[3:10]=="natoms=":
            self.natom=int(line[10:])

        if line[3:11] == "nbf(AO)=":
            nmo = int(line.split('=')[1])
            self.set_attribute('nbasis', nmo)
            self.set_attribute('nmo', nmo)

        # Extract the version number and optionally the build number.
        searchstr = ": TURBOMOLE"
        index = line.find(searchstr)
        if index > -1:
            line = line[index + len(searchstr):]
            tokens = line.split()
            self.metadata["package_version"] = tokens[0][1:].replace("-", ".")
            # Don't add revision information to the main package version for now.
            if tokens[1] == "(":
                revision = tokens[2]

        ## Atomic coordinates in job.last:
        #              +--------------------------------------------------+
        #              | Atomic coordinate, charge and isotop information |
        #              +--------------------------------------------------+
        #
        #
        #              atomic coordinates              atom shells charge pseudo isotop
        #    -2.69176330   -0.00007129   -0.44712612    c      3    6.000    0     0
        #    -1.69851645   -0.00007332    2.06488947    c      3    6.000    0     0
        #     0.92683848   -0.00007460    2.49592179    c      3    6.000    0     0
        #     2.69176331   -0.00007127    0.44712612    c      3    6.000    0     0
        #     1.69851645   -0.00007331   -2.06488947    c      3    6.000    0     0
        #...
        #    -7.04373606    0.00092244    2.74543891    h      1    1.000    0     0
        #    -9.36352819    0.00017229    0.07445322    h      1    1.000    0     0
        #    -0.92683849   -0.00007461   -2.49592179    c      3    6.000    0     0
        #    -1.65164853   -0.00009927   -4.45456858    h      1    1.000    0     0
        if 'Atomic coordinate, charge and isotop information' in line:
            while 'atomic coordinates' not in line:
                line = next(inputfile)

            atomcoords = []
            atomnos = []
            line = next(inputfile)
            while len(line) > 2:
                atomnos.append(self.periodic_table.number[line.split()[3].upper()])
                atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom") 
                                   for x in line.split()[:3]])
                line = next(inputfile)

            self.append_attribute('atomcoords', atomcoords)
            self.set_attribute('atomnos', atomnos)
            self.set_attribute('natom', len(atomcoords))

        # Frequency values in aoforce.out
        #        mode               7        8        9       10       11       12
        #
        #      frequency          53.33    88.32   146.85   171.70   251.75   289.44
        #
        #      symmetry            a        a        a        a        a        a
        #
        #         IR               YES      YES      YES      YES      YES      YES
        # |dDIP/dQ|   (a.u.)     0.0002   0.0000   0.0005   0.0004   0.0000   0.0000
        # intensity (km/mol)       0.05     0.00     0.39     0.28     0.00     0.00
        # intensity (  %   )       0.05     0.00     0.40     0.28     0.00     0.00
        #
        #        RAMAN             YES      YES      YES      YES      YES      YES
        #
        #   1   c           x   0.00000  0.00001  0.00000 -0.01968 -0.04257  0.00001
        #                   y  -0.08246 -0.08792  0.02675 -0.00010  0.00000  0.17930
        #                   z   0.00001  0.00003  0.00004 -0.10350  0.11992 -0.00003
        if 'NORMAL MODES and VIBRATIONAL FREQUENCIES (cm**(-1))' in line:
            vibfreqs, vibsyms, vibirs, vibdisps = [], [], [], []
            while '****  force : all done  ****' not in line:
                if line.strip().startswith('frequency'):
                    freqs = [float(i.replace('i', '-')) for i in line.split()[1:]]
                    vibfreqs.extend(freqs)
                    self.skip_line(inputfile, ['b'])
                    line = next(inputfile)
                    if line.strip().startswith('symmetry'):
                        syms = line.split()[1:]
                        vibsyms.extend(syms)

                    self.skip_lines(inputfile, ['b', 'IR', 'dQIP'])
                    line = next(inputfile)
                    if line.strip().startswith('intensity (km/mol)'):
                        irs = [self.float(f) for f in line.split()[2:]]
                        vibirs.extend(irs)

                    self.skip_lines(inputfile, ['intensity', 'b', 'raman', 'b'])
                    line = next(inputfile)
                    x, y, z = [], [], []
                    while line.split():
                        x.append([float(i) for i in line.split()[3:]])
                        line = next(inputfile)
                        y.append([float(i) for i in line.split()[1:]])
                        line = next(inputfile)
                        z.append([float(i) for i in line.split()[1:]])
                        line = next(inputfile)

                    for j in range(len(x[0])):
                        disps = []
                        for i in range(len(x)):
                            disps.append([x[i][j], y[i][j], z[i][j]])
                        vibdisps.append(disps)

                line = next(inputfile)

            self.set_attribute('vibfreqs', vibfreqs)
            self.set_attribute('vibsyms', vibsyms)
            self.set_attribute('vibirs', vibirs)
            self.set_attribute('vibdisps', vibdisps)

        # In this section we are parsing mocoeffs and moenergies from
        # the files like: mos, alpha and beta.
        # $scfmo    scfconv=6   format(4d20.14)
        # # SCF total energy is     -382.3457535740 a.u.
        # #
        #      1  a      eigenvalue=-.97461484059799D+01   nsaos=60
        # 0.69876828353937D+000.32405121159405D-010.87670894913921D-03-.85232349313288D-07
        # 0.19361534257922D-04-.23841194890166D-01-.81711001390807D-020.13626356942047D-02
        # ...
        # ...
        # $end
        if (line.startswith('$scfmo') or line.startswith('$uhfmo')) and line.find('scfconv') > 0:
            if line.strip().startswith('$uhfmo_alpha'):
                self.unrestricted = True

            # Need to skip the first line to start with lines starting with '#'.
            line = next(inputfile)
            while line.strip().startswith('#') and not line.find('eigenvalue') > 0:
                line = next(inputfile)

            moenergies = []
            mocoeffs = []

            while not line.strip().startswith('$'):
                info = re.match(".*eigenvalue=(?P<moenergy>[0-9D\.+-]{20})\s+nsaos=(?P<count>\d+).*", line)
                eigenvalue = self.float(info.group('moenergy'))
                orbital_energy = utils.convertor(eigenvalue, 'hartree', 'eV')
                moenergies.append(orbital_energy)
                single_coeffs = []
                nsaos = int(info.group('count'))

                while(len(single_coeffs) < nsaos):
                    line = next(inputfile)
                    single_coeffs.extend(Turbomole.split_molines(line))

                mocoeffs.append(single_coeffs)
                line = next(inputfile)

            max_nsaos = max([len(i) for i in mocoeffs])
            for i in mocoeffs:
                while len(i) < max_nsaos:
                    i.append(numpy.nan)

            if not hasattr(self, 'mocoeffs'):
                self.mocoeffs = []

            if not hasattr(self, 'moenergies'):
                self.moenergies = []

            self.mocoeffs.append(mocoeffs)
            self.moenergies.append(moenergies)

        # Parsing the scfenergies, scfvalues and scftargets from job.last file.
        # scf convergence criterion : increment of total energy < .1000000D-05
        #                  and increment of one-electron energy < .1000000D-02
        #
        # ...
        # ...
        #                                              current damping :  0.700
        # ITERATION  ENERGY          1e-ENERGY        2e-ENERGY     NORM[dD(SAO)]  TOL
        #   1  -382.34543727790    -1396.8009423     570.56292464    0.000D+00 0.556D-09
        #                            Exc =   -57.835278090846     N = 69.997494722
        #          max. resid. norm for Fia-block=  2.782D-05 for orbital     33a
        # ...
        # ...
        #                                              current damping :  0.750
        # ITERATION  ENERGY          1e-ENERGY        2e-ENERGY     NORM[dD(SAO)]  TOL
        #   3  -382.34575357399    -1396.8009739     570.56263988    0.117D-03 0.319D-09
        #                            Exc =   -57.835593208072     N = 69.999813370
        #          max. resid. norm for Fia-block=  7.932D-06 for orbital     33a
        #          max. resid. fock norm         =  8.105D-06 for orbital     33a
        #
        # convergence criteria satisfied after  3 iterations
        #
        #
        #                  ------------------------------------------
        #                 |  total energy      =   -382.34575357399  |
        #                  ------------------------------------------
        #                 :  kinetic energy    =    375.67398458525  :
        #                 :  potential energy  =   -758.01973815924  :
        #                 :  virial theorem    =      1.98255043001  :
        #                 :  wavefunction norm =      1.00000000000  :
        #                  ..........................................
        if 'scf convergence criterion' in line:
            total_energy_threshold = self.float(line.split()[-1])
            one_electron_energy_threshold = self.float(next(inputfile).split()[-1])
            scftargets = [total_energy_threshold, one_electron_energy_threshold]
            self.append_attribute('scftargets', scftargets)
            iter_energy = []
            iter_one_elec_energy = []
            while 'convergence criteria satisfied' not in line:
                if 'ITERATION  ENERGY' in line:
                    line = next(inputfile)
                    info = line.split()
                    iter_energy.append(self.float(info[1]))
                    iter_one_elec_energy.append(self.float(info[2]))
                line = next(inputfile)

            assert len(iter_energy) == len(iter_one_elec_energy), \
                'Different number of values found for total energy and one electron energy.'
            scfvalues = [[x - y, a - b] for x, y, a, b in 
                         zip(iter_energy[1:], iter_energy[:-1], iter_one_elec_energy[1:], iter_one_elec_energy[:-1])]
            self.append_attribute('scfvalues', scfvalues)
            while 'total energy' not in line:
                line = next(inputfile)

            scfenergy = utils.convertor(self.float(line.split()[4]), 'hartree', 'eV')
            self.append_attribute('scfenergies', scfenergy)

        #  **********************************************************************
        #  *                                                                    *
        #  *   RHF  energy                             :    -74.9644564256      *
        #  *   MP2 correlation energy (doubles)        :     -0.0365225363      *
        #  *                                                                    *
        #  *   Final MP2 energy                        :    -75.0009789619      *
        # ...
        #  *   Norm of MP1 T2 amplitudes               :      0.0673494687      *
        #  *                                                                    *
        #  **********************************************************************
        # OR
        #  **********************************************************************
        #  *                                                                    *
        #  *   RHF  energy                             :    -74.9644564256      *
        #  *   correlation energy                      :     -0.0507799360      *
        #  *                                                                    *
        #  *   Final CCSD energy                       :    -75.0152363616      *
        #  *                                                                    *
        #  *   D1 diagnostic                           :      0.0132            *
        #  *                                                                    *
        #  **********************************************************************
        if 'C C S D F 1 2   P R O G R A M' in line:
            while 'ccsdf12 : all done' not in line:
                if 'Final MP2 energy' in line:
                    mp2energy = [utils.convertor(self.float(line.split()[5]), 'hartree', 'eV')]
                    self.append_attribute('mpenergies', mp2energy)

                if 'Final CCSD energy' in line:
                    ccenergy = [utils.convertor(self.float(line.split()[5]), 'hartree', 'eV')]
                    self.append_attribute('ccenergies', ccenergy)

                line = next(inputfile)

        #  *****************************************************
        #  *                                                   *
        #  *      SCF-energy   :     -74.49827196840999        *
        #  *      MP2-energy   :      -0.19254365976227        *
        #  *      total        :     -74.69081562817226        *
        #  *                                                   *
        #  *     (MP2-energy evaluated from T2 amplitudes)     *
        #  *                                                   *
        #  *****************************************************
        if 'm p g r a d - program' in line:
            while 'ccsdf12 : all done' not in line:
                if 'MP2-energy' in line:
                    line = next(inputfile)
                    if 'total' in line:
                        mp2energy = [utils.convertor(self.float(line.split()[3]), 'hartree', 'eV')]
                        self.append_attribute('mpenergies', mp2energy)
                line = next(inputfile)

    def deleting_modes(self, vibfreqs, vibdisps, vibirs):
        """Deleting frequencies relating to translations or rotations"""
        i = 0
        while i < len(vibfreqs):
            if vibfreqs[i] == 0.0:
                # Deleting frequencies that have value 0 since they
                # do not correspond to vibrations.
                del vibfreqs[i], vibdisps[i], vibirs[i]
                i -= 1
            i += 1

    def after_parsing(self):
        if hasattr(self, 'vibfreqs'):
            self.deleting_modes(self.vibfreqs, self.vibdisps, self.vibirs)


class OldTurbomole(logfileparser.Logfile):
    """A Turbomole output file. Code is outdated and is not being used."""

    def __init__(self, *args):
        # Call the __init__ method of the superclass
        super(Turbomole, self).__init__(logname="Turbomole", *args)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Turbomole output file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Turbomole("%s")' % (self.filename)

    def atlist(self, atstr):
        # turn atstr from atoms section into array

        fields=atstr.split(',')
        list=[]
        for f in fields:
            if(f.find('-')!=-1):
                rangefields=f.split('-')
                start=int(rangefields[0])
                end=int(rangefields[1])
                
                for j in range(start, end+1, 1):
                    list.append(j-1)
            else:
                list.append(int(f)-1)
        return(list)

    def normalisesym(self, label):
        """Normalise the symmetries used by Turbomole."""
        return ans

    def before_parsing(self):
        self.geoopt = False # Is this a GeoOpt? Needed for SCF targets/values.

    def split_molines(self, inline):
        line=inline.replace("D", "E")
        f1=line[0:20]
        f2=line[20:40]
        f3=line[40:60]
        f4=line[60:80]

        if(len(f4)>1):
            return( (float(f1), float(f2), float(f3), float(f4)) )
        if(len(f3)>1):
            return( (float(f1), float(f2), float(f3)) )
        if(len(f2)>1):
            return( (float(f1), float(f2)) )
        if(len(f1)>1):
            return([float(f1)])
        return
    
    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        if line[3:11]=="nbf(AO)=":
            nmo=int(line[11:])
            self.nbasis=nmo
            self.nmo=nmo
        if line[3:9]=="nshell":
            temp=line.split('=')
            homos=int(temp[1])

        if line[0:6] == "$basis":
            print("Found basis")
            self.basis_lib=[]
            line = inputfile.next()
            line = inputfile.next()

            while line[0] != '*' and line[0] != '$':
                temp=line.split()
                line = inputfile.next()
                while line[0]=="#":
                    line = inputfile.next()
                self.basis_lib.append(AtomBasis(temp[0], temp[1], inputfile))
                line = inputfile.next()
        if line == "$ecp\n":
            self.ecp_lib=[]
            
            line = inputfile.next()
            line = inputfile.next()
            
            while line[0] != '*' and line[0] != '$':
                fields=line.split()
                atname=fields[0]
                ecpname=fields[1]
                line = inputfile.next()
                line = inputfile.next()
                fields=line.split()
                ncore = int(fields[2])

                while line[0] != '*':
                    line = inputfile.next()
                self.ecp_lib.append([atname, ecpname, ncore])
        
        if line[0:6] == "$coord":
            if line[0:11] == "$coordinate":
#                print "Breaking"
                return

#            print "Found coords"
            self.atomcoords = []
            self.atomnos = []
            atomcoords = []
            atomnos = []

            line = inputfile.next()
            if line[0:5] == "$user":
#                print "Breaking"
                return

            while line[0] != "$":
                temp = line.split()
                atsym=temp[3].capitalize()
                atomnos.append(self.table.number[atsym])
                atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom")
                                   for x in temp[0:3]])
                line = inputfile.next()
            self.atomcoords.append(atomcoords)
            self.atomnos = numpy.array(atomnos, "i")

        if line[14:32] == "atomic coordinates":
            atomcoords = []
            atomnos = []

            line = inputfile.next()
           
            while len(line) > 2:
                temp = line.split()
                atsym = temp[3].capitalize()
                atomnos.append(self.table.number[atsym])
                atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom")
                                    for x in temp[0:3]])
                line = inputfile.next()

            if not hasattr(self,"atomcoords"):
                self.atomcoords = []

            self.atomcoords.append(atomcoords)
            self.atomnos = numpy.array(atomnos, "i")

        if line[0:6] == "$atoms":
            print("parsing atoms")
            line = inputfile.next()
            self.atomlist=[]
            while line[0]!="$":
                temp=line.split()
                at=temp[0]
                atnosstr=temp[1]
                while atnosstr[-1] == ",":
                    line = inputfile.next()
                    temp=line.split()
                    atnosstr=atnosstr+temp[0]
#                print "Debug:", atnosstr
                atlist=self.atlist(atnosstr)

                line = inputfile.next()

                temp=line.split()
#                print "Debug basisname (temp):",temp
                basisname=temp[2]
                ecpname=''
                line = inputfile.next()
                while(line.find('jbas')!=-1 or line.find('ecp')!=-1 or
                      line.find('jkbas')!=-1):
                    if line.find('ecp')!=-1:
                        temp=line.split()
                        ecpname=temp[2]
                    line = inputfile.next()

                self.atomlist.append( (at, basisname, ecpname, atlist))

# I have no idea what this does, so "comment" out
        if line[3:10]=="natoms=":
#        if 0:

            self.natom=int(line[10:])

            basistable=[]

            for i in range(0, self.natom, 1):
                for j in range(0, len(self.atomlist), 1):
                    for k in range(0, len(self.atomlist[j][3]), 1):
                        if self.atomlist[j][3][k]==i:
                            basistable.append((self.atomlist[j][0],
                                                   self.atomlist[j][1],
                                               self.atomlist[j][2]))
            self.aonames=[]
            counter=1
            for a, b, c in basistable:
                ncore=0
                if len(c) > 0:
                    for i in range(0, len(self.ecp_lib), 1):
                        if self.ecp_lib[i][0]==a and \
                           self.ecp_lib[i][1]==c:
                            ncore=self.ecp_lib[i][2]
                           
                for i in range(0, len(self.basis_lib), 1):
                    if self.basis_lib[i].atname==a and self.basis_lib[i].basis_name==b:
                        pa=a.capitalize()
                        basis=self.basis_lib[i]

                        s_counter=1
                        p_counter=2
                        d_counter=3
                        f_counter=4
                        g_counter=5
# this is a really ugly piece of code to assign the right labels to
# basis functions on atoms with an ecp
                        if ncore == 2:
                            s_counter=2
                        elif ncore == 10:
                            s_counter=3
                            p_counter=3
                        elif ncore == 18:
                            s_counter=4
                            p_counter=4
                        elif ncore == 28:
                            s_counter=4
                            p_counter=4
                            d_counter=4
                        elif ncore == 36:
                            s_counter=5
                            p_counter=5
                            d_counter=5
                        elif ncore == 46:
                            s_counter=5
                            p_counter=5
                            d_counter=6
                            
                        for j in range(0, len(basis.symmetries), 1):
                            if basis.symmetries[j]=='s':
                                self.aonames.append("%s%d_%d%s" % \
                                              (pa, counter, s_counter, "S"))
                                s_counter=s_counter+1
                            elif basis.symmetries[j]=='p':
                                self.aonames.append("%s%d_%d%s" % \
                                              (pa, counter, p_counter, "PX"))
                                self.aonames.append("%s%d_%d%s" % \
                                              (pa, counter, p_counter, "PY"))
                                self.aonames.append("%s%d_%d%s" % \
                                              (pa, counter, p_counter, "PZ"))
                                p_counter=p_counter+1
                            elif basis.symmetries[j]=='d':
                                self.aonames.append("%s%d_%d%s" % \
                                         (pa, counter, d_counter, "D 0"))
                                self.aonames.append("%s%d_%d%s" % \
                                         (pa, counter, d_counter, "D+1"))
                                self.aonames.append("%s%d_%d%s" % \
                                         (pa, counter, d_counter, "D-1"))
                                self.aonames.append("%s%d_%d%s" % \
                                         (pa, counter, d_counter, "D+2"))
                                self.aonames.append("%s%d_%d%s" % \
                                         (pa, counter, d_counter, "D-2"))
                                d_counter=d_counter+1
                            elif basis.symmetries[j]=='f':
                                 self.aonames.append("%s%d_%d%s" % \
                                       (pa, counter, f_counter, "F 0"))
                                 self.aonames.append("%s%d_%d%s" % \
                                       (pa, counter, f_counter, "F+1"))
                                 self.aonames.append("%s%d_%d%s" % \
                                       (pa, counter, f_counter, "F-1"))
                                 self.aonames.append("%s%d_%d%s" % \
                                       (pa, counter, f_counter, "F+2"))
                                 self.aonames.append("%s%d_%d%s" % \
                                       (pa, counter, f_counter, "F-2"))
                                 self.aonames.append("%s%d_%d%s" % \
                                       (pa, counter, f_counter, "F+3"))
                                 self.aonames.append("%s%d_%d%s" % \
                                        (pa, counter, f_counter, "F-3"))
                            elif basis.symmetries[j]=='g':
                                self.aonames.append("%s%d_%d%s" % \
                                       (pa, counter, f_counter, "G 0"))
                                self.aonames.append("%s%d_%d%s" % \
                                       (pa, counter, f_counter, "G+1"))
                                self.aonames.append("%s%d_%d%s" % \
                                       (pa, counter, f_counter, "G-1"))
                                self.aonames.append("%s%d_%d%s" % \
                                        (pa, counter, g_counter, "G+2"))
                                self.aonames.append("%s%d_%d%s" % \
                                         (pa, counter, g_counter, "G-2"))
                                self.aonames.append("%s%d_%d%s" % \
                                         (pa, counter, g_counter, "G+3"))
                                self.aonames.append("%s%d_%d%s" % \
                                          (pa, counter, g_counter, "G-3"))
                                self.aonames.append("%s%d_%d%s" % \
                                          (pa, counter, g_counter, "G+4"))
                                self.aonames.append("%s%d_%d%s" % \
                                          (pa, counter, g_counter, "G-4"))
                        break
                counter=counter+1
                
        if line=="$closed shells\n":
            line = inputfile.next()
            temp = line.split()
            occs = int(temp[1][2:])
            self.homos = numpy.array([occs-1], "i")

        if line == "$alpha shells\n":
            line = inputfile.next()
            temp = line.split()
            occ_a = int(temp[1][2:])
            line = inputfile.next() # should be $beta shells
            line = inputfile.next() # the beta occs
            temp = line.split()
            occ_b = int(temp[1][2:])
            self.homos = numpy.array([occ_a-1,occ_b-1], "i")

        if line[12:24]=="OVERLAP(CAO)":
            line = inputfile.next()
            line = inputfile.next()
            overlaparray=[]
            self.aooverlaps=numpy.zeros( (self.nbasis, self.nbasis), "d")
            while line != "       ----------------------\n":
                temp=line.split()
                overlaparray.extend(map(float, temp))
                line = inputfile.next()
            counter=0

            for i in range(0, self.nbasis, 1):
                for j in range(0, i+1, 1):
                    self.aooverlaps[i][j]=overlaparray[counter]
                    self.aooverlaps[j][i]=overlaparray[counter]
                    counter=counter+1

        if ( line[0:6] == "$scfmo" or line[0:12] == "$uhfmo_alpha" ) and line.find("scf") > 0:
            temp = line.split()

            if temp[1][0:7] == "scfdump":
#                self.logger.warning("SCF not converged?")
                print("SCF not converged?!")

            if line[0:12] == "$uhfmo_alpha": # if unrestricted, create flag saying so
                unrestricted = 1
            else:
                unrestricted = 0

            self.moenergies=[]
            self.mocoeffs=[]

            for spin in range(unrestricted + 1): # make sure we cover all instances
                title = inputfile.next()
                while(title[0] == "#"):
                    title = inputfile.next()

#                mocoeffs = numpy.zeros((self.nbasis, self.nbasis), "d")
                moenergies = []
                moarray=[]

                if spin == 1 and title[0:11] == "$uhfmo_beta":
                    title = inputfile.next()
                    while title[0] == "#":
                        title = inputfile.next()

                while(title[0] != '$'):
                    temp=title.split()

                    orb_symm=temp[1]

                    try:
                        energy = float(temp[2][11:].replace("D", "E"))
                    except ValueError:
                        print(spin, ": ", title)

                    orb_en = utils.convertor(energy,"hartree","eV")

                    moenergies.append(orb_en)
                    single_mo = []
                    
                    while(len(single_mo)<self.nbasis):
                        self.updateprogress(inputfile, "Coefficients", self.cupdate)
                        title = inputfile.next()
                        lines_coeffs=self.split_molines(title)
                        single_mo.extend(lines_coeffs)
                        
                    moarray.append(single_mo)
                    title = inputfile.next()

#                for i in range(0, len(moarray), 1):
#                    for j in range(0, self.nbasis, 1):
#                        try:
#                            mocoeffs[i][j]=moarray[i][j]
#                        except IndexError:
#                            print "Index Error in mocoeffs.", spin, i, j
#                            break

                mocoeffs = numpy.array(moarray,"d")
                self.mocoeffs.append(mocoeffs)
                self.moenergies.append(moenergies)

        if line[26:49] == "a o f o r c e - program":
            self.vibirs = []
            self.vibfreqs = []
            self.vibsyms = []
            self.vibdisps = []

#            while line[3:31] != "****  force : all done  ****":

        if line[12:26] == "ATOMIC WEIGHTS":
#begin parsing atomic weights
           self.vibmasses=[]
           line=inputfile.next() # lines =======
           line=inputfile.next() # notes
           line=inputfile.next() # start reading
           temp=line.split()
           while(len(temp) > 0):
                self.vibmasses.append(float(temp[2]))
                line=inputfile.next()
                temp=line.split()

        if line[5:14] == "frequency":
            if not hasattr(self,"vibfreqs"):
                self.vibfreqs = []
                self.vibfreqs = []
                self.vibsyms = []
                self.vibdisps = []
                self.vibirs = []

            temp=line.replace("i","-").split()

            freqs = [self.float(f) for f in temp[1:]]
            self.vibfreqs.extend(freqs)
                    
            line=inputfile.next()
            line=inputfile.next()

            syms=line.split()
            self.vibsyms.extend(syms[1:])

            line=inputfile.next()
            line=inputfile.next()
            line=inputfile.next()
            line=inputfile.next()

            temp=line.split()
            irs = [self.float(f) for f in temp[2:]]
            self.vibirs.extend(irs)

            line=inputfile.next()
            line=inputfile.next()
            line=inputfile.next()
            line=inputfile.next()

            x=[]
            y=[]
            z=[]

            line=inputfile.next()
            while len(line) > 1:
                temp=line.split()
                x.append(map(float, temp[3:]))

                line=inputfile.next()
                temp=line.split()
                y.append(map(float, temp[1:]))

                line=inputfile.next()
                temp=line.split()
                z.append(map(float, temp[1:]))
                line=inputfile.next()

# build xyz vectors for each mode

            for i in range(0, len(x[0]), 1):
                disp=[]
                for j in range(0, len(x), 1):
                    disp.append( [x[j][i], y[j][i], z[j][i]])
                self.vibdisps.append(disp)

#        line=inputfile.next()

    def after_parsing(self):

        # delete all frequencies that correspond to translations or rotations
        
        if hasattr(self,"vibfreqs"):
            i = 0
            while i < len(self.vibfreqs):
                if self.vibfreqs[i]==0.0:
                    del self.vibfreqs[i]
                    del self.vibdisps[i]
                    del self.vibirs[i]
                    del self.vibsyms[i]
                    i -= 1
                i += 1


