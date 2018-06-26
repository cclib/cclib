"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 635 $"

import re

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

import utils
import logfileparser

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
    """A Turbomole output file"""

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

#EOF

