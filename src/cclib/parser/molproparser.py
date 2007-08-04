"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 661 $"


import re

# If numpy is not installed, try to import Numeric instead.
try:
    import numpy
except ImportError:
    import Numeric as numpy

import utils
import logfileparser


class Molpro(logfileparser.Logfile):
    """Molpro file parser"""

    def __init__(self,*args):
        # Call the __init__ method of the superclass
        super(Molpro, self).__init__(logname="Molpro",*args)

    def __str__(self):
        """Return a string representation of the object."""
        return "Molpro log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Molpro("%s")' % (self.filename)

    def normalisesym(self, label):
        """Normalise the symmetries used by Molpro."""
        ans = label.replace("`", "'").replace("``", "''")
        return ans

    def before_parsing(self):
        
        self.electronorbitals = ""
        self.insidescf = False

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        if line[1:19] == "ATOMIC COORDINATES":
            
            if not hasattr(self,"atomcoords"):
                self.atomcoords = []
                self.atomnos = []
            line = inputfile.next()
            line = inputfile.next()
            line = inputfile.next()
            atomcoords = []
            atomnos = []
            
            line = inputfile.next()
            while line.strip():
                temp = line.strip().split()
                atomcoords.append([utils.convertor(float(x),"bohr","Angstrom") for x in temp[3:6]]) #bohrs to angs
                atomnos.append(int(round(float(temp[2]))))
                line = inputfile.next()
                
            self.atomnos = numpy.array(atomnos, "i")
            self.atomcoords.append(atomcoords)
            self.natom = len(self.atomnos)
        
        # Use BASIS DATA to parse input for aonames and atombasis.
        # This is always the first place this information is printed, so no attribute check is needed.
        if line[1:11] == "BASIS DATA":
            
            blank = inputfile.next()
            header = inputfile.next()
            blank = inputfile.next()
            self.aonames = []
            self.atombasis = []
            self.gbasis = []
            for i in range(self.natom):
                self.atombasis.append([])
                self.gbasis.append([])
            
            line = "dummy"
            while line.strip() != "":
                line = inputfile.next()
                funcnr = line[1:6]
                funcsym = line[7:9]
                funcatom_ = line[11:14]
                functype_ = line[16:22]
                funcexp = line[25:38]
                funccoeffs = line[38:]

                # If a new function type is printed or the BASIS DATA block ends,
                #   then the previous functions can be added to gbasis.
                # When translating the Molpro function type name into a gbasis code,
                #   note that Molpro prints all components, and we want to add
                #   only one to gbasis, with the proper code (S,P,D,F,G).
                # Warning! The function types differ for cartesian/spherical functions.
                # Skip the first printed function type, however.
                if (functype_.strip() and line[3] != '1') or line.strip() == "":
                    funcbasis = None
                    if functype in ['1s', 's']:
                        funcbasis = 'S'
                    if functype in ['x', '2px']:
                        funcbasis = 'P'
                    if functype in ['xx', '3d0']:
                        funcbasis = 'D'
                    if functype in ['xxx', '4f0']:
                        funcbasis = 'F'
                    if functype in ['xxxx', '5g0']:
                        funcbasis = 'G'
                    if funcbasis:
                        # The function is split into as many columns as there are.
                        for i in range(len(coefficients[0])):
                            func = (funcbasis, [])
                            for j in range(len(exponents)):
                                func[1].append((exponents[j],coefficients[j][i]))
                            self.gbasis[funcatom-1].append(func)

                # If it is a new type, set up the variables for the next shell(s).
                if functype_.strip():
                    exponents = []
                    coefficients = []
                    functype = functype_.strip()
                    funcatom = int(funcatom_.strip())

                # Add exponents nad coefficients to lists.
                if line.strip():
                    funcexp = float(funcexp)
                    funccoeffs = [float(s) for s in funccoeffs.split()]
                    exponents.append(funcexp)
                    coefficients.append(funccoeffs)

                # If the function number is there, add to atombasis and aonames.
                if funcnr.strip():
                    funcnr = int(funcnr.split('.')[0])
                    self.atombasis[funcatom-1].append(funcnr-1)
                    element = self.table.element[self.atomnos[funcatom-1]]
                    aoname = "%s%i_%s" %(element, funcatom, functype)
                    self.aonames.append(aoname)

        if line[1:23] == "NUMBER OF CONTRACTIONS":
            
            nbasis = int(line.split()[3])
            if hasattr(self, "nbasis"):
                assert nbasis == self.nbasis
            else:
                self.nbasis = nbasis

        # This is used to signalize whether we are inside an SCF calculation.
        if line[1:8] == "PROGRAM" and line[14:18] == "-SCF":

            self.insidescf = True

        # Use this information instead of 'SETTING ...', in case the defaults are standard.
        # Note that this is sometimes printed in each geometry optimization step.
        if line[1:20] == "NUMBER OF ELECTRONS":
            
            spinup = int(line.split()[3][:-1])
            spindown = int(line.split()[4][:-1])
            # Nuclear charges (atomnos) should be parsed by now.
            nuclear = numpy.sum(self.atomnos)
            charge = nuclear - spinup - spindown
            mult = spinup - spindown + 1
            
            # Copy charge, or assert for exceptions if already exists.
            if not hasattr(self, "charge"):
                self.charge = charge
            else:
                assert self.charge == charge
            
            # Copy multiplicity, or assert for exceptions if already exists.
            if not hasattr(self, "mult"):
                self.mult = mult
            else:
                assert self.mult == mult
        
        # Convergenve thresholds for SCF cycle, should be contained in a line such as:
        #   CONVERGENCE THRESHOLDS:    1.00E-05 (Density)    1.40E-07 (Energy)
        if self.insidescf and line[1:24] == "CONVERGENCE THRESHOLDS:":

            if not hasattr(self, "scftargets"):
                self.scftargets = []

            scftargets = map(float, line.split()[2::2])
            self.scftargets.append(scftargets)
            # Usually two criteria, but save the names this just in case.
            self.scftargetnames = line.split()[3::2]

        # Read in the print out of the SCF cycle - for scfvalues. For RHF looks like:
        # ITERATION    DDIFF          GRAD             ENERGY        2-EL.EN.            DIPOLE MOMENTS         DIIS
        #     1      0.000D+00      0.000D+00      -379.71523700   1159.621171   0.000000   0.000000   0.000000    0
        #     2      0.000D+00      0.898D-02      -379.74469736   1162.389787   0.000000   0.000000   0.000000    1
        #     3      0.817D-02      0.144D-02      -379.74635529   1162.041033   0.000000   0.000000   0.000000    2
        #     4      0.213D-02      0.571D-03      -379.74658063   1162.159929   0.000000   0.000000   0.000000    3
        #     5      0.799D-03      0.166D-03      -379.74660889   1162.144256   0.000000   0.000000   0.000000    4
        if self.insidescf and line[1:10] == "ITERATION":
        
            if not hasattr(self, "scfvalues"):
                self.scfvalues = []
        
            line = inputfile.next()
            energy = 0.0
            scfvalues = []
            while line.strip() != "":
                if line.split()[0].isdigit():
                
                    ddiff = float(line.split()[1].replace('D','E'))
                    newenergy = float(line.split()[3])
                    ediff = newenergy - energy
                    energy = newenergy

                    # The convergence thresholds must have been read above.
                    # Presently, we recognize MAX DENSITY and MAX ENERGY thresholds.
                    numtargets = len(self.scftargetnames)
                    values = [numpy.nan]*numtargets
                    for n,name in zip(range(numtargets),self.scftargetnames):
                        if "ENERGY" in name.upper():
                            values[n] = ediff
                        elif "DENSITY" in name.upper():
                            values[n] = ddiff
                    scfvalues.append(values)

                line = inputfile.next()
            self.scfvalues.append(numpy.array(scfvalues))

        # SCF result - RHF/UHF and DFT (RKS) energies.
        if line[1:5] in ["!RHF", "!UHF", "!RKS"] and line[16:22] == "ENERGY":
            
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            scfenergy = float(line.split()[4])
            self.scfenergies.append(utils.convertor(scfenergy, "hartree", "eV"))
            
            # We are now done with SCF cycle (after a few lines).
            self.insidescf = False

        # MP2 energies.
        if line[1:5] == "!MP2":
        
            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []
            mp2energy = float(line.split()[-1])
            mp2energy = utils.convertor(mp2energy, "hartree", "eV")
            self.mpenergies.append([mp2energy])
            
        # MP2 energies if MP3 or MP4 is also calculated.
        if line[1:5] == "MP2:":
        
            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []
            mp2energy = float(line.split()[2])
            mp2energy = utils.convertor(mp2energy, "hartree", "eV")
            self.mpenergies.append([mp2energy])
            
        # MP3 (D) and MP4 (DQ or SDQ) energies.
        if line[1:8] == "MP3(D):":
        
            mp3energy = float(line.split()[2])
            mp2energy = utils.convertor(mp3energy, "hartree", "eV")
            line = inputfile.next()
            self.mpenergies[-1].append(mp2energy)
            if line[1:9] == "MP4(DQ):":
                mp4energy = float(line.split()[2])
                line = inputfile.next()
                if line[1:10] == "MP4(SDQ):":
                    mp4energy = float(line.split()[2])
                mp4energy = utils.convertor(mp4energy, "hartree", "eV")
                self.mpenergies[-1].append(mp4energy)

        # The CCSD program operates all closed-shel coupled cluster runs.
        if line[1:15] == "PROGRAM * CCSD":
        
            if not hasattr(self, "ccenergies"):
                self.ccenergies = []
            while line[1:20] != "Program statistics:":
                # The last energy (most exact) will be read last and thus saved.
                if line[1:5] == "!CCD" or line[1:6] == "!CCSD" or line[1:9] == "!CCSD(T)":
                    ccenergy = float(line.split()[-1])
                    ccenergy = utils.convertor(ccenergy, "hartree", "eV")
                line = inputfile.next()
            self.ccenergies.append(ccenergy)

        # Read the occupancy (index of HOMO s).
        # For restricted calculations, there is one line here. For unrestricted, two:
        #   Final alpha occupancy:  ...
        #   Final beta  occupancy:  ...
        if line[1:17] == "Final occupancy:":
            self.homos = [int(line.split()[-1])-1]
        if line[1:23] == "Final alpha occupancy:":
            self.homos = [int(line.split()[-1])-1]
            line = inputfile.next()
            self.homos.append(int(line.split()[-1])-1)

        # From this block atombasis, moenergies, and mocoeffs can be parsed.
        # Note that Molpro does not print this by default, you must add this in the input:
        #   GPRINT,ORBITALS
        # What's more, this prints only the occupied orbitals. To get virtuals, add also:
        #   ORBPTIN,NVIRT
        #   where NVIRT is how many to print (can be some large number, like 99999, to print all).
        # The block is in general flipped when compared to other programs (GAMESS, Gaussian), and
        #   MOs in the rows. Also, it does not cut the table into parts, rather each MO row has
        #   as many lines as it takes to print all the coefficients, as shown below:
        #
        # ELECTRON ORBITALS
        # =================
        #
        #
        #   Orb  Occ    Energy  Couls-En    Coefficients
        #
        #                                   1 1s      1 1s      1 2px     1 2py     1 2pz     2 1s   (...)
        #                                   3 1s      3 1s      3 2px     3 2py     3 2pz     4 1s   (...)
        # (...)
        #
        #   1.1   2   -11.0351  -43.4915  0.701460  0.025696 -0.000365 -0.000006  0.000000  0.006922 (...)
        #                                -0.006450  0.004742 -0.001028 -0.002955  0.000000 -0.701460 (...)
        # (...)
        #
        # For unrestricted calcualtions, ELECTRON ORBITALS is followed on the same line
        #   by FOR POSITIVE SPIN or FOR NEGATIVE SPIN.
        # For examples, see data/Molpro/basicMolpro2006/dvb_sp*.
        if line[1:18] == "ELECTRON ORBITALS" or self.electronorbitals:
            
            # Detect if we are reading beta (negative spin) orbitals.
            spin = 0
            if line[19:36] == "FOR NEGATIVE SPIN" or self.electronorbitals[19:36] == "FOR NEGATIVE SPIN":
                spin = 1
            
            if not self.electronorbitals:
                dashes = inputfile.next()
            blank = inputfile.next()
            blank = inputfile.next()
            headers = inputfile.next()
            blank = inputfile.next()
            
            # Parse the list of atomic orbitals if atombasis or aonames is missing.
            line = inputfile.next()
            if not hasattr(self, "atombasis") or not hasattr(self, "aonames"):
                self.atombasis = []
                for i in range(self.natom):
                    self.atombasis.append([])
                self.aonames = []
                aonum = 0
                while line.strip() != "":
                    for s in line.split():
                        if s.isdigit():
                            atomno = int(s)
                            self.atombasis[atomno-1].append(aonum)
                            aonum += 1
                        else:
                            functype = s
                            element = self.table.element[self.atomnos[atomno-1]]
                            aoname = "%s%i_%s" %(element, atomno, functype)
                            self.aonames.append(aoname)
                    line = inputfile.next()
            else:
                while line.strip() != "":
                    line = inputfile.next()
            line = inputfile.next()
            
            # The energies and coefficients for each MO follow.
            if not hasattr(self, "moenergies"):
                self.moenergies = [[]]
            else:
                self.moenergies.append([])
            if not hasattr(self, "mocoeffs"):
                self.mocoeffs = [[]]
            else:
                self.mocoeffs.append([])
            line = inputfile.next()
            while len(line.split()) >= 10:
                coeffs = []
                while line.strip() != "":
                    if len(line.split()) > 10:
                        moenergy = float(line.split()[2])
                        moenergy = utils.convertor(moenergy, "hartree", "eV")
                        self.moenergies[spin].append(moenergy)
                    coeffs += [float(s) for s in line[31:].split()]
                    line = inputfile.next()
                self.mocoeffs[spin].append(coeffs)
                line = inputfile.next()
            
            # Check if last line begins the next ELECTRON ORBITALS section.
            if line[1:18] == "ELECTRON ORBITALS":
                self.electronorbitals = line
            else:
                self.electronorbitals = ""

        # Is the MATROP program was called appropriately,
        #   the atomic obital overlap matrix S is printed.
        # The matrix is printed straight-out, ten elements in each row, both halves.
        # Note that is the entire matrix is not printed, then aooverlaps
        #   will not have dimensions nbasis x nbasis.
        if line[1:9] == "MATRIX S":
        
            blank = inputfile.next()
            symblocklabel = inputfile.next()
            if not hasattr(self, "aooverlaps"):
                self.aooverlaps = [[]]
            line = inputfile.next()
            while line.strip() != "":
                elements = [float(s) for s in line.split()]
                if len(self.aooverlaps[-1]) + len(elements) <= self.nbasis:
                    self.aooverlaps[-1] += elements
                else:
                    n = len(self.aooverlaps[-1]) + len(elements) - self.nbasis
                    self.aooverlaps[-1] += elements[:-n]
                    self.aooverlaps.append([])
                    self.aooverlaps[-1] += elements[-n:]
                line = inputfile.next()

        if line[1:13] == "Normal Modes" and not hasattr(self,"vibfreqs"):
            
            self.vibfreqs = []
            self.vibirs = []
            self.vibsyms = []
            line = inputfile.next()
            line = inputfile.next()
            
            for i in line.split():
                try: map(int, i)
                except: self.vibsyms.append(i)
                
            while line.strip():
                if line[1:12] == "Wavenumbers":
                    freqs = line.strip().split()[2:]
                    self.vibfreqs.extend(map(float,freqs))
                if line[1:21] == "Intensities [km/mol]":
                    irIntensity = line.strip().split()[2:]
                    self.vibirs.extend(map(float,irIntensity))
                line = inputfile.next()
                if not line.strip():
                    line = inputfile.next()
                    if not line.strip(): break
                    else: self.vibsyms += self.normalisesym(line.split()[1::2])
                    
        if line[1:20] == "Normal Modes of low" and hasattr(self,"vibfreqs"):
            
            line = inputfile.next()
            line = inputfile.next()
            self.lowfreqs = []
            self.lowirs = []
            
            while line.strip():
                if line[1:12] == "Wavenumbers":
                    freqs = line.strip().split()[2:]
                    self.lowfreqs.extend(map(float,freqs))
                if line[1:21] == "Intensities [km/mol]":
                    irIntensity = line.strip().split()[2:]
                    self.lowirs.extend(map(float,irIntensity))
                line = inputfile.next()
                if not line.strip():
                    line = inputfile.next()
                    if not line.strip(): break
                    else: continue
            self.vibfreqs = self.lowfreqs + self.vibfreqs
            self.vibirs = self.lowirs + self.vibirs

        if line[1:16] == "Force Constants":
            
            self.logger.info("Creating attribute hessian")
            self.hessian = []
            line = inputfile.next()
            hess = []
            tmp = []
            
            while line.strip():
                try: map(float, line.strip().split()[2:])
                except: 
                    line = inputfile.next()
                line.strip().split()[1:]
                hess.extend([map(float,line.strip().split()[1:])])
                line = inputfile.next()
            lig = 0
            
            while (lig==0) or (len(hess[0]) > 1):
                tmp.append(hess.pop(0))
                lig += 1
            k = 5
            
            while len(hess) != 0:
                tmp[k] += hess.pop(0)
                k += 1
                if (len(tmp[k-1]) == lig): break
                if k >= lig: k = len(tmp[-1])
            for l in tmp: self.hessian += l
            
        if line[1:14] == "Atomic Masses" and hasattr(self,"hessian"):
            
            line = inputfile.next()
            self.amass = map(float, line.strip().split()[2:])
            
            while line.strip():
                line = inputfile.next()
                self.amass += map(float, line.strip().split()[2:])


if __name__ == "__main__":
    import doctest, molproparser
    doctest.testmod(molproparser, verbose=False)
