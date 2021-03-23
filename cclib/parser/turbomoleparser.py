# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Turbomole output files."""

import re

import numpy

from cclib.parser import logfileparser
from cclib.parser import utils
from cclib.parser import data

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
        """Normalise the symmetries used by Turbomole.

        The labels are standardized except for the first character being lowercase.
        """
        # TODO more work could be required, but we don't have any logfiles
        # with non-C1 symmetry.
        return label.capitalize()

    def before_parsing(self):
        self.geoopt = False # Is this a GeoOpt? Needed for SCF targets/values.
        self.periodic_table = utils.PeriodicTable()
        self.new_module()

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
        
    def new_module(self):
        """
        This method is called when we start parsing a new module (Turbomole subprogram).
        """
        self.DFT = False
        self.metadata['success'] = False

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
            #self.set_attribute('nmo', nmo)
            
        # The DFT functional.
        # This information is printed by dscf but not in an easily parsable format, so we'll take it from the control file instead...
        # Additionally, turbomole stores functional names in lower case. This looks odd, so we'll convert to uppercase (?)
        if line[3:13] == "functional":
            # $dft
            #    functional b-p
            self.metadata['functional'] = line.split()[1].upper()
            self.DFT = True
        
        if "density functional" in line:
            self.DFT = True
        
        # Extract the version number and optionally the build number.
        searchstr = ": TURBOMOLE"
        index = line.find(searchstr)
        if index > -1:
            # We have entered a new module.
            self.new_module()
            
            line = line[index + len(searchstr):]
            tokens = line.split()
            # This line could look like any of the following:
            # - TURBOMOLE rev. V7.4.1 (2bfdd732)
            # - TURBOMOLE V7.2 ( 21471 ) 11 Oct 2017 at 17:04:51
            # - TURBOMOLE V5-9-0 29 Nov 2006 at 22:06:41
            
            # Determine where the version string is.
            if tokens[0] == "rev.":
                version_index = 1
            else:
                version_index = 0

            package_version = tokens[version_index][1:].replace("-", ".")
            self.metadata["package_version"] = package_version
            self.metadata["legacy_package_version"] = package_version
            
            # The revision may or may not contain whitespace inside the brackets.
            if "(" == tokens[version_index +1][0] and ")" == tokens[version_index +1][-1]:
                revision = tokens[version_index +1][1:-1]
                self.metadata["package_version"] = "{}.r{}".format(package_version, revision)
                
            elif tokens[version_index +1] == "(":
                revision = tokens[version_index +2]
                self.metadata["package_version"] = "{}.r{}".format(package_version, revision)
                
            
        if ": all done  ****" in line:
            # End of module, set success flag.
            self.metadata['success'] = True
            
        ## Basis set info from dscf.
        #               +--------------------------------------------------+
        #               |               basis set information              |
        #               +--------------------------------------------------+
        # 
        #               we will work with the 1s 3p 5d 7f 9g ... basis set
        #               ...i.e. with spherical basis functions...
        # 
        #    type   atoms  prim   cont   basis
        #    ---------------------------------------------------------------------------
        #     o        1     15      5   sto-3g hondo  [2s1p|6s3p]
        #     h        2      3      1   sto-3g hondo  [1s|3s]
        #    ---------------------------------------------------------------------------
        if "type   atoms  prim   cont   basis" in line:
            line = next(inputfile)
            line = next(inputfile)
            basis_sets = []
            while "---------------------------------------------------------------------------" not in line:
                basis_sets.append(" ".join(line.split()[4:-1]))
                line = next(inputfile)
            
            # Turbomole gives us the basis set for each atom, but we're only interested if the same basis set is used throughout (for now).
            if len(set(basis_sets)) == 1:
                self.metadata["basis_set"] = list(set(basis_sets))[0]
        
        # Molecular charge info from dscf.
        #               nuc           elec       ->  total
        #  ------------------------------------------------------------------------------
        #                           charge      
        #  ------------------------------------------------------------------------------
        #           43.000000     -42.000000       1.000000
        if "nuc           elec       ->  total" in line:
            line = next(inputfile)
            line = next(inputfile)
            if "charge" in line:
                line = next(inputfile)
                line = next(inputfile)
                
                total_charge = float(line.split()[2])
                total_charge_int = round(total_charge)
                
                # Check we wont loose information converting to int.
                if total_charge != total_charge_int:
                    self.logger.warning("Converting non integer total charge '{}' to integer".format(total_charge))
                
                # Set regardless.
                self.set_attribute("charge", total_charge_int)
                
        
        ## Orbital occupation info from dscf.
        #  orbitals $scfmo  will be written to file mos
        # 
        #     irrep                  1a          2a          3a          4a          5a   
        #  eigenvalues H        -20.25992    -1.24314    -0.57053    -0.46144    -0.39295
        #             eV        -551.3047    -33.8279    -15.5250    -12.5564    -10.6929
        #  occupation              2.0000      2.0000      2.0000      2.0000      2.0000
        # 
        #     irrep                  6a          7a   
        #  eigenvalues H          0.55091     0.64409
        #             eV          14.9910     17.5268
        ## Or
        #  orbitals $uhfmo_beta  will be written to file beta
        # 
        #  orbitals $uhfmo_alpha  will be written to file alpha
        #  
        #  alpha: 
        # 
        #     irrep                 31a         32a         33a         34a         35a   
        #  eigenvalues H         -0.47570    -0.46573    -0.40741    -0.39213    -0.35411
        #             eV         -12.9446    -12.6733    -11.0862    -10.6705     -9.6358
        #  occupation              1.0000      1.0000      1.0000      1.0000      1.0000 
        # 
        #     irrep                 36a         37a         38a         39a         40a   
        #  eigenvalues H         -0.18634    -0.10035    -0.09666    -0.02740     0.06072
        #             eV          -5.0705     -2.7306     -2.6303     -0.7455      1.6522
        #  
        #  beta:  
        # 
        #     irrep                 30a         31a         32a         33a         34a   
        #  eigenvalues H         -0.49118    -0.47348    -0.44470    -0.39020    -0.37919
        #             eV         -13.3658    -12.8842    -12.1009    -10.6181    -10.3184
        #  occupation              1.0000      1.0000      1.0000      1.0000      1.0000 
        # 
        #     irrep                 35a         36a         37a         38a         39a   
        #  eigenvalues H         -0.28091    -0.15088    -0.09343    -0.07531    -0.00688
        #             eV          -7.6440     -4.1058     -2.5424     -2.0493     -0.1873
        if "will be written to file mos" in line and hasattr(self, "mosyms"):
            orbitals, line = self.parse_dscf_orbitals(inputfile, line)
            self.set_attribute("homos", [self.determine_homo(self.mosyms[0], orbitals)])
        
        if "alpha:" in line and hasattr(self, "mosyms"):
            orbitals, line = self.parse_dscf_orbitals(inputfile, line)
            homo = self.determine_homo(self.mosyms[0], orbitals)
            if not hasattr(self, "homos"):
                self.set_attribute('homos', [homo])
            else:
                self.homos[0] = homo
            
        if "beta:" in line and hasattr(self, "mosyms"):
            orbitals, line = self.parse_dscf_orbitals(inputfile, line)
            homo = self.determine_homo(self.mosyms[1], orbitals)
            if not hasattr(self, "homos"):
                self.set_attribute('homos', [homo])
            elif len(self.homos) == 1:
                self.homos.append(homo)
            else:
                self.homos[1] = homo
        
        # Coordinates and gradients from statpt.
        #   *************************************************************************
        #  ATOM                      CARTESIAN COORDINATES
        #   1 c      -2.67642286424381      0.00038527796998     -0.44566112589039
        #   2 c      -1.69183504162310      0.00075591605173      2.05352357416443
        #   3 c       0.92311977253458      0.00052122921704      2.48381308900370
        #   4 c       2.67642286424491      0.00038528140768      0.44566112589851
        #   5 c       1.69183504161656      0.00075592009652     -2.05352357415432
        #   6 h      -2.98983238801665      0.00149374762706      3.67174843088014
        #   7 h       1.64279730150941      0.00038458181684      4.43158185858240
        #   8 h       2.98983238800216      0.00149376082387     -3.67174843087059
        #   9 c       5.44975417206469     -0.00039526372012      1.01184691031725
        #  10 c       7.34299179214000     -0.00071986894317     -0.68271530533475
        #  11 h       7.01332300460381     -0.00029648805084     -2.72785997302809
        #  12 h       5.91422629512709     -0.00052260839470      3.03917498747826
        #  13 h       9.32456490822245     -0.00139866158379     -0.07769877084777
        #  14 c      -5.44975417205833     -0.00039526324514     -1.01184691032156
        #  15 h      -5.91422629512079     -0.00052260600240     -3.03917498748360
        #  16 c      -7.34299179213657     -0.00071986698481      0.68271530532303
        #  17 h      -7.01332300460557     -0.00029648750734      2.72785997301777
        #  18 h      -9.32456490821458     -0.00139865699894      0.07769877082649
        #  19 c      -0.92311977253569      0.00052122610063     -2.48381308899233
        #  20 h      -1.64279730151053      0.00038457006746     -4.43158185856859
        # *************************************************************************
        #  ATOM                      CARTESIAN GRADIENTS  
        #   1 c       0.00001122250995      0.00000667056223      0.00001828623085
        #   2 c       0.00000193381649     -0.00001184358472      0.00001337145426
        #   3 c      -0.00000578603102      0.00000187474419      0.00000987215044
        #   4 c      -0.00001122250794      0.00000667041002     -0.00001828623169
        #   5 c      -0.00000193381900     -0.00001184350870     -0.00001337145668
        #   6 h      -0.00001402399615      0.00001821929428     -0.00000611193875
        #   7 h       0.00001742260704     -0.00000462667818     -0.00000195847724
        #   8 h       0.00001402399610      0.00001821957100      0.00000611193919
        #   9 c      -0.00006557131494     -0.00002136001589      0.00002475454295
        #  10 c       0.00006109842779     -0.00000007565809     -0.00004223963891
        #  11 h      -0.00002879535900      0.00000574108155     -0.00000337279221
        #  12 h       0.00004754086240      0.00000863480152     -0.00000658955451
        #  13 h      -0.00000914653785     -0.00000320753859      0.00003066647865
        #  14 c       0.00006557131503     -0.00002135978825     -0.00002475454202
        #  15 h      -0.00004754086227      0.00000863481272      0.00000658955373
        #  16 c      -0.00006109842855     -0.00000007566316      0.00004223963942
        #  17 h       0.00002879535887      0.00000574109116      0.00000337279251
        #  18 h       0.00000914653868     -0.00000320756462     -0.00003066647871
        #  19 c       0.00000578603127      0.00000187466275     -0.00000987214898
        #  20 h      -0.00001742260682     -0.00000462696679      0.00000195847629
        # *************************************************************************
        #
        # In older versions (5.9)
        # ATOM               CARTESIAN GRADIENTS
        #  1 c      -0.00006249573129     -0.00001452041025     -0.00002380577094
        #  2 c       0.00007390978833     -0.00000084052225      0.00000584172986
        #  3 c      -0.00000362984325      0.00000466314052     -0.00002733252443
        #  4 c       0.00006249783071     -0.00001452158676      0.00002380177091
        #  5 c      -0.00007390342509     -0.00000084270316     -0.00000584066725
        #  6 h      -0.00005950589946      0.00000095087455      0.00005137853319
        #  7 h      -0.00000645073463     -0.00000179899463      0.00003328307243
        #  8 h       0.00005950123260      0.00000095090879     -0.00005137678599
        #  9 c      -0.00004885046058     -0.00002593955699     -0.00004887642413
        # 10 c       0.00000015222157      0.00004027444288      0.00003204866033
        # 11 h       0.00001239310045      0.00001308636889     -0.00000092060218
        # 12 h       0.00001561988506     -0.00001608132438      0.00002717833564
        # 13 h      -0.00000616830304      0.00000016851704     -0.00003103539353
        # 14 c       0.00004884621752     -0.00002593547764      0.00004887866674
        # 15 h      -0.00001561780942     -0.00001607888272     -0.00002718198410
        # 16 c      -0.00000015329549      0.00004026946593     -0.00003204885735
        # 17 h      -0.00001239289077      0.00001308531162      0.00000092210161
        # 18 h       0.00000616838338      0.00000016858902      0.00003103523777
        # 19 c       0.00000362187006      0.00000466096401      0.00002733645543
        # 20 h       0.00000645786975     -0.00000179869042     -0.00003328556910
        #     Optint: norm of internal gradient       0.00029710
        #
        if "CARTESIAN GRADIENTS" in line:
            atomcoords = []
            atomnos = []
            line = next(inputfile)
            while len(line.split()) == 5:
                atomnos.append(self.periodic_table.number[line.split()[-4].capitalize()])
                atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom") 
                                   for x in line.split()[-3:]])
                line = next(inputfile)
            
            self.set_attribute('atomnos', atomnos)
            self.set_attribute('natom', len(atomcoords))
            self.append_attribute("grads", atomcoords)
                
        
        
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
        #
        # During a jobex driven optimisation, the above geometry section will be printed at the start of each module (subprogram).
        # We only want one coord entry per step, so we need to be careful not to add the same geometry multiple times.
        if 'Atomic coordinate, charge and isotop information' in line:
            while 'atomic coordinates' not in line:
                line = next(inputfile)

            atomcoords = []
            atomnos = []
            line = next(inputfile)
            while len(line.split()) >= 3:
                atomnos.append(self.periodic_table.number[line.split()[3].capitalize()])
                atomcoords.append([utils.convertor(float(x), "bohr", "Angstrom") 
                                   for x in line.split()[:3]])
                line = next(inputfile)
            
            # Check the coordinates we just parsed are not already in atomcoords.
            if not hasattr(self, 'atomcoords') or self.atomcoords[-1] != atomcoords:
                self.append_attribute('atomcoords', atomcoords)
                self.set_attribute('atomnos', atomnos)
                self.set_attribute('natom', len(atomcoords))
        
        # Flag that indicates we are doing an opt.
        if "OPTIMIZATION CYCLE" in line:
            self.append_attribute("optstatus", data.ccData.OPT_UNKNOWN)
            
            if "OPTIMIZATION CYCLE 1" in line:
                # This is the start of the opt.
                self.optstatus[-1] += data.ccData.OPT_NEW
        
        # Optimisation convergence criteria using statpt.
        #
        # ****************************************************************** 
        #                     CONVERGENCE INFORMATION
        # 
        #                          Converged?     Value      Criterion
        #        Energy change         no       0.0000011   0.0000010
        #        RMS of displacement   yes      0.0001152   0.0005000
        #        RMS of gradient       yes      0.0000548   0.0005000
        #        MAX displacement      yes      0.0001409   0.0010000
        #        MAX gradient          yes      0.0000670   0.0010000
        # ****************************************************************** 
        if "CONVERGENCE INFORMATION" in line:
            # This is an optimisation.
            # Skip lines.
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            
            convergence = []
            geovalues = []
            geotargets = []
            
            # There are a variable number of criteria.
            while len(line.split()) > 3:
                parts = line.split()
                # lower() is for (unnecessary?) future proofing...
                converged = parts[-3].lower() == "yes"
                value = float(parts[-2])
                criterion = float(parts[-1])
                
                # TODO: possibly require some unit conversions?
                convergence.append(converged)
                geovalues.append(value)
                geotargets.append(criterion)
                
                # Next.                
                line = next(inputfile)
                
            self.set_attribute("geotargets", geotargets)
            self.append_attribute("geovalues", geovalues)
            
            if all(convergence):
                # This iteration has converged.
                self.append_attribute("optdone", len(self.geovalues) -1)
                self.optstatus[-1] += data.ccData.OPT_DONE
            else:
                # Not converged.
                if not hasattr(self, 'optdone'):
                    self.set_attribute("optdone", [])
                #self.optstatus[-1] += data.ccData.OPT_UNCONVERGED
                

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
        # ....
        #
        # reduced mass(g/mol)     3.315    2.518    2.061    3.358    3.191    2.323

        if 'NORMAL MODES and VIBRATIONAL FREQUENCIES (cm**(-1))' in line:
            has_raman = False
            while line[7:11] != 'mode':
                line = next(inputfile)
                if line.startswith(" differential RAMAN cross sections"):
                    has_raman = True
            vibfreqs, vibsyms, vibirs, vibdisps, vibrmasses = [], [], [], [], []
            while 'all done  ****' not in line:

                if line.strip().startswith('mode'):
                    self.skip_line(inputfile, 'b')
                    line = next(inputfile)
                    assert line.strip().startswith('frequency')
                    freqs = [float(i.replace('i', '-')) for i in line.split()[1:]]
                    vibfreqs.extend(freqs)
                    self.skip_lines(inputfile, ['b'])
                    line = next(inputfile)
                    assert line.strip().startswith('symmetry')
                    syms = line.split()[1:]
                    vibsyms.extend(syms)

                    self.skip_lines(inputfile, ['b', 'IR', 'dDIP/dQ'])
                    line = next(inputfile)
                    assert line.strip().startswith('intensity (km/mol)')
                    irs = [utils.float(f) for f in line.split()[2:]]
                    vibirs.extend(irs)

                    self.skip_lines(inputfile, ['intensity %', 'b', 'RAMAN'])
                    if has_raman:
                        self.skip_lines(
                            inputfile,
                            ['(par,par)', '(ort,ort)', '(ort,unpol)', 'depol. ratio']
                        )
                    line = next(inputfile)
                    assert not line.strip()
                    line = next(inputfile)
                    x, y, z = [], [], []
                    atomcounter = 0
                    while atomcounter < self.natom:
                        x.append([float(i) for i in line.split()[3:]])
                        line = next(inputfile)
                        y.append([float(i) for i in line.split()[1:]])
                        line = next(inputfile)
                        z.append([float(i) for i in line.split()[1:]])
                        line = next(inputfile)
                        atomcounter += 1

                    for j in range(len(x[0])):
                        disps = []
                        for i in range(len(x)):
                            disps.append([x[i][j], y[i][j], z[i][j]])
                        vibdisps.append(disps)

                    line = next(inputfile)
                    assert line.startswith('reduced mass(g/mol)')
                    rmasses = [utils.float(f) for f in line.split()[2:]]
                    vibrmasses.extend(rmasses)

                line = next(inputfile)

            self.set_attribute('vibfreqs', vibfreqs)
            self.set_attribute('vibsyms', vibsyms)
            self.set_attribute('vibirs', vibirs)
            self.set_attribute('vibdisps', vibdisps)
            self.set_attribute('vibrmasses', vibrmasses)

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
            mosyms = []

            while not line.strip().startswith('$'):
                number, sym = line.split()[:2]
                info = re.match(r".*eigenvalue=(?P<moenergy>[0-9D\.+-]{20})\s+nsaos=(?P<count>\d+).*", line)
                eigenvalue = utils.float(info.group('moenergy'))
                orbital_energy = utils.convertor(eigenvalue, 'hartree', 'eV')
                moenergies.append(orbital_energy)
                mosyms.append(self.normalisesym(sym))
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
            
            # We now need to sort our orbitals (because Turbomole groups them by symm).
            mos = list(zip(moenergies, mocoeffs, mosyms))
            mos.sort(key = lambda mo: mo[0])
            moenergies, mocoeffs, mosyms = zip(*mos)
            
            self.append_attribute("moenergies", moenergies)
            self.append_attribute("mocoeffs", mocoeffs)
            self.append_attribute("mosyms", mosyms)
            self.set_attribute("nmo", len(moenergies))

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
            total_energy_threshold = utils.float(line.split()[-1])
            one_electron_energy_threshold = utils.float(next(inputfile).split()[-1])
            scftargets = [total_energy_threshold, one_electron_energy_threshold]
            self.append_attribute('scftargets', scftargets)
            iter_energy = []
            iter_one_elec_energy = []
            while 'convergence criteria satisfied' not in line:
                # nbasis
                if  "number of basis functions" in line:
                    self.set_attribute("nbasis", int(line.split()[-1]))
                
                if 'ITERATION  ENERGY' in line:
                    line = next(inputfile)
                    info = line.split()
                    iter_energy.append(utils.float(info[1]))
                    iter_one_elec_energy.append(utils.float(info[2]))
                line = next(inputfile)

            assert len(iter_energy) == len(iter_one_elec_energy), \
                'Different number of values found for total energy and one electron energy.'
            scfvalues = [[x - y, a - b] for x, y, a, b in 
                         zip(iter_energy[1:], iter_energy[:-1], iter_one_elec_energy[1:], iter_one_elec_energy[:-1])]
            self.append_attribute('scfvalues', scfvalues)
            while 'total energy' not in line:
                line = next(inputfile)

            scfenergy = utils.convertor(utils.float(line.split()[4]), 'hartree', 'eV')
            self.append_attribute('scfenergies', scfenergy)
            
            # We need to determine whether this is a HF or DFT energy for metadata.
            if self.DFT:
                self.metadata['methods'].append("DFT")
            else:
                self.metadata['methods'].append("HF")

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
        # Look for MP energies.
        for mp_level in range(2,6):
            if "Final MP{} energy".format(mp_level) in line:
                mpenergy = utils.convertor(utils.float(line.split()[5]), 'hartree', 'eV')
                if mp_level == 2:
                    self.append_attribute('mpenergies', [mpenergy])
                else:
                    self.mpenergies[-1].append(mpenergy)
                self.metadata['methods'].append("MP{}".format(mp_level))

        if 'Final CCSD energy' in line or 'Final CC2 energy' in line:
            self.append_attribute(
                'ccenergies',
                utils.convertor(utils.float(line.split()[5]), 'hartree', 'eV')
            )
            self.metadata['methods'].append("CCSD")

                #line = next(inputfile)

        #  *****************************************************
        #  *                                                   *
        #  *      SCF-energy   :     -74.49827196840999        *
        #  *      MP2-energy   :      -0.19254365976227        *
        #  *      total        :     -74.69081562817226        *
        #  *                                                   *
        #  *     (MP2-energy evaluated from T2 amplitudes)     *
        #  *                                                   *
        #  *****************************************************
        if 'MP2-energy' in line:
            line = next(inputfile)
            if 'total' in line:
                mp2energy = [utils.convertor(utils.float(line.split()[3]), 'hartree', 'eV')]
                self.append_attribute('mpenergies', mp2energy)
                self.metadata['methods'].append("MP2")
                
        # Support for (the now outdated?) rimp2
        # ------------------------------------------------
        #     Method          :  MP2     
        #     Total Energy    :    -75.0009789796
        # ------------------------------------------------
        if "Method          :  MP2" in line:
            line = next(inputfile)
            mp2energy = [utils.convertor(utils.float(line.split()[3]), 'hartree', 'eV')]
            self.append_attribute('mpenergies', mp2energy)
            self.metadata['methods'].append("MP2")
    
    def split_irrep(self, irrep):
        """
        Split a Turbomole irrep into number and symmetry.
        """
        rematch = re.match(r"^([0-9]+)(.+)$", irrep)
        number = int(rematch.group(1))
        sym = self.normalisesym(rematch.group(2))
        return (number, sym)
    
    def parse_dscf_orbitals(self, inputfile, line):
        """
        Extract orbital occupation and energies from a dscf logfile.
        
        Returns
        -------
        tuple
            a two membered tuple where the first element is a list of dictionaries of the the orbitals parsed, while the second is the line on which parsing should continue.
        """
        ## Orbital occupation info from dscf.
        #  orbitals $scfmo  will be written to file mos
        # 
        #     irrep                  1a          2a          3a          4a          5a   
        #  eigenvalues H        -20.25992    -1.24314    -0.57053    -0.46144    -0.39295
        #             eV        -551.3047    -33.8279    -15.5250    -12.5564    -10.6929
        #  occupation              2.0000      2.0000      2.0000      2.0000      2.0000
        # ...
        #     irrep                  6a          7a   
        #  eigenvalues H          0.55091     0.64409
        #             eV          14.9910     17.5268
        ## Or
        #  orbitals $uhfmo_beta  will be written to file beta
        # 
        #  orbitals $uhfmo_alpha  will be written to file alpha
        #  
        #  alpha: 
        # 
        #     irrep                 31a         32a         33a         34a         35a   
        #  eigenvalues H         -0.47570    -0.46573    -0.40741    -0.39213    -0.35411
        #             eV         -12.9446    -12.6733    -11.0862    -10.6705     -9.6358
        #  occupation              1.0000      1.0000      1.0000      1.0000      1.0000 
        # ...
        #     irrep                 36a         37a         38a         39a         40a   
        #  eigenvalues H         -0.18634    -0.10035    -0.09666    -0.02740     0.06072
        #             eV          -5.0705     -2.7306     -2.6303     -0.7455      1.6522
        #  
        #  beta:  
        # 
        #     irrep                 30a         31a         32a         33a         34a   
        #  eigenvalues H         -0.49118    -0.47348    -0.44470    -0.39020    -0.37919
        #             eV         -13.3658    -12.8842    -12.1009    -10.6181    -10.3184
        #  occupation              1.0000      1.0000      1.0000      1.0000      1.0000 
        # ...
        #     irrep                 35a         36a         37a         38a         39a   
        #  eigenvalues H         -0.28091    -0.15088    -0.09343    -0.07531    -0.00688
        #             eV          -7.6440     -4.1058     -2.5424     -2.0493     -0.1873
        # Skip blank line.
        line = next(inputfile)
        
        orbitals = []
        while True: 
            irreps = []
            energies_hartree = []
            energies_eV = []
            occupations = []
            
            # MO index
            line = next(inputfile)
            
            # Check we're still in the right section.
            if "irrep" not in line:
                # All done.
                break
            else:
                # Turbomole lists orbitals of different symmetry separately.
                irreps = line.split()[1:]
            
            # Energy in H.
            line = next(inputfile)
            energies_hartree = [float(energy) for energy in line.split()[2:]]
            
            # Energy in eV.
            line = next(inputfile)
            energies_eV = [float(energy) for energy in line.split()[1:]]
            
            # Occupation.
            # This line will be missing if the orbitals are virtual (unoccupied).
            line = next(inputfile)
            if "occupation" in line:
                occupations = [float(occupation) for occupation in line.split()[1:]]
                line = next(inputfile)
            
            # If we have any missing occupations, fill with 0
            occupations.extend([0.0] * (len(irreps) - len(occupations)))
                
            # Add to list.
            orbitals.extend([
                 {'irrep': irrep, 'energy_H': energy_H, 'energy_eV': energy_eV, 'occupancy': occupation}
                 for irrep, energy_H, energy_eV, occupation
                 in zip(irreps, energies_hartree, energies_eV, occupations)
            ])
            
        return orbitals, line
    
    
    def determine_homo(self, mosyms, dscf_mos):
        """
        Determine the highest occupied molecular orbital.
        
        Returns
        -------
        int
            the index of the HOMO.
        """
        # First, we need a mapping between each irrep and its index in mosyms etc.
        symm_count = {}
        irreps = []
        
        for symm in mosyms:
            try:
                symm_count[symm] += 1
            except KeyError:
                symm_count[symm] = 1
            
            irreps.append((symm_count[symm], symm))
            
        # We also have occupancy info from dscf (but probably only for those orbitals close to the HOMO/LUMO gap).
        # We now need to determine which orbital with occupancy is highest.
        homo = 0
        for mo in dscf_mos:
            if mo['occupancy'] > 0:
                # This orbital has electrons, determine its position out of all orbitals.
                index = irreps.index(self.split_irrep(mo['irrep']))
                if index > homo:
                    homo = index
                    
        return homo
                

    def deleting_modes(self, vibfreqs, vibdisps, vibirs, vibrmasses):
        """Deleting frequencies relating to translations or rotations"""
        i = 0
        while i < len(vibfreqs):
            if vibfreqs[i] == 0.0:
                # Deleting frequencies that have value 0 since they
                # do not correspond to vibrations.
                del vibfreqs[i], vibdisps[i], vibirs[i], vibrmasses[i]
                i -= 1
            i += 1

    def after_parsing(self):
        if hasattr(self, 'vibfreqs'):
            self.deleting_modes(self.vibfreqs, self.vibdisps, self.vibirs, self.vibrmasses)
            
        # Try and determine our multiplicity from our orbitals.
        if hasattr(self, "homos"):
            if len(self.homos) == 1:
                # If we are restricted (len(homos) == 1); assume we have to be singlet.
                self.set_attribute("mult", 1)
            else:
                # Unrestricted, the difference in homos should tell us the no. of unpaired e-.
                unpaired_e = abs(self.homos[0] - self.homos[1])
                self.set_attribute("mult", unpaired_e +1)
                
        # Set a flag if we stopped part way through an opt.
        if hasattr(self, "optstatus") and self.optstatus[-1] != data.ccData.OPT_DONE:
            self.optstatus[-1] += data.ccData.OPT_UNCONVERGED


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

            freqs = [utils.float(f) for f in temp[1:]]
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
            irs = [utils.float(f) for f in temp[2:]]
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


