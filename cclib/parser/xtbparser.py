import re
from itertools import groupby

from cclib.parser import logfileparser, utils

import numpy


class XTB(logfileparser.Logfile):
    def __init__(self, *args, **kwargs):
        # Call the __init__ method of the superclass
        super(XTB, self).__init__(logname="XTB", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return "XTB log file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'XTB("%s")' % (self.filename)

    def normalisesym(self, label):
        """xTB does not require normalizing symmetry labels."""
        return label

    def before_parsing(self):
        self.atomprop = {}
        self.bondprop = {}

    def after_parsing(self):
        if self.atomprop == {}:
            delattr(self, "atomprop")
        if self.bondprop == {}:
            delattr(self, "bondprop")

    def extract(self, inputfile, line):
        # Extract xtb version
        if "* xtb version" == line.strip()[:13]:
            version = line.split()[3]
            self.metadata["legacy_package_version"] = version
            self.set_attribute("xtbversion", version)

        #   -------------------------------------------------
        #  |                Calculation Setup                |
        #    -------------------------------------------------
        #
        #   program call               : <command>
        #   hostname                   : <hostname>
        #   coordinate file            : <filename>.<<filetype>
        #   omp threads                :                    12
        #   number of atoms            :                    18
        #   number of electrons        :                    66
        #   charge                     :                     0
        #   spin                       :                   0.0
        #   first test random number   :      0.87181443679343

        if "coordinate file" == line.strip()[:15]:
            self.metadata["coord_type"] = line.split()[3].split(".")[-1].lower()

        # Grab total charge
        if "charge" == line.strip()[:6]:
            charge = int(line.split()[2])
            self.set_attribute("charge", charge)

        # Multiplicity = Spin + 1
        if "spin" == line.strip()[:4]:
            spin = float(line.split()[2]) + 1
            self.set_attribute("mult", spin)

        # Grabbing SCF energies from geometry optimization steps
        #
        # ........................................................................
        # .............................. CYCLE    1 ..............................
        # ........................................................................
        #    1    -37.9590419 -0.379590E+02  0.116E-04    2.97       0.0  T
        #    2    -37.9590419 -0.674305E-11  0.982E-05    2.97     240.0  T
        #    3    -37.9590419 -0.129340E-09  0.425E-05    2.97     554.7  T
        #      SCC iter.                  ...        0 min,  0.004 sec
        #      gradient                   ...        0 min,  0.002 sec
        #  * total energy  :   -37.3373741 Eh     change       -0.1278963E-08 Eh
        #    gradient norm :     0.0004128 Eh/α   predicted     0.0000000E+00 (-100.00%)
        #    displ. norm   :     0.0058476 α      lambda       -0.2852929E-06
        #    maximum displ.:     0.0054052 α      in ANC's #1, #7, #3, ...
        #
        #   *** GEOMETRY OPTIMIZATION CONVERGED AFTER 1 ITERATIONS ***

        if "CYCLE" == line.replace(".", "").strip()[:5]:
            scfenergies = []
            while (
                line.strip()[4:41] != "GEOMETRY OPTIMIZATION CONVERGED AFTER"
                or line.strip()[4:44] != "FAILED TO CONVERGE GEOMETRY OPTIMIZATION"
            ):
                line = next(inputfile)
                if line.strip()[2:14] == "total energy":
                    energy = float(line.split()[4])
                    scfenergies.append(energy)
                if line.strip()[4:41] == "GEOMETRY OPTIMIZATION CONVERGED AFTER":
                    break
            self.set_attribute("scfenergies", scfenergies)

        # Grab the optimized geometry
        # xtb only gives the final optimized geometry
        #
        # For xyz input file - gives xyz output coords
        #
        # ================
        # final structure:
        # ================
        # 18
        # xtb: 6.3.0 (007a174)
        # C        -1.06986995496938    0.00558244661493   -0.31957205560407
        # O        -0.17255304749241   -0.97627397755570   -0.37573083147008
        # ...
        # H         1.66147402211077    1.54957926818092    0.61042104923243
        # H         5.25469642314449   -2.38603632440806    0.38612072928186
        # H         3.38728187535945   -4.19429919486669   -0.44692022446750
        #
        # For sdf\mol input file - gives sdf\mol output coords
        #
        # ================
        # final structure:
        # ================
        # energy: -37.337374171651 gnorm: 0.000298918814 xtb: 6.4.0 (d4b70c2)
        # xtb     01212213163D
        # xtb: 6.4.1 (unknown)
        # 18 19  0     0  0            999 V2000
        #   -1.0699    0.0056   -0.3196 C   0  0  0  0  0  0  0  0  0  0  0  0
        #   -0.1726   -0.9763   -0.3757 O   0  0  0  0  0  0  0  0  0  0  0  0
        #   -0.4790    1.1791    0.0854 C   0  0  0  0  0  0  0  0  0  0  0  0
        # ...
        #    1.6615    1.5496    0.6104 H   0  0  0  0  0  0  0  0  0  0  0  0
        #    5.2547   -2.3860    0.3861 H   0  0  0  0  0  0  0  0  0  0  0  0
        #    3.3873   -4.1943   -0.4469 H   0  0  0  0  0  0  0  0  0  0  0  0
        #  1  3  2  0  0  0  0
        #  2  1  1  0  0  0  0
        #  2  7  1  0  0  0  0
        # ...
        #  14 13  2  0  0  0  0
        #  15  1  1  0  0  0  0
        #  18 14  1  0  0  0  0
        # M  END
        if "final structure" == line.strip()[:15]:
            self.skip_line(inputfile, "=")

            if self.metadata["coord_type"] == "xyz":
                atomnos = []
                atomcoords = []
                for line in inputfile:
                    # Ending criteria for xyz is a blank line at the end of the coords block
                    if line == " \n":
                        break
                    if line[0].isupper():
                        atom, x, y, z = line.split()
                        atomnos.append(self.table.number[atom])
                        atomcoords.append([float(x), float(y), float(z)])
                self.set_attribute("natom", len(atomnos))
                self.set_attribute("atomnos", atomnos)
                self.set_attribute("atomcoords", atomcoords)
                # print(lines)

            elif self.metadata["coord_type"] == "sdf" or self.metadata["coord_type"] == "mol":
                atomnos = []
                atomcoords = []
                # Ending criteria for sdf\mol is the END at the end of the coord block
                while "END" != line.strip()[-3:]:
                    # Atoms block start with 3 blank spaces, bonds block starts with 1
                    if line[:3] == "   ":
                        x, y, z, atom = line.split()[:4]
                        atomnos.append(self.table.number[atom])
                        atomcoords.append([float(x), float(y), float(z)])
                    line = next(inputfile)
                self.set_attribute("natom", len(atomnos))
                self.set_attribute("atomnos", atomnos)
                self.set_attribute("atomcoords", atomcoords)

            else:
                pass

        # Get Molecular Orbitals energies and HOMO index
        # xTB trunctaes the MO list so we need to take care of that.
        # Unkown energies will be given NaN as a value
        #
        # * Orbital Energies and Occupations

        #      #    Occupation            Energy/Eh            Energy/eV
        #   -------------------------------------------------------------
        #      1        2.0000           -0.7817342             -21.2721
        #    ...           ...                  ...                  ...
        #     21        2.0000           -0.5177364             -14.0883
        #     22        2.0000           -0.5133906             -13.9701
        #     23        2.0000           -0.5119411             -13.9306
        #     24        2.0000           -0.5103339             -13.8869
        #     25        2.0000           -0.5064217             -13.7804
        #     26        2.0000           -0.4793904             -13.0449
        #     27        2.0000           -0.4762317             -12.9589
        #     28        2.0000           -0.4705819             -12.8052
        #     29        2.0000           -0.4558376             -12.4040
        #     30        2.0000           -0.4505134             -12.2591
        #     31        2.0000           -0.4390552             -11.9473
        #     32        2.0000           -0.4371482             -11.8954
        #     33        2.0000           -0.4083272             -11.1111 (HOMO)
        #     34                         -0.2990289              -8.1370 (LUMO)
        #     35                         -0.2703399              -7.3563
        #     36                         -0.2376187              -6.4659
        #     37                         -0.2246900              -6.1141
        #     38                         -0.2213822              -6.0241
        #     39                         -0.2016539              -5.4873
        #     40                         -0.1317437              -3.5849
        #     41                         -0.1173862              -3.1942
        #     42                          0.0207011               0.5633
        #     43                          0.0378419               1.0297
        #     44                          0.0843351               2.2949
        #    ...                                ...                  ...
        #     60                          1.1799189              32.1072
        #   -------------------------------------------------------------
        #               HL-Gap            0.1092983 Eh            2.9742 eV
        #          Fermi-level           -0.3536781 Eh           -9.6241 eV
        #

        if "* Orbital Energies and Occupations" == line.strip():
            # Skip 4 lines to get to the table
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)

            mooccnos = [[]]
            moenergies = [[]]
            monumbers = []

            # Ending criteria is the dashed line at the end of the MO block.
            while line.strip() != "-------------------------------------------------------------":
                line_split = line.split()
                monumbers.append(line_split[0])

                # Parsing the lines before the LUMO line and the HOMO line itself.
                # All MOs are occupied
                if (len(line_split) == 4 and line_split[-1] != "(LUMO)") or (len(line_split) == 5):
                    mooccnos[0].append(line_split[1])
                    moenergies[0].append(line_split[3])

                # For the LUMO line and after we assume 0 electrons
                # Since they are not explicit
                if len(line_split) == 3 or line_split[-1] == "(LUMO)":
                    mooccnos[0].append(0.0)
                    moenergies[0].append(line_split[2])

                # xTB gives the index of the HOMO
                # Occupation can be not an integer which complicates the way
                # to parse the HOMOs for unrestricted and openshell calculations.
                # Keeping it this way for now.
                if line_split[-1] == "(HOMO)" and not hasattr(self, "homos"):
                    self.set_attribute("homos", [int(line_split[0]) - 1])

                line = next(inputfile)

            # Find in index of the "..."
            # And fixing the type of the rest of the values
            fill_in_idx = []
            for idx, monumber in enumerate(monumbers):
                if monumber == "...":
                    fill_in_idx.append(idx)
                else:
                    monumbers[idx] = int(monumbers[idx])
                    moenergies[0][idx] = float(moenergies[0][idx])
                    mooccnos[0][idx] = float(mooccnos[0][idx])

            # Filling in the gaps since xTB truncates the list to include only
            # the first and last MOs, and some number of MOs before and
            # after the HOMO and LUMO.
            #
            # Electron occupency is assumed to be 2 for the missing values
            # before the HOMO, and 0 for the missing values after the LUMO.
            #
            # NaN is put for the missing MO energies since we have no way
            # to extrapolate them.

            missing_mos_num = 0
            for idx in fill_in_idx:
                fixed_idx = idx + missing_mos_num
                first_mo = int(monumbers[fixed_idx - 1])
                last_mo = int(monumbers[fixed_idx + 1])
                missing_mos_num = last_mo - first_mo - 1
                na_list = [numpy.nan] * missing_mos_num

                if mooccnos[0][fixed_idx - 1] == 2:
                    twos_list = [2.0] * missing_mos_num
                    mooccnos[0] = mooccnos[0][:fixed_idx] + twos_list + mooccnos[0][fixed_idx + 1 :]
                else:
                    zeros_list = [0.0] * missing_mos_num
                    mooccnos[0] = (
                        mooccnos[0][:fixed_idx] + zeros_list + mooccnos[0][fixed_idx + 1 :]
                    )
                monumbers = (
                    monumbers[:fixed_idx]
                    + list(numpy.arange(first_mo + 1, last_mo))
                    + monumbers[fixed_idx + 1 :]
                )
                moenergies[0] = moenergies[0][:fixed_idx] + na_list + moenergies[0][fixed_idx + 1 :]

                missing_mos_num -= 1

            self.set_attribute("moenergies", moenergies)

        # Grabbing atomic properties: Coordination number CN, Atomic partial charge q, Dispersion coefficient C6, Polarizability alpha:
        #
        #  #   Z          covCN         q      C6AA      α(0)
        #  1   6 C        3.056     0.087    26.011     8.364
        #  2   8 O        1.729    -0.155    16.635     5.507
        # ...
        # 17   1 H        0.927     0.059     2.210     2.325
        # 18   1 H        0.926     0.066     2.131     2.283
        #
        if "#   Z" == line.strip()[:5]:
            line = next(inputfile)
            atom_convcn = []
            atom_q = []
            atom_c6aa = []
            atom_alpha = []
            while line.strip() != "":
                line_split = line.strip().split()
                atom_convcn.append(float(line_split[3]))
                atom_q.append(float(line_split[4]))
                atom_c6aa.append(float(line_split[5]))
                atom_alpha.append(float(line_split[6]))

                line = next(inputfile)

            self.atomprop["convcn"] = atom_convcn
            self.atomprop["q"] = atom_q
            self.atomprop["c6aa"] = atom_c6aa
            self.atomprop["alpha"] = atom_alpha

        # Grabbing the Wiberg bond orders (WBO)
        #
        # Wiberg/Mayer (AO) data.
        # largest (>0.10) Wiberg bond orders for each atom
        #
        #  ---------------------------------------------------------------------------
        #      #   Z sym  total        # sym  WBO       # sym  WBO       # sym  WBO
        #  ---------------------------------------------------------------------------
        #      1   6 C    3.893 --     3 C    1.482     2 O    1.174    15 H    0.964
        #                              7 C    0.116
        #      2   8 O    2.471 --     1 C    1.174     7 C    1.096
        #      3   6 C    3.982 --     1 C    1.482     6 C    1.254     4 C    1.103
        #      4   6 C    3.951 --     5 N    2.811     3 C    1.103
        #      5   7 N    3.010 --     4 C    2.811
        #      6   6 C    3.982 --     7 C    1.505     3 C    1.254    16 H    0.970
        #      7   6 C    3.913 --     6 C    1.505     8 C    1.116     2 O    1.096
        #                              1 C    0.116
        #      8   6 C    3.919 --    10 C    1.399     9 O    1.131     7 C    1.116
        #                             14 C    0.115
        #      9   8 O    2.468 --    14 C    1.133     8 C    1.131
        #     10   6 C    3.979 --     8 C    1.399    13 C    1.240    11 C    1.119
        #     11   6 C    3.952 --    12 N    2.789    10 C    1.119
        #     12   7 N    2.997 --    11 C    2.789
        #     13   6 C    3.982 --    14 C    1.611    10 C    1.240    17 H    0.973
        #     14   6 C    3.895 --    13 C    1.611     9 O    1.133    18 H    0.964
        #                              8 C    0.115
        #     15   1 H    0.994 --     1 C    0.964
        #     16   1 H    0.996 --     6 C    0.970
        #     17   1 H    0.996 --    13 C    0.973
        #     18   1 H    0.994 --    14 C    0.964
        #  ---------------------------------------------------------------------------
        #

        if "Wiberg/Mayer (AO) data." == line.strip():
            # Skip 6 lines to get to the first line of data
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)

            wbo = []
            # Iterating over lines until the end dashes line
            while (
                line.strip()
                != "---------------------------------------------------------------------------"
            ):
                if line[5] != " ":
                    line_split = line.strip().split()
                    new_atom_index = int(line_split[0]) - 1
                    wbo_total = float(line_split[3])

                    wbo.append([wbo_total])
                    for i in range(0, len(line_split) - 5, 3):
                        wbo_next = line_split[i + 5 : i + 8]
                        wbo_idx = int(wbo_next[0]) - 1
                        wbo_order = float(wbo_next[2])

                        wbo[-1].append([wbo_idx, wbo_order])
                else:
                    line_split = line.strip().split()

                    for i in range(0, len(line_split), 3):
                        wbo_next = line_split[i : i + 3]
                        wbo_idx = int(wbo_next[0]) - 1
                        wbo_order = float(wbo_next[2])
                        wbo[-1].append([wbo_idx, wbo_order])
                line = next(inputfile)

            self.bondprop["wbo"] = wbo

        # Get Fukui Indecis for each atom
        #    #        f(+)     f(-)     f(0)
        #  1N      -0.075   -0.084   -0.080
        #  2C      -0.021   -0.014   -0.018
        #  3N      -0.084   -0.092   -0.088
        # ....................
        # 14H      -0.059   -0.057   -0.058
        # 15H      -0.070   -0.065   -0.068
        # 16H      -0.067   -0.064   -0.065
        # Creates a list of lists the size of (n,3) where n is the number of atoms.
        # Each atom had the following info, in that order:
        #   - f(+)
        #   - f(-)
        #   - f(0)

        if "#        f(+)     f(-)     f(0)" == line.strip():
            line = next(inputfile)
            atom_fp = []
            atom_fn = []
            atom_fz = []
            for i in range(self.natom):
                try:
                    atom_fp.append(float(line[9:19]))
                except:
                    atom_fp.append(-1000000.0)
                try:
                    atom_fn.append(float(line[19:28]))
                except:
                    atom_fn.append(-1000000.0)
                try:
                    atom_fz.append(float(line[28:]))
                except:
                    atom_fz.append(1000000.0)
                line = next(inputfile)

            self.atomprop["fukui"] = [atom_fp, atom_fn, atom_fz]

        # Get LMO data
        #
        #     LMO Fii/eV  ncent    charge center   contributions...
        #     1 sigma -21.33   1.84  12.23319  -5.99934  -6.73018   13O :  0.64   12C :  0.37
        #     2 sigma -21.05   1.82  11.03031  -5.52370  -7.69423   13O :  0.66   31H :  0.34
        #     3 sigma -20.30   1.94  13.15771  -3.57423  -2.57537    7N :  0.58    6C :  0.42
        #     4 sigma -20.29   1.94  13.85375  -4.46664   4.16798   11N :  0.56   10C :  0.45
        # .......................
        #    37 sigma -17.29   1.97  10.58062   1.21250  -3.64843    4C :  0.52   21H :  0.48
        #    38 sigma -17.22   1.98  12.43345   1.40582  -2.52759    4C :  0.52   22H :  0.48
        #    39 sigma -17.22   1.98  13.80996   3.97178  -5.25134    3C :  0.53   19H :  0.48
        #    40 LP    -15.87   1.05  14.66884  -4.88820   5.98676   11N :  0.98
        #    41 LP    -15.58   1.02  13.59986   8.28526  -4.78561    1N :  0.99
        #
        # Creates a list of lists size of (n,4) where n is the number of atoms.
        # Each atom has the following info, in this order:
        #   - Higest priority LMO type encoding (0 for LP, 1 for pi, 2 for delpi, 3 for sigma) - where 0 (LP) is the highest priority.
        #   - Average contribution of the highest priority LMO
        #   - Average Fii/eV of the highest priority LMO
        #   - Average ncent of the highest priority LMO

        if line.startswith(" LMO Fii/eV"):  # and 'donescf' in self.attributes.keys():
            line = next(inputfile)
            lmo_list = []
            while line[:5].strip().isnumeric():
                split = [x for x in re.split("\s+|:", line) if x != ""]
                lmo_num = split[0]
                lmo_type = split[1]
                try:
                    lmo_fii = float(split[2])
                except:
                    lmo_fii = split[2]
                try:
                    lmo_ncent = float(split[3])
                except:
                    lmo_ncent = split[3]

                lmo_cont = split[7:]
                # if lmo_type in ['pi','LP']:
                for i in range(int(len(lmo_cont) / 2)):
                    if (
                        (lmo_type == "pi" and float(lmo_cont[2 * i - 1]) > 0.3)
                        or (lmo_type == "LP" and float(lmo_cont[2 * i - 1]) > 0.7)
                        or (lmo_type in ["sigma", "delpi"])
                    ):
                        count_atoms = re.findall("(\d+|\D+)", lmo_cont[2 * i - 2])
                        lmo_list.append(
                            {
                                "AtomIdx": int(count_atoms[0]),
                                "Contribution": float(lmo_cont[2 * i - 1]),
                                "LMO Num": lmo_num,
                                "LMO Type": lmo_type,
                                "Fii/eV": lmo_fii,
                                "ncent": lmo_ncent,
                            }
                        )
                line = next(inputfile)

            LMO_ORDER = {"LP": 0, "pi": 1, "delpi": 2, "sigma": 3}
            lmo_list = sorted(lmo_list, key=lambda x: (x["AtomIdx"], LMO_ORDER[x["LMO Type"]]))
            lmo_list_cleaned = []
            keys_list = []
            for (key, key_len), group in groupby(
                lmo_list, key=lambda x: (x["AtomIdx"], LMO_ORDER[x["LMO Type"]])
            ):
                if key not in keys_list:
                    keys_list.append(key)
                    temp_list = []
                    for k in group:
                        temp_list.append(k)
                    # atom_cont = max(temp_list, key = lambda x: x['Contribution'])['Contribution']
                    atom_cont = sum(d["Contribution"] for d in temp_list) / len(temp_list)
                    # atom_fii = max(temp_list, key = lambda x: x['Fii/eV'])['Fii/eV']
                    atom_fii = sum(d["Fii/eV"] for d in temp_list) / len(temp_list)
                    # atom_ncent = max(temp_list, key = lambda x: x['ncent'])['ncent']
                    atom_ncent = sum(d["ncent"] for d in temp_list) / len(temp_list)
                    lmo_code = LMO_ORDER[temp_list[0]["LMO Type"]]
                    lmo_list_cleaned.append([lmo_code, atom_cont, atom_fii, atom_ncent])

            self.atomprop["lmo"] = lmo_list_cleaned

        # find if ended successfuly
        if "* finished run on" == line.strip()[:17]:
            self.metadata["success"] = True
        return
