# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for NBO output files"""

from cclib.parser import logfileparser


class NBO(logfileparser.Logfile):
    """A NBO log file"""

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

        """ NATURAL POPULATIONS:  Natural atomic orbital occupancies

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
            24    H  3  pz     Ryd( 2p)     0.00117       2.28805


            Summary of Natural Population Analysis:"""

        if "NAO Atom No lang   Type(AO)    Occupancy      Energy" in line:
            line = next(inputfile)
            line = next(inputfile)

            naos, atoms, nos, langs, types, occupancies, energies = [], [], [], [], [], [], []

            if not hasattr(self, "populations"):
                self.populations = dict()

            # Skip empty lines
            while "Summary of Natural Population Analysis:" not in line:
                if len(line.strip()) <= 0:
                    line = next(inputfile)
                    continue

                nao = int(line[0:5].strip())
                atom = line[8:9]
                no = int(line[9:14].strip())
                lang = line[14:21].strip()
                type_ao = line[21:29]
                occupancy = float(line[33:42])
                energy = float(line[47:56])

                naos.append(nao)
                atoms.append(atom)
                nos.append(no)
                langs.append(lang)
                types.append(type_ao)
                occupancies.append(occupancy)
                energies.append(energy)

                line = next(inputfile)

            npa_dict = {
                "nao": naos,
                "atom": atoms,
                "no": nos,
                "lang": langs,
                "type": types,
                "occupancy": occupancies,
                "energy": energies,
            }

            self.populations["npa"] = npa_dict

        """ Summary of Natural Population Analysis:

                                            Natural Population
                    Natural    ---------------------------------------------
        Atom No    Charge        Core      Valence    Rydberg      Total
        --------------------------------------------------------------------
            O  1   -0.92878      1.99998     6.91655    0.01225     8.92878
            H  2    0.46439      0.00000     0.52980    0.00581     0.53561
            H  3    0.46439      0.00000     0.52980    0.00581     0.53561
        ====================================================================
        * Total *  0.00000      1.99998     7.97616    0.02386    10.00000"""

        if "  Atom No    Charge" in line:
            if not hasattr(self, "atomcharges"):
                self.atomcharges = dict()

            line = next(inputfile)
            line = next(inputfile)

            charges = []

            while "==============" not in line:
                population_analysis = line.split()

                atom = population_analysis[0]
                no = int(population_analysis[1])
                natural_charge = float(population_analysis[2])
                core = float(population_analysis[3])
                valence = float(population_analysis[4])
                rydberg = float(population_analysis[5])  # noqa: F841
                total = float(population_analysis[6])  # noqa: F841

                # TODO append to attibutes
                charges.append(natural_charge)

                line = next(inputfile)

            self.atomcharges["nbo"] = charges

            if not hasattr(self, "natom"):
                self.set_attribute("natom", len(self.atomcharges["nbo"]))

        #                                  Natural Population
        #  ---------------------------------------------------------
        #    Core                       1.99998 ( 99.9990% of    2)
        #    Valence                    7.97616 ( 99.7019% of    8)
        #    Natural Minimal Basis      9.97614 ( 99.7614% of   10)
        #    Natural Rydberg Basis      0.02386 (  0.2386% of   10)
        #  ---------------------------------------------------------

        if line[33:51] == "Natural Population":
            line = next(inputfile)
            line = next(inputfile)

            core = float(line.split()[1])  # noqa: F841
            # TODO append to attibutes

            line = next(inputfile)

            valence = float(line.split()[1])  # noqa: F841
            # TODO append to attibutes

            line = next(inputfile)

            natural_minimal_basis = float(line.split()[3])  # noqa: F841
            # TODO append to attibutes

            line = next(inputfile)

            natural_rydberg_basis = float(line.split()[3])  # noqa: F841
            # TODO append to attibutes

        #     Atom No         Natural Electron Configuration
        #  ----------------------------------------------------------------------------
        #       O  1      [core]2s( 1.75)2p( 5.17)3d( 0.01)
        #       H  2            1s( 0.53)
        #       H  3            1s( 0.53)

        if "Natural Electron Configuration" in line:
            line = next(inputfile)
            line = next(inputfile)

            while len(line.strip()):
                configuration_line = line.split()
                atom = configuration_line[0]
                configuration = "".join(configuration_line[2:])  # noqa: F841

                # TODO append to attibutes

                line = next(inputfile)

        # NATURAL BOND ORBITAL ANALYSIS:

        #                             Occupancies       Lewis Structure    Low   High
        #          Max    Occ     -------------------  -----------------   occ   occ
        #   Cycle  Ctr   Thresh    Lewis   non-Lewis     CR  BD  nC  LP    (L)   (NL)
        #  ============================================================================
        #     1     2     1.90     9.99255   0.00745      1   2   0   2     0      0
        #  ----------------------------------------------------------------------------

        if "NATURAL BOND ORBITAL ANALYSIS" in line:
            # Skip to the values
            for _ in range(6):
                line = next(inputfile)

            nbo_analysis = line.split()

            cycle = int(nbo_analysis[0])  # noqa: F841
            max_ctr = int(nbo_analysis[1])  # noqa: F841
            occ_thresh = float(nbo_analysis[2])  # noqa: F841
            occ_lewis = float(nbo_analysis[3])  # noqa: F841
            occ_non_lewis = float(nbo_analysis[4])  # noqa: F841
            lewis_cr = float(nbo_analysis[5])  # noqa: F841
            lewis_bd = int(nbo_analysis[6])  # noqa: F841
            lewis_nc = int(nbo_analysis[7])  # noqa: F841
            lewis_lp = int(nbo_analysis[8])  # noqa: F841
            low_occ = int(nbo_analysis[9])  # noqa: F841
            high_occ = int(nbo_analysis[10])  # noqa: F841

            # TODO append to attibutes

        #         NATURAL BOND ORBITALS (Summary):

        #                                                      Principal Delocalizations
        #            NBO                 Occupancy    Energy   (geminal,vicinal,remote)
        #  ===============================================================================
        #  Molecular unit  1  (H2O)
        #  ------ Lewis --------------------------------------
        #     1. CR ( 1) O  1             1.99998   -20.54367
        #     2. LP ( 1) O  1             1.99764    -0.49214  18(v),22(v)
        #     3. LP ( 2) O  1             1.99711    -0.76980  17(v),21(v)
        #     4. BD ( 1) O  1- H  2       1.99891    -0.89832  23(v)
        #     5. BD ( 1) O  1- H  3       1.99891    -0.89832  19(v)
        #  ------ non-Lewis ----------------------------------
        #     6. BD*( 1) O  1- H  2       0.00038     0.72889
        #     7. BD*( 1) O  1- H  3       0.00038     0.72889
        #     8. RY ( 1) O  1             0.00001     1.39675
        #     9. RY ( 2) O  1             0.00000     1.92591
        #    10. RY ( 3) O  1             0.00000     2.97200
        #    11. RY ( 4) O  1             0.00000     3.08380
        #    12. RY ( 5) O  1             0.00000     3.11162
        #    13. RY ( 6) O  1             0.00000     3.15738
        #    14. RY ( 7) O  1             0.00000     3.29348
        #    15. RY ( 8) O  1             0.00000     1.20441
        #    16. RY ( 9) O  1             0.00000     1.51316
        #    17. RY ( 1) H  2             0.00152     0.74247
        #    18. RY ( 2) H  2             0.00118     1.94640
        #    19. RY ( 3) H  2             0.00064     1.74999
        #    20. RY ( 4) H  2             0.00000     2.82233
        #    21. RY ( 1) H  3             0.00152     0.74247
        #    22. RY ( 2) H  3             0.00118     1.94640
        #    23. RY ( 3) H  3             0.00064     1.74999
        #    24. RY ( 4) H  3             0.00000     2.82233
        #           -------------------------------
        #                  Total Lewis    9.99255  ( 99.9255%)
        #            Valence non-Lewis    0.00075  (  0.0075%)
        #            Rydberg non-Lewis    0.00669  (  0.0669%)
        #           -------------------------------
        #                Total unit  1   10.00000  (100.0000%)
        #               Charge unit  1    0.00000

        if "Principal Delocalizations" in line:
            if not hasattr(self, "nbo"):
                self.nbo = []

            while " Lewis " not in line:
                line = next(inputfile)

            line = next(inputfile)

            while "   ---" not in line:
                if "-----" in line or len(line[7:28].strip()) < 1:
                    line = next(inputfile)
                    continue

                nao = line[7:28].strip()
                occupancy = float(line[30:40].strip())
                energy = float(line[40:52].strip())

                # TODO
                # geminal, vicinal, remote = None, None, None

                # if len(line) > 52: geminal = line[53:58]
                # if len(line) > 58: vicinal = line[59:64]
                # if len(line) > 62: remote  = line[65:70]

                nbo_dict = {
                    "nao": nao,
                    "occupancy": occupancy,
                    "energy": energy,
                    # TODO
                    # 'delocalizations': {
                    #     'geminal': geminal,
                    #     'vicinal': vicinal,
                    #     'remote' : remote
                    # }
                }

                self.append_attribute("nbo", nbo_dict)

                line = next(inputfile)
