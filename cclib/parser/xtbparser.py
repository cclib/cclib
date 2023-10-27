from datetime import timedelta
from typing import List, Literal, Optional, Tuple

from cclib.parser import logfileparser, utils
from cclib.parser.logfilewrapper import FileWrapper

import numpy as np


class XTB(logfileparser.Logfile):
    """An output parser for the xTB code"""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="xTB", *args, **kwargs)

    def __str__(self) -> str:
        """Return a string representation of the object."""
        return f"xTB log file {self.filename}"

    def __repr__(self) -> str:
        """Return a representation of the object."""
        return f'xTB("{self.filename}")'

    def normalisesym(self, label):
        """xTB does not require normalizing symmetry labels."""
        return label

    def before_parsing(self) -> None:
        """Set attributes before parsing"""
        self.atomprop = {}
        self.bondprop = {}

    def after_parsing(self) -> None:
        """Delete empty attributes after parsing"""
        if not self.atomprop:
            delattr(self, "atomprop")
        if not self.bondprop:
            delattr(self, "bondprop")

    def _extract_version(self, line: str) -> Optional[str]:
        """
        Extract xtb version from the following:

        -----------------------------------------------------------
        |                   =====================                   |
        |                           x T B                           |
        |                   =====================                   |
        |                         S. Grimme                         |
        |          Mulliken Center for Theoretical Chemistry        |
        |                    University of Bonn                     |
        -----------------------------------------------------------

        * xtb version 6.6.1 (8d0f1dd) compiled by 'conda@1efc2f54142f' on 2023-08-01
        """

        return line.split()[3] if "xtb version" in line else None

    def _extract_coord_file(self, line: str) -> Optional[str]:
        """
        Extract the coordinate filename, from which we can strip out the type

        -------------------------------------------------
        |                Calculation Setup                |
        -------------------------------------------------

        program call               : xtb coord.xyz --opt
        coordinate file            : coord.xyz
        omp threads                :                    20
        """
        return line.split()[-1] if "coordinate file" in line else None

    def _extract_charge(self, line: str) -> Optional[int]:
        """
        Extract the total charge. It can always be found in the
        summary formatted as a float.

        Format for GFN-xTB:

         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         ::                     SUMMARY                     ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         :: total energy              -5.070544440612 Eh    ::
         :: gradient norm              0.000057326562 Eh/a0 ::
         :: HOMO-LUMO gap             14.391809984508 eV    ::
         ::.................................................::
         :: SCC energy                -5.104920280363 Eh    ::
         :: -> isotropic ES            0.031458595179 Eh    ::
         :: -> anisotropic ES          0.000396760551 Eh    ::
         :: -> anisotropic XC         -0.000881430881 Eh    ::
         :: -> dispersion             -0.000141085082 Eh    ::
         :: repulsion energy           0.034375839725 Eh    ::
         :: add. restraining           0.000000000000 Eh    ::
         :: total charge              -0.000000000000 e     ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::

         Format for GFN-FF:
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         ::                     SUMMARY                     ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         :: total energy              -0.327405667208 Eh    ::
         :: gradient norm              0.021455947754 Eh/a0 ::
         ::.................................................::
         :: bond energy               -0.269813650323 Eh    ::
         :: angle energy               0.001167514553 Eh    ::
         :: torsion energy             0.000000000000 Eh    ::
         :: repulsion energy           0.030786304642 Eh    ::
         :: electrostat energy        -0.089406523183 Eh    ::
         :: dispersion energy         -0.000139312897 Eh    ::
         :: HB energy                  0.000000000000 Eh    ::
         :: XB energy                  0.000000000000 Eh    ::
         :: bonded atm energy          0.000000000000 Eh    ::
         :: external energy            0.000000000000 Eh    ::
         :: add. restraining           0.000000000000 Eh    ::
         :: total charge              -0.000000000000 e     ::
         ::.................................................::
         :: atomisation energy         0.000000000000 Eh    ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
        """
        return round(float(line.split()[-3])) if "total charge" in line else None

    def _extract_dispersion(self, line: str) -> Optional[float]:
        """
        Extraction the dispersion energy. Refer to the summary tables above.
        """
        return (
            utils.convertor(float(line.split()[-3]), "hartree", "eV")
            if "-> dispersion" in line or "dispersion energy" in line
            else None
        )

    def _extract_final_energy(self, line: str) -> Optional[float]:
        """
        Extract the final total energy from the result table.

           -------------------------------------------------
          | TOTAL ENERGY               -5.070544323569 Eh   |
          | GRADIENT NORM               0.000458081396 Eh/α |
          | HOMO-LUMO GAP              14.381252816459 eV   |
           -------------------------------------------------
        """
        return (
            utils.convertor(float(line.split()[3]), "hartree", "eV")
            if "TOTAL ENERGY" in line
            else None
        )

    def _extract_geom_energy(self, line: str) -> Optional[float]:
        """
                Extract the energies for a geometry step.

        ........................................................................
        .............................. CYCLE    1 ..............................
        ........................................................................

                 iter      E             dE          RMSdq      gap      omega  full diag
                   1     -5.1048382 -0.510484E+01  0.417E-06   14.38       0.0  T
                   2     -5.1048382  0.000000E+00  0.234E-06   14.38   24706.5  T
                   3     -5.1048382  0.000000E+00  0.437E-07   14.38  100000.0  T
                     SCC iter.                  ...        0 min,  0.005 sec
                     gradient                   ...        0 min,  0.010 sec
                 * total energy  :    -5.0705443 Eh     change       -0.4369838E-12 Eh
                   gradient norm :     0.0004582 Eh/α   predicted     0.0000000E+00 (-100.00%)
                   displ. norm   :     0.0005728 α      lambda       -0.1688374E-06
                   maximum displ.:     0.0005029 α      in ANC's #2, #3, #1, ...

                   *** GEOMETRY OPTIMIZATION CONVERGED AFTER 1 ITERATIONS ***
        """

        return (
            utils.convertor(float(line.split()[4]), "hartree", "eV")
            if "* total energy" in line
            else None
        )

    def _extract_symbol_coords(
        self, line: str, mode: Literal["xyz", "mol", "sdf"]
    ) -> Optional[Tuple[str, List[float, float, float]]]:
        """
        Extract the symbol and X, Y, Z coordinates.

        For XYZ files:

        ================
        final structure:
        ================
        3
        xtb: 6.6.1 (8d0f1dd)
        O            1.07015391331798       -0.01769828395654        0.04981203402603
        H            2.02952514441169       -0.00813780275851        0.03338237327689
        H            0.76845094227033        0.44031608671506       -0.73728440730292

        For SDF/mol files:
        ================
        final structure:
        ================
        energy: -37.337374171651 gnorm: 0.000298918814 xtb: 6.4.0 (d4b70c2)
        xtb     01212213163D
        xtb: 6.4.1 (unknown)
        18 19  0     0  0            999 V2000
          -1.0699    0.0056   -0.3196 C   0  0  0  0  0  0  0  0  0  0  0  0
          -0.1726   -0.9763   -0.3757 O   0  0  0  0  0  0  0  0  0  0  0  0
          -0.4790    1.1791    0.0854 C   0  0  0  0  0  0  0  0  0  0  0  0
        ...
           1.6615    1.5496    0.6104 H   0  0  0  0  0  0  0  0  0  0  0  0
           5.2547   -2.3860    0.3861 H   0  0  0  0  0  0  0  0  0  0  0  0
           3.3873   -4.1943   -0.4469 H   0  0  0  0  0  0  0  0  0  0  0  0
         1  3  2  0  0  0  0
         2  1  1  0  0  0  0
         2  7  1  0  0  0  0
        ...
         14 13  2  0  0  0  0
         15  1  1  0  0  0  0
         18 14  1  0  0  0  0
        M  END
        """
        line_split = line.split()
        if mode == "xyz":
            if line_split.isupper():
                return line_split[0], [float(coord) for coord in line_split[1:]]
        elif mode in {"mol", "sdf"}:
            if line_split[3].isupper():
                return line_split[3], [float(coord) for coord in line_split[:3]]
        else:
            raise ValueError(f"Unsupported coordinate file type: {mode}")

    def _extract_orbitals(self, line: str) -> Optional[Tuple[int, float, float, bool]]:
        """
        Extract the orbital index, occupation, energy, and if it's the HOMO.

        * Orbital Energies and Occupations

             #    Occupation            Energy/Eh            Energy/eV
          -------------------------------------------------------------
             1        2.0000           -0.7817342             -21.2721
           ...           ...                  ...                  ...
            21        2.0000           -0.5177364             -14.0883
            22        2.0000           -0.5133906             -13.9701
            23        2.0000           -0.5119411             -13.9306
            24        2.0000           -0.5103339             -13.8869
            25        2.0000           -0.5064217             -13.7804
            26        2.0000           -0.4793904             -13.0449
            27        2.0000           -0.4762317             -12.9589
            28        2.0000           -0.4705819             -12.8052
            29        2.0000           -0.4558376             -12.4040
            30        2.0000           -0.4505134             -12.2591
            31        2.0000           -0.4390552             -11.9473
            32        2.0000           -0.4371482             -11.8954
            33        2.0000           -0.4083272             -11.1111 (HOMO)
            34                         -0.2990289              -8.1370 (LUMO)
            35                         -0.2703399              -7.3563
            36                         -0.2376187              -6.4659
            37                         -0.2246900              -6.1141
            38                         -0.2213822              -6.0241
            39                         -0.2016539              -5.4873
            40                         -0.1317437              -3.5849
            41                         -0.1173862              -3.1942
            42                          0.0207011               0.5633
            43                          0.0378419               1.0297
            44                          0.0843351               2.2949
           ...                                ...                  ...
            60                          1.1799189              32.1072
          -------------------------------------------------------------
                      HL-Gap            0.1092983 Eh            2.9742 eV
                 Fermi-level           -0.3536781 Eh           -9.6241 eV
        """
        line_split = line.split()
        if "(HOMO)" in line or (len(line_split) == 4 and "MO" not in line):
            return (
                int(line_split[0]),
                float(line_split[1]),
                float(line_split[3]),
                bool("(HOMO) in line"),
            )
        if "(LUMO)" in line or (len(line_split) == 3 and "MO" not in line):
            return int(line_split[0]), 0.0, float(line_split[3]), False

    def _extract_mulliken_charge(self, line: str) -> Optional[float]:
        """
        Extract Mulliken charge.

        #   Z          covCN         q      C6AA      α(0)
        1   8 O        1.611    -0.565    24.356     6.661
        2   1 H        0.805     0.282     0.777     1.384
        3   1 H        0.805     0.282     0.777     1.384
        """
        line_split = line.split()
        return float(line_split[4]) if len(line) == 6 else None

    def _extract_wall_time(self, line: str) -> Optional[List[timedelta]]:
        """
        Extract the wall time.

         * wall-time:     0 d,  0 h,  0 min,  0.091 sec
        """
        line_split = line.split()
        return (
            [
                timedelta(
                    days=float(line_split[2]),
                    hours=float(line_split[4]),
                    minutes=float(line_split[6]),
                    seconds=float(line_split[8]),
                )
            ]
            if "*" in line and "wall-time" in line
            else None
        )

    def _extract_cpu_time(self, line: str) -> Optional[List[timedelta]]:
        """
        Extract the CPU time.

         *  cpu-time:     0 d,  0 h,  0 min,  1.788 sec
        """
        line_split = line.split()
        return (
            [
                timedelta(
                    days=float(line_split[2]),
                    hours=float(line_split[4]),
                    minutes=float(line_split[6]),
                    seconds=float(line_split[8]),
                )
            ]
            if "*" in line and "cpu-time" in line
            else None
        )

    # TODO: Get other headers.
    def _extract_method(self, line: str) -> Optional[str]:
        """
        Extract the method.

        Example 1:
           -------------------------------------------------
          |                 G F N 2 - x T B                 |
           -------------------------------------------------

        Example 2:
           -------------------------------------------------
          |                   G F N - F F                   |
          |          A general generic force-field          |
          |                  Version 1.0.3                  |
           -------------------------------------------------
        """
        return (
            "GFN2-xTB" if "G F N 2 - x T B" in line else "GFN-FF" if "G F N - F F" in line else None
        )

    def _is_cycle_line(self, line: str) -> bool:
        """
        Determine if the line indicates it is a geometry optimization cycle.

        ........................................................................
        .............................. CYCLE    1 ..............................
        ........................................................................
        """

        return "CYCLE" in line

    def _is_geom_end_line(self, line: str) -> bool:
        """
        Determine if the line indicates the optimization is over.

        *** GEOMETRY OPTIMIZATION CONVERGED AFTER 1 ITERATIONS ***

        or

        *** FAILED TO CONVERGE GEOMETRY OPTIMIZATION ***
        """
        return "***" in line and "GEOMETRY OPTIMIZATION" in line

    def _is_success(self, line: str) -> bool:
        """
        Determine if the job finished.

        ------------------------------------------------------------------------
        * finished run on 2023/10/26 at 13:18:28.705
        ------------------------------------------------------------------------
        """
        return "* finished run" in line

    def _is_end_of_structure_block(self, line: str, mode: Literal["xyz", "mol", "sdf"]) -> bool:
        """
        Determine if the line indicates the end of a structure block.
        Refer to _extract_symbol_coords for examples of structure blocks.
        """
        if mode == "xyz":
            return line == "\n"
        elif mode in {"mol", "sdf"}:
            return "M" in line and "END" in line
        else:
            raise ValueError(f"Unsupported coordinate file type: {mode}")

    def _is_orbitals_line(self, line: str) -> bool:
        """
        Determine if the line indicates the start of a orbital block.

        * Orbital Energies and Occupations

             #    Occupation            Energy/Eh            Energy/eV
          -------------------------------------------------------------
             1        2.0000           -0.7817342             -21.2721
           ...           ...                  ...                  ...
        """
        return "* Orbital Energies and Occupations" in line

    def _is_atomwise_properties(self, line: str) -> bool:
        """
        Determine if there are atom-wise properties.

        #   Z          covCN         q      C6AA      α(0)
        """
        return line.split()[0:4] == ["#", "Z", "covCN", "q"]

    def _is_geom_opt_converged(self, line: str) -> bool:
        """
        Determine if the geometry optimization is converged.

        *** GEOMETRY OPTIMIZATION CONVERGED AFTER 1 ITERATIONS ***
        """
        return "***" in line and "GEOMETRY OPTIMIZATION CONVERGED" in line

    def extract(self, inputfile: FileWrapper, line: str) -> None:
        # Initialize as False. Will be overwritten to True if/when appropriate.
        self.metadata["success"] = False
        self.optstatus = False
        self.optdone = False

        if version := self._extract_version(line):
            self.metadata["package_version"] = version

        if method := self._extract_method(line):
            self.metadata["methods"] = [method]

        if coord_filename := self._extract_coord_file(line):
            self.metadata["coord_type"] = coord_filename.split(".")[-1]

        if charge := self._extract_charge(line):
            self.charge = charge

        # Cycle through the geometry steps to get the energies
        scf_energies = []
        dispersion_energies = []
        if self._is_cycle_line(line):
            while not self._is_geom_end_line(line):
                if scf_energy := self._extract_geom_energy(line):
                    scf_energies.append(scf_energy)
                if dispersion_energy := self._extract_dispersion(line):
                    dispersion_energies.append(dispersion_energy)
                line = next(inputfile)

            # Check if the geometry optimization was successful
            if self._is_geom_opt_converged(line):
                self.optdone = True

        # Get the final geometry
        if "final structure:" in line:
            atomnos = []
            atomcoords = []
            coord_type = self.metadata.get("coord_type")

            while not self._is_end_of_structure_block(line, coord_type):
                if symbol_coords := self._extract_symbol_coords(line, coord_type):
                    symbol, coords = symbol_coords
                    atomnos.append(self.table.number[symbol])
                    atomcoords.append(coords)
                line = next(inputfile)

            if atomnos:
                self.natom = len(atomnos)
                self.atomnos = atomnos
            if atomcoords:
                self.atomcoords = atomcoords

        # TODO: natom and atomnos are not defined above for static calculations
        # but are often reported in the log file.

        if self._is_orbitals_line(line):
            # Skip 4 lines to get to the table
            next(inputfile)
            next(inputfile)
            next(inputfile)
            line = next(inputfile)

            homos = []
            mooccnos = []
            moenergies = []

            i = 0
            while "------" not in line:
                orbital_info = self._extract_orbitals(line)

                if orbital_info:
                    orbital_idx = orbital_info[0] - 1
                    orbital_occ = orbital_info[1]
                    orbital_energy = orbital_info[2]
                    is_homo = orbital_info[3]

                    # Account for "..." offset
                    offset = orbital_idx - len(moenergies)
                    if offset != 0:
                        moenergies.extend([np.nan] * offset)
                        mooccnos.extend([np.nan] * offset)

                    moenergies.append(orbital_energy)
                    mooccnos.append(orbital_occ)
                    if is_homo:
                        homos.append(orbital_idx)

                i += 1
                line = next(inputfile)

            # TODO: Unrestricted calculations
            if moenergies:
                self.moenergies = [np.array(moenergies)]
            if mooccnos:
                self.mooccnos = np.array(mooccnos)
            if homos:
                self.homos = np.array(homos)

        # TODO: Is "q" mulliken in xTB? Not specified in output...
        mullikens = []
        if self._is_atomwise_properties(line):
            line = next(inputfile)
            while line != "\n":
                q = self._extract_mulliken_charge(line)
                if q:
                    mullikens.append(q)
                line = next(inputfile)

        if mullikens:
            self.atomcharges = {"mulliken": np.array(mullikens)}

        if final_energy := self._extract_final_energy(line):
            if scf_energies:
                # Patch the final total energy to be the last SCF energy
                # since it is higher precision and also always available
                scf_energies[-1] = final_energy
            else:
                # We only have the final energy to store
                scf_energies = [final_energy]

        if scf_energies:
            self.set_attribute("scfenergies", np.array(scf_energies))

        if dispersion_energies:
            self.set_attribute("dispersionenergies", np.array(dispersion_energies))

        if wall_time := self._extract_wall_time(line):
            self.metadata["wall_time"] = wall_time

        if cpu_time := self._extract_cpu_time(line):
            self.metadata["cpu_time"] = cpu_time

        if self._is_success(line):
            self.metadata["success"] = True

        # TODO: vibfreqs etc.
        # TODO: enthalpy, entropy, freeenergy
